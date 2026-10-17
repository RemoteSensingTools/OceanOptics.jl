# =============================================================================
# Tabulated Greek-coefficient phase function
# =============================================================================
#
# Neutral input route for phase (matrix) expansions computed elsewhere —
# published benchmark tables, Mie codes, fits to measured VSFs — given in
# this package's convention: p(μ) = Σ β_ℓ P_ℓ(μ) with β₀ = 1, the
# Sanghavi (2014) B-matrix layout and Hovenier-positive γ (see
# `phase_matrix_moments`).
#
# Coefficients in other conventions (e.g. de Rooij–van der Stap /
# Mishchenko α₁…α₄, β₁, β₂) must be mapped before construction; the sign
# flips known for VLIDORT-format data are tabulated in vSmartMOM's
# docs/src/pages/conventions.md. Angular Mueller tables (F₁₁(Θ), F₁₂(Θ), …)
# need a reviewed generalized-spherical-function expansion first; this
# type deliberately accepts coefficients only.
# =============================================================================

"""
$(TYPEDEF)

Phase function, or full phase matrix, given directly by its Greek
expansion coefficients.

# Fields
$(TYPEDFIELDS)

# Construction

    TabulatedGreekPhase(β; α, γ, δ, ϵ, ζ, source = "", normalize = false)
    TabulatedGreekPhase(λ, β; α, γ, δ, ϵ, ζ, source = "", normalize = false)

The first form is wavelength independent: `β` (and optionally the other
five) are vectors indexed `[ℓ+1]`. The second takes a strictly increasing
wavelength grid `λ` `[nm]` and `(L+1, nλ)` matrices indexed `[ℓ+1, iλ]`;
queries at any other wavelength throw — no interpolation is performed.

Pass none of `α, γ, δ, ϵ, ζ` (scalar-only, `is_polarizable(pf) == false`)
or all five. Every column must satisfy `β₀ = 1` within `√eps(FT)`,
`|β_ℓ| ≤ 2ℓ+1`, and give a backscatter fraction in `[0, 1]` (a table
truncated well before convergence can fail this); `α, γ, ϵ, ζ` must
vanish at `ℓ < 2`. Published tables rounded to a few digits can be rescaled
with `normalize = true`, which divides all six coefficients of each column
by its `β₀`; this is an explicit choice, never applied silently.

Requests for `ℓ_max < L` truncate the expansion (as for every phase
function in this package); requests for `ℓ_max > L` pad with zeros,
which is exact for the tabulated expansion.

# Example

```julia
# Rayleigh scattering without depolarization, as a full phase matrix
β = [1.0, 0.0, 0.5];  α = [0.0, 0.0, 3.0];  γ = [0.0, 0.0, sqrt(1.5)]
δ = [0.0, 1.5, 0.0];  ϵ = zeros(3);         ζ = zeros(3)
pf = TabulatedGreekPhase(β; α, γ, δ, ϵ, ζ, source = "Hansen & Travis (1974)")
phase_matrix_moments(pf, 8).γ[3]   # √(3/2)
```
"""
struct TabulatedGreekPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Wavelength grid `[nm]`; empty for a wavelength-independent table"
    λ::Vector{FT}
    "β_ℓ indexed `[ℓ+1, iλ]` (a single column when wavelength independent)"
    β::Matrix{FT}
    "α_ℓ; `0×0` when scalar-only"
    α::Matrix{FT}
    "γ_ℓ (Hovenier sign, positive for Rayleigh); `0×0` when scalar-only"
    γ::Matrix{FT}
    "δ_ℓ; `0×0` when scalar-only"
    δ::Matrix{FT}
    "ϵ_ℓ; `0×0` when scalar-only"
    ϵ::Matrix{FT}
    "ζ_ℓ; `0×0` when scalar-only"
    ζ::Matrix{FT}
    "Provenance: citation, table, and any transformation applied"
    source::String

    function TabulatedGreekPhase{FT}(λ, β, α, γ, δ, ϵ, ζ, source) where {FT<:AbstractFloat}
        pf = new{FT}(λ, β, α, γ, δ, ϵ, ζ, source)
        _validate_tabulated_greek(pf)
        return pf
    end
end

function TabulatedGreekPhase(β::AbstractVector; α = nothing, γ = nothing, δ = nothing,
                             ϵ = nothing, ζ = nothing, source::AbstractString = "",
                             normalize::Bool = false)
    FT  = float(eltype(β))
    col = x -> x === nothing ? nothing : reshape(collect(FT, x), :, 1)
    return _tabulated_greek_phase(FT[], col(β), col(α), col(γ), col(δ), col(ϵ), col(ζ),
                                  source, normalize)
end

function TabulatedGreekPhase(λ::AbstractVector, β::AbstractMatrix; α = nothing,
                             γ = nothing, δ = nothing, ϵ = nothing, ζ = nothing,
                             source::AbstractString = "", normalize::Bool = false)
    FT  = float(eltype(β))
    mat = x -> x === nothing ? nothing : Matrix{FT}(x)
    return _tabulated_greek_phase(collect(FT, λ), mat(β), mat(α), mat(γ), mat(δ),
                                  mat(ϵ), mat(ζ), source, normalize)
end

# Keyword front end: fill omitted polarized coefficients with 0×0 matrices
# and apply the explicit `normalize` option; validation is in the inner
# constructor.
function _tabulated_greek_phase(λ::Vector{FT}, β::Matrix{FT}, α, γ, δ, ϵ, ζ,
                                source, normalize) where {FT}
    greek = map(m -> m === nothing ? Matrix{FT}(undef, 0, 0) : m,
                (α = α, γ = γ, δ = δ, ϵ = ϵ, ζ = ζ))
    if normalize
        β₀ = β[1:1, :]
        all(x -> isfinite(x) && x > zero(FT), β₀) ||
            throw(ArgumentError("TabulatedGreekPhase: normalize = true needs finite, " *
                                "positive β₀"))
        β     = β ./ β₀
        greek = map(m -> isempty(m) ? m : m ./ β₀, greek)
    end
    return TabulatedGreekPhase{FT}(λ, β, greek.α, greek.γ, greek.δ, greek.ϵ, greek.ζ,
                                   String(source))
end

function _validate_tabulated_greek(pf::TabulatedGreekPhase)
    who = "TabulatedGreekPhase"
    (; λ, β) = pf
    size(β, 1) ≥ 1 || throw(ArgumentError("$who: β needs at least the ℓ = 0 row"))
    if isempty(λ)
        size(β, 2) == 1 ||
            throw(DimensionMismatch("$who: a wavelength-independent table has one " *
                                    "column, got $(size(β, 2))"))
    else
        _check_exact_λ_grid(λ, who)
        size(β, 2) == length(λ) ||
            throw(DimensionMismatch("$who: β must have one column per wavelength " *
                                    "($(length(λ))), got $(size(β, 2))"))
    end
    greek = (α = pf.α, γ = pf.γ, δ = pf.δ, ϵ = pf.ϵ, ζ = pf.ζ)
    given = [name for (name, m) in pairs(greek) if !isempty(m)]
    if !isempty(given)
        length(given) == length(greek) ||
            throw(ArgumentError("$who: pass all of α, γ, δ, ϵ, ζ for a phase matrix or " *
                                "none for a scalar phase function; got only " *
                                join(given, ", ")))
        for (name, m) in pairs(greek)
            size(m) == size(β) ||
                throw(DimensionMismatch("$who: $name has size $(size(m)), β has size " *
                                        "$(size(β))"))
            all(isfinite, m) ||
                throw(ArgumentError("$who: $name contains non-finite values"))
        end
        _validate_low_order_greek(greek, who)
    end
    _validate_beta_moments(β, who)
    for j in axes(β, 2)
        B = _backscatter_fraction_from_moments(view(β, :, j))
        zero(B) ≤ B ≤ one(B) ||
            throw(ArgumentError("$who: column $j has backscatter fraction $B ∉ [0, 1]; " *
                                "the expansion is not a valid phase function (table " *
                                "truncated before convergence?)"))
    end
    return nothing
end

is_polarizable(pf::TabulatedGreekPhase)          = !isempty(pf.α)
is_wavelength_dependent(pf::TabulatedGreekPhase) = !isempty(pf.λ)

# Column of the table that applies at λ (or everywhere, if λ-independent).
_greek_column(pf::TabulatedGreekPhase, λ::Real) =
    is_wavelength_dependent(pf) ? _exact_λ_index(pf.λ, λ, "TabulatedGreekPhase") : 1

function _greek_column(pf::TabulatedGreekPhase)
    is_wavelength_dependent(pf) &&
        throw(ArgumentError("TabulatedGreekPhase: this table is wavelength dependent; " *
                            "query it with an explicit wavelength, e.g. " *
                            "`phase_function_moments(pf, ℓ_max, λ)`"))
    return 1
end

# Moments ℓ = 0:ℓ_max of one column: truncated, or zero-padded beyond the table.
function _tabulated_moments(M::Matrix{FT}, j::Integer, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    out = zeros(FT, ℓ_max + 1)
    n   = min(ℓ_max + 1, size(M, 1))
    out[1:n] .= view(M, 1:n, j)
    return out
end

function _tabulated_matrix_moments(pf::TabulatedGreekPhase, j::Integer, ℓ_max::Integer)
    is_polarizable(pf) || throw(_scalar_only_error(pf))
    return (α = _tabulated_moments(pf.α, j, ℓ_max), β = _tabulated_moments(pf.β, j, ℓ_max),
            γ = _tabulated_moments(pf.γ, j, ℓ_max), δ = _tabulated_moments(pf.δ, j, ℓ_max),
            ϵ = _tabulated_moments(pf.ϵ, j, ℓ_max), ζ = _tabulated_moments(pf.ζ, j, ℓ_max))
end

phase_function_moments(pf::TabulatedGreekPhase, ℓ_max::Integer) =
    _tabulated_moments(pf.β, _greek_column(pf), ℓ_max)
phase_function_moments(pf::TabulatedGreekPhase, ℓ_max::Integer, λ::Real) =
    _tabulated_moments(pf.β, _greek_column(pf, λ), ℓ_max)

phase_matrix_moments(pf::TabulatedGreekPhase, ℓ_max::Integer) =
    _tabulated_matrix_moments(pf, _greek_column(pf), ℓ_max)
phase_matrix_moments(pf::TabulatedGreekPhase, ℓ_max::Integer, λ::Real) =
    _tabulated_matrix_moments(pf, _greek_column(pf, λ), ℓ_max)

phase_function_value(pf::TabulatedGreekPhase{FT}, μ::Real) where {FT} =
    _reconstruct_from_moments(view(pf.β, :, _greek_column(pf)), FT(μ))

backscatter_fraction(pf::TabulatedGreekPhase) =
    _backscatter_fraction_from_moments(view(pf.β, :, _greek_column(pf)))
backscatter_fraction(pf::TabulatedGreekPhase, λ::Real) =
    _backscatter_fraction_from_moments(view(pf.β, :, _greek_column(pf, λ)))

function asymmetry_parameter(pf::TabulatedGreekPhase{FT}) where {FT}
    size(pf.β, 1) ≥ 2 || return zero(FT)
    return pf.β[2, _greek_column(pf)] / FT(3)
end

function Base.show(io::IO, pf::TabulatedGreekPhase{FT}) where {FT}
    print(io, "TabulatedGreekPhase{$FT}(ℓ ≤ ", size(pf.β, 1) - 1, ", ",
          is_polarizable(pf) ? "phase matrix" : "scalar",
          is_wavelength_dependent(pf) ? ", $(length(pf.λ)) λ" : "",
          isempty(pf.source) ? "" : ", \"$(pf.source)\"", ")")
end

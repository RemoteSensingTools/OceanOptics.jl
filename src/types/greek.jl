# =============================================================================
# Greek (generalized-spherical-function) expansion coefficients
# =============================================================================
#
# One phase matrix, expanded as in de Haan, Bosma & Hovenier (1987) and
# Sanghavi (2014) — the convention of vSmartMOM's `GreekCoefs`, whose field
# names this type shares:
#
#   B_ℓ = [β  γ   0   0
#          γ  α   0   0
#          0  0   ζ   ϵ
#          0  0  -ϵ   δ]
#
# with p(μ) = Σ_ℓ β_ℓ P_ℓ(μ), β₀ = 1, and the Hovenier sign of γ and ϵ
# (γ₂ = +√(3/2) for Rayleigh). A scalar expansion carries β only.
# =============================================================================

"""
$(TYPEDEF)

Expansion coefficients of one phase function or phase matrix, indexed
`[ℓ + 1]` for `ℓ = 0, 1, …, ℓ_max`.

# Fields
$(TYPEDFIELDS)

# Convention

`p(μ) = Σ_ℓ β_ℓ P_ℓ(μ)` with the `(2ℓ+1)` factor included, so `β₀ = 1`.
The phase-matrix layout and signs are those of vSmartMOM's `GreekCoefs`
(Sanghavi 2014 Eq. 16, Hovenier sign of `γ` and `ϵ`: `γ₂ = +√(3/2)` for
Rayleigh scattering). [`evaluate_scattering_matrix`](@ref) converts back to the
standard angular scattering matrix.

A scalar expansion has empty `α, γ, δ, ϵ, ζ`; a polarized one has all six
of equal length.

# Construction

    GreekCoefficients(β; α, γ, δ, ϵ, ζ)

Pass none of `α, γ, δ, ϵ, ζ` (scalar) or all five (polarized). The
constructor checks finiteness, `β₀ = 1` within `√eps`, `|β_ℓ| ≤ 2ℓ+1`
(true for every non-negative phase function), and that `α, γ, ϵ, ζ` vanish
for `ℓ < 2`, where their generalized spherical functions are identically
zero. These are necessary conditions, not a proof of a correct import:
compare `asymmetry_parameter` with published values as well.
"""
struct GreekCoefficients{FT<:Real}
    "α_ℓ: B[2,2] (empty when scalar)"
    α::Vector{FT}
    "β_ℓ: B[1,1], the phase-function expansion"
    β::Vector{FT}
    "γ_ℓ: B[1,2] = B[2,1] (empty when scalar)"
    γ::Vector{FT}
    "δ_ℓ: B[4,4] (empty when scalar)"
    δ::Vector{FT}
    "ϵ_ℓ: B[3,4] = −B[4,3] (empty when scalar)"
    ϵ::Vector{FT}
    "ζ_ℓ: B[3,3] (empty when scalar)"
    ζ::Vector{FT}

    function GreekCoefficients{FT}(α, β, γ, δ, ϵ, ζ) where {FT<:Real}
        g = new{FT}(α, β, γ, δ, ϵ, ζ)
        _validate_greek(g)
        return g
    end
end

function GreekCoefficients(β::AbstractVector; α = nothing, γ = nothing, δ = nothing,
                           ϵ = nothing, ζ = nothing)
    FT     = float(eltype(β))
    given  = (α = α, γ = γ, δ = δ, ϵ = ϵ, ζ = ζ)
    names  = [k for (k, v) in pairs(given) if v !== nothing]
    isempty(names) || length(names) == 5 ||
        throw(ArgumentError("GreekCoefficients: pass all of α, γ, δ, ϵ, ζ for a phase " *
                            "matrix or none for a phase function; got only " *
                            join(names, ", ")))
    vec(x) = x === nothing ? FT[] : collect(FT, x)
    return GreekCoefficients{FT}(vec(α), collect(FT, β), vec(γ), vec(δ), vec(ϵ), vec(ζ))
end

function _validate_greek(g::GreekCoefficients{FT}) where {FT}
    who = "GreekCoefficients"
    (; β) = g
    isempty(β) && throw(ArgumentError("$who: β needs at least the ℓ = 0 term"))
    all(isfinite, β) || throw(ArgumentError("$who: β contains non-finite values"))
    tol = sqrt(eps(float(FT)))
    abs(β[1] - 1) ≤ tol ||
        throw(ArgumentError("$who: β₀ = $(β[1]); phase functions must be normalized so " *
                            "that β₀ = 1 (tolerance $tol)"))
    for ℓ in 1:length(β)-1
        abs(β[ℓ+1]) ≤ (2ℓ + 1) * (1 + tol) ||
            throw(ArgumentError("$who: |β_$ℓ| = $(abs(β[ℓ+1])) > 2ℓ+1 = $(2ℓ+1); no " *
                                "non-negative phase function allows this (was the (2ℓ+1) " *
                                "factor applied twice?)"))
    end
    pol = (α = g.α, γ = g.γ, δ = g.δ, ϵ = g.ϵ, ζ = g.ζ)
    nonempty = [k for (k, v) in pairs(pol) if !isempty(v)]
    isempty(nonempty) && return nothing
    length(nonempty) == 5 ||
        throw(ArgumentError("$who: polarized coefficients need all of α, γ, δ, ϵ, ζ; " *
                            "got only " * join(nonempty, ", ")))
    for (k, v) in pairs(pol)
        length(v) == length(β) ||
            throw(DimensionMismatch("$who: $k has $(length(v)) terms, β has $(length(β))"))
        all(isfinite, v) || throw(ArgumentError("$who: $k contains non-finite values"))
    end
    for k in (:α, :γ, :ϵ, :ζ), ℓ in 0:min(1, length(β) - 1)
        abs(getfield(g, k)[ℓ+1]) ≤ tol ||
            throw(ArgumentError("$who: $(k)_$ℓ = $(getfield(g, k)[ℓ+1]) must be zero (its " *
                                "generalized spherical function vanishes for ℓ < 2); is " *
                                "the table shifted by a row?"))
    end
    return nothing
end

"Number of expansion terms, `ℓ_max + 1`."
n_moments(g::GreekCoefficients) = length(g.β)

"Whether `g` carries all six phase-matrix coefficients."
is_polarized(g::GreekCoefficients) = !isempty(g.α)

Base.:(==)(a::GreekCoefficients, b::GreekCoefficients) =
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(GreekCoefficients))
Base.isapprox(a::GreekCoefficients, b::GreekCoefficients; kwargs...) =
    all(f -> isapprox(getfield(a, f), getfield(b, f); kwargs...), fieldnames(GreekCoefficients))

"""
    resize_moments(g::GreekCoefficients, ℓ_max) -> GreekCoefficients

Expansion with terms `ℓ = 0:ℓ_max`: truncated, or padded with zeros (exact
for the expansion as given). Truncation applies no δ-M or other
correction; that is the solver's choice.
"""
function resize_moments(g::GreekCoefficients{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    n   = ℓ_max + 1
    fit = v -> isempty(v) ? v : (length(v) ≥ n ? v[1:n] : [v; zeros(FT, n - length(v))])
    return GreekCoefficients{FT}(fit(g.α), fit(g.β), fit(g.γ), fit(g.δ), fit(g.ϵ), fit(g.ζ))
end

"""
    depolarized(g::GreekCoefficients) -> GreekCoefficients

The *completely depolarizing* phase matrix with phase function `β`:
`α = γ = δ = ϵ = ζ = 0`, i.e. `diag(p(Θ), 0, 0, 0)`. Scattered light is
unpolarized and its intensity does not depend on the incident
polarization. An explicit modelling approximation for scatterers whose
phase matrix is unknown.
"""
function depolarized(g::GreekCoefficients{FT}) where {FT}
    z = zeros(FT, n_moments(g))
    return GreekCoefficients{FT}(copy(z), copy(g.β), copy(z), copy(z), copy(z), z)
end

"""
    mix_greek(weights, gs) -> GreekCoefficients

Weighted mixture `Σ wᵢ gᵢ / Σ wᵢ` of expansions with equal length and
polarization — the mixing law for scatterers in one volume, with weights
proportional to their scattering coefficients. Weights must be
non-negative with a positive sum.
"""
function mix_greek(weights::AbstractVector{<:Real}, gs::AbstractVector{<:GreekCoefficients})
    length(weights) == length(gs) ||
        throw(DimensionMismatch("mix_greek: $(length(weights)) weights for $(length(gs)) expansions"))
    isempty(gs) && throw(ArgumentError("mix_greek: nothing to mix"))
    all(≥(0), weights) || throw(ArgumentError("mix_greek: weights must be non-negative"))
    W = sum(weights)
    W > 0 || throw(ArgumentError("mix_greek: weights sum to zero"))
    g1 = first(gs)
    all(g -> n_moments(g) == n_moments(g1) && is_polarized(g) == is_polarized(g1), gs) ||
        throw(ArgumentError("mix_greek: expansions differ in length or polarization"))
    FT  = promote_type(eltype(g1.β), float(eltype(weights)))
    acc = f -> isempty(getfield(g1, f)) ? FT[] :
               sum((w / W) .* getfield(g, f) for (w, g) in zip(weights, gs))
    return GreekCoefficients{FT}(acc(:α), acc(:β), acc(:γ), acc(:δ), acc(:ϵ), acc(:ζ))
end

"""
    backscatter_fraction(g::GreekCoefficients) -> FT

Backscatter fraction `B = (1/2)∫₋₁⁰ p(μ) dμ` of the finite expansion,
evaluated exactly term by term. For a truncated expansion of a strongly
peaked phase function this differs from the backscatter fraction of the
phase function itself.
"""
backscatter_fraction(g::GreekCoefficients) = _backscatter_fraction_from_moments(g.β)

"Asymmetry parameter `g = ⟨cos Θ⟩ = β₁ / 3` of the expansion."
asymmetry_parameter(g::GreekCoefficients{FT}) where {FT} =
    n_moments(g) ≥ 2 ? g.β[2] / 3 : zero(FT)

function Base.show(io::IO, g::GreekCoefficients{FT}) where {FT}
    print(io, "GreekCoefficients{$FT}(ℓ ≤ $(n_moments(g) - 1), ",
          is_polarized(g) ? "phase matrix" : "phase function", ")")
end

"""
    _backscatter_fraction_from_moments(β) -> FT

`B = (1/2) [β₀ − Σ_{ℓ odd} βₗ Iₗ]` with `Iₗ = ∫₀¹ Pₗ dμ =
(P_{ℓ-1}(0) − P_{ℓ+1}(0)) / (2ℓ+1)`; even orders ℓ ≥ 2 integrate to zero
over a hemisphere.
"""
function _backscatter_fraction_from_moments(β::AbstractVector{FT}) where {FT}
    s = β[1]
    P_prev = one(FT)                  # P_{ℓ-1}(0) for ℓ = 1
    for ℓ in 1:2:length(β)-1
        P_next = -FT(ℓ) / FT(ℓ + 1) * P_prev   # P_{ℓ+1}(0) = -ℓ/(ℓ+1) P_{ℓ-1}(0)
        s -= β[ℓ + 1] * (P_prev - P_next) / FT(2ℓ + 1)
        P_prev = P_next
    end
    return s / 2
end

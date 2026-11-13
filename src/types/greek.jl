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

# Dense storage for an array field of a result struct. A dense vector of
# element type `FT` (a `Vector`, a device array) is stored as it is; one of
# another element type is converted where it lives; anything else (a range,
# a view, a tuple, a generator) is collected on the host.
_dense(::Type{FT}, x::DenseVector{FT}) where {FT} = x
_dense(::Type{FT}, x::DenseVector) where {FT} = map(FT, x)
_dense(::Type{FT}, x) where {FT} = collect(FT, x)

# As `_dense`, but never the array that was passed in.
_dense_copy(::Type{FT}, x::DenseVector) where {FT} = map(FT, x)
_dense_copy(::Type{FT}, x) where {FT} = collect(FT, x)

# Host copy of an array, for code that reads single elements.
_host(x::Array) = x
_host(x::AbstractArray) = Array(x)

# First argument of the inner constructors that store their arguments without
# validation. Only the `Adapt` rules use them: an adapted struct holds the
# values of a validated one, and a backend's kernel adaptor must be able to
# rebuild the struct without running reductions.
struct _Unchecked end

# The `Adapt` rules move arrays and keep the element type `FT`; `x` is one
# adapted array of the struct.
_check_adapt_eltype(to, ::Type{FT}, x) where {FT} =
    eltype(x) === FT ||
        throw(ArgumentError("adapt: $to changes the element type from $FT to $(eltype(x)); " *
                            "build the optics in the wanted precision (for example " *
                            "`column_optics(Float32, …)`) and adapt with an array type " *
                            "without element type, such as `CuArray`"))

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

    GreekCoefficients{FT}(α, β, γ, δ, ϵ, ζ)

All six vectors by position (empty `α, γ, δ, ϵ, ζ` for a scalar
expansion), stored in the array type of `β`.

# Devices

The six vectors share one array type `V`. `Adapt.adapt(to, g)` moves them
to the device of `to` (for example `adapt(CuArray, g)`), keeps the element
type (an adaptor that changes it, such as `CuArray{Float32}` for a `Float64`
expansion, is refused) and does not validate again, so the expansion can be
captured by a device broadcast or kernel. Everything else that reads single
coefficients ([`phase_id`](@ref), [`backscatter_fraction`](@ref),
[`asymmetry_parameter`](@ref), [`mix_greek`](@ref),
[`evaluate_scattering_matrix`](@ref)) is a host operation.
"""
struct GreekCoefficients{FT<:Real, V<:AbstractVector{FT}}
    "α_ℓ: B[2,2] (empty when scalar)"
    α::V
    "β_ℓ: B[1,1], the phase-function expansion"
    β::V
    "γ_ℓ: B[1,2] = B[2,1] (empty when scalar)"
    γ::V
    "δ_ℓ: B[4,4] (empty when scalar)"
    δ::V
    "ϵ_ℓ: B[3,4] = −B[4,3] (empty when scalar)"
    ϵ::V
    "ζ_ℓ: B[3,3] (empty when scalar)"
    ζ::V

    function GreekCoefficients{FT, V}(α, β, γ, δ, ϵ, ζ) where {FT<:Real, V<:AbstractVector{FT}}
        g = new{FT, V}(α, β, γ, δ, ϵ, ζ)
        _validate_greek(g)
        return g
    end

    GreekCoefficients{FT, V}(::_Unchecked, α, β, γ, δ, ϵ, ζ) where {FT<:Real, V<:AbstractVector{FT}} =
        new{FT, V}(α, β, γ, δ, ϵ, ζ)
end

function GreekCoefficients{FT}(α, β, γ, δ, ϵ, ζ) where {FT<:Real}
    b = _dense(FT, β)
    return GreekCoefficients{FT, typeof(b)}(α, b, γ, δ, ϵ, ζ)
end

function Adapt.adapt_structure(to, g::GreekCoefficients{FT}) where {FT}
    β = adapt(to, g.β)
    _check_adapt_eltype(to, FT, β)
    return GreekCoefficients{FT, typeof(β)}(_Unchecked(), adapt(to, g.α), β, adapt(to, g.γ),
                                            adapt(to, g.δ), adapt(to, g.ϵ), adapt(to, g.ζ))
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
    vec(x) = x === nothing ? FT[] : _dense_copy(FT, x)
    return GreekCoefficients{FT}(vec(α), _dense_copy(FT, β), vec(γ), vec(δ), vec(ϵ), vec(ζ))
end

# The conditions on β_ℓ and on the polarized coefficients below ℓ = 2, per
# element, so that the validation is one broadcast and one reduction on the
# device that holds the coefficients.
_admissible_β(x, ℓ, tol) =
    isfinite(x) & (ℓ == 0 ? abs(x - 1) ≤ tol : abs(x) ≤ (2ℓ + 1) * (1 + tol))
_vanishes_below_2(x, ℓ, tol) = isfinite(x) & ((ℓ ≥ 2) | (abs(x) ≤ tol))

function _validate_greek(g::GreekCoefficients{FT}) where {FT}
    (; β) = g
    isempty(β) && throw(ArgumentError("GreekCoefficients: β needs at least the ℓ = 0 term"))
    Base.require_one_based_indexing(β)
    tol = sqrt(eps(float(FT)))
    ℓ   = 0:length(β)-1
    pol = (g.α, g.γ, g.δ, g.ϵ, g.ζ)
    valid = all(_admissible_β.(β, ℓ, tol)) &&
            (all(isempty, pol) ||
             (all(v -> length(v) == length(β), pol) && all(isfinite, g.δ) &&
              all(v -> all(_vanishes_below_2.(v, ℓ, tol)), (g.α, g.γ, g.ϵ, g.ζ))))
    valid || _throw_greek_error(g, tol)
    return nothing
end

# Says which condition of `_validate_greek` failed and where, on host copies.
function _throw_greek_error(g::GreekCoefficients, tol)
    who = "GreekCoefficients"
    β   = _host(g.β)
    all(isfinite, β) || throw(ArgumentError("$who: β contains non-finite values"))
    abs(β[1] - 1) ≤ tol ||
        throw(ArgumentError("$who: β₀ = $(β[1]); phase functions must be normalized so " *
                            "that β₀ = 1 (tolerance $tol)"))
    for ℓ in 1:length(β)-1
        abs(β[ℓ+1]) ≤ (2ℓ + 1) * (1 + tol) ||
            throw(ArgumentError("$who: |β_$ℓ| = $(abs(β[ℓ+1])) > 2ℓ+1 = $(2ℓ+1); no " *
                                "non-negative phase function allows this (was the (2ℓ+1) " *
                                "factor applied twice?)"))
    end
    pol = (α = _host(g.α), γ = _host(g.γ), δ = _host(g.δ), ϵ = _host(g.ϵ), ζ = _host(g.ζ))
    nonempty = [k for (k, v) in pairs(pol) if !isempty(v)]
    length(nonempty) == 5 ||
        throw(ArgumentError("$who: polarized coefficients need all of α, γ, δ, ϵ, ζ; " *
                            "got only " * join(nonempty, ", ")))
    for (k, v) in pairs(pol)
        length(v) == length(β) ||
            throw(DimensionMismatch("$who: $k has $(length(v)) terms, β has $(length(β))"))
        all(isfinite, v) || throw(ArgumentError("$who: $k contains non-finite values"))
    end
    for k in (:α, :γ, :ϵ, :ζ), ℓ in 0:min(1, length(β) - 1)
        abs(pol[k][ℓ+1]) ≤ tol ||
            throw(ArgumentError("$who: $(k)_$ℓ = $(pol[k][ℓ+1]) must be zero (its " *
                                "generalized spherical function vanishes for ℓ < 2); is " *
                                "the table shifted by a row?"))
    end
    throw(ArgumentError("$who: invalid coefficients"))
end

"""
    n_moments(x)

Return the number of expansion terms, `ℓ_max + 1`, of an expansion, a scattering
component, a layer or a column.
"""
n_moments(g::GreekCoefficients) = length(g.β)

"""
    is_polarized(g::GreekCoefficients)

Return whether `g` carries all six phase-matrix coefficients.
"""
is_polarized(g::GreekCoefficients) = !isempty(g.α)

Base.:(==)(a::GreekCoefficients, b::GreekCoefficients) =
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(GreekCoefficients))
# Consistent with `==`, which compares values across element types; `phase_id`
# separates element types and cannot serve here.
Base.hash(g::GreekCoefficients, h::UInt) =
    hash(map(_host, (g.α, g.β, g.γ, g.δ, g.ϵ, g.ζ)), hash(:GreekCoefficients, h))
Base.isapprox(a::GreekCoefficients, b::GreekCoefficients; kwargs...) =
    all(f -> isapprox(getfield(a, f), getfield(b, f); kwargs...), fieldnames(GreekCoefficients))

"""
    resize_moments(g::GreekCoefficients, ℓ_max) -> GreekCoefficients

Return the expansion `g` with terms `ℓ = 0:ℓ_max`, truncated or padded with zeros.

Padding is exact for the expansion as given. Truncation applies no δ-M or other
correction; that is the solver's choice.
"""
function resize_moments(g::GreekCoefficients{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    n   = ℓ_max + 1
    fit = v -> isempty(v) ? v : length(v) ≥ n ? v[1:n] :
                                vcat(v, fill!(similar(v, n - length(v)), zero(FT)))
    return GreekCoefficients{FT}(fit(g.α), fit(g.β), fit(g.γ), fit(g.δ), fit(g.ϵ), fit(g.ζ))
end

"""
    depolarized(g::GreekCoefficients) -> GreekCoefficients

Return the *completely depolarizing* phase matrix with the phase function of `g`.

The result has `α = γ = δ = ϵ = ζ = 0`, i.e. the matrix `diag(p(Θ), 0, 0, 0)`: scattered
light is unpolarized and its intensity does not depend on the incident polarization. It
is an explicit modelling approximation for scatterers whose phase matrix is unknown.
"""
function depolarized(g::GreekCoefficients{FT}) where {FT}
    z = zero(g.β)
    return GreekCoefficients{FT}(copy(z), copy(g.β), copy(z), copy(z), copy(z), z)
end

"""
    zero_padded(g::GreekCoefficients)

Return `g` with all six coefficient vectors of length `ℓ_max + 1`.

A polarized `g` is returned as it is; a scalar expansion gets
`α = γ = δ = ϵ = ζ = 0` (equal to [`depolarized`](@ref)`(g)`). For
consumers that size their arrays by `length(g.α)`; the zeros carry no
polarization information.
"""
zero_padded(g::GreekCoefficients) = is_polarized(g) ? g : depolarized(g)

# FNV-1a (64 bit) over the bytes of the content: defined by the algorithm
# alone, so keys do not depend on the session, the Julia version or the
# identity of the arrays.
const _FNV_OFFSET = 0xcbf29ce484222325
const _FNV_PRIME  = 0x00000100000001b3

_fnv1a(h::UInt64, x::UInt8) = (h ⊻ x) * _FNV_PRIME
function _fnv1a(h::UInt64, x::Unsigned)          # little-endian byte order
    for shift in 0:8:(8 * sizeof(x) - 8)
        h = _fnv1a(h, (x >> shift) % UInt8)
    end
    return h
end
_fnv1a(h::UInt64, s::AbstractString) = foldl(_fnv1a, codeunits(s); init = h)
# `x + zero(x)` maps −0.0 to 0.0, so that values equal under `==` share a key.
_fnv1a(h::UInt64, x::Union{Float16, Float32, Float64}) = _fnv1a(h, reinterpret(Unsigned, x + zero(x)))
_fnv1a(h::UInt64, x::Real) = _fnv1a(h, string(x))

"""
    phase_id(g::GreekCoefficients) -> key

Return the content key of the expansion `g`, a `UInt64`.

The key is the 64-bit FNV-1a hash of the element type name and of the length and
values of `α, β, γ, δ, ϵ, ζ`, in that order.

Expansions that are `==` and have the same element type have the same
key, in every session and Julia version; expansions that differ in a
coefficient, in `ℓ_max`, in element type or in being scalar or polarized
have different keys (up to the 2⁻⁶⁴ collision probability of the hash).
A solver can therefore cache phase operators by `(phase_id, Fourier mode)`.
The key of coefficients held on a device is computed from a host copy and
equals the key of that copy.
"""
function phase_id(g::GreekCoefficients{FT}) where {FT}
    h = _fnv1a(_FNV_OFFSET, string(FT))
    for v in map(_host, (g.α, g.β, g.γ, g.δ, g.ϵ, g.ζ))
        h = _fnv1a(h, UInt64(length(v)))
        for x in v
            h = _fnv1a(h, x)
        end
    end
    return h
end

"""
    mix_greek(weights, gs) -> GreekCoefficients

Return the weighted mixture `Σ wᵢ gᵢ / Σ wᵢ` of expansions of equal length and
polarization.

This is the mixing law for scatterers in one volume, with weights proportional to their
scattering coefficients. Weights must be non-negative with a positive sum.
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

Return the backscatter fraction `B = (1/2)∫₋₁⁰ p(μ) dμ` of the finite expansion `g`.

The integral is evaluated exactly term by term. For a truncated expansion of a strongly
peaked phase function it differs from the backscatter fraction of the phase function
itself.
"""
backscatter_fraction(g::GreekCoefficients) = _backscatter_fraction_from_moments(g.β)

"""
    asymmetry_parameter(g::GreekCoefficients)

Return the asymmetry parameter `⟨cos Θ⟩ = β₁ / 3` of the expansion `g`.
"""
asymmetry_parameter(g::GreekCoefficients{FT}) where {FT} =
    n_moments(g) ≥ 2 ? g.β[2] / 3 : zero(FT)

function Base.show(io::IO, g::GreekCoefficients{FT}) where {FT}
    print(io, "GreekCoefficients{$FT}(ℓ ≤ $(n_moments(g) - 1), ",
          is_polarized(g) ? "phase matrix" : "phase function", ")")
end

"""
    _backscatter_fraction_from_moments(β) -> FT

Return `B = (1/2) [β₀ − Σ_{ℓ odd} βₗ Iₗ]` with `Iₗ = ∫₀¹ Pₗ dμ =
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

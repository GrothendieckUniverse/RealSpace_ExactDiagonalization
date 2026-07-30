const PHASE_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULT_ROOT = joinpath(PHASE_ROOT, "results")
const FIGURE_ROOT = joinpath(PHASE_ROOT, "figures")
const CHECKPOINT_ROOT = joinpath(PHASE_ROOT, "checkpoints")

const TPP_DENOMINATOR = 2 + 2 * sqrt(2)

"Physical model parameters. `t′′` is replaced by `x / TPP_DENOMINATOR` per job."
const BASE_PARAMS = Dict{String,Float64}(
    "t" => -1.0,
    "t′_1" => -1.4 / (2 + sqrt(2)),
    "t′_2" => 1.4 / (2 + sqrt(2)),
    "t′′" => -1.0 / TPP_DENOMINATOR,
    "ϕ_over_2π" => 1 / 8,
    "V1" => 2.0,
    "V2" => 0.45,
    "V3" => 0.2,
    "λ" => 2.8,
)

# The regular extensive grid. Exact diagnostic points are merged into the
# production sweep below so every diagnostic has a matching zero-flux ED scan.
const REGULAR_SWEEP_NUMERATORS = round.(collect(-3.0:0.1:1.5); digits=10)
const STUDY_GEOMETRIES = [(3, 4), (3, 5), (3, 6)]
const CHARGE_GAP_GEOMETRIES = [(3, 3), (3, 4), (3, 5), (3, 6), (3, 7)]
const FCI_PROJECTOR_PLOT_WINDOW = (-0.35, -0.20)

struct PhaseSpec
    name::Symbol
    manifold_size::Int
end

# Phase labels carry only metadata.  There is deliberately no single default
# parameter for a phase: the characteristic-point campaign compares several
# physical t'' values on each side of the suspected transitions.
const PHASE_SPECS = Dict{Symbol,PhaseSpec}(
    :AHC => PhaseSpec(:AHC, 3),
    :FCI => PhaseSpec(:FCI, 3),
    :CDW => PhaseSpec(:CDW, 3),
)

# Physical t'' values, not numerator values.  `CDW` remains a historical
# working label until the positive-side order has been established.
const CHARACTERISTIC_TPP_VALUES = Dict{Symbol,Vector{Float64}}(
    :AHC => [-0.60, -0.55, -0.50, -0.45],
    :FCI => [-0.30, -0.15, 0.00, 0.05, 0.10],
    :CDW => [0.20, 0.30],
)

const DIAGNOSTIC_TPP_VALUES =
    sort!(unique(vcat(values(CHARACTERISTIC_TPP_VALUES)...)))
const DIAGNOSTIC_SWEEP_NUMERATORS = TPP_DENOMINATOR .* DIAGNOSTIC_TPP_VALUES
const SWEEP_NUMERATORS = sort!(unique(vcat(
    REGULAR_SWEEP_NUMERATORS,
    DIAGNOSTIC_SWEEP_NUMERATORS,
)))

# The symmetry slots that form the already-established FCI manifold at its
# reference point.  Each entry is ((k1,k2), level_in_sector), not merely a
# sector label: on 3x6 all three states belong to (0,3).  At another parameter,
# a slot does not by itself guarantee wavefunction continuity through a
# same-sector avoided crossing; that requires overlaps/fidelity tracking.
const FCI_REFERENCE_MANIFOLD = Dict{Tuple{Int,Int},Vector{Tuple{Tuple{Int,Int},Int}}}(
    (3, 4) => [((2, 2), 1), ((1, 2), 1), ((0, 2), 1)],
    (3, 5) => [((0, 0), 1), ((1, 0), 1), ((2, 0), 1)],
    (3, 6) => [((0, 3), 1), ((0, 3), 2), ((0, 3), 3)],
)

actual_tpp(x::Real) = Float64(x) / TPP_DENOMINATOR
numerator_at_tpp(tpp::Real) = Float64(tpp) * TPP_DENOMINATOR
geometry_tag(sample::Tuple{Int,Int}) = "$(sample[1])x$(sample[2])"
geometry_tag(sample::AbstractVector{<:Integer}) = geometry_tag((Int(sample[1]), Int(sample[2])))
tpp_tag(x::Real) = replace(@sprintf("%.4f", Float64(x)), "-" => "m", "." => "p")
phase_point_tag(x::Real) = "tpp_$(tpp_tag(actual_tpp(x)))"

function phase_spec(name)
    key = Symbol(uppercase(String(name)))
    haskey(PHASE_SPECS, key) || error("Unknown phase $name; choose AHC, FCI, or CDW.")
    return PHASE_SPECS[key]
end

function characteristic_tpp_values(name)
    spec = phase_spec(name)
    return copy(CHARACTERISTIC_TPP_VALUES[spec.name])
end

"Production diagnostics: PES is meaningful here only for FCI candidates."
function default_diagnostic_observables(name)
    observables = [:structure, :flow, :pump]
    phase_spec(name).name == :FCI && push!(observables, :pes)
    return observables
end

function require_phase_numerator(name, x::Union{Nothing,Real})
    x !== nothing && return Float64(x)
    values = characteristic_tpp_values(name)
    formatted = join([@sprintf("%.2f", value) for value in values], ", ")
    error("Phase $(phase_spec(name).name) has multiple characteristic t′′ values " *
          "($formatted). Supply --tpp (physical value) or --x (numerator).")
end

"Recommended solver mode by geometry and task."
function mode_for(sample::Tuple{Int,Int}, task::Symbol=:sweep)
    if sample == (4, 6)
        return :matrixfree
    else
        return :matrix
    end
end

function params_at_numerator(x::Real)
    params = copy(BASE_PARAMS)
    params["t′′"] = actual_tpp(x)
    return params
end

function default_particle_number(sample::Tuple{Int,Int})
    ncell = prod(sample)
    ncell % 3 == 0 || error("ν=1/3 band filling requires L1*L2 divisible by 3; got $sample.")
    return ncell ÷ 3
end

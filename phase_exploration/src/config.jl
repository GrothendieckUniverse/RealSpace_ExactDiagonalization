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

# A uniform extensive grid. Extra points may be supplied to each CLI runner.
const SWEEP_NUMERATORS = round.(collect(-3.0:0.1:1.5); digits=10)
const STUDY_GEOMETRIES = [(3, 5), (3, 6), (3, 7), (4, 6)]
const CHARGE_GAP_GEOMETRIES = [(3, 4), (3, 5), (3, 6), (3, 7), (4, 6)]

struct PhaseSpec
    name::Symbol
    numerator::Float64
    manifold_size::Int
end

# Deep, large-gap points inferred from the archived 3x4/3x5 study.
# These are deliberately centralized: edit here before generating the jobs if
# a new preliminary sweep suggests better common points for all geometries.
const PHASE_SPECS = Dict{Symbol,PhaseSpec}(
    :AHC => PhaseSpec(:AHC, -3.0, 3),
    :FCI => PhaseSpec(:FCI, -1.0, 3),
    :CDW => PhaseSpec(:CDW, 1.5, 3),
)

actual_tpp(x::Real) = Float64(x) / TPP_DENOMINATOR
geometry_tag(sample::Tuple{Int,Int}) = "$(sample[1])x$(sample[2])"
geometry_tag(sample::AbstractVector{<:Integer}) = geometry_tag((Int(sample[1]), Int(sample[2])))
tpp_tag(x::Real) = replace(@sprintf("%.4f", Float64(x)), "-" => "m", "." => "p")

function phase_spec(name)
    key = Symbol(uppercase(String(name)))
    haskey(PHASE_SPECS, key) || error("Unknown phase $name; choose AHC, FCI, or CDW.")
    return PHASE_SPECS[key]
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

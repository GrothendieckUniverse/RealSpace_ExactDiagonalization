module CheckerboardPhaseStudy

using CairoMakie
using Distributed
using JLD2
using LinearAlgebra
using Printf
using RealSpace_ExactDiagonalization
using Statistics
using TightBinding

include("config.jl")
include("cli.jl")
include("io.jl")
include("model.jl")
include("sweep.jl")
include("diagnostics.jl")
include("charge_gap.jl")
include("plots.jl")

export BASE_PARAMS,
    TPP_DENOMINATOR,
    SWEEP_NUMERATORS,
    STUDY_GEOMETRIES,
    PHASE_SPECS,
    CHARACTERISTIC_TPP_VALUES,
    FCI_REFERENCE_MANIFOLD,
    PhaseSpec,
    actual_tpp,
    numerator_at_tpp,
    phase_point_tag,
    characteristic_tpp_values,
    geometry_tag,
    parse_cli,
    parse_geometry,
    parse_bool,
    parse_float_list,
    mode_for,
    phase_spec,
    run_sweep_point,
    run_phase_diagnostics,
    run_charge_gap_point,
    plot_sweep_results,
    plot_ed_spectrum_results,
    plot_structure_factor_results,
    plot_diagnostic_results,
    plot_charge_gap_results,
    plot_all_results

end

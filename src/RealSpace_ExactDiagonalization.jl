module RealSpace_ExactDiagonalization

using Distributed
using LinearAlgebra
using SparseArrays, KrylovKit, Arpack
using MLStyle
using CairoMakie
using Test
using Printf

using TightBinding

include("bitwise_operations.jl")
export BitWise_Operations

include("second_quantized_model.jl")

include("symmetry_resolved_ed.jl")

include("../observables/spectrum_flow.jl")
include("../observables/charge_pump.jl")
include("../observables/density_distribution.jl")
include("../observables/static_structure_factor.jl")
include("../observables/off_diagonal_long_range_order.jl")
include("../observables/entanglement_spectrum.jl")
include("../observables/many_body_chern_number.jl")

include("../test/bosonic_fci.jl")
include("../test/fermionic_fci.jl")
include("../test/fermionic_fci_phase_explore.jl")

export encode_configuration_to_bit_mask, decode_bit_mask_to_configuration, decode_bit_mask_to_configuration!,
    Bosonic, Fermionic, Particle_Statistics,
    Real_Space_Second_Quantized_Model,
    Symmetry_Operation, Finite_Symmetry_Group, OneDim_Irrep,
    Symmetry_Orbit_Catalog, Symmetry_Sector_Basis, Symmetry_Resolved_ED_Data,
    CanonicalMap, MatrixFreeHamiltonian,
    build_identity_group, build_translation_group,
    build_identity_irrep_list, build_translation_irrep_list, build_irrep_list,
    build_symmetry_orbit_catalog, update_orbit_stabilizer_phases!,
    build_symmetry_sector_basis, build_ed_data,
    build_ed_Hamiltonian_symmetry_block, build_ed_Hamiltonian_symmetry_block_distributed,
    hamiltonian_linear_operator, hamiltonian_linear_operator_distributed,
    apply_hamiltonian!, apply_hamiltonian_distributed!,
    populate_canonical_map!,
    ed_scan!, ed_scan_at_irrep_matrix!, ed_scan_at_irrep_matrixfree!,
    diagonalize_block_dense, diagonalize_block_arpack, diagonalize_block_matrixfree,
    full_ed, print_spectrum, plot_spectrum, plot_ed_scan_res,
    save_checkpoint, load_checkpoint,
    ed_scan_checkpoint_filename,
    update_second_quantized_model_with_twisted_phases!,
    flux_spectrum_flow,
    many_body_position_phases, flux_charge_pump,
    vertices_occupation_distribution_full_ed,
    static_structure_factor, structure_factor_allowed_momenta,
    static_structure_factor_manifold_average, structure_factor_manifold_allowed_momenta,
    compute_structure_factor_map, compute_structure_factor_manifold_average_map,
    plot_structure_factor_map, plot_structure_factor_map_panels,
    plot_structure_factor_allowed_momenta, plot_structure_factor_allowed_momenta_panels,
    off_diagonal_long_range_order, compute_odlro_map, plot_odlro_map, plot_odlro_map_panels,
    entanglement_spectrum, plot_entanglement_spectrum,
    particle_entanglement_spectrum, plot_particle_entanglement_spectrum,
    many_body_chern_number,
    build_zero_flux_bosonic_fci_second_quantized_model, default_fci_sectors,
    test_bosonic_fci_spectrum_flow, test_bosonic_fci_charge_pump,
    test_bosonic_fci_odlro_demo, test_bosonic_fci_static_structure_factor_demo,
    params_Sun_Gu_Katsura_Sarma,
    build_zero_flux_fermionic_fci_second_quantized_model,
    default_fci_sectors_fermionic,
    test_fermionic_fci_full_ed,
    test_fermionic_fci_spectrum_flow, test_fermionic_fci_charge_pump,
    my_optimal_param,
    fermionic_phase_explore_full_ed, fermionic_phase_explore_spectrum_flow,
    fermionic_phase_explore_charge_pump, fermionic_phase_explore_structure_factor,
    fermionic_phase_explore_allowed_structure_factor,
    fermionic_phase_explore_demo



end # module RealSpace_ExactDiagonalization

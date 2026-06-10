module RealSpace_ExactDiagonalization

using Distributed
using LinearAlgebra
using SparseArrays, KrylovKit, Arpack
using MLStyle
using CairoMakie
using Test

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

include("../test/bosonic_fci.jl")
include("../test/fermionic_fci.jl")

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
    static_structure_factor, compute_structure_factor_map, plot_structure_factor_map,
    off_diagonal_long_range_order, compute_odlro_map, plot_odlro_map, plot_odlro_map_panels,
    build_zero_flux_bosonic_fci_second_quantized_model, default_fci_sectors,
    test_bosonic_fci_spectrum_flow, test_bosonic_fci_charge_pump,
    test_bosonic_fci_odlro_demo,
    params_Sun_Gu_Katsura_Sarma,
    build_zero_flux_fermionic_fci_second_quantized_model,
    default_fci_sectors_fermionic,
    test_fermionic_fci_full_ed,
    test_fermionic_fci_spectrum_flow, test_fermionic_fci_charge_pump



end # module RealSpace_ExactDiagonalization

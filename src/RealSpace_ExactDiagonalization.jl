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

export Bosonic, Fermionic, Statistics,
    Real_Space_Second_Quantized_Model,
    Symmetry_Operation, Finite_Symmetry_Group, OneDim_Irrep,
    Symmetry_Orbit_Catalog, Symmetry_Sector_Basis, Symmetry_Resolved_ED_Data,
    CanonicalMap, MatrixFreeHamiltonian,
    build_identity_group, build_translation_group,
    build_identity_irrep_list, build_translation_irrep_list, build_irrep_list,
    build_symmetry_orbit_catalog, build_symmetry_sector_basis, build_ed_data,
    build_ed_Hamiltonian_symmetry_block, build_ed_Hamiltonian_symmetry_block_distributed,
    hamiltonian_linear_operator, hamiltonian_linear_operator_distributed,
    apply_hamiltonian!, apply_hamiltonian_distributed!,
    populate_canonical_map!,
    ed_scan!, ed_scan_at_irrep_matrix!, ed_scan_at_irrep_matrixfree!,
    diagonalize_block_dense, diagonalize_block_arpack, diagonalize_block_matrixfree,
    full_ed, print_spectrum, plot_spectrum, plot_ed_scan_res,
    save_checkpoint, load_checkpoint



end # module RealSpace_ExactDiagonalization

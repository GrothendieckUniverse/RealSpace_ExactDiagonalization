using Test: @test, @testset

using RealSpace_ExactDiagonalization

include("../src/test.jl")
const ModelChecks = Main.Test

@testset "Haldane hard-core boson 2x3 FCI" begin
    vals, labels = ModelChecks.sanity_check_bose_hubbard_2x3()
    @test abs(vals[1] + 7.16380536) < 1e-6
    @test abs(vals[2] + 7.16337536) < 1e-6
    @test labels[1][2] == (0, 0)
    @test labels[2][2] == (1, 0)
end

@testset "Matrix and matrix-free agreement" begin
    maxdiff, _, _ = ModelChecks.validate_matrix_matrixfree_agreement(; sample_size=[2, 3], n_filled=3, nev=3)
    @test maxdiff < 1e-10
end

@testset "Spinful Fermi Hubbard 3x4 reference" begin
    E0, _ = ModelChecks.sanity_check_spinful_fermi_hubbard_3x4()
    @test abs(E0 + 4.913259209075605) < 1e-7
end

@testset "Generic spinless fermion statistics" begin
    E0, expected = ModelChecks.sanity_check_spinless_fermion_chain()
    @test abs(E0 - expected) < 1e-10
end

using Test, LinearAlgebra, SparseArrays, RealSpace_ExactDiagonalization
const ED = RealSpace_ExactDiagonalization
include(joinpath(@__DIR__, "..", "phase_exploration", "src", "CheckerboardPhaseStudy.jl"))
const PE = CheckerboardPhaseStudy
@testset "Flux projection into full Fock space" begin
    for stats in (Bosonic(), Fermionic())
        model,_,_,filling = PE.build_checkerboard_problem((3,2),0.0)
        model.particle_statistics=stats
        for flux in ([0.,0.],[.23,.17],[1.,0.])
            update_second_quantized_model_with_twisted_phases!(model;twisted_phases_over_2π=flux)
            d=build_ed_data(model;filling_fraction=filling,symmetry_group=build_translation_group(model.lattice,flux))
            full=build_ed_data(model;filling_fraction=filling,symmetry_group=build_identity_group(model.lattice.n_site))
            function block(data,idx)
                b=build_symmetry_sector_basis(data.orbit_catalog,data.irrep_list[idx])
                c=ED.CanonicalMap(data.symmetry_group,stats,data.orbit_catalog)
                H=ED.build_ed_Hamiltonian_symmetry_block(b,model.bilinear_terms,model.density_density_terms,c)
                return b,c,Matrix(H)
            end
            bf,_,Hf=block(full,1); allvals=Float64[]
            for idx in eachindex(d.irrep_list)
                b,c,H=block(d,idx); n=size(H,1); B=zeros(ComplexF64,size(Hf,1),n)
                for j in 1:n
                    v=zeros(ComplexF64,n); v[j]=1
                    for (mask,amp) in ED._expand_sector_state_to_fock_amplitudes(v,b,stats)
                        B[bf.representative_mask_to_mask_idx_map[mask],j]=amp
                    end
                end
                @test H ≈ H' atol=1e-11
                @test B'B ≈ I atol=1e-11
                @test H ≈ B'Hf*B atol=1e-11
                @test Hf*B ≈ B*H atol=1e-11
                mf=ED.MatrixFreeHamiltonian(b,model.bilinear_terms,model.density_density_terms,c)
                x=ComplexF64.(1:n) .+ im
                @test mf(x) ≈ H*x atol=1e-10
                append!(allvals,eigvals(Hermitian(H)))
                for mask in bf.representative_mask_list
                    a=ED.project_to_sector(mask,b,c)
                    z=ED.project_to_unnormalized_sector(mask,b,stats)
                    @test (a===nothing)==(z===nothing)
                    if a!==nothing
                        @test a[1]==z[1]
                        @test a[2]≈z[2] atol=1e-11
                    end
                end
            end
            @test sort(allvals) ≈ eigvals(Hermitian(Hf)) atol=1e-10
        end
    end
end
@testset "Flux observable checkpoint and full-space regressions" begin
    mktempdir() do dir
        model,_,_,filling=PE.build_checkerboard_problem((3,2),0.0)
        fluxes=[0.,.23,1.]
        a=flux_spectrum_flow(model,[(0,0)];filling_fraction=filling,nev=4,twisted_phases_over_2π_list=fluxes,checkpoint_dir=joinpath(dir,"sectors"))
        b=flux_spectrum_flow(model,:identity;filling_fraction=filling,nev=4,twisted_phases_over_2π_list=fluxes,checkpoint_dir=joinpath(dir,"sectors"))
        @test all(isfinite,b.energies)
        @test a.global_energies ≈ b.global_energies atol=1e-10
        @test length(load_checkpoint(a.checkpoint_paths[2]).ed_scan_res)==6
        @test a.global_energies[1,:]≈a.global_energies[end,:] atol=1e-10
        resumed=flux_spectrum_flow(model,[(0,0)];filling_fraction=filling,nev=4,twisted_phases_over_2π_list=fluxes,checkpoint_dir=joinpath(dir,"sectors"))
        @test resumed.global_energies≈a.global_energies atol=1e-12
    end
    @test_throws ErrorException ED._unit_det_link([Dict(ED.Mask(1)=>1. +0im)],[Dict(ED.Mask(2)=>1. +0im)])
end

# test/test_assembly_simple.jl - Simplified unit tests for Assembly modules
# Focus on basic functionality to validate the module extraction

using Test
using BridgeFEM
using SparseArrays
using LinearAlgebra

@testset "Assembly Module Basic Tests" begin

    @testset "Matrix Assembly Basic Tests" begin
        
        @testset "assemble_matrices function - basic" begin
            # Create simple bridge configuration for testing  
            E_T = [0.0 200e9; 100.0 160e9]  # Temperature-dependent Young's modulus
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(3, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Test at reference temperature
            T = 20.0
            M, K = assemble_matrices(bridge_opts, T)
            
            # Verify types and dimensions
            @test M isa SparseMatrixCSC{Float64, Int64}
            @test K isa SparseMatrixCSC{Float64, Int64}
            @test size(M) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            @test size(K) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            
            # Verify matrix properties
            @test issymmetric(M)
            @test issymmetric(K)
            @test all(diag(M) .> 0)  # Mass matrix positive definite
            @test all(diag(K) .>= 0) # Stiffness matrix positive semi-definite
            
            println("✅ Basic assembly_matrices test passed")
        end
        
        @testset "assemble_stiffness! function - basic" begin
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(3, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Pre-allocate stiffness matrix
            K = spzeros(bridge_opts.n_dofs, bridge_opts.n_dofs)
            
            # Test assembly with known properties
            EA = 20e6  # N
            EI = 200.0  # N⋅m²
            K_result = assemble_stiffness!(K, bridge_opts, EA, EI)
            
            # Verify return value and properties
            @test K_result === K
            @test issymmetric(K_result)
            @test nnz(K_result) > 0
            @test all(diag(K_result) .>= 0)
            
            println("✅ Basic assemble_stiffness! test passed")
        end
        
        @testset "Temperature dependence" begin
            E_T = [0.0 210e9; 50.0 200e9; 100.0 180e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 7850.0, 0.02, 0.0008, E_T, 100.0)
            
            # Test at different temperatures
            M1, K1 = assemble_matrices(bridge_opts, 0.0)   # Low temperature
            M2, K2 = assemble_matrices(bridge_opts, 100.0) # High temperature
            
            # Mass should be identical (temperature independent)
            @test norm(M1 - M2) < 1e-10
            
            # Stiffness should differ due to temperature-dependent E
            @test norm(K2) < norm(K1)  # Higher temp = lower E = lower stiffness
            
            println("✅ Temperature dependence test passed")
        end
    end
    
    @testset "DOF Mapping Basic Tests" begin
        
        @testset "DOF mapping with simple support" begin
            E_T = [0.0 200e9; 50.0 200e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Create test support using proper constructor
            E_T_support = [0.0 200e9; 50.0 200e9]
            support = SupportElement(
                2,      # connection_node
                [1, 2], # connection_dofs (x, y)
                0.0,    # angle
                2,      # n_elem
                0.1,    # A
                0.001,  # I
                E_T_support, # E_T matrix
                5.0,    # L
                [1, 2] # bc_bottom
            )
            
            support_dof_maps, total_dofs = create_support_dof_mapping(bridge_opts, [support])
            
            # Basic validation
            @test length(support_dof_maps) == 1
            @test total_dofs > bridge_opts.n_dofs
            @test all(support_dof_maps[1] .> 0)  # All DOFs assigned
            
            println("✅ Basic DOF mapping test passed")
        end
    end
    
    @testset "Helper Function Basic Tests" begin
        
        @testset "create_expanded_transformation" begin
            angle = π/4  # 45 degrees
            n_nodes = 2
            
            T_exp = create_expanded_transformation(angle, n_nodes)
            
            # Verify dimensions and properties
            expected_size = 3 * n_nodes
            @test size(T_exp) == (expected_size, expected_size)
            @test norm(T_exp' * T_exp - I) < 1e-10  # Orthogonal
            
            println("✅ Transformation matrix test passed")
        end
    end
    
    @testset "Integration Tests" begin
        
        @testset "Assembly functions work with Core types" begin
            # Test realistic bridge configuration
            E_T = [0.0 210e9; 50.0 200e9; 100.0 180e9]  
            bc = BridgeBC([[1, "all"], [5, "all"]])  # Fixed ends
            bridge_opts = BridgeOptions(8, bc, 20.0, 7850.0, 0.02, 0.0008, E_T, 100.0)
            
            # Should work without errors
            M, K = assemble_matrices(bridge_opts, 25.0)
            @test size(M) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            @test size(K) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            @test issymmetric(M)
            @test issymmetric(K)
            
            println("✅ Integration test passed")
        end
        
        @testset "Numerical precision" begin
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Test repeatability
            M1, K1 = assemble_matrices(bridge_opts, 25.0)
            M2, K2 = assemble_matrices(bridge_opts, 25.0)
            
            # Results should be identical to machine precision
            @test norm(M1 - M2) < 1e-15
            @test norm(K1 - K2) < 1e-15
            
            println("✅ Numerical precision test passed")
        end
    end

end

println("🎉 All Assembly module tests completed!")
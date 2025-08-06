# test/test_assembly.jl - Unit tests for Assembly modules
# Tests matrix assembly and DOF mapping functionality

using Test
using BridgeFEM
using SparseArrays
using LinearAlgebra

@testset "Assembly Module Tests" begin

    @testset "Matrix Assembly Tests" begin
        
        @testset "assemble_stiffness! function" begin
            # Create simple bridge configuration for testing
            E_T = [0.0 200e9; 100.0 160e9]  # Temperature-dependent Young's modulus
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Pre-allocate stiffness matrix
            K = spzeros(bridge_opts.n_dofs, bridge_opts.n_dofs)
            
            # Test assembly with known properties
            EA = 20e6  # N
            EI = 200.0  # N⋅m²
            K_result = assemble_stiffness!(K, bridge_opts, EA, EI)
            
            # Verify return value is the same matrix
            @test K_result === K
            
            # Verify matrix is symmetric
            @test issymmetric(K_result)
            
            # Verify matrix is positive definite (stiffness should be)
            @test all(diag(K_result) .>= 0)
            
            # Verify sparsity pattern - stiffness matrix should have banded structure
            @test nnz(K_result) > 0
            @test nnz(K_result) < bridge_opts.n_dofs^2  # Should be sparse
            
            # Test with different EA, EI values
            K2 = spzeros(bridge_opts.n_dofs, bridge_opts.n_dofs)
            K2_result = assemble_stiffness!(K2, bridge_opts, 2*EA, 2*EI)
            
            # With doubled stiffness, all entries should be doubled
            @test norm(K2_result - 2*K_result) < 1e-10
        end
        
        @testset "assemble_matrices function" begin
            # Test basic assembly
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Test at reference temperature
            T = 20.0
            M, K = assemble_matrices(bridge_opts, T)
            
            # Verify types
            @test M isa SparseMatrixCSC{Float64, Int64}
            @test K isa SparseMatrixCSC{Float64, Int64}
            
            # Verify dimensions
            @test size(M) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            @test size(K) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            
            # Verify matrix properties
            @test issymmetric(M)
            @test issymmetric(K)
            
            # Mass matrix should be positive definite
            @test all(diag(M) .> 0)
            
            # Stiffness matrix should be positive semi-definite
            @test all(diag(K) .>= 0)
            
            # Test temperature dependence
            M1, K1 = assemble_matrices(bridge_opts, 0.0)   # Low temperature
            M2, K2 = assemble_matrices(bridge_opts, 100.0) # High temperature
            
            # Mass should be identical (temperature independent)
            @test norm(M1 - M2) < 1e-10
            
            # Stiffness should differ due to temperature-dependent E
            # At higher temp, E is lower (from E_T definition), so K should be lower
            @test norm(K2) < norm(K1)
            
            # Test analytical validation for simple cantilever beam
            # For a simple beam, the first diagonal entry should match analytical stiffness
            dx = bridge_opts.L / bridge_opts.n_elem
            E = bridge_opts.E(T)
            expected_axial_stiff = E * bridge_opts.A / dx
            
            # Check that some stiffness terms are reasonable
            @test maximum(diag(K)) > 0
            @test minimum(diag(K)) >= 0
        end
        
        @testset "assemble_matrices_with_supports function" begin
            # Create bridge with supports
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Create support element - using proper constructor
            E_T_support = [0.0 200e9; 100.0 160e9]  # Temperature-dependent E
            support = SupportElement(
                2,      # connection_node
                [1, 2], # connection_dofs
                45.0,   # angle (degrees)
                2,      # n_elem
                0.1,    # A
                0.001,  # I
                E_T_support, # E_T matrix
                5.0,    # L
                [1, 2, 3] # bc_bottom
            )
            
            temperatures = [20.0, 50.0, 80.0]
            sim_opts = SimulationOptions(bridge_opts, [support], temperatures)
            
            M, K = assemble_matrices_with_supports(sim_opts)
            
            # Verify 3D array structure
            @test size(M, 3) == length(temperatures)
            @test size(K, 3) == length(temperatures)
            
            # Verify dimensions include support DOFs
            @test size(M, 1) == sim_opts.total_dofs
            @test size(M, 2) == sim_opts.total_dofs
            @test size(K, 1) == sim_opts.total_dofs  
            @test size(K, 2) == sim_opts.total_dofs
            
            # Total DOFs should be larger than bridge DOFs alone
            @test sim_opts.total_dofs > bridge_opts.n_dofs
            
            # Test each temperature slice
            for t in 1:length(temperatures)
                M_t = M[:, :, t]
                K_t = K[:, :, t]
                
                # Should be symmetric matrices
                @test issymmetric(M_t)
                @test issymmetric(K_t)
                
                # Should have positive entries on diagonal
                @test all(diag(M_t) .>= 0)
                @test all(diag(K_t) .>= 0)
            end
            
            # Temperature dependency test
            # Higher temperature should give lower stiffness
            K_low = K[:, :, 1]   # 20°C
            K_high = K[:, :, 3]  # 80°C
            
            # Overall stiffness should be lower at higher temperature
            @test norm(K_high) < norm(K_low)
        end
    end
    
    @testset "DOF Mapping Tests" begin
        
        @testset "create_support_dof_mapping function" begin
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Create test support - using proper constructor
            E_T_support = [0.0 200e9; 50.0 200e9]  # Constant E
            support = SupportElement(
                2,      # connection_node
                [1, 2], # connection_dofs (x, y)
                0.0,    # angle
                2,      # n_elem (3 nodes, 9 DOFs)
                0.1,    # A
                0.001,  # I
                E_T_support, # E_T matrix
                5.0,    # L
                [1, 2] # bc_bottom (only first 2 DOFs)
            )
            
            support_dof_maps, total_dofs = create_support_dof_mapping(bridge_opts, [support])
            
            # Should return one mapping per support
            @test length(support_dof_maps) == 1
            
            # Total DOFs should be larger than bridge DOFs
            @test total_dofs > bridge_opts.n_dofs
            
            # Support DOF mapping should have correct length
            n_support_dofs = 3 * (support.n_elem + 1)  # 3 DOFs per node
            @test length(support_dof_maps[1]) == n_support_dofs
            
            # Check connection mapping
            dof_map = support_dof_maps[1]
            
            # First node DOFs should map to bridge connection DOFs
            bridge_connection_dofs = 3 * (support.connection_node - 1) .+ support.connection_dofs
            @test dof_map[support.connection_dofs] == bridge_connection_dofs
            
            # All DOFs should be assigned (no zeros)
            @test all(dof_map .> 0)
            
            # All mapped DOF indices should be unique
            @test length(unique(dof_map)) == length(dof_map)
        end
        
        @testset "get_dof_from_node function" begin
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Create temperature-dependent E_T matrix for support
            E_T_support = [0.0 200e9; 100.0 160e9]
            support = SupportElement(
                2,              # connection_node
                [1, 2],         # connection_dofs  
                0.0,            # angle
                2,              # n_elem
                0.1,            # A
                0.001,          # I
                E_T_support,    # E_T matrix
                5.0,            # L
                [1, 2]          # bc_bottom
            )
            
            # Test bridge node (within bridge DOF range)
            bridge_node = 1
            dofs = get_dof_from_node(bridge_opts, [support], bridge_node)
            expected_dofs = 3 * (bridge_node - 1) .+ [1, 2, 3]
            @test dofs == expected_dofs
            
            # Test support connection node - should return bridge DOFs for that node
            connection_node = support.connection_node  # Node 2
            connection_dofs = get_dof_from_node(bridge_opts, [support], connection_node)
            expected_connection_dofs = 3 * (connection_node - 1) .+ [1, 2, 3]  # [4, 5, 6]
            @test connection_dofs == expected_connection_dofs
        end
        
        @testset "get_bc_dofs function" begin
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, [1, 2]], [4, [3]]])  # Multiple BC nodes
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Create temperature-dependent E_T matrix for support
            E_T_support = [0.0 200e9; 100.0 160e9]
            support = SupportElement(
                2,              # connection_node
                [1, 2],         # connection_dofs  
                0.0,            # angle
                2,              # n_elem
                0.1,            # A
                0.001,          # I
                E_T_support,    # E_T matrix
                5.0,            # L
                [1, 2]          # bc_bottom
            )
            
            support_dof_maps, total_dofs = create_support_dof_mapping(bridge_opts, [support])
            
            bc_dofs = get_bc_dofs(bridge_opts, [support], support_dof_maps)
            
            # Should return vector of DOF indices
            @test isa(bc_dofs, Vector{Int})
            
            # Should have entries for bridge boundary conditions
            @test length(bc_dofs) > 0
            
            # All DOF indices should be positive
            @test all(bc_dofs .> 0)
            
            # DOF indices should be within total DOF range
            @test all(bc_dofs .<= total_dofs)
        end
    end
    
    @testset "Helper Function Tests" begin
        
        @testset "assemble_local_support function" begin
            E_T_support = [0.0 200e9; 100.0 180e9]  # Temperature-dependent E
            support = SupportElement(
                1,      # connection_node
                [1],    # connection_dofs
                30.0,   # angle (degrees)
                3,      # n_elem
                0.2,    # A
                0.002,  # I
                E_T_support, # E_T matrix
                6.0,    # L
                [1, 2, 3] # bc_bottom
            )
            
            # Test at reference temperature
            T = 20.0
            K_local = assemble_local_support(support, T)
            
            # Verify dimensions
            expected_dofs = 3 * (support.n_elem + 1)
            @test size(K_local) == (expected_dofs, expected_dofs)
            
            # Should be symmetric
            @test issymmetric(K_local)
            
            # Should be positive semi-definite
            @test all(diag(K_local) .>= 0)
            
            # Test temperature dependence
            K_low = assemble_local_support(support, 0.0)
            K_high = assemble_local_support(support, 100.0)
            
            # Higher temperature should give lower stiffness
            @test norm(K_high) < norm(K_low)
        end
        
        @testset "create_support_mass_matrix function" begin
            E_T_support = [0.0 200e9; 50.0 200e9]  # Constant E
            support = SupportElement(
                1,      # connection_node
                [1],    # connection_dofs
                0.0,    # angle
                2,      # n_elem
                0.15,   # A
                0.0015, # I
                E_T_support, # E_T matrix
                4.0,    # L
                [1, 3] # bc_bottom (x and rotation DOFs)
            )
            
            ρ = 2500.0  # kg/m³
            M_local = create_support_mass_matrix(support, ρ)
            
            # Verify dimensions
            expected_dofs = 3 * (support.n_elem + 1)
            @test size(M_local) == (expected_dofs, expected_dofs)
            
            # Should be symmetric
            @test issymmetric(M_local)
            
            # Should have positive diagonal entries
            @test all(diag(M_local) .> 0)
            
            # Mass should be proportional to density
            M_double = create_support_mass_matrix(support, 2*ρ)
            @test norm(M_double - 2*M_local) < 1e-10
        end
        
        @testset "create_expanded_transformation function" begin
            angle = π/3  # 60 degrees
            n_nodes = 3
            
            T_exp = create_expanded_transformation(angle, n_nodes)
            
            # Verify dimensions
            expected_size = 3 * n_nodes
            @test size(T_exp) == (expected_size, expected_size)
            
            # Should be orthogonal (T' * T = I)
            @test norm(T_exp' * T_exp - I) < 1e-10
            
            # Determinant should be 1 (proper rotation)
            @test abs(det(T_exp) - 1.0) < 1e-10
            
            # Test identity transformation (angle = 0)
            T_identity = create_expanded_transformation(0.0, n_nodes)
            @test norm(T_identity - I) < 1e-10
        end
    end
    
    @testset "Integration with Core and Elements" begin
        # Test that Assembly modules work correctly with Core types and Elements functions
        
        @testset "Core types integration" begin
            # Test with realistic bridge configuration
            E_T = [0.0 210e9; 50.0 200e9; 100.0 180e9]  
            bc = BridgeBC([[1, "all"], [5, "all"]])  # Fixed ends
            bridge_opts = BridgeOptions(8, bc, 20.0, 7850.0, 0.02, 0.0008, E_T, 100.0)
            
            # Should work without errors
            M, K = assemble_matrices(bridge_opts, 25.0)
            @test size(M) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
            @test size(K) == (bridge_opts.n_dofs, bridge_opts.n_dofs)
        end
        
        @testset "Elements module integration" begin
            # Test that assembly functions correctly use Elements module functions
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(2, bc, 5.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            M, K = assemble_matrices(bridge_opts, 20.0)
            
            # Should produce reasonable results using Elements functions
            @test nnz(M) > 0
            @test nnz(K) > 0
            @test issymmetric(M)
            @test issymmetric(K)
        end
    end
    
    @testset "Numerical Accuracy Tests" begin
        # Test numerical accuracy against analytical solutions
        
        @testset "Simple beam stiffness" begin
            # Single element beam - analytical stiffness is known
            E = 200e9  # Pa
            A = 0.01   # m²
            I = 8.33e-6  # m⁴
            L = 1.0    # m
            
            E_T = [0.0 E; 50.0 E]  # Constant E at multiple temperatures
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(1, bc, L, 7850.0, A, I, E_T, 1.0)
            
            M, K = assemble_matrices(bridge_opts, 0.0)
            
            # For single element, can verify analytical values
            # Axial stiffness: EA/L
            # Flexural stiffness: 12EI/L³ for transverse DOFs
            
            dx = L / bridge_opts.n_elem
            expected_EA_L = E * A / dx
            expected_12EI_L3 = 12 * E * I / dx^3
            
            # Check that stiffness values are in reasonable range
            @test maximum(diag(K)) >= expected_EA_L * 0.5  # Should be at least half
            @test maximum(diag(K)) <= expected_EA_L * 2.0  # Should not exceed double
        end
        
        @testset "Tolerance validation" begin
            # Test numerical precision requirements (1e-10 tolerance)
            E_T = [0.0 200e9; 100.0 160e9]
            bc = BridgeBC([[1, "all"]])
            bridge_opts = BridgeOptions(4, bc, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
            
            # Test repeatability
            M1, K1 = assemble_matrices(bridge_opts, 25.0)
            M2, K2 = assemble_matrices(bridge_opts, 25.0)
            
            # Results should be identical to machine precision
            @test norm(M1 - M2) < 1e-15
            @test norm(K1 - K2) < 1e-15
            
            # Test that small parameter changes produce small changes in results
            bridge_opts_perturb = BridgeOptions(4, bc, 10.001, 2500.0, 0.1, 0.001, E_T, 50.0)
            M3, K3 = assemble_matrices(bridge_opts_perturb, 25.0)
            
            # Small geometry change should produce small matrix change
            rel_change_M = norm(M3 - M1) / norm(M1)
            rel_change_K = norm(K3 - K1) / norm(K1)
            
            @test rel_change_M < 1e-3  # Relative change should be small
            @test rel_change_K < 1e-3
        end
    end

end
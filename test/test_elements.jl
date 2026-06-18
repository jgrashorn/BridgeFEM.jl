# Unit tests for Elements module - Finite Element Computations
using Test
using BridgeFEM
using LinearAlgebra

@testset "Elements Module Tests" begin
    
    @testset "Frame Element Stiffness Matrix Tests" begin
        # Test 1: Basic functionality with known values
        @testset "Basic Stiffness Matrix" begin
            EA = 200e9 * 0.01  # E * A = 2e9 N
            EI = 200e9 * 1e-6  # E * I = 200 N⋅m²
            L = 2.0            # Length = 2 m
            
            K = frame_elem_stiffness(EA, EI, L)
            
            # Check matrix dimensions
            @test size(K) == (6, 6)
            
            # Check symmetry
            @test isapprox(K, K', atol=1e-10)
            
            # Check expected coefficient values for this specific case
            c1 = EA / L          # = 1e9
            c2 = 12.0 * EI / L^3 # = 12 * 200 / 8 = 300
            c3 = 6.0 * EI / L^2  # = 6 * 200 / 4 = 300
            c4 = 4.0 * EI / L    # = 4 * 200 / 2 = 400
            c5 = 2.0 * EI / L    # = 2 * 200 / 2 = 200
            
            # Test key matrix entries
            @test K[1,1] ≈ c1 atol=1e-6      # Axial stiffness
            @test K[2,2] ≈ c2 atol=1e-6      # Transverse stiffness
            @test K[3,3] ≈ c4 atol=1e-6      # Rotational stiffness
            @test K[1,4] ≈ -c1 atol=1e-6     # Axial coupling
            @test K[2,5] ≈ -c2 atol=1e-6     # Transverse coupling
            
            # Test zeros where expected
            @test K[1,2] ≈ 0.0 atol=1e-10    # No axial-transverse coupling
            @test K[1,3] ≈ 0.0 atol=1e-10    # No axial-rotational coupling
        end
        
        @testset "Unit Beam Stiffness (Analytical Comparison)" begin
            # Standard unit beam: EA = EI = L = 1
            EA = EI = L = 1.0
            K = frame_elem_stiffness(EA, EI, L)
            
            # Expected matrix for unit beam
            K_expected = [
                1.0   0.0   0.0   -1.0   0.0   0.0;
                0.0  12.0   6.0    0.0 -12.0   6.0;
                0.0   6.0   4.0    0.0  -6.0   2.0;
               -1.0   0.0   0.0    1.0   0.0   0.0;
                0.0 -12.0  -6.0    0.0  12.0  -6.0;
                0.0   6.0   2.0    0.0  -6.0   4.0
            ]
            
            @test isapprox(K, K_expected, atol=1e-10)
        end
        
        @testset "Positive Definiteness" begin
            # Stiffness matrix should be positive semi-definite for physical elements
            EA = 210e9 * 0.005  # Steel beam
            EI = 210e9 * 2e-5
            L = 3.0
            
            K = frame_elem_stiffness(EA, EI, L)
            
            # Check that eigenvalues are non-negative (allowing for rigid body modes)
            eigenvals = real(eigvals(K))
            @test all(eigenvals .>= -1e-7)  # Account for numerical precision in unconstrained elements
        end
    end
    
    @testset "Frame Element Mass Matrix Tests" begin
        @testset "Basic Mass Matrix" begin
            ρ = 7850.0    # Steel density kg/m³
            A = 0.01      # Cross-sectional area m²
            L = 2.0       # Length m
            
            M = frame_elem_mass(ρ, A, L)
            
            # Check matrix dimensions
            @test size(M) == (6, 6)
            
            # Check symmetry
            @test isapprox(M, M', atol=1e-10)
            
            # Check positive definiteness (mass matrix should be positive definite)
            eigenvals = real(eigvals(M))
            @test all(eigenvals .> 1e-10)
        end
        
        @testset "Unit Element Mass (Analytical Comparison)" begin
            # Unit element: ρ = A = L = 1
            ρ = A = L = 1.0
            M = frame_elem_mass(ρ, A, L)
            
            mass_factor = ρ * A * L / 420.0  # = 1/420
            
            # Expected matrix for unit element
            M_expected = mass_factor * [
                140   0      0     70    0      0    ;
                0     156    22    0     54    -13   ;
                0     22     4     0     13    -3    ;
                70    0      0     140   0      0    ;
                0     54     13    0     156   -22   ;
                0    -13    -3     0    -22     4
            ]
            
            @test isapprox(M, M_expected, atol=1e-12)
        end
        
        @testset "Mass Conservation" begin
            # Total mass should equal ρ * A * L
            ρ = 2500.0  # Concrete density
            A = 0.02    # Area
            L = 4.0     # Length
            
            M = frame_elem_mass(ρ, A, L)
            
            # For consistent mass matrix, sum of translational entries gives total mass
            total_mass_x = M[1,1] + M[1,4] + M[4,1] + M[4,4]  # Axial DOFs
            total_mass_y = M[2,2] + M[2,5] + M[5,2] + M[5,5]  # Transverse DOFs
            
            expected_total_mass = ρ * A * L
            
            @test isapprox(total_mass_x, expected_total_mass, rtol=1e-10)
            @test isapprox(total_mass_y, expected_total_mass, rtol=1e-10)
        end
    end
    
    @testset "Transformation Matrix Tests" begin
        @testset "Identity Transformation (0 degrees)" begin
            T = transformation_matrix(0.0)
            
            # Should be identity matrix for 0 rotation
            @test size(T) == (6, 6)
            @test isapprox(T, Matrix{Float64}(I, 6, 6), atol=1e-12)
        end
        
        @testset "90 Degree Rotation" begin
            T = transformation_matrix(90.0)
            
            # Check specific values for 90° rotation
            @test T[1,1] ≈ 0.0 atol=1e-12    # cos(90°) = 0
            @test T[1,2] ≈ -1.0 atol=1e-12   # -sin(90°) = -1
            @test T[2,1] ≈ 1.0 atol=1e-12    # sin(90°) = 1
            @test T[2,2] ≈ 0.0 atol=1e-12    # cos(90°) = 0
            
            # Same for nodes 2 (DOFs 4,5)
            @test T[4,4] ≈ 0.0 atol=1e-12
            @test T[4,5] ≈ -1.0 atol=1e-12
            @test T[5,4] ≈ 1.0 atol=1e-12
            @test T[5,5] ≈ 0.0 atol=1e-12
            
            # Rotational DOFs (3,6) should remain unchanged
            @test T[3,3] ≈ 1.0 atol=1e-12
            @test T[6,6] ≈ 1.0 atol=1e-12
        end
        
        @testset "45 Degree Rotation" begin
            T = transformation_matrix(45.0)
            sqrt2_half = √2 / 2
            
            # Check trigonometric values (corrected rotation direction)
            @test T[1,1] ≈ sqrt2_half atol=1e-12   # cos(45°)
            @test T[1,2] ≈ -sqrt2_half atol=1e-12  # -sin(45°)
            @test T[2,1] ≈ sqrt2_half atol=1e-12   # sin(45°)
            @test T[2,2] ≈ sqrt2_half atol=1e-12   # cos(45°)
        end
        
        @testset "Orthogonality Properties" begin
            angles = [0.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0, 270.0]
            
            for θ in angles
                T = transformation_matrix(θ)
                
                # Transformation matrix should be orthogonal: T' * T = I
                @test isapprox(T' * T, Matrix{Float64}(I, 6, 6), atol=1e-10)
                
                # Determinant should be 1 (proper rotation, no reflection)
                @test isapprox(det(T), 1.0, atol=1e-10)
            end
        end
        
        @testset "Vector Transformation" begin
            # Test transformation of a unit vector along x-axis
            θ = 30.0
            T = transformation_matrix(θ)
            
            # Unit vector in local x-direction at node 1: [1, 0, 0, 0, 0, 0]
            v_local = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            v_global = T * v_local
            
            # Should give [cos(30°), sin(30°), 0, 0, 0, 0] in global coordinates
            # (unit vector along local x-axis rotated 30° counterclockwise)
            @test v_global[1] ≈ cosd(30.0) atol=1e-12
            @test v_global[2] ≈ sind(30.0) atol=1e-12
            @test v_global[3] ≈ 0.0 atol=1e-12
            @test v_global[4] ≈ 0.0 atol=1e-12
            @test v_global[5] ≈ 0.0 atol=1e-12
            @test v_global[6] ≈ 0.0 atol=1e-12
        end
    end
    
    @testset "Function Interface Consistency" begin
        # Test that function signatures match exactly with original implementation
        @testset "Parameter Types and Return Types" begin
            # Test frame_elem_stiffness
            EA, EI, L_e = 1e6, 1e3, 2.0
            K = frame_elem_stiffness(EA, EI, L_e)
            @test isa(K, Matrix{Float64})
            @test size(K) == (6, 6)
            
            # Test frame_elem_mass with exact type annotations
            ρ, A, L = 7850.0, 0.01, 2.0
            M = frame_elem_mass(ρ, A, L)
            @test isa(M, Matrix{Float64})
            @test size(M) == (6, 6)
            
            # Test transformation_matrix
            θ = 45.0
            T = transformation_matrix(θ)
            @test isa(T, Matrix{Float64})
            @test size(T) == (6, 6)
        end
    end
end
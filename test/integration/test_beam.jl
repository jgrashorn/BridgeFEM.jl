# Simply Supported Beam Static Analysis Test  
# Converted from tests/simply_supported_beam.jl to Test.jl framework
# Updated to use proper modular structure

using Test
using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

# Use proper BridgeFEM module import
using BridgeFEM

@testset "Simply Supported Beam Static Analysis" begin

    @testset "Setup and Model Configuration" begin
        # Beam and material parameters
        L = 10.0               # Beam length (m)
        n_elem = 12           # Number of finite elements
        n_node = n_elem + 1   # Number of nodes
        ρ = 7800.0            # Density (kg/m^3)
        A = .1                # Cross-section area (m^2)
        I = .001              # Moment of inertia (m^4)
        E0 = 207e9            # Base Young's modulus (Pa)
        α = -1e5              # E-temperature slope (Pa/K)
        cutoff_freq = 500.0   # Cutoff frequency for modes (Hz)

        # Simply supported boundary conditions (different from cantilever)
        bc = BridgeBC([
            [1, "trans"],      # Node 1: translational DOFs constrained
            [n_node, "y"]      # Last node: vertical DOF constrained
        ])

        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)
        supports = SupportElement[]  # No additional supports
        Ts = [10.0, 20.0]

        sim_opts = SimulationOptions(bo, supports, collect(Ts), damping_ratio=0.02)

        @test isa(bo, BridgeOptions)
        @test isa(sim_opts, SimulationOptions)
        @test length(Ts) == 2
        @test n_elem == 12
        @test n_node == 13
        @test length(supports) == 0  # No additional supports for simply supported beam
    end

    @testset "Matrix Assembly and System Setup" begin
        # Setup (repeated for test isolation)
        L = 10.0; n_elem = 12; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .001; E0 = 207e9; α = -1e5; cutoff_freq = 500.0
        
        bc = BridgeBC([
            [1, "trans"],
            [n_node, "y"]
        ])
        
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)
        supports = SupportElement[]
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, supports, collect(Ts), damping_ratio=0.02)

        M_T, K_T = setup_physical(sim_opts)

        @test isa(M_T, Function) || isa(M_T, Interpolations.AbstractInterpolation)
        @test isa(K_T, Function) || isa(K_T, Interpolations.AbstractInterpolation)
        
        # Test matrix evaluation at specific temperature
        M_test = M_T(20.0)
        K_test = K_T(20.0)
        
        @test isa(M_test, Array)
        @test isa(K_test, Array)
        @test size(M_test) == size(K_test)
        @test size(M_test)[1] == sim_opts.total_dofs
    end

    @testset "Force Application and Boundary Conditions" begin
        # Setup
        L = 10.0; n_elem = 12; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .001; E0 = 207e9; α = -1e5; cutoff_freq = 500.0
        
        bc = BridgeBC([
            [1, "trans"],
            [n_node, "y"]
        ])
        
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)
        supports = SupportElement[]
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, supports, collect(Ts), damping_ratio=0.02)

        M_T, K_T = setup_physical(sim_opts)

        # Apply force at center of beam (different from cantilever)
        f0 = -10000.0
        f = zeros(n_node*3)
        f[(n_node ÷ 2)*3 + 2] = f0  # Center node, vertical DOF
        
        M_, K_ = apply_bc(M_T(20.0), K_T(20.0), sim_opts)

        @test f[(n_node ÷ 2)*3 + 2] == f0
        @test sum(abs.(f)) == abs(f0)  # Only one non-zero force component
        @test isa(M_, Matrix)
        @test isa(K_, Matrix)
        
        # Test that boundary conditions maintain matrix size (penalty method approach)
        @test size(M_) == size(M_T(20.0))
        @test size(K_) == size(K_T(20.0))
    end

    @testset "FEM Solution vs Analytical Simply Supported Beam Theory" begin
        # Complete setup and solution
        L = 10.0; n_elem = 12; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .001; E0 = 207e9; α = -1e5; cutoff_freq = 500.0
        
        bc = BridgeBC([
            [1, "trans"],
            [n_node, "y"]
        ])
        
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)
        supports = SupportElement[]
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, supports, collect(Ts), damping_ratio=0.02)

        M_T, K_T = setup_physical(sim_opts)

        # Force at center
        f0 = -10000.0
        f = zeros(n_node*3)
        f[(n_node ÷ 2)*3 + 2] = f0

        M_, K_ = apply_bc(M_T(20.0), K_T(20.0), sim_opts)
        u = K_ \ f

        # Analytical solution for simply supported beam with center loading
        x = range(0, L, length=n_node)
        u_analytical = zeros(n_node)
        for (i, x_) in enumerate(x)
            if x_ <= L/2
                u_analytical[i] = f0 * x_ / (48 * E0 * I) * (3 * L^2 - 4.0 * x_^2)
            else
                u_analytical[i] = f0 * (L-x_) / (48 * E0 * I) * (3 * L^2 - 4.0 * (L-x_)^2)
            end
        end

        # Test FEM solution against analytical simply supported beam theory
        u_fem = u[2:3:end]  # Extract vertical displacements
        
        @test isapprox(u_analytical, u_fem, atol=1e-4)
        @test maximum(abs.(u_analytical - u_fem)) < 1e-4
        
        # Test maximum deflection occurs at center
        center_idx = n_node ÷ 2 + 1
        @test abs(u_analytical[center_idx]) == maximum(abs.(u_analytical))
        @test abs(u_fem[center_idx]) ≈ maximum(abs.(u_fem)) atol=1e-10
        
        # Test symmetry of deflection
        for i in 1:n_node÷2
            @test abs(u_analytical[i] - u_analytical[end-i+1]) < 1e-12
            @test abs(u_fem[i] - u_fem[end-i+1]) < 1e-10
        end
    end

    @testset "Beam Properties and Boundary Condition Validation" begin
        # Setup
        L = 10.0; n_elem = 12; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .001; E0 = 207e9; α = -1e5; cutoff_freq = 500.0
        
        bc = BridgeBC([
            [1, "trans"],
            [n_node, "y"]
        ])
        
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)
        supports = SupportElement[]
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, supports, collect(Ts), damping_ratio=0.02)

        M_T, K_T = setup_physical(sim_opts)

        f0 = -10000.0
        f = zeros(n_node*3)
        f[(n_node ÷ 2)*3 + 2] = f0

        M_, K_ = apply_bc(M_T(20.0), K_T(20.0), sim_opts)
        u = K_ \ f

        u_y = u[2:3:end]  # Extract vertical displacements

        # Test boundary conditions are satisfied
        # For simply supported beam: vertical displacement at ends should be zero
        @test abs(u_y[1]) < 1e-12    # First node vertical displacement ≈ 0
        @test abs(u_y[end]) < 1e-12  # Last node vertical displacement ≈ 0
        
        # Test that deflection has correct sign (downward for negative force)
        @test all(u_y .<= 1e-12)  # All deflections should be ≤ 0 (downward)
        
        # Test maximum deflection at center
        center_idx = n_node ÷ 2 + 1
        @test abs(u_y[center_idx]) == maximum(abs.(u_y))
        
        # Test analytical maximum deflection formula
        # Maximum deflection for simply supported beam with center load: δ_max = -PL³/(48EI)
        delta_max_analytical = abs(f0 * L^3 / (48 * E0 * I))
        delta_max_fem = abs(u_y[center_idx])
        
        @test abs(delta_max_fem - delta_max_analytical) / delta_max_analytical < 0.01  # Within 1%
    end
end 
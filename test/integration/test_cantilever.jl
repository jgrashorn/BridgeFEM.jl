# Fixed Cantilever Static Analysis Test
# Converted from tests/fixed_cantilever.jl to Test.jl framework
# Maintains current include patterns during transition period

using Test
using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

# Include current monolithic structure
include("../../src/bridge_model.jl")
include("../../src/model_reduction.jl")
include("../../src/utils.jl")
include("../../src/dynamic_simulation.jl")

@testset "Fixed Cantilever Static Analysis" begin
    
    @testset "Setup and Model Configuration" begin
        # Beam and material parameters
        L = 10.0               # Beam length (m)
        n_elem = 11           # Number of finite elements
        n_node = n_elem + 1   # Number of nodes
        ρ = 7800.0            # Density (kg/m^3)
        A = .1                # Cross-section area (m^2)
        I = .0001             # Moment of inertia (m^4)
        E0 = 207e9            # Base Young's modulus (Pa)
        α = -1e5              # E-temperature slope (Pa/K)
        cutoff_freq = 100.0   # Cutoff frequency for modes (Hz)

        bc = BridgeBC([  # Node 1: both translational and rotational DOFs fixed
            [1, "all"]
        ])

        E_bridge = [
            10.0 E0;
            20.0  E0
        ]

        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)

        # Support element parameters with temperature dependence
        Ts = [10.0, 20.0]

        # Create comprehensive simulation options
        sim_opts = SimulationOptions(
            bo, collect(Ts), damping_ratio=0.02
        )

        @test isa(bo, BridgeOptions)
        @test isa(sim_opts, SimulationOptions)
        @test length(Ts) == 2
        @test n_elem == 11
        @test n_node == 12
    end

    @testset "Matrix Assembly and Decomposition" begin
        # Beam and material parameters (repeated for test isolation)
        L = 10.0; n_elem = 11; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .0001; E0 = 207e9; α = -1e5; cutoff_freq = 100.0
        
        bc = BridgeBC([[1, "all"]])
        E_bridge = [10.0 E0; 20.0 E0]
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, collect(Ts), damping_ratio=0.02)

        M, K = assemble_matrices_with_supports(sim_opts)
        M_, K_ = apply_bc(M, K, sim_opts)

        @test isa(M, Array)
        @test isa(K, Array)
        @test size(M) == size(K)
        @test size(M_) == size(K_)
        
        # Test that boundary conditions maintain matrix size (penalty method approach)
        @test size(M_) == size(M)
        @test size(K_) == size(K)
    end

    @testset "Modal Analysis and Force Application" begin
        # Setup (repeated for test isolation)
        L = 10.0; n_elem = 11; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .0001; E0 = 207e9; α = -1e5; cutoff_freq = 100.0
        
        bc = BridgeBC([[1, "all"]])
        E_bridge = [10.0 E0; 20.0 E0]
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, collect(Ts), damping_ratio=0.02)

        M, K = assemble_matrices_with_supports(sim_opts)
        M_, K_ = apply_bc(M, K, sim_opts)
        
        _, _, λs, vectors, vectors_unnormalized = assemble_and_decompose(sim_opts)

        @test isa(λs, Array)
        @test isa(vectors, Array)
        @test size(vectors, 3) == length(Ts)
        @test all(λs[:, 1] .> 0)  # All eigenvalues should be positive
        
        # Force at free end
        f0 = -10000.0
        f = zeros(n_node*3)
        f[end-1] = f0

        @test f[end-1] == f0
        @test sum(abs.(f)) == abs(f0)  # Only one non-zero force component
    end

    @testset "FEM Solution vs Analytical Beam Theory" begin
        # Complete setup and solution
        L = 10.0; n_elem = 11; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .0001; E0 = 207e9; α = -1e5; cutoff_freq = 100.0
        
        bc = BridgeBC([[1, "all"]])
        E_bridge = [10.0 E0; 20.0 E0]
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, collect(Ts), damping_ratio=0.02)

        M, K = assemble_matrices_with_supports(sim_opts)
        M_, K_ = apply_bc(M, K, sim_opts)
        _, _, λs, vectors, vectors_unnormalized = assemble_and_decompose(sim_opts)

        # Force application
        f0 = -10000.0
        f = zeros(n_node*3)
        f[end-1] = f0

        # Modal force
        f_m = vectors[:, :, 1]' * f

        # Solve for displacement
        u = K_[:,:,1] \ f

        # Solve for modal displacement
        q = f_m ./ (2π .* λs[:,1]).^2
        u_m = vectors[:,:,1] * q

        # Analytical solutions for cantilever beam
        x = range(0, L, length=n_node)
        u_analytical = f0 .* x.^2 ./ (6 .*E0.*I) .* (3 .*L .- x)
        ϕ_analytical = f0 .* x ./ (2 .*E0.*I) .* (2 .*L .- x)

        # Test FEM solution against analytical beam theory
        @test isapprox(u_analytical, u[2:3:end], atol=1e-4)
        @test isapprox(ϕ_analytical, u[3:3:end], atol=1e-4)
        
        # Test modal solution against analytical beam theory  
        @test isapprox(u_analytical, u_m[2:3:end], atol=1e-4)
        @test isapprox(ϕ_analytical, u_m[3:3:end], atol=1e-4)
        
        # Additional validation tests
        @test maximum(abs.(u_analytical - u[2:3:end])) < 1e-4
        @test maximum(abs.(ϕ_analytical - u[3:3:end])) < 1e-4
        @test maximum(abs.(u_analytical - u_m[2:3:end])) < 1e-4
        @test maximum(abs.(ϕ_analytical - u_m[3:3:end])) < 1e-4
    end

    @testset "Solution Properties Validation" begin
        # Verify solution properties
        L = 10.0; n_elem = 11; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .0001; E0 = 207e9; α = -1e5; cutoff_freq = 100.0
        
        bc = BridgeBC([[1, "all"]])
        E_bridge = [10.0 E0; 20.0 E0]
        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, collect(Ts), damping_ratio=0.02)

        M, K = assemble_matrices_with_supports(sim_opts)
        M_, K_ = apply_bc(M, K, sim_opts)
        _, _, λs, vectors, vectors_unnormalized = assemble_and_decompose(sim_opts)

        f0 = -10000.0
        f = zeros(n_node*3)
        f[end-1] = f0
        u = K_[:,:,1] \ f

        # Test boundary conditions are satisfied (fixed end)
        @test abs(u[2]) < 1e-12  # Vertical displacement at fixed end should be ~0
        @test abs(u[3]) < 1e-12  # Rotation at fixed end should be ~0
        
        # Test that deflection increases along beam length
        u_y = u[2:3:end]  # Extract vertical displacements
        @test all(diff(u_y) .< 0)  # Downward deflection increases (becomes more negative)
        
        # Test maximum deflection is at free end
        @test abs(u_y[end]) == maximum(abs.(u_y))
    end
end 
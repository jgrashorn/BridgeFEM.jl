# Simply Supported Beam Dynamic Analysis Test
# Converted from tests/simply_supported_beam_dynamic.jl to Test.jl framework
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

@testset "Simply Supported Beam Dynamic Analysis" begin

    @testset "Dynamic Model Setup and Configuration" begin
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

        # Simply supported boundary conditions
        bc = BridgeBC([
            [1, "trans"],      # Node 1: translational DOFs constrained  
            [n_node, "y"]      # Last node: vertical DOF constrained
        ])

        bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)
        supports = SupportElement[]
        Ts = [10.0, 20.0]
        sim_opts = SimulationOptions(bo, supports, collect(Ts), damping_ratio=0.02)

        @test isa(bo, BridgeOptions)
        @test isa(sim_opts, SimulationOptions)
        @test cutoff_freq > 100.0  # Higher cutoff for dynamic analysis
        @test length(Ts) == 2
        @test sim_opts.damping_ratio == 0.02
        @test length(supports) == 0  # No additional supports
    end

    @testset "Dynamic System Setup and Modal Analysis" begin
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
        λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

        @test isa(M_T, Function) || isa(M_T, Interpolations.AbstractInterpolation)
        @test isa(K_T, Function) || isa(K_T, Interpolations.AbstractInterpolation)
        @test isa(λ_T, Function) || isa(λ_T, Interpolations.AbstractInterpolation)
        @test isa(Φ_T, Function) || isa(Φ_T, Interpolations.AbstractInterpolation)
        @test n_modes > 0
        @test n_modes <= sim_opts.total_dofs
        
        # Test that we have valid eigenfrequencies
        λ_test = λ_T(20.0)  # Get eigenvalues at 20°C
        @test all(λ_test[:, 1] .> 0)  # All eigenvalues should be positive
        @test n_modes <= length(λ_test[:, 1])  # Number of modes should not exceed available eigenvalues
    end

    @testset "Analytical Solution and Force Setup" begin
        # Setup
        L = 10.0; n_elem = 12; n_node = n_elem + 1
        ρ = 7800.0; A = .1; I = .001; E0 = 207e9; α = -1e5; cutoff_freq = 500.0
        
        # Force parameters
        force_dof = (n_node ÷ 2)*3 + 2  # Center node, vertical DOF
        f0 = -10000.0

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

        @test force_dof == (n_node ÷ 2)*3 + 2  # Force at center
        @test length(u_analytical) == n_node
        @test u_analytical[1] ≈ 0.0 atol=1e-12  # Zero at left support
        @test u_analytical[end] ≈ 0.0 atol=1e-12  # Zero at right support
        
        # Test symmetry of analytical solution
        center_idx = n_node ÷ 2 + 1
        @test abs(u_analytical[center_idx]) == maximum(abs.(u_analytical))
        
        # Test symmetry
        for i in 1:n_node÷2
            @test abs(u_analytical[i] - u_analytical[end-i+1]) < 1e-12
        end
    end

    @testset "Dynamic Response Parameters and Initial Conditions" begin
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
        λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

        # Force and initial conditions
        force_dof = (n_node ÷ 2)*3 + 2
        f0 = -10000.0

        temp = t -> 20.0  # Constant temperature function
        load_vector = (t, dof) -> begin 
            f = zeros(Float64, length(dof))
            f[force_dof] = f0
            return f
        end

        # Initial conditions in state space
        u0 = zeros(sim_opts.total_dofs * 2)
        
        u0_modal_ = Φ_T(temp(0.0))' * M_T(temp(0.0)) * u0[1:sim_opts.total_dofs]
        du0_modal_ = Φ_T(temp(0.0))' * u0[sim_opts.total_dofs+1:end]
        u0_modal = [u0_modal_; du0_modal_]

        @test length(u0) == 2 * sim_opts.total_dofs
        @test all(u0 .== 0.0)  # Zero initial conditions
        @test length(u0_modal) == 2 * n_modes
        @test all(u0_modal .== 0.0)  # Zero modal initial conditions
        
        # Test load vector
        test_dofs = collect(1:sim_opts.total_dofs)
        f_test = load_vector(0.0, test_dofs)
        @test length(f_test) == sim_opts.total_dofs
        @test f_test[force_dof] == f0
        @test sum(abs.(f_test)) == abs(f0)  # Only one non-zero force component
    end

    @testset "Modal vs Physical Space Dynamic Response" begin
        # Complete setup
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
        λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

        # Force and parameters
        force_dof = (n_node ÷ 2)*3 + 2
        f0 = -10000.0

        temp = t -> 20.0
        load_vector = (t, dof) -> begin 
            f = zeros(Float64, length(dof))
            f[force_dof] = f0
            return f
        end

        u0 = zeros(sim_opts.total_dofs * 2)
        u0_modal_ = Φ_T(temp(0.0))' * M_T(temp(0.0)) * u0[1:sim_opts.total_dofs]
        du0_modal_ = Φ_T(temp(0.0))' * u0[sim_opts.total_dofs+1:end]
        u0_modal = [u0_modal_; du0_modal_]

        # Rayleigh damping
        α, β = 0.001, 0.001
        C = α * M_T(temp(20.0)) + β * K_T(temp(20.0))
        C_modal = Φ_T(20.0)' * C * Φ_T(20.0)
        ζ = diag(C_modal)

        tspan = (0.0, 1.0)  # Shorter simulation for test

        # Modal space problem
        prob_modal = ODEProblem(beam_modal_ode!, u0_modal, tspan,
                            (; 
                                n_modes=n_modes,
                                n_dofs=sim_opts.total_dofs,
                                T_func = temp,
                                λ_interp = λ_T,
                                Φ_interp = Φ_T,
                                ζ = ζ,
                                load_vector = load_vector,
                            ))

        # Physical space problem
        prob_physical = ODEProblem(beam_physical_ode!, u0, tspan,
                            (; 
                                n_dofs=sim_opts.total_dofs,
                                bc_dofs=sim_opts.bc_dofs,
                                T_func = temp,
                                M_interp = M_T,
                                K_interp = K_T,
                                α = α,
                                β = β,
                                load_vector = load_vector,
                            ))

        @test isa(prob_modal, ODEProblem)
        @test isa(prob_physical, ODEProblem)
        
        # Solve both problems
        sol_modal = solve(prob_modal, saveat=0.01)
        sol_physical = solve(prob_physical, saveat=0.01)

        @test sol_modal.retcode == :Success
        @test sol_physical.retcode == :Success
        @test length(sol_modal.t) == length(sol_physical.t)
        
        # Reconstruct physical response from modal solution
        q = reduce(hcat, sol_modal.u)
        u_modal, du_modal = reconstruct_physical(sim_opts, q, Φ_T, temp, sol_modal.t)
        
        # Extract physical response
        u_physical_raw = reduce(hcat, sol_physical.u)
        u_physical = u_physical_raw[1:sim_opts.total_dofs, :]
        du_physical = u_physical_raw[sim_opts.total_dofs+1:end, :]

        @test size(u_modal) == (sim_opts.total_dofs, length(sol_modal.t))
        @test size(u_physical) == (sim_opts.total_dofs, length(sol_physical.t))
        
        # Compare responses at center (force_dof)
        u_modal_center = u_modal[force_dof, :]
        u_physical_center = u_physical[force_dof, :]
        
        # Test that modal and physical responses are reasonably close
        # (allowing for modal truncation effects)
        relative_error = maximum(abs.(u_modal_center .- u_physical_center)) / maximum(abs.(u_physical_center))
        @test relative_error < 0.15  # Within 15% (modal approximation for simply supported beam)
        
        # Test that both responses show dynamic behavior (non-zero displacement)
        @test maximum(abs.(u_modal_center)) > 1e-6
        @test maximum(abs.(u_physical_center)) > 1e-6
    end

    @testset "Eigenfrequency Validation and Simply Supported Beam Properties" begin
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
        λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

        # Test eigenfrequencies at different temperatures
        for temp_val in [10.0, 20.0]
            λ = λ_T(temp_val)
            frequencies = sqrt.(λ[:, 1]) / (2*π)  # Convert to Hz
            
            @test all(frequencies .> 0)  # All frequencies should be positive
            @test issorted(frequencies)  # Frequencies should be sorted
            @test frequencies[1] < cutoff_freq  # At least first mode below cutoff
            
            # Test analytical simply supported beam first frequency (approximate)
            # First natural frequency of simply supported beam: f1 = (π²/2) * sqrt(EI/(ρAL^4))
            E_test = E0  # Use base modulus for analytical comparison
            f1_analytical = (π^2/2) * sqrt(E_test * I / (ρ * A * L^4)) / (2*π)
            
            # Allow for numerical differences (finite element vs analytical)
            @test abs(frequencies[1] - f1_analytical) / f1_analytical < 0.2  # Within 20% (simply supported beam)
        end
        
        # Test that mode shapes are orthogonal
        Φ = Φ_T(20.0)
        M = M_T(20.0)
        modal_mass = Φ' * M * Φ
        
        # Modal mass matrix should be approximately diagonal
        off_diagonal_norm = norm(modal_mass - diagm(diag(modal_mass)))
        diagonal_norm = norm(diagm(diag(modal_mass)))
        @test off_diagonal_norm / diagonal_norm < 0.01  # Off-diagonal terms < 1% of diagonal
    end

    @testset "Steady State Response Validation" begin
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
        λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

        # Force parameters
        force_dof = (n_node ÷ 2)*3 + 2
        f0 = -10000.0

        temp = t -> 20.0
        load_vector = (t, dof) -> begin 
            f = zeros(Float64, length(dof))
            f[force_dof] = f0
            return f
        end

        u0 = zeros(sim_opts.total_dofs * 2)
        u0_modal_ = Φ_T(temp(0.0))' * M_T(temp(0.0)) * u0[1:sim_opts.total_dofs]
        du0_modal_ = Φ_T(temp(0.0))' * u0[sim_opts.total_dofs+1:end]
        u0_modal = [u0_modal_; du0_modal_]

        α, β = 0.001, 0.001
        C = α * M_T(temp(20.0)) + β * K_T(temp(20.0))
        C_modal = Φ_T(20.0)' * C * Φ_T(20.0)
        ζ = diag(C_modal)

        tspan = (0.0, 2.0)  # Longer simulation to reach steady state

        prob_modal = ODEProblem(beam_modal_ode!, u0_modal, tspan,
                            (; 
                                n_modes=n_modes,
                                n_dofs=sim_opts.total_dofs,
                                T_func = temp,
                                λ_interp = λ_T,
                                Φ_interp = Φ_T,
                                ζ = ζ,
                                load_vector = load_vector,
                            ))

        sol_modal = solve(prob_modal, saveat=0.01)
        q = reduce(hcat, sol_modal.u)
        u_modal, du_modal = reconstruct_physical(sim_opts, q, Φ_T, temp, sol_modal.t)

        # Analytical solution for steady state
        x = range(0, L, length=n_node)
        u_analytical = zeros(n_node)
        for (i, x_) in enumerate(x)
            if x_ <= L/2
                u_analytical[i] = f0 * x_ / (48 * E0 * I) * (3 * L^2 - 4.0 * x_^2)
            else
                u_analytical[i] = f0 * (L-x_) / (48 * E0 * I) * (3 * L^2 - 4.0 * (L-x_)^2)
            end
        end

        # Compare final modal solution with analytical steady state
        u_modal_final = u_modal[2:3:end, end]  # Vertical displacements at final time
        
        @test isapprox(u_analytical, u_modal_final, atol=1e-4)
        @test maximum(abs.(u_analytical - u_modal_final)) < 1e-4
        
        # Test that system reaches steady state (derivative approaches zero)
        du_modal_final = du_modal[2:3:end, end]  # Final velocities
        @test maximum(abs.(du_modal_final)) < 1e-3  # Small final velocities indicate steady state
    end
end 
# Dynamics/simulation.jl - Dynamic simulation and ODE solving functions
# Extracted from dynamic_simulation.jl

using LinearAlgebra
using DifferentialEquations

# Import types from Core module - assumes this is included within BridgeFEM module
# where Core types are already available

"""
    beam_modal_ode!(du, u, p, t)

Modal space ODE function for beam dynamic analysis.
"""
function beam_modal_ode!(du, u, p, t)
    T = p.T_func(t)
    # q     = u[1:p.n_modes]        # modal displacements
    # qdot  = u[p.n_modes+1:end]    # modal velocities
    # # Interpolate natural frequencies at current T (convert Hz to rad/s)
    # ω = 2π .* p.λ_interp(p.T_func(t))
    # # Damping ratios (could be constant, interpolated, or Rayleigh-like)
    # ζ = p.ζ                      # vector of damping ratios per mode
    # # Interpolate mode shapes at T
    # Φ = p.Φ_interp(p.T_func(t))
    # # Assemble global load vector
    # f = p.load_vector(t,1:p.n_dofs)
    # # Project force onto each mode
    # fhat = p.Φ_interp(p.T_func(t))' * p.load_vector(t,1:p.n_dofs)
    # Modal accelerations
    # qddot = -p.ζ .* u[p.n_modes+1:end] .- (2π .* p.λ_interp(p.T_func(t)) .^2) .* u[1:p.n_modes] .+ p.Φ_interp(p.T_func(t))' * p.load_vector(t,1:p.n_dofs)
    # Fill derivative vector
    du[1:p.n_modes] .= u[p.n_modes+1:end]
    du[p.n_modes+1:end] .= -p.ζ .* u[p.n_modes+1:end] .- ((2π .* p.λ_interp(T)) .^2) .* u[1:p.n_modes] .+ p.Φ_interp(T)' * p.load_vector(t,1:p.n_dofs)
end

"""
    beam_physical_ode!(du, u, p, t)

Physical space ODE function for beam dynamic analysis.
"""
function beam_physical_ode!(du, u, p, t)
    # Unpack state vector
    du .= 0.0  # Reset derivative vector
    u_ = u[1:p.n_dofs]        # modal displacements
    udot_ = u[p.n_dofs+1:end] # modal velocities

    # Interpolate natural frequencies and mode shapes
    M = p.M_interp(p.T_func(t))
    K = p.K_interp(p.T_func(t))

    M_ff, K_ff, retained, removed = remove_fixed_dofs(M, K, p.bc_dofs, p.n_dofs)

    M_fc = M[retained, removed]
    K_fc = K[retained, removed]

    f_add_fc = K_fc * u_[removed]

    u_dofs = [retained; p.n_dofs .+ retained]
    n_retained = length(retained)

    # Assemble global load vector
    f = p.load_vector(t, 1:p.n_dofs)
    f = f[retained]
    f .+= f_add_fc

    # Construct system matrix blocks
    Z = zeros(n_retained, n_retained)
    I = Matrix{Float64}(LinearAlgebra.I, n_retained, n_retained)

    Minv = pinv(M_ff)
    # Rayleigh-like modal damping (if needed, adjust as appropriate)
    D = p.α .* M_ff + p.β .* K_ff

    A = [Z I; -Minv*K_ff -Minv*D]

    b = [zeros(n_retained); Minv*f]
    du_ff = A * u[u_dofs] + b

    du[u_dofs] .= du_ff
end

function solve_modal_simulation(sim_opts::SimulationOptions, u0::Vector{Float64}, tspan::Tuple{Float64, Float64}, temp_func::Function, load_vector::Function; saveat=0.01, alg = Rodas4(autodiff=AutoFiniteDiff()))

    λ_T, Φ_T, n_modes = setup_ROM(sim_opts)
    M_T, K_T = setup_physical(sim_opts)

    u0_modal_ = Φ_T(temp_func(0))' * M_T(temp_func(0)) * u0[1:sim_opts.total_dofs]
    du0_modal_ = Φ_T(temp_func(0))' * M_T(temp_func(0)) * u0[sim_opts.total_dofs+1:end]
    u0_modal = [u0_modal_; du0_modal_]

    α = sim_opts.damping_ratio
    β = 0.1*α

    C = α * M_T(temp_func(0)) + β * K_T(temp_func(0))  # Rayleigh damping matrix
    C_modal = Φ_T(temp_func(0))' * C * Φ_T(temp_func(0))  # Project to modal space
    ζ = diag(C_modal)

    prob_modal = ODEProblem(beam_modal_ode!, u0_modal, tspan,
                    (; 
                        n_modes=n_modes,
                        n_dofs=sim_opts.total_dofs,
                        T_func = temp_func,
                        λ_interp = λ_T,
                        ζ = ζ,  # Constant damping
                        Φ_interp = Φ_T,
                        load_vector = load_vector,
                    ))

    @info "Solving dynamic response with $(n_modes) modes in modal space..."
    @time u_ = solve(prob_modal, alg, saveat=saveat);
    q = reduce(hcat, u_.u)
    u, du = reconstruct_physical(sim_opts, q, Φ_T, temp_func, u_.t)

    return u, du, u_.t, temp_func.(u_.t)

end

function solve_physical_simulation(sim_opts::SimulationOptions, u0::Vector{Float64}, tspan::Tuple{Float64, Float64}, temp_func::Function, load_vector::Function; saveat=0.01)

    M_T, K_T = setup_physical(sim_opts)

    α = sim_opts.damping_ratio
    β = 0.1*α

    prob_physical = ODEProblem(beam_physical_ode!, u0, tspan,
                    (; 
                        n_dofs=sim_opts.total_dofs,
                        bc_dofs=sim_opts.bc_dofs,
                        T_func = temp_func,
                        M_interp = M_T,
                        K_interp = K_T,
                        α = α,
                        β = β,
                        load_vector = load_vector,
                    ))

    @info "Solving dynamic response with $(sim_opts.total_dofs) DOFs in physical space..."
    @time u_ = solve(prob_physical, saveat=saveat);
    u = u_[1:sim_opts.total_dofs, :]
    du = u_[sim_opts.total_dofs+1:end, :]

    return u, du, u_.t, temp_func.(u_.t)

end

"""
    solve_modal_parareal(sim_opts::SimulationOptions, u0::Vector{Float64}, 
                        tspan::Tuple{Float64, Float64}, temp_func::Function, 
                        load_vector::Function; 
                        dt_fine=0.001, dt_coarse=0.01, n_iterations=5, n_processors=4)

Parareal time-parallel integration for modal space dynamics.

# Arguments
- `sim_opts`: Simulation options
- `u0`: Initial state vector [q; qdot]
- `tspan`: Time span (t0, tf)
- `temp_func`: Temperature function of time
- `load_vector`: External load vector function
- `dt_fine`: Fine solver timestep (default: 0.001)
- `dt_coarse`: Coarse solver timestep (default: 0.01)
- `n_iterations`: Number of Parareal iterations (default: 5)
- `n_processors`: Number of time slices for parallelization (default: 4)

# Returns
- `u`: Physical displacement history
- `du`: Physical velocity history
- `t`: Time vector
- `T`: Temperature history
"""
function solve_modal_parareal(sim_opts::SimulationOptions, u0::Vector{Float64}, 
                              tspan::Tuple{Float64, Float64}, temp_func::Function, 
                              load_vector::Function; 
                              dt_fine=0.001, dt_coarse=0.01, n_iterations=5, n_processors=4,
                              coarse_alg=Rodas4(autodiff=AutoFiniteDiff()), fine_alg=Tsit5())
    
    λ_T, Φ_T, n_modes = setup_ROM(sim_opts)
    M_T, K_T = setup_physical(sim_opts)

    # Convert initial conditions to modal space
    u0_modal_ = Φ_T(temp_func(0))' * M_T(temp_func(0)) * u0[1:sim_opts.total_dofs]
    du0_modal_ = Φ_T(temp_func(0))' * M_T(temp_func(0)) * u0[sim_opts.total_dofs+1:end]
    u0_modal = [u0_modal_; du0_modal_]

    # Setup damping
    α = sim_opts.damping_ratio
    β = 0.1*α
    C = α * M_T(temp_func(0)) + β * K_T(temp_func(0))
    C_modal = Φ_T(temp_func(0))' * C * Φ_T(temp_func(0))
    ζ = diag(C_modal)

    # Pre-allocate cache arrays for ODE function
    params = (
        n_modes = n_modes,
        n_dofs = sim_opts.total_dofs,
        T_func = temp_func,
        λ_interp = λ_T,
        Φ_interp = Φ_T,
        ζ = ζ,
        load_vector = load_vector,
        # Cache arrays to reduce allocations
        ω_cache = zeros(n_modes),
        Φ_cache = zeros(sim_opts.total_dofs, n_modes),
        f_cache = zeros(sim_opts.total_dofs),
        fhat_cache = zeros(n_modes),
    )

    # Divide time domain into slices
    t_total = tspan[2] - tspan[1]
    Δt_slice = t_total / n_processors
    t_boundaries = range(tspan[1], tspan[2], length=n_processors+1)
    
    @info "Parareal setup: $(n_processors) time slices, $(n_iterations) iterations"
    @info "Coarse solver: $(coarse_alg), Fine solver: $(fine_alg)"

    # Initialize solution at time slice boundaries
    U = Vector{Vector{Float64}}(undef, n_processors+1)
    U[1] = copy(u0_modal)
    
    # Step 1: Initial coarse sequential solve to get boundary values
    @info "Initial coarse sequential solve..."
    @time for k in 1:n_processors
        t_span_k = (t_boundaries[k], t_boundaries[k+1])
        prob_coarse = ODEProblem(beam_modal_ode!, U[k], t_span_k, params)
        sol_coarse = solve(prob_coarse, coarse_alg, dt=dt_coarse, 
                          adaptive=true, save_everystep=false)
        U[k+1] = sol_coarse.u[end]
    end

    # Store previous iteration values
    U_old = deepcopy(U)
    
    # Parareal iterations
    for iter in 1:n_iterations
        @info "Parareal iteration $(iter)/$(n_iterations)..."
        
        # Parallel fine solve on each time slice
        U_fine = Vector{Vector{Float64}}(undef, n_processors+1)
        U_fine[1] = U[1]
        
        @time begin
            fine_solutions = Vector{Any}(undef, n_processors)
            
            Threads.@threads for k in 1:n_processors
                t_span_k = (t_boundaries[k], t_boundaries[k+1])
                prob_fine = ODEProblem(beam_modal_ode!, U_old[k], t_span_k, params)
                fine_solutions[k] = solve(prob_fine, fine_alg, dt=dt_fine,
                                         adaptive=true, save_everystep=false)
            end
            
            # Extract final values
            for k in 1:n_processors
                U_fine[k+1] = fine_solutions[k].u[end]
            end
        end
        
        # Sequential coarse correction
        for k in 1:n_processors
            t_span_k = (t_boundaries[k], t_boundaries[k+1])
            
            # Coarse solve with new initial condition
            prob_coarse_new = ODEProblem(beam_modal_ode!, U[k], t_span_k, params)
            sol_coarse_new = solve(prob_coarse_new, coarse_alg, dt=dt_coarse,
                                   adaptive=true, save_everystep=false)
            
            # Coarse solve with old initial condition
            prob_coarse_old = ODEProblem(beam_modal_ode!, U_old[k], t_span_k, params)
            sol_coarse_old = solve(prob_coarse_old, coarse_alg, dt=dt_coarse,
                                   adaptive=true, save_everystep=false)
            
            # Parareal update: U_new = G(U_new) + F(U_old) - G(U_old)
            U[k+1] = sol_coarse_new.u[end] .+ U_fine[k+1] .- sol_coarse_old.u[end]
        end
        
        # Check convergence
        max_diff = maximum(norm(U[k] .- U_old[k]) for k in 1:n_processors+1)
        @info "  Max difference from previous iteration: $(max_diff)"
        
        if max_diff < 1e-6
            @info "  Converged after $(iter) iterations"
            break
        end
        
        U_old = deepcopy(U)
    end
    
    # Fine solve on all slices with converged initial conditions to get full history
    @info "Final fine solve for full time history..."
    all_solutions = Vector{Any}(undef, n_processors)
    
    @time Threads.@threads for k in 1:n_processors
        t_span_k = (t_boundaries[k], t_boundaries[k+1])
        prob_final = ODEProblem(beam_modal_ode!, U[k], t_span_k, params)
        all_solutions[k] = solve(prob_final, fine_alg, saveat=dt_fine, adaptive=true)
    end
    
    # Concatenate solutions properly
    t_all = Float64[]
    u_all_list = Vector{Float64}[]
    
    for k in 1:n_processors
        sol = all_solutions[k]
        if k < n_processors
            # Exclude last point to avoid duplication
            append!(t_all, sol.t[1:end-1])
            append!(u_all_list, sol.u[1:end-1])
        else
            # Include all points for last slice
            append!(t_all, sol.t)
            append!(u_all_list, sol.u)
        end
    end
    
    # Convert to matrix
    u_all = reduce(hcat, u_all_list)
    
    # Reconstruct physical space solution
    u, du = reconstruct_physical(sim_opts, u_all, Φ_T, temp_func, t_all)
    
    return u, du, t_all, temp_func.(t_all)
end

# End of Dynamics/simulation.jl

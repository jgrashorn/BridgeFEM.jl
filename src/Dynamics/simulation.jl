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
    q     = u[1:p.n_modes]        # modal displacements
    qdot  = u[p.n_modes+1:end]    # modal velocities
    # Interpolate natural frequencies at current T (convert Hz to rad/s)
    ω = 2π .* p.λ_interp(T)
    # Damping ratios (could be constant, interpolated, or Rayleigh-like)
    ζ = p.ζ                      # vector of damping ratios per mode
    # Interpolate mode shapes at T
    Φ = p.Φ_interp(T)
    # Assemble global load vector
    f = p.load_vector(t,1:p.n_dofs)
    # Project force onto each mode
    fhat = Φ' * f
    # Modal accelerations
    qddot = -ζ .* qdot .- (ω .^2) .* q .+ fhat
    # Fill derivative vector
    du[1:p.n_modes] .= qdot
    du[p.n_modes+1:end] .= qddot
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

function solve_modal_simulation(sim_opts::SimulationOptions, u0::Vector{Float64}, tspan::Tuple{Float64, Float64}, temp_func::Function, load_vector::Function; saveat=0.01)

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
    @time u_ = solve(prob_modal, saveat=saveat);
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

# End of Dynamics/simulation.jl

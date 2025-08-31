using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

# Use proper BridgeFEM module import
using BridgeFEM

# Beam and material parameters
L = 10.0               # Beam length (m)
n_elem = 12           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = .1              # Cross-section area (m^2)
I = .001              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
α = -1e5              # E-temperature slope (Pa/K)
cutoff_freq = 500.0  # Cutoff frequency for modes (Hz)

# analytical solution
force_dof = (n_node ÷ 2)*3 + 2
f0 = -10000.0

x = range(0, L, length=n_node)
u_analytical = zeros(n_node)
for (i,x_) in enumerate(x)
    if x_ <= L/2
        u_analytical[i] = f0 .* x_ ./ (48 .*E0.*I) .* (3 .*L.^2 .- 4.0 .* x_.^2)
    else
        u_analytical[i] = f0 .* (L-x_) ./ (48 .*E0.*I) .* (3 .*L.^2 .- 4.0 .* (L-x_).^2)
    end
end

bc = BridgeBC([
    [1, "trans"],
    [n_node, "y"]
])

bo = BridgeOptions(n_elem, bc, L, ρ, A, I, [-100 E0; 70 E0], cutoff_freq)

# Support element parameters with temperature dependence
E_bridge = [
    20.0  E0
]

supports = SupportElement[]

Ts = [10.0, 20.0]

# Create comprehensive simulation options
sim_opts = SimulationOptions(
    bo, supports, collect(Ts), damping_ratio=0.02
)

M_T, K_T = setup_physical(sim_opts)
λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

temp = t -> 20.0  # Constant temperature for this test

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

α, β = 0.001, 0.001  # Rayleigh damping coefficients

C = α * M_T(temp(20.0)) + β * K_T(temp(20.0))  # Rayleigh damping matrix
C_modal = Φ_T(20.0)' * C * Φ_T(20.0)  # Project to modal space
ζ = diag(C_modal)

tspan = (0.0, 10.0)  # Simulation time span

prob_modal = ODEProblem(beam_modal_ode!, u0_modal, tspan,
                    (; 
                        n_modes=n_modes,
                        n_dofs=sim_opts.total_dofs,
                        T_func = temp,
                        λ_interp = λ_T,
                        Φ_interp = Φ_T,
                        ζ = ζ,  # Constant damping
                        load_vector = load_vector,
                    ))

prob_physical = ODEProblem(beam_physical_ode!, u0, tspan,
                    (; 
                        n_dofs=sim_opts.total_dofs,
                        bc_dofs=sim_opts.bc_dofs,
                        T_func = temp,
                        M_interp = M_T,
                        K_interp = K_T,
                        # ζ = fill(damping_ratio, total_dofs),  # Constant damping
                        α = α,
                        β = β,
                        load_vector = load_vector,
                    ))

@info "Solving dynamic response with $(n_modes) modes in modal space..."
@time sol_modal = solve(prob_modal, saveat=0.01);
q = reduce(hcat, sol_modal.u)
u_modal, du_modal = reconstruct_physical(sim_opts, q, Φ_T, temp, sol_modal.t)

@info "Solving dynamic response with $(sim_opts.total_dofs) DOFs in physical space..."
@time sol_physical = solve(prob_physical, saveat=0.01)
u_ = reduce(hcat, sol_physical.u)
u_physical = u_[1:sim_opts.total_dofs, :]
du_physical = u_[sim_opts.total_dofs+1:end, :]

plot(sol_modal.t, u_modal[force_dof, :], label="Modal Displacement at Free End")
plot!(sol_physical.t, u_[force_dof, :], label="Physical Displacement at Free End", linestyle=:dash)

plot(x, u_analytical, label="Analytical Solution", linestyle=:dot, color=:red, linewidth=3)
plot!(x, u_modal[2:3:end, end], label="Modal Solution at end", linestyle=:dashdot, color=:green, linewidth=2)
plot!(x, u_physical[2:3:end, end], label="Physical Solution at end", linestyle=:dashdot, color=:blue)
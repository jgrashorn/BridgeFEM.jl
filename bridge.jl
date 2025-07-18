using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

include("src/bridge_model.jl")
include("src/model_reduction.jl")
include("src/utils.jl")
include("src/dynamic_simulation.jl")

# Beam and material parameters
L = 20.0               # Beam length (m)
n_elem = 12           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = 10.25              # Cross-section area (m^2)
I = 0.00271              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
α = -1e5              # E-temperature slope (Pa/K)
cutoff_freq = 100.0  # Cutoff frequency for modes (Hz)

bc = BridgeBC([  # Node 1: both translational and rotational DOFs fixed
    [1, "trans"],
    [n_node, "y"]
])

E_bridge = [
    -10.0  E0;   # E at -10°C
     20.0  E0;   # E at 20°C
     50.0  E0    # E at 50°C (thermal expansion reduces stiffness)
]

bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)

E_support = E_bridge

A_support = 0.01        # Cross-section area (m^2)
I_support = 0.001       # Moment of inertia (m^4)
L_support = 5.0       # Length of support element (m)

se = [SupportElement(
        n_elem ÷ 3,          
        [1, 2, 3],               # Connect x,y,ϕ DOFs
        -150.0,                   # angle (degrees)
        5,                       # 5 elements in support
        A_support,               # Cross-sectional area
        I_support,               # Moment of inertia
        E_support,          # Temperature-dependent Young's modulus
        L_support,               # Length
        BCTypes["trans"]           # Fix all DOFs at bottom
    ),
    SupportElement(
        n_elem - (n_elem ÷ 3),   
        [1, 2, 3],               # Connect x,y,ϕ DOFs
        -30.0,                   # angle (degrees)
        5,                       # 5 elements in support
        A_support,               # Cross-sectional area
        I_support,               # Moment of inertia
        E_support,          # Temperature-dependent Young's modulus
        L_support,               # Length
        BCTypes["trans"]           # Fix all DOFs at bottom
    )]

# Example workflow with supports:
supports = se  # Your support element
# supports = SupportElement[]

Ts = [-10.0, 20.0, 50.0]  # Temperatures in degrees Celsius

# Create comprehensive simulation options
sim_opts = SimulationOptions(
    bo, supports, collect(Ts), damping_ratio=0.02
)
save_simulation_options(sim_opts, "data/bridge_simulation_config.json")

M_T, K_T = setup_physical(sim_opts)
λ_T, Φ_T, n_modes = setup_ROM(sim_opts)

force_node = (3*(n_node ÷ 2)) + 2

tspan = (0.0, 10.0)

load_vector = (t, dof) -> begin 
    f = zeros(Float64, length(dof))
    f[force_node] = -1000.0  # Apply force at force_node
    return f
end

α, β = 0.01, 0.01  # Rayleigh damping coefficients

C = α * M_T(20.0) + β * K_T(20.0)  # Rayleigh damping matrix
C_modal = Φ_T(20.0)' * C * Φ_T(20.0)  # Project to modal space
ζ = diag(C_modal)

Φ = Φ_T(20.0)  # Mode shapes at T=20°C

u0_physical = zeros(sim_opts.total_dofs * 2)

u0_modal_ = Φ_T(20.0)' * M_T(20.0) * u0_physical[1:sim_opts.total_dofs]
du0_modal_ = Φ_T(20.0)' * u0_physical[sim_opts.total_dofs+1:end]
u0_modal = [u0_modal_; du0_modal_]

T_func = (t) -> 20.0

prob_modal = ODEProblem(beam_modal_ode!, u0_modal, tspan,
                    (; 
                        n_modes=n_modes,
                        n_dofs=sim_opts.total_dofs,
                        T_func = T_func,
                        λ_interp = λ_T,
                        ζ = ζ,  # Constant damping
                        Φ_interp = Φ_T,
                        load_vector = load_vector,
                    ))

prob_physical = ODEProblem(beam_physical_ode!, u0_physical, tspan,
                    (; 
                        n_dofs=sim_opts.total_dofs,
                        bc_dofs=sim_opts.bc_dofs,
                        T_func = T_func,
                        M_interp = M_T,
                        K_interp = K_T,
                        α = α,
                        β = β,
                        load_vector = load_vector,
                    ))

@info "Solving dynamic response with $(n_modes) modes in modal space..."
@time sol_modal = solve(prob_modal, saveat=0.01);
q = reduce(hcat, sol_modal.u)
u_modal, du_modal = reconstruct_physical(sim_opts, q, Φ_T, T_func, sol_modal.t)

@info "Solving dynamic response with $(sim_opts.total_dofs) DOFs in physical space..."
@time sol_physical = solve(prob_physical, saveat=0.01);
u_ = reduce(hcat, sol_physical.u)
u_physical = u_[1:sim_opts.total_dofs, :]
du_physical = u_[sim_opts.total_dofs+1:end, :]

p1 = plot(size=(800,600), dpi=300)
plot!(p1,u_physical[5:6:end,:]',label="physical")
plot!(p1,u_modal[5:6:end,:]',linestyle=:dash,label="modal")
display(p1)
savefig(p1, "bridge_dynamic_response.png")
# # 1. Plot structure only
# plot_bridge_with_supports(bo, supports)

time_subsample = sol_modal.t
u_subsample = u_modal
anim_fast = animate_dynamic_response(bo, supports, u_subsample, time_subsample,
                                   scale_factor=1000.0,
                                   n_frames=400,
                                   fps=24,
                                   filename="bridge_dynamics.gif")

# support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)

# 2. Plot specific mode shape
# mode_num = 3
# mode_shape = Φ_T(50.0)[:,mode_num]  # First mode at first temperature

# plot_mode_shape(bo, supports, mode_shape, mode_num, scale_factor=10000.0)

# # 3. Animate a mode
# mode_num = 1
# mode_shape = vectors_unnormalized[:, mode_num, 1]  # Second mode
# anim = animate_mode(bo, supports, mode_shape, mode_num, scale_factor=10.0, fsize=(800, 600))
# gif(anim, "mode2_animation.gif", fps=15)
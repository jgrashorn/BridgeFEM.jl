using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

include("../src/bridge_model.jl")
include("../src/model_reduction.jl")
include("../src/utils.jl")
include("../src/dynamic_simulation.jl")

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

bc = BridgeBC([  # Node 1: both translational and rotational DOFs fixed
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

# M, K = assemble_matrices_with_supports(sim_opts)
M_T, K_T = setup_physical(sim_opts)

f0 = -10000.0
f = zeros(n_node*3)
f[(n_node ÷ 2)*3 + 2] = f0

M_, K_ = apply_bc(M_T(20.0), K_T(20.0), sim_opts)

u = K_ \ f

x = range(0, L, length=n_node)

u_analytical = zeros(n_node)
for (i,x_) in enumerate(x)
    if x_ <= L/2
        u_analytical[i] = f0 .* x_ ./ (48 .*E0.*I) .* (3 .*L.^2 .- 4.0 .* x_.^2)
    else
        u_analytical[i] = f0 .* (L-x_) ./ (48 .*E0.*I) .* (3 .*L.^2 .- 4.0 .* (L-x_).^2)
    end
end

plot(x, u[2:3:end], label="y")
plot!(x, u_analytical, label="y analytical", linestyle=:dash, xlabel="x (m)", ylabel="y (m)")

using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using JSON

include("../src/bridge_model.jl")
include("../src/model_reduction.jl")
include("../src/utils.jl")
include("../src/dynamic_simulation.jl")

# Beam and material parameters
L = 10.0               # Beam length (m)
n_elem = 11           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = .1              # Cross-section area (m^2)
I = .0001              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
α = -1e5              # E-temperature slope (Pa/K)
cutoff_freq = 100.0  # Cutoff frequency for modes (Hz)

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

M, K = assemble_matrices_with_supports(sim_opts)
M_, K_ = apply_bc(M, K, sim_opts)

_, _, λs, vectors, vectors_unnormalized = assemble_and_decompose(sim_opts)

# force at free end
f0 = -10000.0
f = zeros(n_node*3)
f[end-1] = f0

# modal force
f_m = vectors[:, :, 1]' * f

# solve for displacement
u = K_[:,:,1] \ f

# solve for modal displacement
q = f_m ./ (2π .* λs[:,1]).^2
u_m = vectors[:,:,1] * q

x = range(0, L, length=n_node)

u_analytical = zeros(n_node*3)
u_analytical = f0 .* x.^2 ./ (6 .*E0.*I) .* (3 .*L .- x)
ϕ_analytical = f0 .* x ./ (2 .*E0.*I) .* (2 .*L .- x)

isapprox(u_analytical, u[2:3:end], atol=1e-4)
isapprox(ϕ_analytical, u[3:3:end], atol=1e-4)
isapprox(u_analytical, u_m[2:3:end], atol=1e-4)
isapprox(ϕ_analytical, u_m[3:3:end], atol=1e-4)
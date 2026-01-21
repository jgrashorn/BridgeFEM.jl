Example calculation:

```julia
using Pkg
Pkg.add(url="https://github.com/jgrashorn/BridgeFEM")
using LinearAlgebra, Plots
using Interpolations
using JSON

# Use proper BridgeFEM module import
using BridgeFEM

# Beam and material parameters
L = 20.0               # Beam length (m)
n_elem = 100           # Number of finite elements
n_node = n_elem + 1   # Number of nodes
ρ = 7800.0            # Density (kg/m^3)
A = 10.25              # Cross-section area (m^2)
I = 0.00271              # Moment of inertia (m^4)
E0 = 207e9             # Base Young's modulus (Pa)
cutoff_freq = 100.0  # Cutoff frequency for modes (Hz)

bc = BridgeBC([  # Node 1: both tranlational DOF fixed
    [1, "trans"],
    [n_node, "y"]
])

E_bridge = [
    -10.0  E0;   # E at -10°C
     20.0  E0*0.9;   # E at 20°C
     50.0  E0*0.8    # E at 50°C (thermal expansion reduces stiffness)
]

bo = BridgeOptions(n_elem, bc, L, ρ, A, I, E_bridge, cutoff_freq)

E_support = E_bridge

A_support = 0.01        # Cross-section area (m^2)
I_support = 0.001       # Moment of inertia (m^4)
L_support = 5.0       # Length of support element (m)

supports = SupportElement[] # currently not supported bc of bug in boundary conditions

Ts = [-10.0, 20.0, 50.0]  # Temperatures in degrees Celsius

# Create comprehensive simulation options
sim_opts = SimulationOptions(
    bo, supports, collect(Ts), damping_ratio=0.02
)

save_simulation_options(sim_opts, "bridge_simulation_config.json")

force_node = (3*(n_node ÷ 2)) + 2

chunk_duration = 60.0

load_vector_f = (t, dof) -> begin 
    f = zeros(Float64, length(dof))
    f[force_node] = -10000.0 * sin(t)  # Apply force at force_node
    return f
end

load_vector = (t_start) -> begin
    return (t, dof) -> load_vector_f(t + t_start, dof)
end

u0 = zeros(sim_opts.total_dofs * 2)

nChunk = 3

u_ = zeros(sim_opts.total_dofs*2, 0)
f_ = zeros(sim_opts.total_dofs, 0)
T_ = zeros(0)

T_func = (t) -> 20.0 * exp(-0.001*t)  # Constant temperature function

for i in 1:nChunk

    @info "Starting chunk $i of $nChunk"

    t_curr = (i-1)*chunk_duration

    tspan = (t_curr, t_curr + chunk_duration)

    u_modal = solve_modal_simulation(sim_opts, u0, tspan, T_func, load_vector_f)
    for t_ in tspan[1]:0.1:tspan[2]
        f_ = hcat(f_, load_vector_f(t_, 1:sim_opts.total_dofs))
        T_ = vcat(T_, T_func(t_))
    end
    u0 = [u_modal[1][:,end]; u_modal[2][:,end]]
    u_ = hcat(u_, [u_modal[1]; u_modal[2]])
end

plot(u_[5:6,:]')
```

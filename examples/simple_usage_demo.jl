# Simple Usage Demonstration for BridgeFEM.jl
# Shows the new modular import patterns after package reorganization

# Import the BridgeFEM package (new pattern)
using BridgeFEM

# Import external dependencies as needed
using LinearAlgebra

println("🌉 BridgeFEM.jl - Simple Usage Demonstration")
println("=" ^ 50)

# 1. Basic Bridge Configuration
println("\n1. Creating Bridge Configuration...")

# Temperature-dependent Young's modulus [Temperature(°C), E(Pa)]
E_T = [0.0 200e9; 50.0 180e9; 100.0 160e9]

# Boundary conditions - fix first node completely
bc = BridgeBC([[1, "all"]])

# Bridge parameters
bridge = BridgeOptions(
    15,          # n_elem: Number of elements
    bc,          # boundary conditions
    20.0,        # L: Length (m)  
    2500.0,      # ρ: Density (kg/m³)
    0.3,         # A: Cross-sectional area (m²)
    0.01,        # I: Moment of inertia (m⁴)
    E_T,         # Temperature-dependent Young's modulus
    100.0        # max_frequency: Maximum frequency for modal analysis (Hz)
)

println("✅ Bridge created: $(bridge.n_elem) elements, length $(bridge.L)m")

# 2. Support Element Configuration
println("\n2. Creating Support Elements...")

# Create a simple support at mid-span
mid_node = div(bridge.n_elem, 2) + 1  # Connect to middle node
support = SupportElement(
    mid_node,    # connection_node: Node to connect to
    [1, 2],      # connection_dofs: Connect u,v DOFs
    90.0,        # angle: Vertical support (90 degrees)
    3,           # n_elem: Elements in support
    0.1,         # A: Support cross-sectional area (m²)  
    0.001,       # I: Support moment of inertia (m⁴)
    E_T,         # Temperature-dependent Young's modulus
    5.0,         # L: Support length (m)
    [1, 2, 3]    # bc_bottom: Fix bottom of support completely
)

supports = [support]
println("✅ Support created at node $(mid_node) with $(support.n_elem) elements")

# 3. Temperature Configuration
println("\n3. Setting up Temperature Analysis...")

temperatures = [20.0, 40.0, 60.0]  # Analysis temperatures (°C)
println("✅ Temperature range: $(temperatures[1])°C to $(temperatures[end])°C")

# 4. Basic Matrix Assembly
println("\n4. Demonstrating Basic Matrix Assembly...")

# Assemble matrices at reference temperature
M, K = assemble_matrices(bridge, temperatures[1])
println("✅ Global matrices assembled:")
println("   - Mass matrix: $(size(M,1))×$(size(M,2))")
println("   - Stiffness matrix: $(size(K,1))×$(size(K,2))")

# 5. Boundary Condition Application
println("\n5. Applying Boundary Conditions...")

# Create simulation options for boundary condition demonstration
sim_opts = SimulationOptions(bridge, supports, temperatures)

# Apply boundary conditions (this would normally be done internally)
println("✅ Boundary conditions configured for $(length(sim_opts.bc_dofs)) constrained DOFs")

# 6. JSON Configuration Persistence
println("\n6. Demonstrating Configuration Persistence...")

# Save configuration to JSON
config_file = "simple_demo_config.json"
save_simulation_options(sim_opts, config_file)
println("✅ Configuration saved to: $config_file")

# Load configuration back
loaded_opts = load_simulation_options(config_file)
println("✅ Configuration loaded successfully")

# Verify round-trip compatibility
if loaded_opts.bridge.L ≈ sim_opts.bridge.L
    println("✅ Round-trip verification: PASSED")
else
    println("❌ Round-trip verification: FAILED")
end

# 7. Summary
println("\n" * "=" ^ 50)
println("🎉 BridgeFEM.jl Module Usage Demonstration Complete!")
println("\nKey Features Demonstrated:")
println("  ✓ Modular package import with 'using BridgeFEM'")
println("  ✓ Bridge and support element configuration")
println("  ✓ Temperature-dependent material modeling")
println("  ✓ Matrix assembly operations")
println("  ✓ Boundary condition setup")
println("  ✓ JSON configuration persistence")
println("\nThe package is now properly organized with:")
println("  • Core types and constants")
println("  • Element computation functions")
println("  • Global matrix assembly")
println("  • Boundary condition application")
println("  • I/O serialization")
println("  • Dynamic simulation capabilities")
println("  • Modal analysis functions")
println("  • Visualization utilities")

# Clean up demo file
rm(config_file, force=true)
println("\n🧹 Demo cleanup completed.")

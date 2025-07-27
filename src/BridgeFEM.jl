"""
# BridgeFEM.jl

A Julia package for finite element analysis of bridge structures with temperature-dependent 
material properties, dynamic analysis capabilities, and modular architecture.

## Features
- Temperature-dependent material modeling
- Static and dynamic structural analysis  
- Support element modeling
- JSON configuration persistence
- Modular architecture for extensibility

## Modules
- **Core**: Fundamental data structures and constants
- **Elements**: Element stiffness and mass matrices (planned)
- **Assembly**: Global matrix assembly (planned)
- **BoundaryConditions**: Constraint application (planned)
- **Dynamics**: Dynamic simulation (planned)
- **IO**: Configuration and results persistence (planned)

## Quick Start
```julia
using BridgeFEM

# Create bridge configuration
E_T = [0.0 200e9; 100.0 160e9]  # Temperature-dependent Young's modulus
bc = BridgeBC([[1, "all"]])     # Fixed boundary condition
bridge = BridgeOptions(10, bc, 30.0, 2500.0, 0.1, 0.001, E_T, 50.0)

# Save/load configuration
save_simulation_options(SimulationOptions(bridge, [20.0, 30.0]), "config.json")
sim_opts = load_simulation_options("config.json")
```
"""
module BridgeFEM

# Core module exports - fundamental types and functionality
include("Core/constants.jl")
include("Core/types.jl")

# Export Core types and constants
export BCTypes

# Export Core data structures
export BridgeBC, BridgeOptions, SupportElement, SimulationOptions

# Export JSON serialization functions
export bridge_options_to_dict, dict_to_bridge_options
export support_element_to_dict, dict_to_support_element  
export simulation_options_to_dict, load_simulation_options, save_simulation_options

# Future module exports (will be uncommented as modules are implemented)
# Elements module
# export frame_elem_stiffness, frame_elem_mass

# Assembly module  
# export assemble_matrices, create_support_dof_mapping

# BoundaryConditions module
# export apply_bc, get_bc_dofs

# IO module
# export save_results, load_results

# Dynamics module
# export simulate_dynamic_response

# ModelReduction module
# export compute_modal_properties

# Visualization module  
# export plot_mode_shape, plot_response

end # module BridgeFEM 
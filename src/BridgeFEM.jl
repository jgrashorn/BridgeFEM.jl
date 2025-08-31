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
- **Elements**: Element stiffness and mass matrices
- **Assembly**: Global matrix assembly and DOF management
- **BoundaryConditions**: Constraint application and DOF management
- **IO**: Configuration and results persistence
- **Dynamics**: Dynamic simulation and ODE solving
- **ModelReduction**: Modal analysis and eigenvalue computation

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

# Elements module - finite element computations
include("Elements/finite_elements.jl")

# Assembly modules - global matrix assembly and DOF management
include("Assembly/matrices.jl")
include("Assembly/dof_mapping.jl")

# BoundaryConditions module - boundary condition application
include("BoundaryConditions/application.jl")

# IO module - JSON serialization and configuration persistence
include("IO/serialization.jl")

# Dynamics module - dynamic simulation and ODE solving
include("Dynamics/simulation.jl")

# ModelReduction module - modal analysis and eigenvalue computation
include("ModelReduction/modal.jl")

# Visualization module - plotting and visualization functions
include("Visualization/plotting.jl")

# Export Core types and constants
export BCTypes

# Export Core data structures
export BridgeBC, BridgeOptions, SupportElement, SimulationOptions



# Elements module exports - finite element computations
export frame_elem_stiffness, frame_elem_mass, transformation_matrix

# Assembly module exports - global matrix assembly and DOF management
export assemble_matrices, assemble_stiffness!, assemble_matrices_with_supports
export create_support_dof_mapping, get_dof_from_node, get_bc_dofs
export assemble_local_support, create_support_mass_matrix, create_expanded_transformation
export interpolate_matrix, setup_matrix_interpolation, setup_physical

# BoundaryConditions module exports - boundary condition application
export apply_bc, remove_fixed_dofs

# IO module exports - JSON serialization and configuration persistence
export bridge_options_to_dict, dict_to_bridge_options
export support_element_to_dict, dict_to_support_element  
export simulation_options_to_dict, load_simulation_options, save_simulation_options

# Dynamics module exports - dynamic simulation and ODE solving
export beam_modal_ode!, beam_physical_ode!

# ModelReduction module exports - modal analysis and eigenvalue computation
export decompose_matrices, assemble_and_decompose, setup_interpolation, 
       interpolate_modes, reconstruct_physical, setup_ROM

# Visualization module exports - plotting and visualization functions
export plot_bridge_with_supports, plot_mode_shape, animate_dynamic_response, 
       animate_modal_response, animate_dof_response

end # module BridgeFEM 
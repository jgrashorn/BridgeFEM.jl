"""
# Core Type Definitions for BridgeFEM.jl

This module provides fundamental data structures and constants for finite element 
bridge analysis, including boundary conditions, material properties, support 
elements, and simulation configuration.

## Exported Types
- `BridgeBC`: Boundary condition specifications
- `BridgeOptions`: Analysis configuration with material and geometric properties  
- `SupportElement`: Support structure modeling
- `SimulationOptions`: Dynamic analysis control parameters

## JSON Serialization
All types support JSON serialization for configuration persistence with exact
format compatibility for existing research workflows.

## Temperature Dependencies
Both `BridgeOptions` and `SupportElement` support temperature-dependent Young's
modulus through Interpolations.jl integration.
"""

# Core type definitions for BridgeFEM analysis

using JSON
using Interpolations
using Dates

# Import constants from the same module
include("constants.jl")

"""
    BridgeBC(conds::Vector{Vector{Any}})

Boundary condition specification for finite element analysis.

# Fields
- `conds::Vector{Vector{Any}}`: Constraint specifications where each entry contains
  [node_id, constraint_type]. Constraint types can be strings (mapped via BCTypes)
  or direct DOF vectors.

# Examples
```julia
# String constraint types (recommended)
bc = BridgeBC([[1, "all"], [5, "trans"]])

# Direct DOF specification  
bc = BridgeBC([[1, [1,2,3]], [5, [1,2]]])
```
"""
struct BridgeBC
    conds::Vector{Vector{Any}}

    function BridgeBC(conds::Vector{Vector{Any}})
        # If the second entry is a String, map it; otherwise, use as is
        c_ = [
            (isa(c[2], String) ? [c[1], BCTypes[c[2]]] : c)
            for c in conds
        ]
        return new(c_)
    end
end

"""
    BridgeOptions

Mutable configuration structure for bridge analysis containing material properties,
geometric specifications, and boundary conditions.

# Fields
- `n_elem::Int`: Number of finite elements
- `n_nodes::Int`: Number of nodes (automatically computed as n_elem + 1)
- `n_dofs::Int`: Total degrees of freedom (3 per node: u, v, θ)
- `bc_nodes::BridgeBC`: Boundary condition specifications
- `L::Float64`: Beam length (m)
- `ρ::Float64`: Material density (kg/m³)
- `A::Float64`: Cross-sectional area (m²)
- `I::Float64`: Moment of inertia (m⁴)
- `E_T::Matrix{Float64}`: Temperature-dependent Young's modulus data [T, E]
- `E::Function`: Young's modulus interpolation function E(T)
- `cutoff_freq::Float64`: Modal analysis cutoff frequency (Hz)

# Constructors
```julia
# Temperature-dependent Young's modulus
E_T = [0.0 200e9; 50.0 180e9; 100.0 160e9]
bc = BridgeBC([[1, "all"]])
bridge = BridgeOptions(10, bc, 30.0, 2500.0, 0.1, 0.001, E_T, 50.0)
```
"""
mutable struct BridgeOptions
    n_elem::Int # Number of finite elements
    n_nodes::Int # Number of nodes
    n_dofs::Int # Number of degrees of freedom (3 per node)
    bc_nodes::BridgeBC # Boundary condition nodes (default: empty)
    L::Float64 # Beam length (m)
    ρ::Float64 # Density (kg/m^3)
    A::Float64 # Cross-section area (m^2)
    I::Float64 # Moment of inertia (m^4)
    E_T::Matrix{Float64} # Young's modulus at a specific temperature (Pa)
    E::Function # Young's modulus as a function of temperature (Pa)
    cutoff_freq::Float64 # Cutoff frequency for modes (Hz)

    function BridgeOptions(n_elem::Int, bc_nodes::BridgeBC, L::Float64, ρ::Float64, A::Float64, I::Float64, E_T::Matrix{Float64}, cutoff_freq::Float64)
        E_interp = interpolate((E_T[:,1],), E_T[:,2], Gridded(Linear()))
        E = T -> E_interp(T)
        return BridgeOptions(n_elem, bc_nodes, L, ρ, A, I, E_T, E, cutoff_freq)
    end

    function BridgeOptions(n_elem::Int, bc_nodes::BridgeBC, L::Float64, ρ::Float64, A::Float64, I::Float64, E_T::Matrix{Float64}, E::Function, cutoff_freq::Float64)
        n_nodes  = n_elem + 1
        n_dofs   = 3 * n_nodes  # 3 DOFs per node (u, v, theta)
        return new(n_elem, n_nodes, n_dofs, bc_nodes, L, ρ, A, I, E_T, E, cutoff_freq)
    end
end

# JSON serialization functions moved to IO module (src/IO/serialization.jl)
# - bridge_options_to_dict(bo::BridgeOptions) 
# - dict_to_bridge_options(dict::Dict)
# - support_element_to_dict(se::SupportElement)
# - dict_to_support_element(dict::Dict)
# - simulation_options_to_dict(opts::SimulationOptions)
# - save_simulation_options(opts, filename)
# - load_simulation_options(filename)
# These functions are now available from the main BridgeFEM module

"""
    SupportElement

Mutable structure for modeling support elements connected to the main bridge structure.

# Fields
- `connection_node::Int`: Node on main bridge this support connects to
- `connection_dofs::Vector{Int}`: Specific DOFs to connect (e.g., [1,2] for x,y only)
- `angle::Float64`: Local to global rotation angle (degrees)
- `n_elem::Int`: Number of elements in the support
- `A::Float64`: Cross-sectional area (m²)
- `I::Float64`: Moment of inertia (m⁴)
- `E_T::Matrix{Float64}`: Temperature-dependent Young's modulus data [T, E]
- `E::Function`: Young's modulus interpolation function E(T)
- `L::Float64`: Support length (m)
- `bc_bottom::Vector{Int}`: Bottom boundary condition DOF types

# Constructors
```julia
# Temperature-dependent Young's modulus
E_T = [0.0 200e9; 100.0 160e9]
support = SupportElement(5, [1,2], 45.0, 3, 0.05, 0.0005, E_T, 10.0, [1,2,3])

# Constant Young's modulus
support = SupportElement(3, [1], 0.0, 2, 150e9, 0.02, 0.0001, 5.0, [1,2])
```
"""
mutable struct SupportElement
    connection_node::Int        # Node on main bridge this support connects to
    connection_dofs::Vector{Int} # Specific DOFs to connect (e.g., [1,2] for x,y only)
    angle::Float64             # Local to global rotation angle in degrees
    n_elem::Int               # Number of elements in the support
    A::Float64                # Cross-sectional area (constant)
    I::Float64                # Moment of inertia (constant)
    E_T::Matrix{Float64}      # Young's modulus vs temperature data
    E::Function               # Young's modulus as function of temperature
    L::Float64
    bc_bottom::Vector{Int}    # Boundary conditions at bottom of support (DOF types)
end

"""
    SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, L, bc_bottom)

Constructor with temperature-dependent Young's modulus.

Uses Interpolations.jl to create temperature-dependent Young's modulus function from discrete data points.
"""
function SupportElement(connection_node::Int, connection_dofs::Vector{Int}, angle::Float64, 
                       n_elem::Int, A::Float64, I::Float64, E_T::Matrix{Float64}, L::Float64, 
                       bc_bottom::Vector{Int})
    E_interp = interpolate((E_T[:,1],), E_T[:,2], Gridded(Linear()))
    E = T -> E_interp(T)
    return SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, E, L, bc_bottom)
end

"""
    SupportElement(connection_node, connection_dofs, angle, n_elem, E_const, A, I, L, bc_bottom)

Constructor with constant Young's modulus.

Creates temperature-independent support by generating dummy temperature data with constant E value.
"""
function SupportElement(connection_node::Int, connection_dofs::Vector{Int}, angle::Float64, 
                       n_elem::Int, E_const::Float64, A::Float64, I::Float64, L::Float64, 
                       bc_bottom::Vector{Int})
    # Create dummy temperature data for constant E
    E_T = [-100.0 E_const; 100.0 E_const]

    return SupportElement(connection_node, connection_dofs, angle, n_elem, A, I, E_T, L, bc_bottom)
end

"""
    SimulationOptions

Configuration structure for dynamic analysis containing bridge configuration,
support elements, environmental conditions, and system properties.

# Fields
- `bridge::BridgeOptions`: Main bridge configuration
- `supports::Vector{SupportElement}`: Support element configurations
- `temperatures::Vector{Float64}`: Temperature profile for analysis
- `damping_ratio::Float64`: System damping ratio
- `total_dofs::Int`: Total system degrees of freedom
- `total_elements::Int`: Total number of elements in system
- `support_dof_mapping::Vector{Vector{Int}}`: Maps support DOFs to global DOFs
- `bc_dofs::Vector{Int}`: All constrained DOFs in the system
- `created_at::String`: Timestamp of creation

# Note
Current implementation uses simplified DOF mapping computation to avoid dependencies
on functions that will be moved to other modules in future stories.
"""
struct SimulationOptions
    bridge::BridgeOptions
    supports::Vector{SupportElement}
    temperatures::Vector{Float64}
    damping_ratio::Float64
    total_dofs::Int
    total_elements::Int
    support_dof_mapping::Vector{Vector{Int}} # Maps support DOFs to global DOFs
    bc_dofs::Vector{Int} # All constrained DOFs in the system
    created_at::String
end

# Simplified constructor for SimulationOptions without immediate DOF mapping computation
# Note: This avoids dependency on create_support_dof_mapping and get_bc_dofs functions
# which will be moved to other modules in future stories
function SimulationOptions(bridge::BridgeOptions, supports::Vector{SupportElement}, temperatures::Vector{Float64}; 
                          damping_ratio=0.02)
    # Temporary implementation - compute basic properties
    total_elements = bridge.n_elem + (isempty(supports) ? 0 : sum(s.n_elem for s in supports))
    # For now, assume simple DOF mapping - will be replaced when Assembly/BC modules are created
    total_dofs = bridge.n_dofs + (isempty(supports) ? 0 : sum(s.n_elem + 1 for s in supports) * 3)
    support_dof_mapping = Vector{Vector{Int}}()  # Empty for now
    bc_dofs = Vector{Int}()  # Empty for now
    
    return SimulationOptions(
        bridge, supports, temperatures, damping_ratio, total_dofs, total_elements, 
        support_dof_mapping, bc_dofs, string(now())
    )
end

function SimulationOptions(bridge::BridgeOptions, temperatures::Vector{Float64}; damping_ratio=0.02)
    supports = SupportElement[]  # No supports by default
    return SimulationOptions(bridge, supports, temperatures; damping_ratio=damping_ratio)
end

 
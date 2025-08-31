"""
# IO/Serialization Functions

JSON serialization functions for BridgeFEM.jl.

Contains functions for saving and loading simulation configurations to/from JSON files,
maintaining exact compatibility with existing configuration file formats.

## Functions
- `bridge_options_to_dict(bo)`: Convert BridgeOptions to dictionary for JSON serialization
- `dict_to_bridge_options(dict)`: Convert dictionary to BridgeOptions from JSON data
- `support_element_to_dict(se)`: Convert SupportElement to dictionary for JSON serialization
- `dict_to_support_element(dict)`: Convert dictionary to SupportElement from JSON data
- `simulation_options_to_dict(opts)`: Convert SimulationOptions to dictionary for JSON serialization
- `save_simulation_options(opts, filename)`: Save SimulationOptions to JSON file
- `load_simulation_options(filename)`: Load SimulationOptions from JSON file

## Dependencies
- JSON.jl: JSON parsing and generation
- Core module: All data structure types (BridgeOptions, SupportElement, etc.)
"""

using JSON

# Import Core types  
using ..BridgeFEM: BridgeOptions, SupportElement, SimulationOptions

"""
    bridge_options_to_dict(bo) -> Dict

Convert BridgeOptions struct to dictionary for JSON serialization.

Preserves exact field structure and data types for compatibility with existing
configuration files.

# Arguments
- `bo`: Bridge configuration options

# Returns
- `Dict`: Dictionary representation suitable for JSON serialization

# Notes
Temperature-Young's modulus data is converted to vector of vectors format
for JSON compatibility while preserving numerical precision.
"""
function bridge_options_to_dict(bo)::Dict
    return Dict(
        "n_elem" => bo.n_elem,
        "bc_nodes" => [[c[1], c[2]] for c in bo.bc_nodes.conds],
        "L" => bo.L,
        "ρ" => bo.ρ,
        "A" => bo.A,
        "I" => bo.I,
        "E_T" => [bo.E_T[i, :] for i in 1:size(bo.E_T, 1)],
        "cutoff_freq" => bo.cutoff_freq
    )
end

"""
    dict_to_bridge_options(dict::Dict) -> BridgeOptions

Convert dictionary from JSON data to BridgeOptions struct.

Reconstructs BridgeOptions with proper type conversions and boundary condition
setup, maintaining exact compatibility with existing configuration files.

# Arguments
- `dict::Dict`: Dictionary from JSON deserialization

# Returns
- Bridge configuration object

# Notes
Handles temperature-Young's modulus matrix reconstruction and boundary condition
conversion while preserving all original data structure relationships.
"""
function dict_to_bridge_options(dict::Dict)
    E_T_mat = Float64.(reduce(vcat, [row' for row in dict["E_T"]]))
    bconds = BridgeBC([c for c in dict["bc_nodes"]])
    return BridgeOptions(
        dict["n_elem"],
        bconds,
        dict["L"],
        dict["ρ"],
        dict["A"],
        dict["I"],
        E_T_mat,
        dict["cutoff_freq"],
    )
end

"""
    support_element_to_dict(se) -> Dict

Convert SupportElement struct to dictionary for JSON serialization.

Preserves support element configuration including temperature-dependent material
properties and boundary condition specifications.

# Arguments
- `se`: Support element configuration

# Returns
- `Dict`: Dictionary representation suitable for JSON serialization

# Notes
Temperature-Young's modulus data and DOF constraint vectors are properly
formatted for JSON while maintaining numerical precision and type information.
"""
function support_element_to_dict(se)::Dict
    return Dict(
        "connection_node" => se.connection_node,
        "connection_dofs" => se.connection_dofs,
        "angle" => se.angle,
        "n_elem" => se.n_elem,
        "A" => se.A,
        "I" => se.I,
        "E_T" => [se.E_T[i, :] for i in 1:size(se.E_T, 1)],  # Save temperature-E data
        "L" => se.L,
        "bc_bottom" => se.bc_bottom
    )
end

"""
    dict_to_support_element(dict::Dict) -> SupportElement

Convert dictionary from JSON data to SupportElement struct.

Reconstructs SupportElement with temperature-dependent Young's modulus
interpolation and proper type conversions for all fields.

# Arguments
- `dict::Dict`: Dictionary from JSON deserialization

# Returns
- Support element configuration object

# Notes
Uses SupportElement constructor to automatically set up temperature-dependent
Young's modulus interpolation function from discrete data points.
"""
function dict_to_support_element(dict::Dict)
    E_T_mat = Float64.(reduce(vcat, [row' for row in dict["E_T"]]))
    return SupportElement(
        dict["connection_node"],
        Vector{Int}(dict["connection_dofs"]),
        dict["angle"],
        dict["n_elem"],
        dict["A"],
        dict["I"],
        E_T_mat,  # Use temperature-dependent data
        dict["L"],
        Vector{Int}(dict["bc_bottom"])
    )
end

"""
    simulation_options_to_dict(opts) -> Dict

Convert SimulationOptions struct to dictionary for complete simulation configuration JSON.

Creates comprehensive dictionary containing bridge configuration, support elements,
temperature ranges, and simulation metadata for complete system persistence.

# Arguments
- `opts`: Complete simulation configuration

# Returns
- `Dict`: Dictionary representation suitable for JSON serialization

# Notes
Recursively converts all nested structures (BridgeOptions, SupportElement arrays)
while preserving creation timestamps and computational parameters.
"""
function simulation_options_to_dict(opts)::Dict
    return Dict(
        "bridge" => bridge_options_to_dict(opts.bridge),
        "supports" => [support_element_to_dict(s) for s in opts.supports],
        "temperatures" => opts.temperatures,
        "damping_ratio" => opts.damping_ratio,
        "total_dofs" => opts.total_dofs,
        "total_elements" => opts.total_elements,
        "created_at" => opts.created_at,
    )
end

"""
    load_simulation_options(filename::String) -> SimulationOptions

Load complete simulation configuration from JSON file.

Reads and reconstructs SimulationOptions with all nested structures,
maintaining exact compatibility with existing configuration file formats.

# Arguments
- `filename::String`: Path to JSON configuration file

# Returns
- Complete simulation configuration object

# Throws
- `SystemError`: If file cannot be read
- `JSON.ParseError`: If file contains invalid JSON
- `KeyError`: If required fields are missing

# Notes
Automatically handles version compatibility and type conversions while
preserving all numerical precision and structural relationships.
"""
function load_simulation_options(filename::String)
    dict_data = JSON.parsefile(filename)
    
    # Reconstruct bridge options
    bridge = dict_to_bridge_options(dict_data["bridge"])
    
    # Reconstruct support elements
    supports = [dict_to_support_element(s) for s in dict_data["supports"]]
    
    # Extract other parameters
    temperatures = Vector{Float64}(dict_data["temperatures"])
    damping_ratio = dict_data["damping_ratio"]
    total_dofs = dict_data["total_dofs"]
    total_elements = dict_data["total_elements"]
    created_at = dict_data["created_at"]
    
    return SimulationOptions(bridge, supports, temperatures, damping_ratio, 
                           total_dofs, total_elements, Vector{Vector{Int}}(), Vector{Int}(), created_at)
end

"""
    save_simulation_options(opts::SimulationOptions, filename::String)

Save complete simulation configuration to JSON file with pretty formatting.

Creates human-readable JSON file with proper indentation and field ordering
for configuration persistence and sharing.

# Arguments
- `opts`: Complete simulation configuration to save
- `filename::String`: Output file path

# Throws
- `SystemError`: If file cannot be created or written

# Notes
Uses 2-space indentation for readable formatting while maintaining exact
numerical precision and complete structural information.
"""
function save_simulation_options(opts, filename::String)
    dict_data = simulation_options_to_dict(opts)
    open(filename, "w") do f
        JSON.print(f, dict_data, 2)  # Pretty print with 2-space indentation
    end
    @info "Simulation options saved to: $filename"
end

# Functions are now part of the main BridgeFEM module

using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using SparseArrays
using JSON
using Arpack
using Dates

# Import Core types and constants
using ..BridgeFEM: BCTypes, BridgeBC, BridgeOptions, SupportElement, SimulationOptions

# Import BoundaryConditions module functions
using ..BridgeFEM: apply_bc, remove_fixed_dofs

# Import IO module functions  
using ..BridgeFEM: bridge_options_to_dict, dict_to_bridge_options
using ..BridgeFEM: support_element_to_dict, dict_to_support_element
using ..BridgeFEM: simulation_options_to_dict, load_simulation_options, save_simulation_options

# Import Elements module functions
using ..BridgeFEM: frame_elem_stiffness, frame_elem_mass, transformation_matrix

# Import Assembly module functions
using ..BridgeFEM: assemble_matrices, assemble_stiffness!, assemble_matrices_with_supports
using ..BridgeFEM: create_support_dof_mapping, get_dof_from_node, get_bc_dofs  
using ..BridgeFEM: assemble_local_support, create_support_mass_matrix, create_expanded_transformation


# Finite element functions moved to Elements module (src/Elements/finite_elements.jl)
# - frame_elem_stiffness(EA, EI, L_e) 
# - frame_elem_mass(ρ, A, L)
# - transformation_matrix(θ)
# These functions are now imported from the main BridgeFEM module above

# Assembly functions moved to Assembly module (src/Assembly/matrices.jl)
# - assemble_stiffness!(K, bo::BridgeOptions, EA, EI)
# - assemble_matrices(bo::BridgeOptions, T::Float64=20.0)
# These functions are now imported from the main BridgeFEM module above

# Boundary condition application functions moved to BoundaryConditions module
# - apply_bc(M, K, so) - for both 2D and 3D arrays
# - remove_fixed_dofs(M, K, bc_dofs, total_dofs)
# These functions are now imported from the main BridgeFEM module above

# DOF mapping and support assembly functions moved to Assembly module
# - create_support_dof_mapping(bo::BridgeOptions, supports::Vector{SupportElement})
# - get_dof_from_node(bridge::BridgeOptions, supports::Vector{SupportElement}, node::Int)
# - get_bc_dofs(bridge::BridgeOptions, supports::Vector{SupportElement}, support_dof_mapping::Vector{Vector{Int}})
# - assemble_local_support(support::SupportElement, T::Float64=20.0)
# - assemble_matrices_with_supports(so::SimulationOptions)
# These functions are now imported from the main BridgeFEM module above



# Support matrix creation functions moved to Assembly module
# - create_expanded_transformation(angle::Float64, n_nodes::Int)
# - create_support_mass_matrix(support::SupportElement, ρ::Float64)
# These functions are now imported from the main BridgeFEM module above

function interpolate_matrix(M::Array{Float64,3}, Ts::Vector{Float64})
    M_interp = interpolate((1:size(M,1), 1:size(M,2), Ts), M, Gridded(Linear()))
    M_T = t -> M_interp(1:size(M,1), 1:size(M,2), t)
    return M_T
end

function setup_interpolation(M::Array{Float64,3}, K::Array{Float64,3}, Ts::Vector{Float64})
    # Create interpolation function for each temperature slice
    M_T = interpolate_matrix(M, Ts)
    K_T = interpolate_matrix(K, Ts)

    # Return a function that evaluates the matrices at a given temperature
    return M_T, K_T
end

function setup_physical(so::SimulationOptions)

    M, K = assemble_matrices_with_supports(so)
    M_T, K_T = setup_interpolation(M, K, so.temperatures)

    return M_T, K_T

end
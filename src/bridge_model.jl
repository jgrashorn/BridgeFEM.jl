using LinearAlgebra, DifferentialEquations, Plots
using Interpolations
using SparseArrays
using JSON
using Arpack
using Dates

# Import Core types and constants
using ..BridgeFEM: BCTypes, BridgeBC, BridgeOptions, SupportElement, SimulationOptions
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

function apply_bc(M::Matrix{Float64}, K::Matrix{Float64}, so::SimulationOptions)

    bc_dofs = so.bc_dofs

    n_dofs = size(K, 1)
    
    for bc in bc_dofs
        node = bc[1]
        # dof_types = bc[2]  # DOF type(s)
        # dof_indices = 3 * (node - 1) .+ dof_types  # Convert to global DOF indices

        dof_indices = node
        
        for d_ in dof_indices
            # @info "Applying boundary condition at DOF $d_"
            if d_ <= n_dofs  # Check bounds
                K[:, d_] .= 0.0
                K[d_, :] .= 0.0
                K[d_, d_] = 1.0
                M[d_, :] .= 0.0
                M[:, d_] .= 0.0
                M[d_, d_] = 0.0
            end
        end
    end
    return M, K
end

function apply_bc(M::Array{Float64,3}, K::Array{Float64,3}, so::SimulationOptions)

    M_, K_ = zeros(size(M)), zeros(size(K))

    for i in axes(M,3)
        M_[:,:,i], K_[:,:,i] = apply_bc(M[:,:,i], K[:,:,i], so)
    end

    return M_, K_
end

# DOF mapping and support assembly functions moved to Assembly module
# - create_support_dof_mapping(bo::BridgeOptions, supports::Vector{SupportElement})
# - get_dof_from_node(bridge::BridgeOptions, supports::Vector{SupportElement}, node::Int)
# - get_bc_dofs(bridge::BridgeOptions, supports::Vector{SupportElement}, support_dof_mapping::Vector{Vector{Int}})
# - assemble_local_support(support::SupportElement, T::Float64=20.0)
# - assemble_matrices_with_supports(so::SimulationOptions)
# These functions are now imported from the main BridgeFEM module above

function remove_fixed_dofs(M, K, bc_dofs::Vector{Int}, total_dofs::Int)
    # Remove fixed DOFs from mass and stiffness matrices
    retained_dofs = setdiff(1:total_dofs, bc_dofs)
    removed_dofs = bc_dofs
    
    if ndims(M) == 2
        # 2D case (single temperature)
        M_ = M[retained_dofs, retained_dofs]
        K_ = K[retained_dofs, retained_dofs]
        return M_, K_, retained_dofs, removed_dofs
    end
    M_ = M[retained_dofs, retained_dofs, :]
    K_ = K[retained_dofs, retained_dofs, :]
    return M_, K_, retained_dofs, removed_dofs
end

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
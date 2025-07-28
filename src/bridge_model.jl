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
using ..BridgeFEM: frame_elem_stiffness, transformation_matrix, rotated_stiffness
using ..BridgeFEM: frame_elem_mass, create_support_mass_matrix
using ..BridgeFEM: assemble_local_support, create_expanded_transformation










# Assemble global stiffness matrix
function assemble_stiffness!(K, bo::BridgeOptions, EA, EI)
    dx = bo.L / bo.n_elem
    for e = 1:bo.n_elem
        ke = frame_elem_stiffness(EA, EI, dx)
        # DOFs for element e: nodes e and e+1, each with 3 DOFs
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        K[dofs, dofs] .+= ke
    end
    return K
end

function assemble_matrices(bo::BridgeOptions, T::Float64=20.0)

    # Discretization
    dx = bo.L / bo.n_elem
    n_dof = bo.n_dofs    # 3 DOFs per node (u, v, theta)

    M = spzeros(n_dof, n_dof)
    K = spzeros(n_dof, n_dof)

    # Mass matrix (lumped for simplicity: translational and rotational inertia per node)
    m_trans = bo.ρ * bo.A * dx / 2  # each node shares half element mass
    m_rot   = bo.ρ * bo.A * dx^3 / 24  # rotational inertia for slender beam

    for i in 1:bo.n_nodes
        # Translational DOFs (u and v)
        M[3*(i-1)+1, 3*(i-1)+1] = m_trans  # u direction
        M[3*(i-1)+2, 3*(i-1)+2] = m_trans  # v direction
        # Rotational DOF (theta)
        M[3*(i-1)+3, 3*(i-1)+3] = m_rot    # rotation
    end

    E = bo.E(T)  # Young's modulus at temperature T
    assemble_stiffness!(K, bo, E * bo.A, E * bo.I)

    # apply_bc!(M, K, bo.bc_nodes.conds)

    return M, K
end

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



function create_support_dof_mapping(bo::BridgeOptions, supports::Vector{SupportElement})
    
    global_dof_offset = bo.n_dofs
    support_dof_maps = Vector{Vector{Int}}()
    
    for (i, support) in enumerate(supports)
        n_support_nodes = support.n_elem + 1
        n_support_dofs = 3 * n_support_nodes
        
        # Initialize mapping vector
        local_to_global = zeros(Int, n_support_dofs)
        
        # SIMPLIFIED: First node (index 1) connects to bridge
        bridge_connection_node = support.connection_node
        bridge_connection_dofs = 3 * (bridge_connection_node - 1) .+ support.connection_dofs
        
        # First node of support connects to bridge
        connection_node_dofs = [1, 2, 3]  # First node DOFs
        
        # Map connected DOFs from first node to bridge
        for (j, connection_dof) in enumerate(support.connection_dofs)
            local_support_dof = connection_node_dofs[connection_dof]
            local_to_global[local_support_dof] = bridge_connection_dofs[j]
        end
        
        # Assign new DOFs to all non-connected DOFs
        for dof_idx in 1:n_support_dofs
            if local_to_global[dof_idx] == 0  # Not yet assigned
                global_dof_offset += 1
                local_to_global[dof_idx] = global_dof_offset
            end
        end
        
        push!(support_dof_maps, local_to_global)
    end
    
    total_dofs = global_dof_offset
    return support_dof_maps, total_dofs
end

function get_node_from_dof(so::SimulationOptions, dof::Int)

    node = dof <= so.bridge.n_dofs ? dof ÷ 3 : begin
        for support_map in so.support_dof_mapping
            node = findfirst(==(dof),support_map)
            if !isempty(node)
                return so.bridge.n_elem + 1 + node
            end
        end
    end

    return node

end

function get_dof_from_node(so::SimulationOptions,node::Int)
    # Get the DOF mapping for the given node
    dof_map = node < so.bridge.n_dofs ? 3 * (node - 1) .+ [1, 2, 3] : begin
        
        # Check if node has support connections
        for support in so.supports
            if support.connection_node == node
                dof_map = support.connection_dofs
                return dof_map
            end
        end
    end
end

function get_dof_from_node(bridge::BridgeOptions, supports::Vector{SupportElement}, node::Int)
    # Get the DOF mapping for the given node
    dof_map = node < bridge.n_dofs ? 3 * (node - 1) .+ [1, 2, 3] : begin
        
        # Check if node has support connections
        for support in supports
            if support.connection_node == node
                dof_map = support.connection_dofs
                return dof_map
            end
        end
    end
end

function get_bc_dofs(bridge::BridgeOptions, supports::Vector{SupportElement}, support_dof_mapping::Vector{Vector{Int}})

    bc_dofs = Vector{Int}()

    for bc_node in bridge.bc_nodes.conds
        node = bc_node[1]
        dofs = get_dof_from_node(bridge, supports, node)
        push!(bc_dofs, dofs[bc_node[2]]...)
    end

    for (i, support) in enumerate(supports)
        end_dofs = support_dof_mapping[i][end-2:end]
        fixed_dofs = end_dofs[support.bc_bottom]
        push!(bc_dofs,fixed_dofs...)
    end

    return bc_dofs

end

function assemble_matrices_with_supports(so::SimulationOptions)
    # Get DOF mappings
    support_dof_maps, total_dofs = create_support_dof_mapping(so.bridge, so.supports)
    nTemps = length(so.temperatures)
    
    # Initialize expanded matrices
    M = zeros(so.total_dofs, so.total_dofs, nTemps)
    K = zeros(so.total_dofs, so.total_dofs, nTemps)
    
    for (t, T) in enumerate(so.temperatures)
        # Assemble main bridge (already temperature-dependent)
        M_bridge, K_bridge = assemble_matrices(so.bridge, T)

        M_ = zeros(so.total_dofs, so.total_dofs)
        K_ = zeros(so.total_dofs, so.total_dofs)

        M_[1:so.bridge.n_dofs, 1:so.bridge.n_dofs] .= M_bridge
        K_[1:so.bridge.n_dofs, 1:so.bridge.n_dofs] .= K_bridge
        
        # Assemble each support (now temperature-dependent)
        for (i, support) in enumerate(so.supports)
            # Get local support matrices at temperature T
            K_local = assemble_local_support(support, T)  # FIXED: Pass temperature
            M_local = create_support_mass_matrix(support, so.bridge.ρ)  # Mass not temperature dependent
            
            # Rotate BOTH matrices to global coordinates
            n_support_nodes = support.n_elem + 1
            T_expanded = create_expanded_transformation(support.angle, n_support_nodes)
            
            K_rotated = T_expanded' * K_local * T_expanded
            M_rotated = T_expanded' * M_local * T_expanded
            
            # Map to global DOFs
            dof_map = so.support_dof_mapping[i]
            K_[dof_map, dof_map] .+= K_rotated
            M_[dof_map, dof_map] .+= M_rotated
            
        end
        M[:,:,t] = M_
        K[:,:,t] = K_
    end

    return M, K

end

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
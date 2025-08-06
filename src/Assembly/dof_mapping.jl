# Assembly/dof_mapping.jl - DOF management and connectivity utilities  
# Extracted from bridge_model.jl lines 127-199

# Import Core types
using ..BridgeFEM: BridgeOptions, SupportElement, SimulationOptions

"""
    create_support_dof_mapping(bo::BridgeOptions, supports::Vector{SupportElement}) -> Tuple{Vector{Vector{Int}}, Int}

Create DOF mapping between support elements and global DOFs.

# Arguments
- `bo::BridgeOptions`: Bridge configuration options
- `supports::Vector{SupportElement}`: Vector of support elements

# Returns  
- `(support_dof_maps, total_dofs)`: Tuple containing vector of DOF mappings for each support and total DOF count
"""
function create_support_dof_mapping(bo::BridgeOptions, supports::Vector{SupportElement})::Tuple{Vector{Vector{Int}}, Int}
    
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

"""
    get_dof_from_node(bridge::BridgeOptions, supports::Vector{SupportElement}, node::Int) -> Vector{Int}

Get DOF indices for a specific node, considering support connections.

# Arguments
- `bridge::BridgeOptions`: Bridge configuration options
- `supports::Vector{SupportElement}`: Vector of support elements  
- `node::Int`: Node number to get DOFs for

# Returns
- Vector of DOF indices for the specified node
"""
function get_dof_from_node(bridge::BridgeOptions, supports::Vector{SupportElement}, node::Int)::Vector{Int}
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
    return dof_map
end

"""
    get_bc_dofs(bridge::BridgeOptions, supports::Vector{SupportElement}, support_dof_mapping::Vector{Vector{Int}}) -> Vector{Int}

Extract all boundary condition DOF indices from bridge and support constraints.

# Arguments
- `bridge::BridgeOptions`: Bridge configuration options
- `supports::Vector{SupportElement}`: Vector of support elements
- `support_dof_mapping::Vector{Vector{Int}}`: DOF mapping for support elements

# Returns  
- Vector of all constrained DOF indices
"""
function get_bc_dofs(bridge::BridgeOptions, supports::Vector{SupportElement}, support_dof_mapping::Vector{Vector{Int}})::Vector{Int}

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

"""
    constraint_dof_indices(boundary_conditions) -> Vector{Int}

Identify DOF indices that have constraints applied.

# Arguments
- `boundary_conditions`: Boundary condition specification

# Returns
- Vector of constrained DOF indices
"""
function constraint_dof_indices(boundary_conditions)::Vector{Int}
    # Implementation would extract constraint DOF indices from boundary conditions
    # This is a placeholder - actual implementation depends on BC data structure
    return Vector{Int}()
end

"""
    get_element_dofs(element_id::Int, connectivity) -> Vector{Int}

Extract DOF indices for a specific element based on connectivity.

# Arguments
- `element_id::Int`: Element identifier
- `connectivity`: Element connectivity data structure

# Returns  
- Vector of DOF indices for the specified element
"""
function get_element_dofs(element_id::Int, connectivity)::Vector{Int}
    # Implementation would extract element DOF indices from connectivity
    # This is a placeholder - actual implementation depends on connectivity data structure
    return Vector{Int}()
end
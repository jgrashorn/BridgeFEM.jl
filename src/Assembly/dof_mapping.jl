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
- `boundary_conditions`: Boundary condition specification (BridgeBC object)

# Returns
- Vector of constrained DOF indices
"""
function constraint_dof_indices(boundary_conditions)::Vector{Int}
    if isa(boundary_conditions, BridgeBC)
        # Extract constrained DOF indices from BridgeBC object
        constrained_dofs = Vector{Int}()
        
        for bc_condition in boundary_conditions.conds
            node = bc_condition[1]
            dof_types = bc_condition[2]
            
            # Convert node and DOF types to global DOF indices
            # DOF numbering: 3*(node-1)+1, 3*(node-1)+2, 3*(node-1)+3
            base_dof = 3 * (node - 1)
            
            if isa(dof_types, Vector)
                # Multiple DOF types specified as vector
                for dof_type in dof_types
                    push!(constrained_dofs, base_dof + dof_type)
                end
            else
                # Single DOF type
                push!(constrained_dofs, base_dof + dof_types)
            end
        end
        
        return sort(unique(constrained_dofs))
    else
        # For other boundary condition formats, return empty for now
        # Can be extended as needed for different BC data structures
        return Vector{Int}()
    end
end

"""
    get_element_dofs(element_id::Int, connectivity) -> Vector{Int}

Extract DOF indices for a specific element based on connectivity.

# Arguments
- `element_id::Int`: Element identifier (1-based indexing)
- `connectivity`: Element connectivity data structure. Can be:
  - `nothing`: Assumes consecutive node connectivity (element e connects nodes e and e+1)
  - `Matrix{Int}`: Each row represents an element, columns are connected node IDs
  - `Vector{Vector{Int}}`: Each entry contains node IDs for that element
  - `Vector{Tuple{Int,Int}}`: Each entry is a tuple of (node1, node2) for that element

# Returns  
- Vector of DOF indices for the specified element (6 DOFs for 2-node frame element)
"""
function get_element_dofs(element_id::Int, connectivity)::Vector{Int}
    # Handle different connectivity data structures
    if connectivity === nothing
        # Default consecutive node connectivity for 1D structures
        # Element e connects nodes e and e+1
        node1 = element_id
        node2 = element_id + 1
        
    elseif isa(connectivity, Matrix{Int})
        # Matrix format: each row is an element, columns are node IDs
        if element_id > size(connectivity, 1)
            throw(BoundsError("Element ID $element_id exceeds connectivity matrix size"))
        end
        node1 = connectivity[element_id, 1]
        node2 = connectivity[element_id, 2]
        
    elseif isa(connectivity, Vector{Vector{Int}})
        # Vector of vectors format
        if element_id > length(connectivity)
            throw(BoundsError("Element ID $element_id exceeds connectivity vector length"))
        end
        element_nodes = connectivity[element_id]
        if length(element_nodes) < 2
            throw(ArgumentError("Element $element_id must have at least 2 nodes"))
        end
        node1 = element_nodes[1]
        node2 = element_nodes[2]
        
    elseif isa(connectivity, Vector{Tuple{Int,Int}})
        # Vector of tuples format
        if element_id > length(connectivity)
            throw(BoundsError("Element ID $element_id exceeds connectivity vector length"))
        end
        node1, node2 = connectivity[element_id]
        
    else
        throw(ArgumentError("Unsupported connectivity data structure: $(typeof(connectivity))"))
    end
    
    # Convert node IDs to DOF indices
    # Each node has 3 DOFs: u, v, theta (numbered 3*(node-1)+1, 3*(node-1)+2, 3*(node-1)+3)
    dofs = [
        3*(node1-1)+1, 3*(node1-1)+2, 3*(node1-1)+3,  # Node 1 DOFs
        3*(node2-1)+1, 3*(node2-1)+2, 3*(node2-1)+3   # Node 2 DOFs
    ]
    
    return dofs
end
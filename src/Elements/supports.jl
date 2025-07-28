"""
Elements/supports.jl

This file contains support element functionality for structural analysis.
Provides support element assembly, DOF mapping, and transformation functions
for modeling flexible supports, springs, and soil-structure interaction.
"""

# File is included in BridgeFEM module scope
# All imports are available from the main module context

"""
    assemble_local_support(support::SupportElement, T::Float64=20.0)

Assemble the local stiffness matrix for a support element.

# Arguments
- `support`: SupportElement struct containing properties
- `T`: Temperature for temperature-dependent material properties (°C)

# Returns
- `K_local`: Local stiffness matrix for the support element

This function assembles the local stiffness matrix for a support element
using temperature-dependent Young's modulus through interpolation.
"""
function assemble_local_support(support::SupportElement, T::Float64=20.0)
    n_nodes = support.n_elem + 1
    n_dofs = 3 * n_nodes
    dx = support.L / support.n_elem
    
    K_local = zeros(n_dofs, n_dofs)
    
    # Get temperature-dependent Young's modulus
    E = support.E(T)
    EA = E * support.A
    EI = E * support.I
    
    # Assemble support elements
    for e = 1:support.n_elem
        ke = frame_elem_stiffness(EA, EI, dx)
        # DOFs for element e: nodes e and e+1
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        K_local[dofs, dofs] .+= ke
    end
    
    return K_local
end

"""
    create_expanded_transformation(angle::Float64, n_nodes::Int)

Create an expanded transformation matrix for multiple nodes.

# Arguments
- `angle`: Rotation angle in degrees
- `n_nodes`: Number of nodes in the element

# Returns
- `T_expanded`: Expanded transformation matrix

This function creates a transformation matrix that applies 2D rotation
to the x,y DOFs of multiple nodes while leaving rotational DOFs unchanged.
The transformation is used for inclined support elements.
"""
function create_expanded_transformation(angle::Float64, n_nodes::Int)
    n_dofs = 3 * n_nodes
    T_expanded = Matrix{Float64}(LinearAlgebra.I, n_dofs, n_dofs)
    
    # Create single node transformation matrix for 2D rotation
    c = cosd(angle)
    s = sind(angle)
    T_single = [c s; -s c]
    
    for node = 1:n_nodes
        dof_start = 3 * (node - 1) + 1
        
        # Apply 2D rotation to x,y DOFs (rotation DOF unchanged)
        T_expanded[dof_start:dof_start+1, dof_start:dof_start+1] = T_single
    end
    
    return T_expanded
end

# End of supports.jl 
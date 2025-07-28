"""
Elements/stiffness.jl

This file contains finite element stiffness matrix computation functions.
Provides frame element stiffness matrices, transformation matrices, and
rotated stiffness computations for structural analysis.
"""

# File is included in BridgeFEM module scope
# All imports are available from the main module context

"""
    frame_elem_stiffness(EA, EI, L_e)

Compute the local stiffness matrix for a 2D frame element.

# Arguments
- `EA`: Axial stiffness (Young's modulus × Cross-sectional area)
- `EI`: Flexural stiffness (Young's modulus × Moment of inertia)
- `L_e`: Element length

# Returns
- `k`: 6×6 local stiffness matrix

The stiffness matrix is organized with DOFs as:
[u₁, v₁, θ₁, u₂, v₂, θ₂] where subscripts denote node numbers
and u, v, θ are axial, transverse, and rotational DOFs respectively.
"""
function frame_elem_stiffness(EA, EI, L_e)
    c1 = EA / L_e
    c2 = 12.0 * EI / L_e^3
    c3 = 6.0 * EI / L_e^2
    c4 = 4.0 * EI / L_e
    c5 = 2.0 * EI / L_e
    
    k = [
        c1    0.0   0.0     -c1   0.0   0.0;
        0.0   c2    c3      0.0   -c2   c3;
        0.0   c3    c4      0.0   -c3   c5;
        -c1   0.0   0.0     c1    0.0   0.0;
        0.0   -c2   -c3     0.0   c2    -c3;
        0.0   c3    c5      0.0   -c3   c4
    ]
    return k
end

"""
    transformation_matrix(θ)

Compute the transformation matrix from local to global coordinates.

# Arguments
- `θ`: Rotation angle in degrees

# Returns
- `T`: 6×6 transformation matrix

The transformation matrix rotates local element coordinates to global
coordinates for a 2D frame element with 3 DOFs per node.
"""
function transformation_matrix(θ)
    c = cosd(θ)
    s = sind(θ)
    T = Matrix{Float64}(LinearAlgebra.I, 6, 6)
    T[1,1] =  c; T[1,2] =  s
    T[2,1] = -s; T[2,2] =  c
    T[4,4] =  c; T[4,5] =  s
    T[5,4] = -s; T[5,5] =  c
    return T
end

"""
    rotated_stiffness(support::SupportElement)

Compute the global stiffness matrix for a rotated support element.

# Arguments
- `support`: SupportElement struct containing geometric and material properties

# Returns
- `ke_global`: 6×6 global stiffness matrix
- `dofs`: DOF connectivity information

This function combines local stiffness computation with coordinate transformation
to produce the global stiffness matrix for an inclined support element.
"""
function rotated_stiffness(support::SupportElement)
    ke_local = frame_elem_stiffness(support.EA, support.EI, support.L)
    Te = transformation_matrix(support.angle)
    ke_global = Te' * ke_local * Te
    return ke_global, support.dofs
end

# End of stiffness.jl 
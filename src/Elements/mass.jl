"""
Elements/mass.jl

This file contains finite element mass matrix computation functions.
Provides frame element mass matrices and support mass functionality
for dynamic structural analysis.
"""

# File is included in BridgeFEM module scope
# All imports are available from the main module context

"""
    frame_elem_mass(ρ::Float64, A::Float64, L::Float64)

Compute the consistent mass matrix for a 2D frame element.

# Arguments
- `ρ`: Material density (kg/m³)
- `A`: Cross-sectional area (m²)
- `L`: Element length (m)

# Returns
- `M`: 6×6 consistent mass matrix

The mass matrix uses the consistent mass formulation for a 2D frame element
with translational and rotational inertia. DOF organization follows:
[u₁, v₁, θ₁, u₂, v₂, θ₂] where subscripts denote node numbers.
"""
function frame_elem_mass(ρ::Float64, A::Float64, L::Float64)
    mass_factor = ρ * A * L / 420.0
    
    M = mass_factor * [
        140   0      0     70    0      0    ;
        0     156    22L   0     54    -13L  ;
        0     22L    4L^2  0     13L   -3L^2 ;
        70    0      0     140   0      0    ;
        0     54     13L   0     156   -22L  ;
        0    -13L   -3L^2  0    -22L    4L^2
    ]
    
    return M
end

"""
    create_support_mass_matrix(support::SupportElement, ρ::Float64)

Create the mass matrix for a support element.

# Arguments
- `support`: SupportElement struct containing geometric properties
- `ρ`: Material density (kg/m³)

# Returns
- `M_local`: Local mass matrix for the support element

This function assembles the mass matrix for a support element using
consistent mass formulation. The mass matrix is not temperature dependent
as density and geometry are assumed constant.
"""
function create_support_mass_matrix(support::SupportElement, ρ::Float64)
    n_nodes = support.n_elem + 1
    n_dofs = 3 * n_nodes
    dx = support.L / support.n_elem
    
    M_local = zeros(n_dofs, n_dofs)
    
    # Mass matrix is not temperature dependent (density and geometry constant)
    # Use reference area for mass calculation
    A_ref = support.A
    
    # Assemble mass matrix for each element
    for e = 1:support.n_elem
        me = frame_elem_mass(ρ, A_ref, dx)
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        M_local[dofs, dofs] .+= me
    end
    
    return M_local
end

# End of mass.jl 
# Elements Module - Finite Element Computations
# 
# This module contains core finite element computation functions for frame elements
# including stiffness matrices, mass matrices, and coordinate transformations.

module FiniteElements

using LinearAlgebra

export frame_elem_stiffness, frame_elem_mass, transformation_matrix

"""
    frame_elem_stiffness(EA, EI, L_e)

Compute the 6x6 element stiffness matrix for a frame element in local coordinates.

# Arguments
- `EA`: Axial stiffness (E × A where E is elastic modulus and A is cross-sectional area)
- `EI`: Flexural stiffness (E × I where I is second moment of area)
- `L_e`: Element length

# Returns
- `k`: 6x6 stiffness matrix in local element coordinates

The stiffness matrix follows standard frame element formulation with degrees of freedom:
[u1, v1, θ1, u2, v2, θ2] where u,v are translations and θ is rotation.
"""
function frame_elem_stiffness(EA::Real, EI::Real, L_e::Real)::Matrix{Float64}
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
    frame_elem_mass(ρ::Float64, A::Float64, L::Float64)

Compute the 6x6 consistent mass matrix for a frame element.

# Arguments
- `ρ`: Material density
- `A`: Cross-sectional area
- `L`: Element length

# Returns
- `M`: 6x6 consistent mass matrix

The mass matrix uses the standard finite element formulation with the same
degrees of freedom as the stiffness matrix: [u1, v1, θ1, u2, v2, θ2].
"""
function frame_elem_mass(ρ::Float64, A::Float64, L::Float64)::Matrix{Float64}
    # Pre-compute mass factor for efficiency
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
    transformation_matrix(θ)

Compute the 6x6 coordinate transformation matrix for a rotated frame element.

# Arguments
- `θ`: Rotation angle in degrees (measured counterclockwise from global x-axis)

# Returns
- `T`: 6x6 transformation matrix to convert from local to global coordinates

The transformation matrix converts element matrices from local element coordinates
to global structural coordinates. For a vector in local coordinates v_local,
the global coordinates are: v_global = T * v_local
"""
function transformation_matrix(θ::Real)::Matrix{Float64}
    c = cosd(θ)
    s = sind(θ)
    T = Matrix{Float64}(LinearAlgebra.I, 6, 6)
    T[1,1] =  c; T[1,2] =  s
    T[2,1] = -s; T[2,2] =  c
    T[4,4] =  c; T[4,5] =  s
    T[5,4] = -s; T[5,5] =  c
    return T
end

end # module FiniteElements
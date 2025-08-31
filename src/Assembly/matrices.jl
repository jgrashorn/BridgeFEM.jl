# Assembly/matrices.jl - Global matrix assembly operations
# Extracted from bridge_model.jl lines 24-63 and related functions

using SparseArrays
using LinearAlgebra
using Interpolations

# Import Core types and functions
using ..BridgeFEM: BridgeOptions, SupportElement, SimulationOptions

# Import Elements module functions  
using ..BridgeFEM: frame_elem_stiffness, frame_elem_mass, transformation_matrix

"""
    assemble_stiffness!(K, bo::BridgeOptions, EA::Real, EI::Real) -> Matrix

Assemble global stiffness matrix K for bridge structure.

# Arguments
- `K`: Pre-allocated global stiffness matrix
- `bo::BridgeOptions`: Bridge configuration options
- `EA::Real`: Axial stiffness (Young's modulus × Area)
- `EI::Real`: Flexural stiffness (Young's modulus × Moment of inertia)

# Returns
- Modified global stiffness matrix K
"""
function assemble_stiffness!(K, bo::BridgeOptions, EA::Real, EI::Real)::typeof(K)
    dx = bo.L / bo.n_elem
    for e = 1:bo.n_elem
        ke = frame_elem_stiffness(EA, EI, dx)
        # DOFs for element e: nodes e and e+1, each with 3 DOFs
        dofs = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*e+1, 3*e+2, 3*e+3]
        K[dofs, dofs] .+= ke
    end
    return K
end

"""
    assemble_matrices(bo::BridgeOptions, T::Float64=20.0) -> Tuple{SparseMatrixCSC, SparseMatrixCSC}

Assemble global mass and stiffness matrices for bridge structure at given temperature.

# Arguments
- `bo::BridgeOptions`: Bridge configuration options
- `T::Float64`: Temperature for temperature-dependent properties (default: 20.0°C)

# Returns
- `(M, K)`: Tuple of global mass matrix M and stiffness matrix K as sparse matrices
"""
function assemble_matrices(bo::BridgeOptions, T::Float64=20.0)::Tuple{SparseMatrixCSC{Float64,Int64}, SparseMatrixCSC{Float64,Int64}}

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

    return M, K
end

"""
    assemble_matrices_with_supports(so::SimulationOptions) -> Tuple{Array{Float64,3}, Array{Float64,3}}

Assemble global mass and stiffness matrices including support elements for multiple temperatures.

# Arguments
- `so::SimulationOptions`: Complete simulation configuration with bridge, supports, and temperatures

# Returns
- `(M, K)`: Tuple of 3D arrays containing mass and stiffness matrices for each temperature
"""
function assemble_matrices_with_supports(so::SimulationOptions)::Tuple{Array{Float64,3}, Array{Float64,3}}
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
            
            # Map to global DOFs - create mapping if not provided
            dof_map = if isempty(so.support_dof_mapping)
                # Fallback: create DOF mapping on-the-fly for this support
                mapping, _ = create_support_dof_mapping(so.bridge, [support])
                mapping[1]  # Get mapping for this single support
            else
                so.support_dof_mapping[i]
            end
            K_[dof_map, dof_map] .+= K_rotated
            M_[dof_map, dof_map] .+= M_rotated
            
        end
        M[:,:,t] = M_
        K[:,:,t] = K_
    end

    return M, K
end

# Helper functions for support matrix assembly (previously in bridge_model.jl)

"""
    assemble_local_support(support::SupportElement, T::Float64=20.0) -> Matrix{Float64}

Assemble local stiffness matrix for a support element at given temperature.

# Arguments  
- `support::SupportElement`: Support element configuration
- `T::Float64`: Temperature for temperature-dependent properties (default: 20.0°C)

# Returns
- Local stiffness matrix for the support element
"""
function assemble_local_support(support::SupportElement, T::Float64=20.0)::Matrix{Float64}
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
    create_support_mass_matrix(support::SupportElement, ρ::Float64) -> Matrix{Float64}

Create mass matrix for a support element (temperature-independent).

# Arguments
- `support::SupportElement`: Support element configuration
- `ρ::Float64`: Material density

# Returns  
- Local mass matrix for the support element
"""
function create_support_mass_matrix(support::SupportElement, ρ::Float64)::Matrix{Float64}
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

"""
    create_expanded_transformation(angle::Float64, n_nodes::Int) -> Matrix{Float64}

Create expanded transformation matrix for rotating support element to global coordinates.

# Arguments
- `angle::Float64`: Rotation angle in radians
- `n_nodes::Int`: Number of nodes in support element

# Returns
- Expanded transformation matrix for all DOFs
"""
function create_expanded_transformation(angle::Float64, n_nodes::Int)::Matrix{Float64}
    n_dofs = 3 * n_nodes
    T_expanded = Matrix{Float64}(LinearAlgebra.I, n_dofs, n_dofs)
    
    T_single = transformation_matrix(angle)
    
    for node = 1:n_nodes
        dof_start = 3 * (node - 1) + 1
        dof_end = dof_start + 2
        
        # Apply 2D rotation to x,y DOFs (rotation DOF unchanged)
        T_expanded[dof_start:dof_start+1, dof_start:dof_start+1] = T_single[1:2, 1:2]
    end
    
    return T_expanded
end

"""
    interpolate_matrix(M::Array{Float64,3}, Ts::Vector{Float64}) -> Function

Create interpolation function for 3D matrix array across temperature dimension.

# Arguments
- `M::Array{Float64,3}`: 3D matrix array with temperature as third dimension
- `Ts::Vector{Float64}`: Temperature values corresponding to third dimension

# Returns
- Function that interpolates matrix values at any temperature
"""
function interpolate_matrix(M::Array{Float64,3}, Ts::Vector{Float64})
    M_interp = interpolate((1:size(M,1), 1:size(M,2), Ts), M, Gridded(Linear()))
    M_T = t -> M_interp(1:size(M,1), 1:size(M,2), t)
    return M_T
end

"""
    setup_matrix_interpolation(M::Array{Float64,3}, K::Array{Float64,3}, Ts::Vector{Float64}) -> Tuple{Function, Function}

Set up interpolation functions for mass and stiffness matrices across temperature.

# Arguments
- `M::Array{Float64,3}`: 3D mass matrix array with temperature as third dimension
- `K::Array{Float64,3}`: 3D stiffness matrix array with temperature as third dimension  
- `Ts::Vector{Float64}`: Temperature values corresponding to third dimension

# Returns
- `(M_T, K_T)`: Tuple of interpolation functions for mass and stiffness matrices
"""
function setup_matrix_interpolation(M::Array{Float64,3}, K::Array{Float64,3}, Ts::Vector{Float64})
    # Create interpolation function for each temperature slice
    M_T = interpolate_matrix(M, Ts)
    K_T = interpolate_matrix(K, Ts)

    # Return a function that evaluates the matrices at a given temperature
    return M_T, K_T
end

"""
    setup_physical(so::SimulationOptions) -> Tuple{Function, Function}

Set up complete physical system with temperature-dependent matrix interpolation.

# Arguments
- `so::SimulationOptions`: Complete simulation configuration with bridge, supports, and temperatures

# Returns
- `(M_T, K_T)`: Tuple of interpolation functions for mass and stiffness matrices
"""
function setup_physical(so::SimulationOptions)
    M, K = assemble_matrices_with_supports(so)
    M_T, K_T = setup_matrix_interpolation(M, K, so.temperatures)
    return M_T, K_T
end
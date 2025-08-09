"""
# BoundaryConditions/Application Functions

Boundary condition application functions for BridgeFEM.jl.

Contains functions for applying boundary conditions to global matrices and
identifying constrained degrees of freedom.

## Functions
- `apply_bc(M, K, so)`: Apply boundary conditions to mass and stiffness matrices
- `remove_fixed_dofs(M, K, bc_dofs, total_dofs)`: Remove fixed DOFs from matrices

## Dependencies
- LinearAlgebra (stdlib): Matrix operations
- SparseArrays (stdlib): Sparse matrix support
- Core module: SimulationOptions and related types
"""

# LinearAlgebra and SparseArrays are already imported in the main module

"""
    apply_bc(M::Matrix{Float64}, K::Matrix{Float64}, so) -> (Matrix{Float64}, Matrix{Float64})

Apply boundary conditions to mass and stiffness matrices for a single temperature case.

Modifies the global mass and stiffness matrices to enforce boundary conditions by:
- Setting constrained DOF rows and columns to zero
- Setting diagonal terms to identity values for stiffness matrix
- Setting diagonal terms to zero for mass matrix (no inertia at fixed DOFs)

# Arguments
- `M::Matrix{Float64}`: Global mass matrix
- `K::Matrix{Float64}`: Global stiffness matrix  
- `so`: Simulation options containing boundary condition information

# Returns
- `(M_modified, K_modified)`: Tuple of modified mass and stiffness matrices

# Notes
Boundary conditions are applied by zeroing rows/columns and setting appropriate diagonal values.
This method preserves matrix symmetry and conditioning for finite element solvers.
"""
function apply_bc(M::Matrix{Float64}, K::Matrix{Float64}, so)::Tuple{Matrix{Float64}, Matrix{Float64}}
    bc_dofs = so.bc_dofs
    n_dofs = size(K, 1)
    +
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

"""
    apply_bc(M::Array{Float64,3}, K::Array{Float64,3}, so) -> (Array{Float64,3}, Array{Float64,3})

Apply boundary conditions to mass and stiffness matrices for multiple temperature cases.

Applies boundary conditions to each temperature slice of the 3D arrays representing
temperature-dependent mass and stiffness matrices.

# Arguments
- `M::Array{Float64,3}`: Global mass matrices (dof × dof × temperature)
- `K::Array{Float64,3}`: Global stiffness matrices (dof × dof × temperature)
- `so`: Simulation options containing boundary condition information

# Returns
- `(M_modified, K_modified)`: Tuple of modified 3D mass and stiffness matrix arrays

# Notes
This function applies boundary conditions consistently across all temperature cases,
preserving the temperature-dependent behavior while enforcing structural constraints.
"""
function apply_bc(M::Array{Float64,3}, K::Array{Float64,3}, so)::Tuple{Array{Float64,3}, Array{Float64,3}}
    M_, K_ = zeros(size(M)), zeros(size(K))

    for i in axes(M,3)
        M_[:,:,i], K_[:,:,i] = apply_bc(M[:,:,i], K[:,:,i], so)
    end

    return M_, K_
end

"""
    remove_fixed_dofs(M, K, bc_dofs::Vector{Int}, total_dofs::Int) -> (Matrix/Array, Matrix/Array, Vector{Int}, Vector{Int})

Remove fixed degrees of freedom from mass and stiffness matrices.

Creates reduced matrices by removing rows and columns corresponding to fixed DOFs,
resulting in smaller systems for computational efficiency.

# Arguments
- `M`: Mass matrix (2D or 3D array)
- `K`: Stiffness matrix (2D or 3D array)
- `bc_dofs::Vector{Int}`: DOF indices to remove (fixed boundary conditions)
- `total_dofs::Int`: Total number of degrees of freedom

# Returns
- `(M_reduced, K_reduced, retained_dofs, removed_dofs)`: Reduced matrices and DOF mappings

# Notes
This function preserves the matrix structure while reducing computational cost by
eliminating known-zero displacement DOFs from the system of equations.
"""
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

# Functions are now part of the main BridgeFEM module

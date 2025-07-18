"""
    decompose_matrices(M, K)

Perform eigenvalue decomposition of mass and stiffness matrices across temperature range.

This function computes natural frequencies and mode shapes for each temperature state,
ensuring mode shape consistency across temperatures using modal correlation.

# Arguments
- `M::Array{Float64,3}`: Mass matrices for each temperature [n_dof × n_dof × n_temps]
- `K::Array{Float64,3}`: Stiffness matrices for each temperature [n_dof × n_dof × n_temps]

# Returns
- `ω_matrix::Matrix{Float64}`: Natural frequencies [n_modes × n_temps] in Hz
- `Φ_tensor::Array{Float64,3}`: Mass-normalized mode shapes [n_dof × n_modes × n_temps]
- `Φ_tensor_unnormalized::Array{Float64,3}`: Unnormalized mode shapes for reference

# Algorithm
1. Solve generalized eigenvalue problem Kφ = λMφ for each temperature
2. Extract natural frequencies: ω = √λ/(2π)
3. Normalize mode shapes to unit length
4. Apply Modal Assurance Criterion (MAC) for mode tracking
5. Mass-normalize final mode shapes: φᵀMφ = I

# Mode Tracking
Mode shapes are tracked across temperatures by comparing with previous temperature
state using dot product correlation. Modes with negative correlation are flipped
to maintain consistent modal coordinate definitions.

# Notes
- Returns both normalized and unnormalized mode shapes for flexibility
- Natural frequencies are converted from rad/s to Hz
- All n_dof modes are retained (no truncation at this stage)

# See Also
- [`assemble_and_decompose`](@ref): Complete assembly and decomposition workflow
- [`Modal Assurance Criterion`](https://en.wikipedia.org/wiki/Modal_assurance_criterion)
"""
function decompose_matrices(M, K)
    # Collect mode shapes at different temperatures
    n_dof = size(M, 1)
    n_temps = size(M, 3)
    n_mode = n_dof # Number of modes to keep

    Φ_tensor = zeros(n_dof, n_mode, n_temps)
    Φ_tensor_unnormalized = zeros(n_dof, n_mode, n_temps)
    ω_matrix = zeros(n_mode, n_temps)

    decomp = [eigen(K[:,:,i], M[:,:,i]) for i in axes(M, 3)]

    λs_all = [d.values for d in decomp]
    vecs_all = [d.vectors for d in decomp]

    for i in axes(M, 3)
        λs = λs_all[i]
        vecs = vecs_all[i]
        # Normalize mode shapes (optional)
        for i in 1:n_mode
            vecs[:,i] ./= norm(vecs[:,i])
        end

        Φ_tensor[:,:,i] .= vecs

        if i > 1
            for dof in 1:n_mode
                dot(Φ_tensor[:,dof,i], Φ_tensor[:,dof,i-1])
                if dot(Φ_tensor[:,dof,i], Φ_tensor[:,dof,i-1]) < 0
                    # @info "Flipping mode shape $i"
                    Φ_tensor[:,dof,i] .*= -1
                end
            end
        end

        ωs = sqrt.(real(λs))./ (2π)  # Convert eigenvalues to natural frequencies (Hz)
        Φ_tensor_unnormalized[:,:,i] .= Φ_tensor[:,:,i]

        # Mass-normalize the mode shapes
        for j in axes(Φ_tensor, 2)
            mi = Φ_tensor[:, j, i]' * M[:,:,i] * Φ_tensor[:, j, i]
            Φ_tensor[:, j, i] ./= sqrt(mi)
        end

        ω_matrix[:,i] .= ωs
    end

    return ω_matrix, Φ_tensor, Φ_tensor_unnormalized  # Return both normalized and unnormalized mode shapes
end

function assemble_and_decompose(so::SimulationOptions)

    nTs = length(so.temperatures)
    @info "Assembling matrices with supports for $nTs temperatures"
    
    M, K = assemble_matrices_with_supports(so)
    M_ = zeros(size(M))
    K_ = zeros(size(K))
    # M_, K_, _ = remove_fixed_dofs(M, K, so.bc_dofs, so.total_dofs)
    # matrices = [apply_bc(M[:,:,i], K[:,:,i], so) for i in axes(M, 3)]

    M_, K_ = apply_bc(M, K, so)

    @info "Decomposing matrices"
    λs, vectors, vectors_unnormalized = decompose_matrices(M_, K_)

    @info "Found $(size(λs, 1)) modes across $nTs temperatures"

    keep_modes = λs[:,1] .< bo.cutoff_freq
    λs = λs[keep_modes, :]

    @info "Keeping $(sum(keep_modes)) modes below cutoff frequency $(bo.cutoff_freq) Hz"

    vectors = vectors[:, keep_modes, :]
    vectors_unnormalized = vectors_unnormalized[:, keep_modes, :]

    return M, K, λs, vectors, vectors_unnormalized
end

function setup_interpolation(λs::Matrix{Float64}, vectors::Array{Float64,3}, Ts::Vector{Float64})

    λ_T = interpolate_modes(λs, Ts)
    Φ_T = interpolate_modes(vectors, Ts)

    return λ_T, Φ_T
end

function interpolate_modes(vec::Matrix{Float64}, Ts)
    n_modes = size(vec, 1)

    λ_interp = interpolate((1:n_modes, Ts), vec, Gridded(Linear()))
    λ_T = t -> λ_interp(1:n_modes, t)

    return λ_T
end

function interpolate_modes(mat::Array{Float64,3}, Ts)
    n_modes = size(mat, 2)
    total_dofs = size(mat, 1)

    Φ_interp = interpolate((1:total_dofs, 1:n_modes, Ts), mat, Gridded(Linear()))
    Φ_T = t -> Φ_interp(1:total_dofs, 1:n_modes, t)

    return Φ_T
end

"""
    reconstruct_physical(so::SimulationOptions, q_full, Φ_interp, T_func, time)

"""
    reconstruct_physical(so::SimulationOptions, q_full, Φ_interp, T_func, time)

function reconstruct_physical(so::SimulationOptions, q_full, Φ_interp, T_func, time)
    
    bo = so.bridge
    # Auto-detect total DOFs from the interpolator dimensions
    total_dofs = so.total_dofs
    
    n_modes_total = size(q_full, 1)
    n_modes = n_modes_total ÷ 2
    n_times = length(time)

    u_full  = zeros(total_dofs, n_times)
    du_full = zeros(total_dofs, n_times)

    for (i, t) in enumerate(time)
        T_now = T_func(t)
        Φ = Φ_interp(T_now)

        q_disp = q_full[1:n_modes, i]
        q_vel  = q_full[n_modes+1:end, i]

        u_full[:, i]  .= Φ * q_disp
        du_full[:, i] .= Φ * q_vel
    end

    return u_full, du_full
end

function setup_ROM(so::SimulationOptions)

    M, K, λs, vectors, vectors_unnormalized = assemble_and_decompose(so)

    n_modes = size(λs, 1)
    λ_T, Φ_T = setup_interpolation(λs, vectors, so.temperatures)

    return λ_T, Φ_T, n_modes

end
"""
# Core Constants for BridgeFEM.jl

This module defines fundamental constants used throughout the finite element analysis,
particularly for boundary condition specification and constraint mapping.
"""

if !@isdefined(BCTypes)
    """
        BCTypes

    Dictionary mapping boundary condition type strings to degrees of freedom (DOF) constraints.

    Used for convenient specification of common boundary condition patterns in structural analysis.

    # Constraint Types
    - `"all"`: All DOFs constrained [1,2,3] (fixed support)
    - `"trans"`: Translational DOFs only [1,2] (pinned support)  
    - `"x"`: X-direction translation only [1]
    - `"y"`: Y-direction translation only [2]
    - `"ϕ"`: Rotation only [3]

    # DOF Numbering
    1. X-direction translation (u)
    2. Y-direction translation (v)  
    3. Rotation about Z-axis (θ)

    # Examples
    ```julia
    # Fixed support at node 1 (all DOFs constrained)
    bc = BridgeBC([[1, "all"]])

    # Pinned support at node 5 (translations constrained, rotation free)
    bc = BridgeBC([[5, "trans"]])

    # Roller support in Y-direction at node 3
    bc = BridgeBC([[3, "y"]])
    ```
    """
    const BCTypes = Dict(
        "all" => [1,2,3], # All DOFs (x, y, rotation)
        "trans" => [1,2], # Both translational DOFs
        "x" => 1,    # X-direction translational DOF
        "y" => 2,    # Y-direction translational DOF
        "ϕ" => 3,    # Rotational DOF
    )
end 
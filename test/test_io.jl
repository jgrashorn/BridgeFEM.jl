"""
Unit tests for IO module.

Tests JSON serialization functionality including:
- Round-trip serialization compatibility with existing configuration files
- Type preservation and reconstruction accuracy
- Numerical precision maintenance
- Compatibility with existing .json configuration file formats

Requires Test.jl framework with @testset organization.
"""

using Test
using JSON
using Interpolations

# Import main BridgeFEM module for types and IO functions
using BridgeFEM
using BridgeFEM: BCTypes, BridgeBC, BridgeOptions, SupportElement, SimulationOptions
using BridgeFEM: bridge_options_to_dict, dict_to_bridge_options
using BridgeFEM: support_element_to_dict, dict_to_support_element
using BridgeFEM: simulation_options_to_dict, load_simulation_options, save_simulation_options

@testset "IO Module Tests" begin
    
    @testset "BridgeOptions Serialization Round-Trip" begin
        # Create test BridgeOptions with temperature-dependent properties
        E_T = [0.0 200e9; 50.0 180e9; 100.0 160e9]  # Multi-point temperature dependence
        bc_nodes = BridgeBC([[1, "all"], [10, "trans"]])  # Mixed boundary conditions (using "trans" instead of "pin")
        
        bridge_orig = BridgeOptions(
            15,                    # n_elem
            bc_nodes,             # bc_nodes
            25.5,                 # L
            2600.0,               # ρ
            0.125,                # A
            0.00083,              # I
            E_T,                  # E_T
            45.5                  # cutoff_freq
        )
        
        # Convert to dictionary
        bridge_dict = bridge_options_to_dict(bridge_orig)
        
        # Test dictionary structure and values
        @test bridge_dict["n_elem"] == 15
        @test bridge_dict["L"] ≈ 25.5 rtol=1e-12
        @test bridge_dict["ρ"] ≈ 2600.0 rtol=1e-12
        @test bridge_dict["A"] ≈ 0.125 rtol=1e-12
        @test bridge_dict["I"] ≈ 0.00083 rtol=1e-12
        @test bridge_dict["cutoff_freq"] ≈ 45.5 rtol=1e-12
        
        # Test boundary condition preservation (checks resolved DOF arrays)
        @test length(bridge_dict["bc_nodes"]) == 2
        @test bridge_dict["bc_nodes"][1] == [1, [1, 2, 3]]  # "all" resolves to [1, 2, 3]
        @test bridge_dict["bc_nodes"][2] == [10, [1, 2]]    # "trans" resolves to [1, 2]
        
        # Test E_T matrix preservation
        @test length(bridge_dict["E_T"]) == 3
        @test bridge_dict["E_T"][1] ≈ [0.0, 200e9] rtol=1e-12
        @test bridge_dict["E_T"][2] ≈ [50.0, 180e9] rtol=1e-12
        @test bridge_dict["E_T"][3] ≈ [100.0, 160e9] rtol=1e-12
        
        # Convert back to BridgeOptions
        bridge_reconstructed = dict_to_bridge_options(bridge_dict)
        
        # Test exact reconstruction
        @test bridge_reconstructed.n_elem == bridge_orig.n_elem
        @test bridge_reconstructed.L ≈ bridge_orig.L rtol=1e-12
        @test bridge_reconstructed.ρ ≈ bridge_orig.ρ rtol=1e-12
        @test bridge_reconstructed.A ≈ bridge_orig.A rtol=1e-12
        @test bridge_reconstructed.I ≈ bridge_orig.I rtol=1e-12
        @test bridge_reconstructed.cutoff_freq ≈ bridge_orig.cutoff_freq rtol=1e-12
        
        # Test boundary condition reconstruction
        @test length(bridge_reconstructed.bc_nodes.conds) == length(bridge_orig.bc_nodes.conds)
        for (orig, recon) in zip(bridge_orig.bc_nodes.conds, bridge_reconstructed.bc_nodes.conds)
            @test orig[1] == recon[1]  # Node number
            @test orig[2] == recon[2]  # Constraint type
        end
        
        # Test E_T matrix reconstruction
        @test size(bridge_reconstructed.E_T) == size(bridge_orig.E_T)
        @test bridge_reconstructed.E_T ≈ bridge_orig.E_T rtol=1e-12
        
        # Test calculated properties consistency
        @test bridge_reconstructed.n_nodes == bridge_orig.n_nodes
        @test bridge_reconstructed.n_dofs == bridge_orig.n_dofs
    end
    
    @testset "SupportElement Serialization Round-Trip" begin
        # Create test SupportElement with complex temperature dependence
        E_T = [-20.0 220e9; 0.0 200e9; 40.0 180e9; 80.0 160e9]  # Multi-point curve
        
        support_orig = SupportElement(
            7,                    # connection_node
            [1, 3],              # connection_dofs (partial DOF connection)
            -15.5,               # angle (negative angle)
            4,                   # n_elem
            0.08,                # A
            0.0006,              # I
            E_T,                 # E_T
            12.75,               # L
            [1, 2, 6]            # bc_bottom (mixed constraint pattern)
        )
        
        # Convert to dictionary
        support_dict = support_element_to_dict(support_orig)
        
        # Test dictionary structure and values
        @test support_dict["connection_node"] == 7
        @test support_dict["connection_dofs"] == [1, 3]
        @test support_dict["angle"] ≈ -15.5 rtol=1e-12
        @test support_dict["n_elem"] == 4
        @test support_dict["A"] ≈ 0.08 rtol=1e-12
        @test support_dict["I"] ≈ 0.0006 rtol=1e-12
        @test support_dict["L"] ≈ 12.75 rtol=1e-12
        @test support_dict["bc_bottom"] == [1, 2, 6]
        
        # Test E_T preservation with negative temperatures
        @test length(support_dict["E_T"]) == 4
        @test support_dict["E_T"][1] ≈ [-20.0, 220e9] rtol=1e-12
        @test support_dict["E_T"][2] ≈ [0.0, 200e9] rtol=1e-12
        @test support_dict["E_T"][3] ≈ [40.0, 180e9] rtol=1e-12
        @test support_dict["E_T"][4] ≈ [80.0, 160e9] rtol=1e-12
        
        # Convert back to SupportElement
        support_reconstructed = dict_to_support_element(support_dict)
        
        # Test exact reconstruction
        @test support_reconstructed.connection_node == support_orig.connection_node
        @test support_reconstructed.connection_dofs == support_orig.connection_dofs
        @test support_reconstructed.angle ≈ support_orig.angle rtol=1e-12
        @test support_reconstructed.n_elem == support_orig.n_elem
        @test support_reconstructed.A ≈ support_orig.A rtol=1e-12
        @test support_reconstructed.I ≈ support_orig.I rtol=1e-12
        @test support_reconstructed.L ≈ support_orig.L rtol=1e-12
        @test support_reconstructed.bc_bottom == support_orig.bc_bottom
        
        # Test E_T matrix reconstruction
        @test size(support_reconstructed.E_T) == size(support_orig.E_T)
        @test support_reconstructed.E_T ≈ support_orig.E_T rtol=1e-12
        
        # Test E function behavior (interpolation consistency)
        test_temps = [-10.0, 20.0, 60.0]
        for T in test_temps
            @test support_reconstructed.E(T) ≈ support_orig.E(T) rtol=1e-10
        end
    end
    
    @testset "SimulationOptions Complete Round-Trip" begin
        # Create comprehensive test configuration
        E_T_bridge = [0.0 200e9; 100.0 160e9]
        bc_nodes = BridgeBC([[1, "all"], [20, "y"]])  # Using "y" instead of "roller"
        bridge = BridgeOptions(20, bc_nodes, 50.0, 2500.0, 0.15, 0.002, E_T_bridge, 60.0)
        
        # Create multiple support elements with different configurations
        E_T_support1 = [0.0 210e9; 50.0 190e9; 100.0 170e9]
        support1 = SupportElement(5, [1, 2], 0.0, 3, 0.1, 0.001, E_T_support1, 8.0, [1, 2, 3])
        
        E_T_support2 = [-50.0 250e9; 0.0 220e9; 100.0 180e9]
        support2 = SupportElement(15, [2, 3], 45.0, 5, 0.12, 0.0015, E_T_support2, 12.0, [1, 3])
        
        supports = [support1, support2]
        temperatures = [10.0, 25.0, 40.0, 55.0, 70.0]
        
        sim_opts_orig = SimulationOptions(
            bridge,
            supports,
            temperatures,
            damping_ratio=0.035
        )
        
        # Convert to dictionary
        sim_dict = simulation_options_to_dict(sim_opts_orig)
        
        # Test top-level structure
        @test haskey(sim_dict, "bridge")
        @test haskey(sim_dict, "supports")
        @test haskey(sim_dict, "temperatures")
        @test haskey(sim_dict, "damping_ratio")
        @test haskey(sim_dict, "total_dofs")
        @test haskey(sim_dict, "total_elements")
        @test haskey(sim_dict, "created_at")
        
        # Test value preservation
        @test sim_dict["temperatures"] ≈ temperatures rtol=1e-12
        @test sim_dict["damping_ratio"] ≈ 0.035 rtol=1e-12
        @test sim_dict["total_dofs"] == sim_opts_orig.total_dofs  # Use actual calculated value
        @test sim_dict["total_elements"] == sim_opts_orig.total_elements  # Use actual calculated value
        @test haskey(sim_dict, "created_at")  # Check that timestamp exists (auto-generated)
        
        # Test nested structure preservation
        @test length(sim_dict["supports"]) == 2
        @test isa(sim_dict["bridge"], Dict)
        @test isa(sim_dict["supports"][1], Dict)
        @test isa(sim_dict["supports"][2], Dict)
        
        # Test JSON serialization compatibility
        json_string = JSON.json(sim_dict, 2)  # Pretty print format
        @test isa(json_string, String)
        @test length(json_string) > 100  # Non-trivial content
        
        # Parse back from JSON
        parsed_dict = JSON.parse(json_string)
        sim_opts_reconstructed = SimulationOptions(
                                    dict_to_bridge_options(parsed_dict["bridge"]),
                                    [dict_to_support_element(s) for s in parsed_dict["supports"]],
                                    Vector{Float64}(parsed_dict["temperatures"]),
                                    damping_ratio=parsed_dict["damping_ratio"]
                                )
        
        # Test complete reconstruction accuracy
        @test sim_opts_reconstructed.bridge.n_elem == sim_opts_orig.bridge.n_elem
        @test sim_opts_reconstructed.bridge.L ≈ sim_opts_orig.bridge.L rtol=1e-12
        @test length(sim_opts_reconstructed.supports) == length(sim_opts_orig.supports)
        @test sim_opts_reconstructed.temperatures ≈ sim_opts_orig.temperatures rtol=1e-12
        @test sim_opts_reconstructed.damping_ratio ≈ sim_opts_orig.damping_ratio rtol=1e-12
        @test sim_opts_reconstructed.total_dofs == sim_opts_orig.total_dofs
        @test sim_opts_reconstructed.total_elements == sim_opts_orig.total_elements
        # Check that both have timestamp format (they will differ slightly due to reconstruction time)
        @test length(sim_opts_reconstructed.created_at) > 15  # Reasonable timestamp length
        @test length(sim_opts_orig.created_at) > 15  # Reasonable timestamp length
    end
    
    @testset "File I/O Operations" begin
        # Create test configuration
        E_T = [0.0 200e9; 100.0 160e9]
        bc_nodes = BridgeBC([[1, "all"]])
        bridge = BridgeOptions(10, bc_nodes, 20.0, 2500.0, 0.1, 0.001, E_T, 50.0)
        
        E_T_support = [0.0 200e9; 100.0 180e9]
        support = SupportElement(5, [1, 2], 0.0, 2, 0.05, 0.0005, E_T_support, 5.0, [1, 2, 3])
        
        sim_opts_orig = SimulationOptions(
            bridge,
            [support],
            [20.0, 40.0, 60.0],
            damping_ratio=0.02
        )
        
        # Test save operation
        test_filename = "test_config_roundtrip.json"
        
        # Clean up any existing test file
        isfile(test_filename) && rm(test_filename)
        
        # Save configuration
        save_simulation_options(sim_opts_orig, test_filename)
        
        # Verify file was created
        @test isfile(test_filename)
        
        # Test file content is valid JSON
        file_content = read(test_filename, String)
        @test length(file_content) > 50  # Non-trivial content
        
        parsed_content = JSON.parse(file_content)
        @test isa(parsed_content, Dict)
        @test haskey(parsed_content, "bridge")
        @test haskey(parsed_content, "supports")
        @test haskey(parsed_content, "temperatures")
        
        # Load configuration
        sim_opts_loaded = load_simulation_options(test_filename)
        
        # Test complete round-trip accuracy
        @test sim_opts_loaded.bridge.n_elem == sim_opts_orig.bridge.n_elem
        @test sim_opts_loaded.bridge.L ≈ sim_opts_orig.bridge.L rtol=1e-12
        @test sim_opts_loaded.bridge.ρ ≈ sim_opts_orig.bridge.ρ rtol=1e-12
        @test sim_opts_loaded.bridge.A ≈ sim_opts_orig.bridge.A rtol=1e-12
        @test sim_opts_loaded.bridge.I ≈ sim_opts_orig.bridge.I rtol=1e-12
        @test sim_opts_loaded.bridge.cutoff_freq ≈ sim_opts_orig.bridge.cutoff_freq rtol=1e-12
        
        # Test E_T matrix preservation
        @test sim_opts_loaded.bridge.E_T ≈ sim_opts_orig.bridge.E_T rtol=1e-12
        
        # Test support element preservation
        @test length(sim_opts_loaded.supports) == length(sim_opts_orig.supports)
        @test sim_opts_loaded.supports[1].connection_node == sim_opts_orig.supports[1].connection_node
        @test sim_opts_loaded.supports[1].connection_dofs == sim_opts_orig.supports[1].connection_dofs
        @test sim_opts_loaded.supports[1].angle ≈ sim_opts_orig.supports[1].angle rtol=1e-12
        @test sim_opts_loaded.supports[1].E_T ≈ sim_opts_orig.supports[1].E_T rtol=1e-12
        
        # Test temperature-dependent function behavior
        test_temps = [0.0, 25.0, 50.0, 75.0, 100.0]
        for T in test_temps
            @test sim_opts_loaded.bridge.E(T) ≈ sim_opts_orig.bridge.E(T) rtol=1e-10
            @test sim_opts_loaded.supports[1].E(T) ≈ sim_opts_orig.supports[1].E(T) rtol=1e-10
        end
        
        # Test other parameters
        @test sim_opts_loaded.temperatures ≈ sim_opts_orig.temperatures rtol=1e-12
        @test sim_opts_loaded.damping_ratio ≈ sim_opts_orig.damping_ratio rtol=1e-12
        @test sim_opts_loaded.total_dofs == sim_opts_orig.total_dofs
        @test sim_opts_loaded.total_elements == sim_opts_orig.total_elements
        @test sim_opts_loaded.created_at == sim_opts_orig.created_at
        
        # Clean up test file
        rm(test_filename)
    end
    
    @testset "Edge Cases and Error Handling" begin
        # Test with minimal configuration
        E_T_minimal = [0.0 200e9; 100.0 200e9]  # Constant E
        bc_minimal = BridgeBC([[1, "all"]])
        bridge_minimal = BridgeOptions(1, bc_minimal, 1.0, 1000.0, 0.01, 0.0001, E_T_minimal, 10.0)
        
        sim_opts_minimal = SimulationOptions(bridge_minimal, SupportElement[], [20.0], damping_ratio=0.01)
        
        # Test minimal configuration round-trip
        dict_minimal = simulation_options_to_dict(sim_opts_minimal)
        
        @test dict_minimal["bridge"]["n_elem"] == 1
        @test dict_minimal["bridge"]["L"] ≈ 1.0 rtol=1e-12
        @test length(dict_minimal["supports"]) == 0
        @test length(dict_minimal["temperatures"]) == 1
        @test dict_minimal["temperatures"][1] ≈ 20.0 rtol=1e-12
        
        # Test with extreme values
        E_T_extreme = [0.0 1e15; 1000.0 1e5]  # Very large/small E values
        bc_extreme = BridgeBC([[1, "all"], [100, "trans"]])  # Large node numbers (using "trans" instead of "pin")
        bridge_extreme = BridgeOptions(100, bc_extreme, 1000.0, 10000.0, 1.0, 1.0, E_T_extreme, 1000.0)
        
        dict_extreme = bridge_options_to_dict(bridge_extreme)
        bridge_reconstructed = dict_to_bridge_options(dict_extreme)
        
        @test bridge_reconstructed.n_elem == 100
        @test bridge_reconstructed.E_T[1, 2] ≈ 1e15 rtol=1e-10
        @test bridge_reconstructed.E_T[2, 2] ≈ 1e5 rtol=1e-10
        @test bridge_reconstructed.bc_nodes.conds[2][1] == 100
        
        # Test precision preservation with small values
        E_T_small = [0.0 1e-6; 100.0 1e-9]
        bridge_small = BridgeOptions(1, bc_minimal, 1e-3, 1e-3, 1e-9, 1e-12, E_T_small, 1e-6)
        
        dict_small = bridge_options_to_dict(bridge_small)
        bridge_small_reconstructed = dict_to_bridge_options(dict_small)
        
        @test bridge_small_reconstructed.L ≈ 1e-3 rtol=1e-12
        @test bridge_small_reconstructed.ρ ≈ 1e-3 rtol=1e-12
        @test bridge_small_reconstructed.A ≈ 1e-9 rtol=1e-12
        @test bridge_small_reconstructed.I ≈ 1e-12 rtol=1e-12
        @test bridge_small_reconstructed.E_T[1, 2] ≈ 1e-6 rtol=1e-12
        @test bridge_small_reconstructed.E_T[2, 2] ≈ 1e-9 rtol=1e-12
    end
    
    @testset "JSON Format Compatibility" begin
        # Test that generated JSON matches expected format structure
        E_T = [0.0 200e9; 100.0 160e9]
        bc_nodes = BridgeBC([[1, "all"]])
        bridge = BridgeOptions(5, bc_nodes, 10.0, 2500.0, 0.1, 0.001, E_T, 50.0)
        
        sim_opts = SimulationOptions(bridge, SupportElement[], [20.0, 40.0], damping_ratio=0.02)
        
        dict_output = simulation_options_to_dict(sim_opts)
        json_output = JSON.json(dict_output, 2)
        
        # Test that JSON contains expected structure
        @test occursin("\"bridge\"", json_output)
        @test occursin("\"supports\"", json_output)
        @test occursin("\"temperatures\"", json_output)
        @test occursin("\"n_elem\"", json_output)
        @test occursin("\"bc_nodes\"", json_output)
        @test occursin("\"E_T\"", json_output)
        
        # Test that numerical values are properly formatted
        @test occursin("2.0e11", json_output)  # 200e9 in scientific notation
        @test occursin("1.6e11", json_output)  # 160e9 in scientific notation
        @test occursin("2500.0", json_output)
        @test occursin("0.1", json_output)
        @test occursin("0.001", json_output)
        
        # Test array structures
        @test occursin("[", json_output)
        @test occursin("]", json_output)
        @test count("[", json_output) >= 3  # Multiple arrays present
    end
    
end  # @testset "IO Module Tests"

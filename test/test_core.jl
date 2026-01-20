using Test
using BridgeFEM
using JSON
using Interpolations

@testset "Core Module Tests" begin
    
    @testset "BCTypes Constants" begin
        @test BCTypes["all"] == [1, 2, 3]
        @test BCTypes["trans"] == [1, 2]
        @test BCTypes["x"] == 1
        @test BCTypes["y"] == 2
        @test BCTypes["ϕ"] == 3
    end
    
    @testset "BridgeBC Type" begin
        # Test with string constraint types
        bc1 = BridgeBC([[1, "all"], [5, "trans"]])
        @test bc1.conds[1] == [1, [1, 2, 3]]
        @test bc1.conds[2] == [5, [1, 2]]
        
        # Test with numeric constraint types
        bc2 = BridgeBC([[1, [1, 2, 3]], [5, [1, 2]]])
        @test bc2.conds[1] == [1, [1, 2, 3]]
        @test bc2.conds[2] == [5, [1, 2]]
        
        # Test mixed constraint types
        bc3 = BridgeBC([[1, "all"], [2, [1]], [3, "y"]])
        @test bc3.conds[1] == [1, [1, 2, 3]]
        @test bc3.conds[2] == [2, [1]]
        @test bc3.conds[3] == [3, 2]
    end
    
    @testset "BridgeOptions Type" begin
        # Test temperature-dependent constructor
        E_T = [0.0 200e9; 50.0 180e9; 100.0 160e9]
        bc = BridgeBC([[1, "all"]])
        bridge = BridgeOptions(10, bc, 30.0, 2500.0, 0.1, 0.001, E_T, 50.0)
        
        @test bridge.n_elem == 10
        @test bridge.n_nodes == 11
        @test bridge.n_dofs == 33  # 3 DOFs per node
        @test bridge.L == 30.0
        @test bridge.ρ == 2500.0
        @test bridge.A == 0.1
        @test bridge.I == 0.001
        @test bridge.cutoff_freq == 50.0
        @test size(bridge.E_T) == (3, 2)
        
        # Test Young's modulus interpolation function
        @test bridge.E(0.0) ≈ 200e9 atol=1e-6
        @test bridge.E(50.0) ≈ 180e9 atol=1e-6
        @test bridge.E(100.0) ≈ 160e9 atol=1e-6
        @test bridge.E(25.0) ≈ 190e9 atol=1e-6  # Interpolated value
    end
    
    @testset "SupportElement Type" begin
        # Test temperature-dependent constructor
        E_T = [0.0 200e9; 100.0 160e9]
        support1 = SupportElement(5, [1, 2], 45.0, 3, 0.05, 0.0005, E_T, 10.0, [1, 2, 3])
        
        @test support1.connection_node == 5
        @test support1.connection_dofs == [1, 2]
        @test support1.angle == 45.0
        @test support1.n_elem == 3
        @test support1.A == 0.05
        @test support1.I == 0.0005
        @test support1.L == 10.0
        @test support1.bc_bottom == [1, 2, 3]
        @test size(support1.E_T) == (2, 2)
        
        # Test Young's modulus interpolation
        @test support1.E(0.0) ≈ 200e9 atol=1e-6
        @test support1.E(100.0) ≈ 160e9 atol=1e-6
        @test support1.E(50.0) ≈ 180e9 atol=1e-6  # Interpolated value
        
        # Test constant Young's modulus constructor
        support2 = SupportElement(3, [1], 0.0, 2, 150e9, 0.02, 0.0001, 5.0, [1, 2])
        
        @test support2.connection_node == 3
        @test support2.connection_dofs == [1]
        @test support2.A == 0.02
        @test support2.I == 0.0001
        @test support2.L == 5.0
        @test support2.bc_bottom == [1, 2]
        
        # Constant E should be same at all temperatures
        @test support2.E(-100.0) ≈ 150e9 atol=1e-6
        @test support2.E(0.0) ≈ 150e9 atol=1e-6
        @test support2.E(100.0) ≈ 150e9 atol=1e-6
    end
    
    @testset "SimulationOptions Type" begin
        # Create components
        E_T = [0.0 200e9; 100.0 160e9]
        bc = BridgeBC([[1, "all"]])
        bridge = BridgeOptions(5, bc, 20.0, 2500.0, 0.1, 0.001, E_T, 100.0)
        support = SupportElement(3, [1, 2], 0.0, 2, 0.05, 0.0005, E_T, 5.0, [1, 2, 3])
        temperatures = [20.0, 30.0, 40.0]
        
        # Test constructor with supports
        sim_opts1 = SimulationOptions(bridge, [support], temperatures, damping_ratio=0.05)
        
        @test sim_opts1.bridge === bridge
        @test length(sim_opts1.supports) == 1
        @test sim_opts1.supports[1] === support
        @test sim_opts1.temperatures == temperatures
        @test sim_opts1.damping_ratio == 0.05
        @test sim_opts1.total_elements == bridge.n_elem + support.n_elem
        @test !isempty(sim_opts1.created_at)
        
        # Test constructor without supports
        sim_opts2 = SimulationOptions(bridge, temperatures)
        
        @test sim_opts2.bridge === bridge
        @test length(sim_opts2.supports) == 0
        @test sim_opts2.temperatures == temperatures
        @test sim_opts2.damping_ratio == 0.02  # Default value
        @test sim_opts2.total_elements == bridge.n_elem
    end
    
    @testset "JSON Serialization" begin
        # Create test objects
        E_T = [0.0 200e9; 50.0 180e9; 100.0 160e9]
        bc = BridgeBC([[1, "all"], [5, "trans"]])
        bridge = BridgeOptions(8, bc, 25.0, 2400.0, 0.08, 0.0008, E_T, 75.0)
        support = SupportElement(4, [1, 2], 30.0, 3, 0.04, 0.0004, E_T, 8.0, [1, 2, 3])
        sim_opts = SimulationOptions(bridge, [support], [15.0, 25.0, 35.0], damping_ratio=0.03)
        
        @testset "BridgeOptions Serialization" begin
            # Test round-trip serialization
            dict = bridge_options_to_dict(bridge)
            bridge_restored = dict_to_bridge_options(dict)
            
            @test bridge_restored.n_elem == bridge.n_elem
            @test bridge_restored.n_nodes == bridge.n_nodes
            @test bridge_restored.n_dofs == bridge.n_dofs
            @test bridge_restored.L == bridge.L
            @test bridge_restored.ρ == bridge.ρ
            @test bridge_restored.A == bridge.A
            @test bridge_restored.I == bridge.I
            @test bridge_restored.cutoff_freq == bridge.cutoff_freq
            @test bridge_restored.E_T ≈ bridge.E_T
            
            # Test that boundary conditions are preserved
            @test length(bridge_restored.bc_nodes.conds) == length(bridge.bc_nodes.conds)
            @test bridge_restored.bc_nodes.conds[1] == bridge.bc_nodes.conds[1]
            @test bridge_restored.bc_nodes.conds[2] == bridge.bc_nodes.conds[2]
            
            # Test Young's modulus interpolation function works
            @test bridge_restored.E(25.0) ≈ bridge.E(25.0) atol=1e-6
        end
        
        @testset "SupportElement Serialization" begin
            # Test round-trip serialization
            dict = support_element_to_dict(support)
            support_restored = dict_to_support_element(dict)
            
            @test support_restored.connection_node == support.connection_node
            @test support_restored.connection_dofs == support.connection_dofs
            @test support_restored.angle == support.angle
            @test support_restored.n_elem == support.n_elem
            @test support_restored.A == support.A
            @test support_restored.I == support.I
            @test support_restored.L == support.L
            @test support_restored.bc_bottom == support.bc_bottom
            @test support_restored.E_T ≈ support.E_T
            
            # Test Young's modulus interpolation function works
            @test support_restored.E(25.0) ≈ support.E(25.0) atol=1e-6
        end
        
        @testset "SimulationOptions Serialization" begin
            # Test serialization to dict
            dict = simulation_options_to_dict(sim_opts)
            
            @test haskey(dict, "bridge")
            @test haskey(dict, "supports")
            @test haskey(dict, "temperatures")
            @test haskey(dict, "damping_ratio")
            @test haskey(dict, "total_dofs")
            @test haskey(dict, "total_elements")
            @test haskey(dict, "created_at")
            
            @test dict["temperatures"] == sim_opts.temperatures
            @test dict["damping_ratio"] == sim_opts.damping_ratio
            @test dict["total_dofs"] == sim_opts.total_dofs
            @test dict["total_elements"] == sim_opts.total_elements
            
            # Test that nested objects are properly serialized
            @test typeof(dict["bridge"]) == Dict{String, Any}
            @test typeof(dict["supports"]) == Vector{Dict{String, Any}}
            @test length(dict["supports"]) == 1
        end
        
        @testset "File I/O" begin
            # Test saving and loading to file
            filename = tempname() * ".json"
            try
                save_simulation_options(sim_opts, filename)
                @test isfile(filename)
                
                sim_opts_loaded = load_simulation_options(filename)
                
                # Verify loaded data matches original
                @test sim_opts_loaded.bridge.n_elem == sim_opts.bridge.n_elem
                @test sim_opts_loaded.bridge.L == sim_opts.bridge.L
                @test length(sim_opts_loaded.supports) == length(sim_opts.supports)
                @test sim_opts_loaded.temperatures == sim_opts.temperatures
                @test sim_opts_loaded.damping_ratio == sim_opts.damping_ratio
                
                # Test that interpolation functions work in loaded data
                @test sim_opts_loaded.bridge.E(25.0) ≈ sim_opts.bridge.E(25.0) atol=1e-6
                @test sim_opts_loaded.supports[1].E(25.0) ≈ sim_opts.supports[1].E(25.0) atol=1e-6
                
            finally
                # Clean up temporary file
                isfile(filename) && rm(filename)
            end
        end
    end
end 
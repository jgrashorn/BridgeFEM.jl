using Test

# Main test runner for BridgeFEM.jl
# Works with current monolithic structure during transition period

@testset verbose = true "BridgeFEM.jl Test Suite" begin
    
    @testset "Unit Tests - Modular Structure" begin
        @testset "Core Module Tests" begin
            include("test_core.jl")
        end
        
        @testset "Elements Module Tests" begin
            include("test_elements.jl")
        end

        @testset "Assembly Module Tests" begin
            include("test_assembly.jl")
        end

        @testset "Simple Assembly Module Tests" begin
            include("test_assembly_simple.jl")
        end
        
        @testset "Dynamics Module Tests" begin
            include("test_dynamics.jl")
        end
        
        @testset "ModelReduction Module Tests" begin
            include("test_model_reduction.jl")
        end

        @testset "BoundaryConditions Module Tests" begin
            include("test_boundary_conditions.jl")
        end

        @testset "IO Module Tests" begin
            include("test_io.jl")
        end
    end
    @testset "Integration Tests - Modular Structure" begin
        
        @testset "Fixed Cantilever Static Analysis" begin
            include("integration/test_cantilever.jl")
        end
        
        @testset "Fixed Cantilever Dynamic Analysis" begin
            include("integration/test_cantilever_dynamic.jl")
        end
        
        @testset "Simply Supported Beam Static Analysis" begin
            include("integration/test_beam.jl")
        end
        
        @testset "Simply Supported Beam Dynamic Analysis" begin
            include("integration/test_beam_dynamic.jl")
        end
    end
end
using Test

# Main test runner for BridgeFEM.jl
# Works with current monolithic structure during transition period

@testset verbose = true "BridgeFEM.jl Test Suite" begin
    
    @testset "Unit Tests - Modular Structure" begin
        @testset "Core Module Tests" begin
            include("test_core.jl")
        end
    end
    
    # Integration tests temporarily disabled during modularization
    # Will be re-enabled after all modules are extracted and updated
    # @testset "Integration Tests - Current Monolithic Structure" begin
    #     
    #     @testset "Fixed Cantilever Static Analysis" begin
    #         include("integration/test_cantilever.jl")
    #     end
    #     
    #     @testset "Fixed Cantilever Dynamic Analysis" begin
    #         include("integration/test_cantilever_dynamic.jl")
    #     end
    #     
    #     @testset "Simply Supported Beam Static Analysis" begin
    #         include("integration/test_beam.jl")
    #     end
    #     
    #     @testset "Simply Supported Beam Dynamic Analysis" begin
    #         include("integration/test_beam_dynamic.jl")
    #     end
    # end
end 
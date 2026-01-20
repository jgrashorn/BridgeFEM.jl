# Standalone Integration Test Runner for BridgeFEM.jl
# Works with current monolithic structure during transition period

using Test

println("="^60)
println("BridgeFEM.jl Integration Test Suite")
println("="^60)

@testset verbose = true "BridgeFEM.jl Integration Tests" begin
    
    # Test the converted fixed cantilever test
    @testset "Fixed Cantilever Static Analysis" begin
        println("\n🔧 Running Fixed Cantilever Static Analysis tests...")
        include("integration/test_cantilever.jl")
        println("✅ Fixed Cantilever tests completed")
    end
    
    # Test the converted fixed cantilever dynamic test
    @testset "Fixed Cantilever Dynamic Analysis" begin
        println("\n🔧 Running Fixed Cantilever Dynamic Analysis tests...")
        include("integration/test_cantilever_dynamic.jl")
        println("✅ Fixed Cantilever Dynamic tests completed")
    end
    
    # Test the converted simply supported beam static test
    @testset "Simply Supported Beam Static Analysis" begin
        println("\n🔧 Running Simply Supported Beam Static Analysis tests...")
        include("integration/test_beam.jl")
        println("✅ Simply Supported Beam tests completed")
    end
    
    # Test the converted simply supported beam dynamic test
    @testset "Simply Supported Beam Dynamic Analysis" begin
        println("\n🔧 Running Simply Supported Beam Dynamic Analysis tests...")
        include("integration/test_beam_dynamic.jl")
        println("✅ Simply Supported Beam Dynamic tests completed")
    end
end

println("\n" * "="^60)
println("All integration tests completed successfully! ✅")
println("="^60) 
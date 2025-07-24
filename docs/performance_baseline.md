# Performance Baseline Documentation

## Overview

This document describes the performance baseline establishment process for BridgeFEM.jl before code reorganization. The baseline provides quantitative metrics to validate that the modularization maintains acceptable performance levels with no more than 5% degradation.

## Baseline Script

**File**: `tests/performance_baseline.jl`  
**Purpose**: Establish baseline performance metrics for key computational operations  
**Target**: Validate ≤5% performance degradation after reorganization

## Baseline Establishment Process

### 1. Model Configuration

The baseline uses a medium-sized bridge model for representative benchmarking:

- **Elements**: 50 finite elements
- **Nodes**: 51 nodes (3 DOFs each = 153 total DOFs)
- **Beam Length**: 25.0 meters
- **Material**: Steel (ρ = 7800 kg/m³, E₀ = 207 GPa)
- **Cross-section**: A = 0.15 m², I = 0.005 m⁴
- **Boundary Conditions**: Simply supported (fixed translations at first node, vertical support at last node)
- **Temperature Range**: -50°C to +100°C with temperature-dependent Young's modulus

This configuration represents a typical medium-complexity structural analysis problem suitable for research workflows.

### 2. Benchmarked Operations

#### Matrix Assembly (Task 2)
- **`assemble_matrices()`**: Complete mass and stiffness matrix assembly
- **`assemble_stiffness!()`**: Stiffness matrix assembly with pre-allocated sparse matrices
- **`frame_elem_stiffness()`**: Element-level stiffness matrix computation
- **Sparse Matrix Analysis**: Memory usage patterns and sparsity characteristics

#### Eigenvalue Computation (Task 3)
- **LinearAlgebra.eigen()**: Standard generalized eigenvalue computation
- **Arpack.eigs()**: Partial eigenvalue computation for modal analysis
- **decompose_matrices()**: Complete modal decomposition workflow across temperatures
- **Eigenvector Processing**: Mass normalization and mode tracking

#### ODE Solving (Task 4)
- **Modal Space Integration**: Using `beam_modal_ode!()` with reduced DOFs
- **Physical Space Integration**: Using `beam_physical_ode!()` with full DOFs
- **Tolerance Scaling**: Performance with different solver tolerances
- **Time Span Scaling**: 1s, 5s, and 10s simulation durations

### 3. System Information Capture

Each baseline includes:
- Julia version
- System architecture and CPU thread count
- Package versions (LinearAlgebra, SparseArrays, DifferentialEquations, Arpack)
- Timestamp and baseline metadata
- Model configuration parameters

## Baseline Metrics Schema

The baseline metrics are stored in `tests/baseline_metrics.json` with the following structure:

```json
{
  "baseline_info": {
    "julia_version": "1.x.x",
    "timestamp": "2024-01-XX...",
    "system_info": {
      "arch": "x86_64-linux-gnu",
      "cpu_threads": 8,
      "word_size": 64
    },
    "package_versions": { ... }
  },
  "model_config": {
    "n_elements": 50,
    "n_nodes": 51,
    "n_dofs": 153,
    "beam_length": 25.0,
    "description": "Medium-sized bridge model..."
  },
  "benchmarks": {
    "matrix_assembly": {
      "description": "Matrix Assembly (assemble_matrices)",
      "mean_time_s": 0.xxx,
      "std_time_s": 0.xxx,
      "min_time_s": 0.xxx,
      "max_time_s": 0.xxx,
      "mean_memory_bytes": xxx,
      "runs": 5,
      "individual_times_s": [...]
    },
    // ... additional benchmarks
  },
  "metadata": {
    "baseline_version": "1.0",
    "purpose": "Performance baseline before code reorganization",
    "validation_target": "≤5% performance degradation requirement"
  }
}
```

## Usage Instructions

### Establishing Initial Baseline

```bash
cd tests
julia performance_baseline.jl
```

This creates `baseline_metrics.json` with comprehensive performance metrics.

### Validation After Reorganization

1. **Re-run baseline script** to generate current metrics:
   ```bash
   julia performance_baseline.jl  # Creates updated baseline_metrics.json
   ```

2. **Load and compare metrics** in Julia:
   ```julia
   include("performance_baseline.jl")
   
   # Load baseline and current metrics
   baseline = load_baseline_metrics("baseline_metrics.json")
   current = load_baseline_metrics("current_metrics.json")
   
   # Compare with 5% tolerance
   results = compare_performance_metrics(current, baseline, 5.0)
   
   # Validate overall performance
   passed = validate_performance_against_baseline()
   ```

3. **Automated validation**:
   ```julia
   # Simple pass/fail check
   performance_ok = validate_performance_against_baseline()
   if !performance_ok
       error("Performance degradation exceeds 5% tolerance")
   end
   ```

## Performance Validation Criteria

### Acceptance Criteria
- **≤5% degradation**: All benchmarked operations must maintain performance within 5% of baseline
- **No functional regression**: All operations must produce mathematically equivalent results
- **Memory efficiency**: Sparse matrix operations should maintain similar memory patterns

### Benchmark Categories
1. **Matrix Assembly**: Critical for FEM setup performance
2. **Eigenvalue Computation**: Essential for modal analysis workflows
3. **ODE Solving**: Core dynamic simulation performance
4. **Overall Workflow**: End-to-end simulation timings

### Expected Performance Characteristics

Based on the medium-sized bridge model (153 DOFs):

- **Matrix Assembly**: ~1-10 ms for typical sparse matrix operations
- **Eigenvalue Computation**: ~10-100 ms for 10-20 modes
- **Modal ODE Integration**: ~50-500 ms for 5-second simulations
- **Physical ODE Integration**: ~200-2000 ms for 5-second simulations

Actual performance depends on system specifications and Julia compilation state.

## Repeatability Considerations

### Ensuring Consistent Results

1. **Warmup Runs**: Script includes 2 warmup runs to ensure Julia compilation
2. **Multiple Runs**: Each benchmark performs 5 runs and reports statistics
3. **Garbage Collection**: Forced GC between runs to minimize memory effects
4. **Fixed Model**: Consistent model parameters across all runs
5. **Random Seed**: Consider setting seeds for any stochastic operations

### Factors Affecting Performance

- **System Load**: Run on dedicated system when possible
- **Julia Version**: Performance may vary between Julia versions
- **Package Updates**: Dependency updates may affect performance
- **Hardware**: CPU architecture and memory speed impact results
- **Operating System**: OS-specific optimizations may vary

### Recommendations

- Run baseline on the same system used for development
- Minimize system load during benchmarking
- Use the same Julia version for baseline and validation
- Document any significant system or environment changes

## Integration with Development Workflow

### Before Reorganization
1. Run baseline script to establish metrics
2. Commit `baseline_metrics.json` to version control
3. Document system configuration and Julia version

### During Reorganization
1. Periodically validate performance against baseline
2. Investigate any significant degradations immediately
3. Update metrics if fundamental algorithmic improvements are made

### After Reorganization
1. Run final validation against original baseline
2. Document any intentional changes in performance characteristics
3. Establish new baseline for future development if needed

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Ensure all packages are available
2. **Memory Issues**: Large matrices may require adequate RAM
3. **Compilation Time**: First run may be slower due to compilation
4. **Numerical Precision**: Small variations in results are expected

### Performance Investigation

If performance degradation is detected:

1. **Identify Specific Benchmarks**: Check which operations are affected
2. **Profile Code**: Use Julia's profiling tools for detailed analysis
3. **Compare Algorithms**: Verify that reorganization preserved algorithms
4. **Memory Analysis**: Check for memory leaks or excessive allocations

## File Locations

- **Baseline Script**: `tests/performance_baseline.jl`
- **Metrics Storage**: `tests/baseline_metrics.json`
- **Comparison Results**: `tests/performance_comparison.json` (generated during validation)
- **Documentation**: `docs/performance_baseline.md` (this file)

## Related Documentation

- **Theory**: `docs/theory.md` - Mathematical background for benchmarked operations
- **Examples**: `docs/examples.md` - Usage examples including performance-critical code
- **API Reference**: `docs/api-reference.md` - Function signatures for benchmarked operations 
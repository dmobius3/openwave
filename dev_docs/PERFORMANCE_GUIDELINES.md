# Performance Guidelines

## Overview

OpenWave simulates subatomic physics at the Planck scale, requiring careful attention to performance optimization. This guide provides best practices for writing high-performance code in the project. It contains general performance principles, Taichi optimization, GPU best practices, and benchmarking guidelines.

## General Performance Principles

### Profile Before Optimizing

- Always measure performance before and after optimization
- Use profiling tools to identify actual bottlenecks
- Focus optimization efforts on hot paths

### Algorithm Selection

- Choose appropriate algorithms for the scale of the problem
- Consider time vs. space complexity tradeoffs
- Use established scientific computing libraries when appropriate

## Python Performance

### NumPy Best Practices

- Vectorize operations instead of using Python loops
- Use appropriate dtypes (float32 vs float64)
- Avoid unnecessary array copies
- Use views and slices when possible

### Memory Management

- Pre-allocate arrays when size is known
- Reuse buffers in iterative calculations
- Be mindful of memory layout (C vs Fortran order)

### Lattice Step Loops: Write Into Fixed Buffers

For research scripts that step a NumPy field on a lattice (a time-step loop, a relaxation, a reproduction run of a shipped package). Every `np.roll` call and every operator in an expression such as `-6.0 * u + np.roll(u, 1, 0)` allocates a new array the size of the whole field and frees it a moment later. On a large 3D grid the loop then spends its time moving memory rather than doing arithmetic, and in a pool every worker competes for the same memory bandwidth.

| Practice | Reason | Check |
| --- | --- | --- |
| Allocate the step's work arrays once, before the loop, and write every result into them with `out=` (`np.multiply(u, dt, out=tmp)`, `np.add(v, tmp, out=v)`) | Removes the per-step allocate-and-free of full-size temporaries | Run the original and the rewrite from the same initial state and compare the fields to a stated tolerance (`1e-12` in float64), not exact equality: the order of additions can change |
| Write a periodic shift as slice adds into the preallocated buffer, not `np.roll` | `np.roll` copies the whole field on every call | The same agreement check, on a field whose size differs along each axis so a wrong axis shows up |

`np.roll` stays the right choice when the copy does not matter: a one-off shift, a small field, or code read far more often than it runs. Rewrite only the loop that runs per step on a large field.

The demonstration is [`utils/demo_preallocated_step.py`](utils/demo_preallocated_step.py): a 3D field stepped both ways from the same random start, the agreement asserted to `1e-12`, and both timings printed for the machine at hand. Its `lap_into` is the slice-add Laplacian to copy from.

Write the fast version as new code next to the original and keep the original as the reference. A script shipped by a model author, or pinned by hash, is reproduced as shipped and never edited for speed.

### Pools: Size by What Bounds the Job

| Job shape | What bounds it | Rule |
| --- | --- | --- |
| Many short jobs, each building a large field | Memory bandwidth: throughput stops rising at a worker count well below the core count | Measure the saturation point and stop there; more workers only add contention |
| A few long single-threaded jobs, each solving one problem | CPU: throughput scales with the core count until memory runs out | Up to one job per core, checking peak memory per job against the machine's RAM |
| One large job that is mostly BLAS / LAPACK (big matrix products, eigensolvers) | Its own threads | Run it alone and let the math library use its threads; do not pin it to 1 |

In every pool, pin each worker's math-library threads to 1 before NumPy is imported: `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and `VECLIB_MAXIMUM_THREADS=1`, set in the environment of the launch command or at the top of the worker before `import numpy`. Otherwise N workers can each start a thread per core and fight over the CPU.

Find the saturation point with a fixed-work batch: the same number of identical jobs at a few worker counts, comparing jobs per minute.

**Compare fixed-work runs only.** A timing compares across machines, settings or code versions only when every run does the same work: a fixed number of steps on a fixed grid. A run that iterates to convergence does a different amount of work each time, so its wall time measures the stopping rule as much as the speed. Keep the numbers with the machine that produced them: the rules above are portable, the timings are not.

## Taichi Optimization

### Kernel Design

- Minimize kernel launches by batching operations
- Use Taichi's parallel for-loops effectively
- Avoid excessive atomic operations
- Structure data for coalesced memory access

### Data Layout

```python
# Good: Structure of Arrays (SoA)
@ti.data_oriented
class Particles:
    def __init__(self, n):
        self.x = ti.field(dtype=ti.f32, shape=n)
        self.y = ti.field(dtype=ti.f32, shape=n)
        self.z = ti.field(dtype=ti.f32, shape=n)

# Consider: Array of Structures (AoS) when accessing all properties together
@ti.data_oriented
class Particles:
    def __init__(self, n):
        self.pos = ti.Vector.field(3, dtype=ti.f32, shape=n)
```

### GPU Optimization

- Ensure sufficient parallelism for GPU utilization
- Minimize CPU-GPU data transfers
- Use appropriate block sizes for GPU kernels
- Avoid divergent branching in GPU code

### Real-World Optimization Example: Spacetime 2D Slider

The Spacetime 2D slider visualization was optimized from a CPU-bound implementation to a GPU-accelerated version with the following techniques:

#### Original Performance Issues

- Serial CPU loops for rendering (O(n²) without parallelization)
- Object recreation on parameter changes (memory allocation overhead)
- Repeated calculations in render loop
- No caching of invariant values

#### Applied Optimizations

1. **GPU-Parallelized Position Calculations**

```python
@ti.kernel
def compute_screen_position(self, offset: ti.f32, size: ti.f32):
    for i, j in ti.ndrange(self.count, self.count):
        universe_pos = self.grid[i, j]
        self.screen_pos[i, j][0] = (universe_pos[0] + offset) / size
        self.screen_pos[i, j][1] = (universe_pos[1] + offset) / size
```

1. **Pre-allocated Memory with Update Methods**

```python
class Lattice:
    def __init__(self, scale_factor):
        self.max_count = 1000  # Pre-allocate max size
        self.grid = ti.Vector.field(2, dtype=ti.f32, shape=(self.max_count, self.max_count))
        
    def update_scale(self, scale_factor):
        # Update existing structure instead of recreating
        self.spacing = 2 * constants.PLANCK_LENGTH * scale_factor * np.e
        self.count = min(int(self.size / self.spacing), self.max_count)
```

1. **Cached Calculations Outside Render Loop**

```python
# Pre-compute constants that don't change during rendering
universe_to_screen_ratio = min(config.SCREEN_RES) / lattice.size
screen_radius = max(granule.radius * universe_to_screen_ratio, 1)
offset = (lattice.size - lattice.spacing * (lattice.count - 1)) / 2

# Only recalculate when scale changes
if scale.value != previous_scale:
    # Update cached values
```

1. **Batch Operations for Data Transfer**

```python
# Get all position at once instead of accessing individually
screen_position = lattice.get_screen_position_numpy()
```

#### Performance Impact

- Reduced position calculation from O(n²) CPU operations to parallel GPU kernel
- Eliminated memory allocation overhead during parameter changes
- Minimized CPU-GPU data transfer to once per frame
- Achieved smooth real-time interaction even with thousands of granules

## Parallel Computing

### Threading Considerations

- Use thread-safe data structures when necessary
- Minimize lock contention
- Consider using thread-local storage for temporary data
- For process pools of NumPy jobs (one thread per worker, sizing by job shape): [§ Pools: Size by What Bounds the Job](#pools-size-by-what-bounds-the-job)

### Vectorization

- Write code that compiler can auto-vectorize
- Use SIMD-friendly data layouts
- Align data structures appropriately

## Caching Strategies

### Computation Caching

- Cache expensive calculations when inputs don't change
- Use memoization for recursive functions
- Implement lazy evaluation where appropriate

### Data Caching

- Keep frequently accessed data in cache-friendly layouts
- Use spatial and temporal locality principles
- Consider cache line sizes (typically 64 bytes)

## Performance Benchmarking

### Benchmark Guidelines

- Create reproducible benchmarks
- Test with realistic data sizes
- Measure both average and worst-case performance
- Document hardware specifications for benchmarks

### Performance Targets

- Granule simulation: >1M particles at 30 FPS
- Wave propagation: Real-time for 2D, near real-time for 3D
- Memory usage: Scale linearly with particle count

## Common Pitfalls to Avoid

### Python-Specific

- Avoid creating objects in inner loops
- Don't use global variables in performance-critical code
- Minimize string operations in numerical code
- Avoid repeated attribute lookups in loops

### Numerical Computing

- Check for numerical stability in iterative methods
- Use appropriate precision for calculations
- Avoid unnecessary type conversions
- Be careful with division by small numbers

## Performance Monitoring

### Metrics to Track

- Frame rate / simulation steps per second
- Memory usage and allocation patterns
- GPU utilization (if applicable)
- Cache hit rates

### Profiling Tools

- Python: cProfile, line_profiler, memory_profiler
- Taichi: Built-in profiler (`ti.profiler`)
- System: perf, Intel VTune, NVIDIA Nsight

## Optimization Workflow

1. **Establish baseline**: Measure current performance
2. **Profile**: Identify bottlenecks
3. **Optimize**: Apply targeted optimizations
4. **Verify**: Ensure correctness is maintained
5. **Measure**: Quantify improvement
6. **Document**: Record optimization rationale and results

## References

- [Taichi Performance Tuning Guide](https://docs.taichi-lang.org/docs/performance)
- [NumPy Performance Tips](https://numpy.org/doc/stable/user/c-info.beyond-basics.html)
- [Python Performance Tips](https://wiki.python.org/moin/PythonSpeed/PerformanceTips)

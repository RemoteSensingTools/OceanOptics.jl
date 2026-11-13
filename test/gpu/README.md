# CUDA device tests

`test/test_device.jl` runs in the default suite on the JLArrays reference backend. This
environment runs the same file on a CUDA device; it is skipped when `CUDA.functional()` is false.

```bash
julia --project=test/gpu -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'   # once
julia --project=test/gpu test/gpu/runtests.jl
```

CUDA uses the first visible device; on a host with several, select one with the environment
variable `CUDA_VISIBLE_DEVICES`. The arrays of the tests are small (tens of elements).

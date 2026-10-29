## Blast Wave Problem

This case uses the shared `MHD3D` class and ideal solver from `../common/`.

- `mhd3d_init_.cpp` defines the initially stationary, high-pressure sphere: `r_0` sets its radius, `ro0/ro1` and `pr0/pr1` set the ambient/interior states, and `b_0`, `bthe`, and `bphi` set the uniform magnetic field.
- `mhd3d_paras.cpp` defines the global domain and boundary flags through `setup_grid()`; the supplied configuration is periodic in all three directions.
- `mymacros.hpp` selects the grid, MPI decomposition, numerical methods, and CFL number.

See [the MPI guide](../README.md) for building, running, restarting, and merging/reading results. When changing the pressure contrast, magnetic-field strength, or solver settings, check the resulting density and pressure for nonphysical values.

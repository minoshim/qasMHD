## Magnetic reconnection

This case uses the shared `DMHD3D` class and dissipative solver from `../common/`, enabled by `DISSIPATION = 1` in its Makefile.

- `dmhd3d_init_.cpp` defines the Harris sheet, velocity/magnetic perturbations, and `setdc()` for viscosity and resistivity profiles.
- `mhd3d_paras.cpp` defines the global domain and boundary flags through `setup_grid()`.
- `mymacros.hpp` selects the grid, MPI decomposition, numerical methods, `RANDOM`, and the Reynolds numbers `REV` and `REM`.

`RANDOM=1` adds a random `vy` perturbation whose X/Z-dependent sequence is shared across Y ranks; `RANDOM=0` disables only this velocity perturbation. In both modes, the magnetic perturbation amplitude `b1amp` is randomized around `b1` with half-width `b1`, using a separate Z-only sequence shared across X/Y ranks. Both generators use `rand_noise_mt()` and `std::mt19937`.

For a reproducible initial state, replace `std::time(nullptr)` in `dmhd3d_init_.cpp` with a fixed unsigned seed and retain the same mesh and MPI decomposition.

The diffusion substep safety coefficient is 9, the 3D counterpart of the 2D
coefficient 6, accounting for the additional direction and fourth-order
differences. The existing ideal/dissipative splitting order is unchanged.
Viscosity and resistivity must be finite and nonnegative; invalid coefficients
or a failed MPI reduction terminate all ranks.

See [the MPI guide](../README.md) for building, running, restarting, and merging/reading results.

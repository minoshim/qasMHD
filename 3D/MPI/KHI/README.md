## Kelvin-Helmholtz instability

This case uses the shared `MHD3D` class and ideal solver from `../common/`.

- `mhd3d_init_.cpp` defines the shear layer, density, magnetic-field orientation, and perturbation amplitudes.
- `mhd3d_paras.cpp` defines the global domain and boundary flags through `setup_grid()`.
- `mymacros.hpp` selects the grid, MPI decomposition, numerical methods, and `RANDOM` flag.

With `RANDOM=1`, the `vy` perturbation uses an X/Z-dependent Mersenne-Twister sequence shared across Y ranks; with `RANDOM=0`, it uses a sinusoid in X. The shear amplitude in `vx` remains randomly perturbed in both modes, using a separate Z-only sequence shared across X/Y ranks. Spatial shear and localization profiles are applied to these amplitudes.

For a reproducible initial state, replace `std::time(nullptr)` in `mhd3d_init_.cpp` with a fixed unsigned seed and retain the same mesh and MPI decomposition. See [the MPI guide](../README.md) for building, running, restarting, and merging/reading results.

## Orszag-Tang Vortex Problem

This case uses the shared `MHD3D` class and ideal solver from `../common/`.

- `mhd3d_init_.cpp` defines the initial fields. Its local `NUM` macro selects the X-Y, Y-Z, or Z-X plane (`0`–`2`), the transposed versions (`3`–`5`), or the 3D configuration (`6`).
- `mhd3d_paras.cpp` defines the periodic global domain and explicitly passes zero coordinate shifts to `setup_grid()`.
- `mymacros.hpp` selects the grid, MPI decomposition, and numerical methods. The default `NUM=0` configuration is independent of Z; choose suitable resolution along each active direction when changing `NUM`.

See [the MPI guide](../README.md) for building, running, restarting, and merging/reading results. The simulation stores CT magnetic fields; the merger provides both cell-center fields and separate CT face files for diagnostics.

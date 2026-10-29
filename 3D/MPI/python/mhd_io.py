"""Readers for CT face fields written by the 3D MPI merge utility."""

from pathlib import Path
import numpy as np


def read_ct_fields(directory, index, nx, ny, nz):
    """Return Bx(nz,ny,nx+1), By(nz,ny+1,nx), Bz(nz+1,ny,nx).

    Upper boundary faces are included; tangential ghost cells are not.
    """
    fields = []
    for component, shape in (("bx", (nz, ny, nx+1)),
                             ("by", (nz, ny+1, nx)),
                             ("bz", (nz+1, ny, nx))):
        path = Path(directory)/f"merge_{component}_face_{index:05d}.dat"
        try:
            raw = path.read_bytes()
        except FileNotFoundError as error:
            raise FileNotFoundError(
                f"{path}: missing CT field; rerun the current merge.out on the raw data"
            ) from error
        expected = shape[0]*shape[1]*shape[2]*np.dtype(np.float32).itemsize
        if len(raw) != expected:
            raise ValueError(f"{path}: expected {expected} bytes, got {len(raw)}")
        fields.append(np.frombuffer(raw, dtype=np.float32).reshape(shape))
    return tuple(fields)

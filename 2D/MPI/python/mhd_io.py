"""Readers for CT face fields written by the 2D MPI merge utility."""

from pathlib import Path
import numpy as np


def read_ct_fields(directory, index, nx, ny):
    """Return Bx(ny,nx+1), By(ny+1,nx), including right/top boundary faces."""
    fields = []
    for component, shape in (("bx", (ny, nx+1)), ("by", (ny+1, nx))):
        path = Path(directory)/f"merge_{component}_face_{index:05d}.dat"
        try:
            raw = path.read_bytes()
        except FileNotFoundError as error:
            raise FileNotFoundError(
                f"{path}: missing CT field; rerun the current merge.out on the raw data"
            ) from error
        expected = shape[0]*shape[1]*np.dtype(np.float32).itemsize
        if len(raw) != expected:
            raise ValueError(f"{path}: expected {expected} bytes, got {len(raw)}")
        fields.append(np.frombuffer(raw, dtype=np.float32).reshape(shape))
    return tuple(fields)

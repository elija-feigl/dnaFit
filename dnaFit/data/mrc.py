#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" collection of scripts to allow manipulation of a cryo-EM map."""
import logging
from pathlib import Path
from typing import List

import MDAnalysis as mda
import mrcfile
import numpy as np
import numpy.typing as npt

logger = logging.getLogger(__name__)


def recenter_mrc(
    path: Path, to_position=np.array([0.0, 0.0, 0.0]), apply=True
) -> npt.NDArray[np.float64]:
    """recenter the mrc file to a specific position. returns required shift"""
    with mrcfile.open(path, mode="r+") as mrc:
        _cell = np.array(mrc.header["cella"])
        half_box: npt.NDArray[np.float64] = 0.5 * np.array([_cell["x"], _cell["y"], _cell["z"]])
        to_origin = to_position - half_box

        _origin = np.array(mrc.header["origin"])
        origin = np.array([_origin["x"], _origin["y"], _origin["z"]])

        shift: npt.NDArray[np.float64] = to_origin - origin
        if apply:
            mrc.header["origin"] = tuple(shift)
        return -shift


def get_mrc_center(path: Path) -> npt.NDArray[np.float64]:
    """return position of the center of an mrc file"""
    with mrcfile.open(path, mode="r+") as mrc:
        _cell = np.array(mrc.header["cella"])
        half_box: npt.NDArray[np.float64] = 0.5 * np.array([_cell["x"], _cell["y"], _cell["z"]])

        _origin = np.array(mrc.header["origin"])
        origin = np.array([_origin["x"], _origin["y"], _origin["z"]])

        return origin + half_box


def get_mrc_box(path: Path) -> List[float]:
    """returns box dimensions of mrc file.
    todo-low: include box angles
    """
    with mrcfile.open(path, mode="r+") as mrc:
        _cell = np.array(mrc.header["cella"])
        box: npt.NDArray[np.float64] = np.array([_cell["x"], _cell["y"], _cell["z"]])
        box_dimensions: List[float] = box.tolist()
        return box_dimensions


def write_mrc_from_atoms(
    path: Path,
    atoms: mda.AtomGroup,
    path_out: Path,
    context: float = 4.0,
    cut_box=True,
    keep_data=False,
):
    """mask a mrc map using an group of atoms"""
    if not atoms:
        logger.warning("Cannot crop empty atom selection. No file written")
        return

    with mrcfile.open(path) as mrc:
        voxel, grid, origin, full_data = _get_mrc_properties(mrc)

    data_mask = _create_voxel_mask(atoms, grid, origin, voxel, context)
    data = full_data * data_mask
    if cut_box:
        data, origin, voxel = _mrc_cutbox(data, full_data, origin, voxel, keep_full=keep_data)

    with mrcfile.new(path_out, overwrite=True) as mrc_out:
        mrc_out.set_data(data.transpose())
        mrc_out.voxel_size = tuple(voxel.tolist())
        mrc_out.header["origin"] = tuple(origin)


def _create_voxel_mask(atoms, grid, origin, voxel_size, context):
    data_mask = np.zeros(grid, dtype=np.float32)
    v_context: npt.NDArray[np.float64] = np.full(3, context / voxel_size).astype(int) + 1
    v_context_x, v_context_y, v_context_z = v_context
    grid_positions = np.rint((atoms.positions - origin) / voxel_size).astype(int)
    for pos in grid_positions:
        x_pos, y_pos, z_pos = pos[0], pos[1], pos[2]  # fast to slow axis
        x_low, x_high = x_pos - v_context_x, x_pos + v_context_x
        y_low, y_high = y_pos - v_context_y, y_pos + v_context_y
        z_low, z_high = z_pos - v_context_z, z_pos + v_context_z
        data_mask[x_low:x_high, y_low:y_high, z_low:z_high] = 1.0
    return data_mask


def _get_mrc_properties(mrc):
    _origin = np.array(mrc.header["origin"])
    origin = np.array([_origin["x"], _origin["y"], _origin["z"]])
    _cell = np.array(mrc.header["cella"])
    box = np.array([_cell["x"], _cell["y"], _cell["z"]])
    grid = np.array([mrc.header["nx"], mrc.header["ny"], mrc.header["nz"]])
    data = mrc.data.transpose()
    voxel_size = box / grid

    return (voxel_size, grid, origin, data)


def _mrc_cutbox(data, full_data, m_origin, m_spacing, keep_full=False):
    idx_data = np.nonzero(data)
    if keep_full:
        data = full_data

    x_min, x_max = np.min(idx_data[0]), np.max(idx_data[0])
    y_min, y_max = np.min(idx_data[1]), np.max(idx_data[1])
    z_min, z_max = np.min(idx_data[2]), np.max(idx_data[2])

    xyz_diff = max(x_max - x_min, y_max - y_min, z_max - z_min)
    x_pad = 0.5 * (xyz_diff + x_min - x_max)
    y_pad = 0.5 * (xyz_diff + y_min - y_max)
    z_pad = 0.5 * (xyz_diff + z_min - z_max)
    x_low = int(x_pad) if (x_pad % 1.0) == 0.0 else int(x_pad) + 1
    y_low = int(y_pad) if (y_pad % 1.0) == 0.0 else int(y_pad) + 1
    z_low = int(z_pad) if (z_pad % 1.0) == 0.0 else int(z_pad) + 1

    data_small: npt.NDArray[np.float64] = np.zeros((xyz_diff, xyz_diff, xyz_diff), dtype=np.float32)
    x_high = -int(x_pad) or None
    y_high = -int(y_pad) or None
    z_high = -int(z_pad) or None
    data_small[
        x_low:x_high,
        y_low:y_high,
        z_low:z_high,
    ] = data[x_min:x_max, y_min:y_max, z_min:z_max]

    # compute new origin
    origin = m_origin + (
        (x_min - x_low) * m_spacing[0],
        (y_min - y_low) * m_spacing[1],
        (z_min - z_low) * m_spacing[2],
    )
    grid = np.shape(data_small)
    cell = grid * m_spacing
    spacing = cell / grid
    return data_small, origin, spacing

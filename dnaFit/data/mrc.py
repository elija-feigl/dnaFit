import logging
import mrcfile
import numpy as np
import MDAnalysis as mda

from pathlib import Path


""" DESCR:
    collection of scripts to allow manipulation of a cryo-EM map.
"""
logger = logging.getLogger(__name__)


def recenter_mrc(path: Path, to_position=np.array([0.0, 0.0, 0.0])) -> np.ndarray:
    with mrcfile.open(path, mode='r+') as mrc:
        c = np.array(mrc.header["cella"])
        half_box = 0.5 * np.array([c["x"], c["y"], c["z"]])
        to_origin = (to_position - half_box)

        h = np.array(mrc.header["origin"])
        origin = np.array([h["x"], h["y"], h["z"]])

        shift = to_origin - origin
        mrc.header["origin"] = tuple(shift)
        return -shift


def write_mrc_from_atoms(path: Path, atoms: mda.AtomGroup,
                         path_out: Path, context: float = 4., cut_box=True):
    if not len(atoms):
        logger.warning(
            f"Cannot crop with empty atom selection. No file written")
        return

    with mrcfile.open(path) as mrc:
        voxel, grid, origin, full_data = _get_mrc_properties(mrc)

    data_mask = _create_voxel_mask(atoms, grid, origin, voxel, context)
    data = full_data * data_mask
    if cut_box:
        data, origin, voxel = _mrc_cutbox(data, origin, voxel)

    with mrcfile.new(path_out, overwrite=True) as mrc_out:
        mrc_out.set_data(data.transpose())
        mrc_out._set_voxel_size(*voxel)
        mrc_out.header["origin"] = tuple(origin)


def _create_voxel_mask(atoms, grid, origin, voxel_size, context):
    data_mask = np.zeros(grid, dtype=np.float32)
    v_context = np.full(3, context / voxel_size).astype(int) + 1
    v_context_x, v_context_y, v_context_z = v_context
    grid_positions = np.rint(
        ((atoms.positions - origin) / voxel_size)).astype(int)
    for pos in grid_positions:
        x, y, z = pos[0], pos[1], pos[2]  # fast to slow axis
        data_mask[x - v_context_x: x + v_context_x,
                  y - v_context_y: y + v_context_y,
                  z - v_context_z: z + v_context_z
                  ] = 1.
    return data_mask


def _get_mrc_properties(mrc):
    o = np.array(mrc.header["origin"])
    origin = np.array([o["x"], o["y"], o["z"]])
    c = np.array(mrc.header["cella"])
    cell = np.array([c["x"], c["y"], c["z"]])
    grid = np.array(
        [mrc.header["nx"], mrc.header["ny"], mrc.header["nz"]])
    data = mrc.data.transpose()
    voxel_size = cell / grid

    return (voxel_size, grid, origin, data)


def _mrc_cutbox(data, m_origin, m_spacing):
    idx_data = np.nonzero(data)

    x_min, x_max = np.min(idx_data[0]), np.max(idx_data[0])
    y_min, y_max = np.min(idx_data[1]), np.max(idx_data[1])
    z_min, z_max = np.min(idx_data[2]), np.max(idx_data[2])

    xyz_diff = max(x_max - x_min, y_max - y_min, z_max - z_min)
    x_pad = 0.5 * (xyz_diff + x_min - x_max)
    y_pad = 0.5 * (xyz_diff + y_min - y_max)
    z_pad = 0.5 * (xyz_diff + z_min - z_max)
    x_low = int(x_pad) if (x_pad % 1.) == 0. else int(x_pad) + 1
    y_low = int(y_pad) if (y_pad % 1.) == 0. else int(y_pad) + 1
    z_low = int(z_pad) if (z_pad % 1.) == 0. else int(z_pad) + 1

    data_small = np.zeros(
        (xyz_diff, xyz_diff, xyz_diff), dtype=np.float32)
    data_small[x_low: -int(x_pad) or None, y_low: -int(y_pad) or None,
               z_low: -int(z_pad) or None] = data[x_min: x_max,
                                                  y_min: y_max,
                                                  z_min: z_max]

    # cumpute new origin
    origin = (m_origin + ((x_min - x_low) * m_spacing[0],
                          (y_min - y_low) * m_spacing[1],
                          (z_min - z_low) * m_spacing[2]
                          )
              )
    grid = np.shape(data_small)
    cell = grid * m_spacing
    spacing = (cell / grid)
    return data_small, origin, spacing

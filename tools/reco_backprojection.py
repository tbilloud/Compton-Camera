# Basic backprojection reconstruction for Compton camera data
# It's slow, ~1 sec per cone -> install cupy if possible

from opengate.logger import global_log
from tools.display_reconstruction import *
import numpy as xp

try:
    import cupy as xp
except ImportError:
    global_log.warning(f"Cupy is not installed. Using numpy instead.")


def reco_bp(cones_df, vpitch, vsize, napari=False, det=False):
    if len(cones_df) > 1:  # avoid logging when used in point source validation
        global_log.info(f'Reconstructing volume with backprojection')

    volume = xp.zeros(vsize, dtype=xp.float32)
    grid_x = xp.linspace(-vsize[0] // 2, vsize[0] // 2, vsize[0]) * vpitch
    grid_y = xp.linspace(-vsize[1] // 2, vsize[1] // 2, vsize[1]) * vpitch
    grid_z = xp.linspace(-vsize[2] // 2, vsize[2] // 2, vsize[2]) * vpitch
    X, Y, Z = xp.meshgrid(grid_x, grid_y, grid_z, indexing='ij')

    for _, c in cones_df.iterrows():
        apex = xp.array([c['Apex_X'], c['Apex_Y'], c['Apex_Z']])
        d = xp.array([c['Direction_X'], c['Direction_Y'], c['Direction_Z']])
        cosT = c['cosT']

        # Compute distance from apex to each voxel
        voxel_vec = xp.stack([X - apex[0], Y - apex[1], Z - apex[2]], axis=-1)
        voxel_distances = xp.linalg.norm(voxel_vec, axis=-1)

        # Compute angle with direction vector
        vox_vec_norm = voxel_vec / xp.expand_dims(voxel_distances, axis=-1)
        dot_products = xp.sum(vox_vec_norm * d, axis=-1)

        # Compute mask of voxels satisfying the Compton cone condition
        tolerance = 0.01  # Adjust tolerance as needed
        cone_mask = xp.abs(dot_products - cosT) < tolerance

        # Accumulate contribution to the volume
        volume[cone_mask] += 1

    volume = xp.swapaxes(volume, 0, 1)
    if napari:
        if xp.__name__ == 'cupy': volume = volume.get()
        display_reconstruction(volume, vsize, vpitch, det)

    return volume

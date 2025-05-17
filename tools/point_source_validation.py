from opengate.logger import global_log
from tools.analysis_cones import *
from tools.reco_backprojection import *
from pathlib import Path
import numpy as xp
import matplotlib.pyplot as plt

try:
    import cupy as xp
except ImportError:
    global_log.warning(f"Cupy is not installed. Using numpy instead.")


# Script to check the precision of Gate9 or Gate10 simulations with a point source
# Can be run as a standalone script or (WIP) as a function in a Gate10 script

# All cones should intersect at the source point
# Possible reasons for bad cones:
# - rayleigh scattering
# - compton scattering with electron not at rest (doppler broadening)
# - particle-induced X-ray emission (fluorescence, Auger)
# - more than 2 coincident events (nSingles)
# - electron/gamma escape
# - time resolution (pile-up, singles with different eventID, true_coinc)
# - energy/spatial resolution

def valid_psource(cones_df, src_pos, vpitch, vsize, plot_seq=False,
                  plot_stk=False):
    stime = time.time()
    global_log.info(f'Offline [source validation]: START')
    if not len(cones_df):
        global_log.error(f"Empty input (no cones in dataframe).")
        global_log.info(f"Offline [source validation]: {get_stop_string(stime)}")
        return
    else:
        global_log.debug(f"Offline [cones tpx]: Input cone dataframe with ({len(cones_df)} entries)")


    # Source position must be in units of voxels in vol
    sp_vox = [int(src_pos[i] / vpitch) + (vsize[i] // 2) for i in range(3)]

    # ######## RECONSTRUCT CONE BY CONE #######################################
    z_slice_stack = xp.zeros((len(cones_df), vsize[0], vsize[1]), dtype=xp.float32)
    nb = 0  # number of bad cones
    cones_df = cones_df.reset_index(drop=True)
    for idx, cone in cones_df.iterrows():
        vol = reco_bp(cone.to_frame().T, vpitch, vsize)
        z_slice = vol[:, :, sp_vox[2]]
        z_slice_stack[idx, :, :] = z_slice
        if z_slice[sp_vox[0], sp_vox[1]] == 0:
            nb += 1  # TODO: sometime cone is bad but z_slice is not 0 (due to error)

        # ##############################################################
        # # Display stack with matplotlib (one by one)
        # ##############################################################
        if plot_seq:
            if xp.__name__ == 'cupy': z_slice = xp.asnumpy(z_slice)
            plt.imshow(z_slice, cmap='gray', origin='lower')
            plt.scatter(sp_vox[0], sp_vox[1], c='r', s=10)
            plt.title(f'EventID: {int(cone["EventID"])}')
            add_secondary_axes(plt.gca(), vpitch)
            plt.colorbar()
            plt.tight_layout()
            plt.show()

    global_log.debug(f"Offline [source validation]: {nb} cones not intersecting source")

    # ##############################################################
    # # Display stack with matplotlib (summed)
    # ##############################################################
    if plot_stk:
        fig, ax = plt.subplots()
        stack = xp.sum(z_slice_stack, axis=0)
        if xp.__name__ == 'cupy': stack = xp.asnumpy(stack)
        ax.imshow(stack, cmap='gray_r', origin='lower')
        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
        add_secondary_axes(ax, vpitch)
        plt.scatter(sp_vox[0], sp_vox[1], c='r', s=10)
        plt.tight_layout()
        plt.show()

    global_log.info(f"Offline [source validation]: {get_stop_string(stime)}")

    return z_slice_stack, vsize

def add_secondary_axes(ax, vpitch):
    Xmm = ax.secondary_xaxis('top')
    Xmm.set_xlabel('X (mm)', color='red')
    Xmm.set_xticks(ax.get_xticks())
    Xmm.set_xticklabels(xp.round(ax.get_xticks() * vpitch, 2), color='red')
    Ymm = ax.secondary_yaxis('right', color='red')
    Ymm.set_ylabel('Y (mm)', color='red')
    Ymm.set_yticks(ax.get_yticks())
    Ymm.set_yticklabels(xp.round(ax.get_yticks() * vpitch, 2), color='red')


if __name__ == "__main__":
    # ###### READING Gate9.2 sequenceCoincidence.root files ##############
    # fname, E0_MeV, vpitch, source_pos = Path('../Gate9/output'), 0.140, 200, [0, 0, 0]
    # nSingles_max = False  # maximum number of singles per coincidence, False to disable
    # true_coinc = False  # filter true coincidences (i.e. avoid singles from different events)
    # nentries = None  # None to read all entries
    # cones_array = seqCoin2ConesArray(fname / 'CC_sequenceCoincidence.root', E0_MeV, vsize, vpitch, er, nSingles_max, true_coinc, nentries)

    # ###### READING Gate9.2 CC_Cones.root files #########################
    # fname, E0_MeV, vpitch, source_pos = Path('../Gate9/output'), 0.140, 200, [0, 0, 0]
    # nSingles_max = False  # maximum number of singles per coincidence, False to disable
    # true_coinc = False  # filter true coincidences (i.e. avoid singles from different events)
    # nentries = None  # None to read all entries
    # cones_array = conesTTree2conesArray(fname / 'CC_Cones.root', E0_MeV, er, nSingles_max, true_coinc, nentries)

    # ###### READING Gate10 hit root files ##############
    fname, E0_MeV, vpitch, source_pos = Path('../output'), 1.0, 200, [0, 0, -50]
    cones_array = gHits2cones_byEvtID(fname / 'CC_Hits.root', E0_MeV)

    # ###### Preprocessing #########
    # print('number of cones with a nan value:', cp.isnan(cones_array).any(axis=1).sum())
    # print(len(cones_array), 'cones before removing nans')
    # cones_array = remove_nans(cones_array)
    # print(len(cones_array), 'cones after removing nans')

    valid_psource(cones_array, vpitch, source_pos)

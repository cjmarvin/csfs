import numpy as np
import sys

vis_file_min = "DATA/fsr_waves_min_vis.txt"
vis_file_max = "DATA/fsr_waves_max_vis.txt"
vis_file_blaze = "DATA/blaze_waves_vis.txt"

nir_file_min = "DATA/fsr_waves_min_nir.txt"
nir_file_max = "DATA/fsr_waves_max_nir.txt"
nir_file_blaze = "DATA/blaze_waves_nir.txt"

try:
    arg = sys.argv[1]
except IndexError:
    arg = False

def nir_wtlist(write=False):
    nir_min = np.loadtxt(nir_file_min)
    nir_max = np.loadtxt(nir_file_max)
    nir_blaze = np.loadtxt(nir_file_blaze)
    waves = np.concatenate([nir_min, nir_max, nir_blaze])
    waves = np.sort(waves) * 1.e7
    print waves
    if write:
        np.savetxt("nir_sim_wtlist.txt", waves)

def vis_wtlist(write=False):
    vis_min = np.loadtxt(vis_file_min)
    vis_max = np.loadtxt(vis_file_max)
    vis_blaze = np.loadtxt(vis_file_blaze)
    waves = np.concatenate([vis_min, vis_max, vis_blaze])
    waves = np.sort(waves) * 1.e7
    print waves
    if write:
        np.savetxt("vis_sim_wtlist.txt", waves)

if __name__ == "__main__":
    nir_wtlist(write=arg)
    vis_wtlist(write=arg)
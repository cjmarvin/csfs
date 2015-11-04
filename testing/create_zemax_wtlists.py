import numpy as np
import sys

nir_file = "zemax_nir_wavelengths.txt"
vis_file = "zemax_vis_wavelengths.txt"

try:
    arg = sys.argv[1]
except IndexError:
    arg = False

def nir_wtlist(write=False):
    m, li, lc, lf = np.loadtxt(nir_file, unpack=True, skiprows=2)
    waves = np.concatenate([li, lc, lf])
    waves = np.sort(waves) * 1.e4
    print waves
    if write:
        np.savetxt("nir_wtlist.txt", waves)

def vis_wtlist(write=False):
    m, li, lc, lf = np.loadtxt(vis_file, unpack=True, skiprows=2)
    waves = np.concatenate([li, lc, lf])
    waves = np.sort(waves) * 1.e4
    print waves
    if write:
        np.savetxt("vis_wtlist.txt", waves)

if __name__ == "__main__":
    nir_wtlist(write=arg)
    vis_wtlist(write=arg)
"""
simgeneral_c.py

Prepares data to be fed into simulation.

"""
import numpy as np
#import matplotlib.pyplot as plt
try:
    import pyfits
except ImportError:
    import astropy.io.fits as pyfits
import time
from physics import n_sell, energy2counts
from printer import Printer
from printtiming import print_timing
#from slitfunctions import *
from slitfunctions import slit_image, point_source
import cfunctions
import wavefuncs


@print_timing
def sim_general(arm, fiber, wavelengths=None, intensities=None,
    fiber_description=None, wavemap=None):
    """
    General simulation function.

    Parameters
    ----------
    arm : SpectralArm object
        Container object of spectrograph and simulation parameters.
    fiber : int {0,1}
        Fiber to simulate. 1 = A and 2 = B.

    """
    if wavelengths is None: wavelengths = arm.wavelengths
    if intensities is None: intensities = arm.intensities
    #if fiber_description is None: fiber_description = arm.fiber_description
    if wavemap is None: wavemap = arm.wavemap

    print "Fiber %s - %s" % (arm.fib_char[fiber], arm.fiber_description)
    t1 = time.time()    # start time

    orders = arm.OSET
    norders = arm.NOSET

    if arm.fib_char[fiber] in ['A', 'B']: #@MZ
        offset = arm.fib_OFFSET[fiber]

    # output arrays
    image = arm.zero_images[fiber] # initialize image CCD array
    if arm.SAMPLING == "mc":
        image = np.array(image, dtype=np.uint)
    if arm.SAMPLING == "grid":
        wavemap_arr = np.zeros(arm.CCD_DIMS, dtype=np.uint)
    elif arm.SAMPLING == "mc":
        wavemap_arr = np.zeros(arm.CCD_DIMS)

        # m_list = np.zeros((0,0), dtype=np.uint)

    # CREATE SLIT ?
    if arm.SAMPLING == "grid":
        ns = arm.ns # input slit samples
        if arm.slit is '0':
            slitx, slity = point_source(arm, offset=offset)
        else:
            slitx, slity = slit_image(arm, offset=offset)
        nslit = slitx.size # effective slit samples
        perturb = arm.slitperturb
    elif arm.SAMPLING == "mc":
        ns = False
        slitx = False
        slity = False
        perturb = False
        if arm.slit[0] in '0 1 2'.split():
            arm.slittyp = 'uniform'
            SLIT_FLAG = 0
            arm.slit_locx = 0.0
            arm.slit_locy = 0.0
            arm.slit_scalex = 0.0
            arm.slit_scaley = 0.0
        elif arm.slit[0] == '3':
            arm.slittyp = 'gaussian'
            SLIT_FLAG = 1
            if len(arm.slit) == 1:
                arm.slit_locx = 0.0
                arm.slit_locy = 0.0
                arm.slit_scalex = 1.0
                arm.slit_scaley = 1.0
            elif len(arm.slit) == 2:
                arm.slit_locx = 0.0
                arm.slit_locy = 0.0
                arm.slit_scalex = float(arm.slit[1])
                arm.slit_scaley = float(arm.slit[1])
            elif len(arm.slit) == 3:
                arm.slit_locx = float(arm.slit[1])
                arm.slit_locy = float(arm.slit[1])
                arm.slit_scalex = float(arm.slit[2])
                arm.slit_scaley = float(arm.slit[2])
            elif len(arm.slit) == 4:
                arm.slit_locx = float(arm.slit[1])
                arm.slit_locy = float(arm.slit[2])
                arm.slit_scalex = float(arm.slit[3])
                arm.slit_scaley = float(arm.slit[3])
            elif len(arm.slit) == 5:
                arm.slit_locx = float(arm.slit[1])
                arm.slit_locy = float(arm.slit[2])
                arm.slit_scalex = float(arm.slit[3])
                arm.slit_scaley = float(arm.slit[4])

    # set simulation flags
    if arm.blaze:
        BLAZE_FLAG = 1
    else:
        BLAZE_FLAG = 0
    LOC_FLAG = 0

    # CDF interpolation flag
    for word in 'flatfield line'.split():
        if word in arm.fiber_description:
            interp_flag = 0 # do not interpolate (for line lists)
        else:
            interp_flag = 3 # interpolate
    if arm.interp_flag:
        interp_flag = arm.interp_flag

    if arm.ARM == "NIR" and arm.gap:
        print "Creating NIR detector gap of %s mm at x=%s mm" % (arm.gap, arm.gapoffset)

    # SOLUTION/COMPUTATION BLOCK
    # =========================================================================
    if arm.SAMPLING == "mc":
        print "Sampling method: Monte Carlo"
        print "Wavelength samping: %s" % arm.sm
        print "Fiber sampling distribution: %s" % arm.slittyp
        counts_tot = np.trapz(energy2counts(wavelengths, intensities), wavelengths) # total counts
        wave_range = 0.0 # wavelength coverage of each spectral order
        for i,m in enumerate(orders):
            waves, weights = wavefuncs.feed_wavelengths(arm, m)
            counts = energy2counts(waves, weights) # flux to photon counts
            wave_range = np.trapz(counts, waves) / counts_tot
            nwaves = waves.size
            nrays = int(arm.nr[fiber] * wave_range)
            output = '  * Calculating order %i (%i of %i), total %.2f%% complete.' % (m, i+1, norders, ((i+1.) * 100.0 / norders))
            if arm.sm == "cdf":
                nrays_org = nrays   # @MZ
                irays = 1e8     # @MZ
                nloop = np.int((nrays_org - 1) / irays) + 1    # @MZ
                for loop in xrange(nloop):  # @MZ
                    nrays = np.int(min(irays, nrays_org-loop*irays))   # @MZ
                    print 'Memory block loop {}/{}, nrays {}/{}'.format(loop+1, nloop, nrays, nrays_org)   # @MZ
                    waves_inp = np.zeros(nrays)
                    print "  * Sampling wavelengths using Cumulative Distribution Function (CDF)..."
                    cfunctions.random_wave_cdf(
                        waves,
                        counts,
                        nwaves,
                        interp_flag,
                        waves_inp,
                        nrays)
                    print output
                    print "  * wavelength coverage = %.2f%%" % (wave_range * 100.0)
                    print "  * nrays = %s" % nrays
                    cfunctions.compute_mc_cdf(
                        arm.ARM_FLAG,
                        BLAZE_FLAG,
                        arm.GAP_FLAG,
                        SLIT_FLAG,
                        nrays,
                        m,
                        arm.gap,
                        arm.gapoffset,
                        arm.XD_0,
                        arm.YD_0,
                        offset,
                        arm.slit_locx,
                        arm.slit_locy,
                        arm.slit_scalex,
                        arm.slit_scaley,
                        np.ascontiguousarray(waves_inp, dtype=np.float64),
                        image,
                        wavemap_arr)
                nrays = nrays_org
            elif arm.sm == "rej":
                print output
                print "  * wavelength coverage = %.2f%%" % (wave_range * 100.0)
                print "  * nrays = %s" % nrays
                cfunctions.compute_mc_rejsamp(
                    arm.ARM_FLAG,
                    BLAZE_FLAG,
                    arm.GAP_FLAG,
                    SLIT_FLAG,
                    nwaves,
                    nrays,
                    m,
                    arm.gap,
                    arm.gapoffset,
                    arm.XD_0,
                    arm.YD_0,
                    weights.max(),
                    offset,
                    arm.slit_locx,
                    arm.slit_locy,
                    arm.slit_scalex,
                    arm.slit_scaley,
                    np.ascontiguousarray(waves, dtype=np.float64),
                    np.ascontiguousarray(weights, dtype=np.float64),
                    image,
                    wavemap_arr)
            arm.fib_nrays[fiber].append(nrays)
        inds = np.where(image != 0.0)
        arm.fib_mean_rays_per_pixel[fiber] = np.mean(image[inds])
        arm.fib_min_rays_per_pixel[fiber] = np.min(image[inds])
        arm.fib_max_rays_per_pixel[fiber] = np.max(image[inds])
        arm.fib_nrays_tot[fiber] = np.sum(image)

    elif arm.SAMPLING == "grid":
        print "Sampling method: Grid"
        print "Fiber location perturbation: %s" % arm.perturb
        for i,m in enumerate(orders):
            output = '  * Calculating order %i (%i of %i), total %.2f%% complete.' % (m, i+1, norders, ((i+1.) * 100.0 / norders))
            print output
            if arm.perturb:
                waves, weights = wavefuncs.feed_wavelengths(arm, m)
                n_g_sell = n_sell(arm.ARM, waves)
                nwaves = waves.size
                cfunctions.compute_grid_perturb(
                    arm.ARM_FLAG,
                    BLAZE_FLAG,
                    arm.GAP_FLAG,
                    nwaves,
                    nslit,
                    m,
                    arm.gap,
                    arm.gapoffset,
                    arm.XD_0,
                    arm.YD_0,
                    arm.slitperturb,
                    np.ascontiguousarray(n_g_sell, dtype=np.float64),
                    np.ascontiguousarray(slitx, dtype=np.float64),
                    np.ascontiguousarray(slity, dtype=np.float64),
                    np.ascontiguousarray(waves, dtype=np.float64),
                    np.ascontiguousarray(weights, dtype=np.float64),
                    image,
                    wavemap_arr)
            else:
                waves, weights = wavefuncs.feed_wavelengths(arm, m)
                n_g_sell = n_sell(arm.ARM, waves)
                nwaves = waves.size
                returnx = np.empty(nwaves*nslit)
                returny = np.empty(nwaves*nslit)
                cfunctions.compute_grid(
                    arm.ARM_FLAG,
                    BLAZE_FLAG,
                    arm.GAP_FLAG,
                    nwaves,
                    nslit,
                    m,
                    arm.gap,
                    arm.gapoffset,
                    arm.XD_0,
                    arm.YD_0,
                    np.ascontiguousarray(n_g_sell, dtype=np.float64),
                    np.ascontiguousarray(slitx, dtype=np.float64),
                    np.ascontiguousarray(slity, dtype=np.float64),
                    np.ascontiguousarray(waves, dtype=np.float64),
                    np.ascontiguousarray(weights, dtype=np.float64),
                    image,
                    wavemap_arr)
        if wavemap:
            inds = np.where(wavemap_arr != 0)
            image[inds] *= 1.0e7 / wavemap_arr[inds] # average each pixel and convert mm to Angstrom
        else:
            image *= (arm.MAX_SNR**2) / image.max()
            print "Normalizing image array to a max of %s SNR." % arm.MAX_SNR

        arm.fib_nslit[fiber] = nslit
    arm.add_image(image)
    arm.fib_sim_time[fiber] = time.time() - t1
    return 0

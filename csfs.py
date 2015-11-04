#!/usr/bin/env python
"""
CARMENES Spectrum Forward Simulator
This program forward simulates the upcoming CARMENES spectrograph which
launches in 2014.
The simulation is purely geometrical and does not take into account any
luminosity effects (blaze efficiencies, aberrations, stray light, cosmic rays,
ghosting, etc).

Authors: Christopher J. Marvin, Mathias Zechmeister
Institution: Institute for Astrophysics, Goettingen
Advisors: Prof. Ansgar Reiners, Dr. Mathias Zechmeister

Tested for Python 2.7.2, numpy 1.6.1, matplotlib 1.1.0

Now under Git version control (BitBucket.org).

Previous History:

Version 3.1.1
  - header in wavelengthtrace.write_region changed to fix loading region error

Version 3.1.0
- Sellmeier equation implemented for wavelength-dependent index of refraction
  in cross-disperser of VIS arm
- updated VIS arm parameters
- Grism tilt implemented for VIS arm
- functions.c added functions that compute xy-coords of spectrograph and also
  bins, reducing overall simulation time
- vis.h added, which contains VIS arm parameters (any future changes must be
  done in both spectralarm.py AND vis.h, see improvements to be made)
- cfunctions.py added to wrap functions.c
- simgeneral_c.py added, equivalent to simgeneral.py but using functions
  written in c
  (simgeneral.py probably deprecated)
- slitfunctions.py now returns 1d numpy arrays, even for point-sources,
  which helps with the general functionality of functions.c
- spectral arm parameters also in header file vis.h (and must be updated
  here as well)
- spectralarm.py: find_blaze_offset() and find_ccd_wavelength_limits() now
  use cfunctions.py
- raytrace can now probably be deprecated
- wavelengthtrace.py utilizes cfunctions.compute_n_bin_wlt, and now loops over
  every spectral order to trace non-FSR wavelengths
- cleardata argument added in order to regenerate offset and wavelength solutions

  Notes:
- due to the n_sellmeier implementation, the NIR arm is not currently functional

Version 2.1.0
- NIR arm now functional
- simgeneral created and implemented
- simulationtest.py implemented with evaluated matrix multiplication equations,
  making the calculations much faster
- spectralarm created as a standalone module, and imports command line
  arguments, and most variables passed through spectralarm object
- SpectralArm.calculate_wavelengths() changed to span entire wavelength range
  of arm
- fabry perot intensities moved to spectralarm.py
- redshift function created in spectralarm.py
- wavetrace.py replaced by wavelengthtrace.py
- writefits now writes fits.gz as default
- detector offset and ccd limits are now saved to disk
- -oset command line argument added
- -solve command line argument added
- other minor fixes (variables generalized, etc.)

Version 2.0.0:
- changes by Mathias Zechmeister (MZ)
- fib A now down fib B up
- obsolete sim_phoenix, sim_spectrum, simlinelist
  cleaned and generalised to sim_fp
- halfmoon.py and pointsource.py merged to slitfunction.py

Version 1.1.0:
- program changed to incorporate command line arguments
- RV shift added
- functions to that can be added or edited (slit function and simulation
functions, etc.) are moved outside of main code for easier access and
maintenance


Improvements to be made:
 - n_sellmeier for infrasil (NIR arm)
 - better naming conventions
 - import parameters into spectralarm.py from vis.h (and nir.h)
 - change wta and wtb to wlt, enable importing a raw FITS image, and tracing
   the solutions? (basically, wavetrace as a standalone entity)
"""
__author__ = 'Christopher J. Marvin, Dr. Mathias Zechmeister'

# Python modules
#import matplotlib.pyplot as plt
import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import os
import subprocess
import time

# CSFS modules
from cmdline import parse_command_line
from physics import *
from printer import Printer
from printtiming import print_timing
from writefits import write_to_fits
from spectralarm import SpectralArm
from simgeneral_c import sim_general
import wavelengthtrace as wlt
import wavefuncs as wf
from noise import add_noise
from version import __version__

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                              MAIN PROGRAM
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@print_timing
def main():

    # COMMAND LINE ARGUMENTS
    args = parse_command_line()

    print "\nCARMENES Spectrum Forward Simulator (CSFS v%s)" % __version__

    # INITIALIZE SPECTRAL ARM
    arm = SpectralArm(args)

    # SIMULATION START TIMES
    t1 = time.time()
    sim_start = time.strftime("%a, %d, %b %Y %H:%M:%S + 0000", time.gmtime())

    if arm.fib_simmode[0] not in arm.SIM_LIST or arm.fib_simmode[1] not in arm.SIM_LIST:
        raise ValueError

    for i in xrange(len(arm.fib_simmode)):
        #fiber = i
        simmode = arm.fib_simmode[i]
        #simmodetype = arm.simmodetypes[i] # seems unnecessary

        if simmode is '0':
            arm.slittyp = "0D point-source"
            arm.add_image(np.zeros(arm.CCD_DIMS, dtype=np.uint16))

        elif simmode is 'B':
            arm.slittyp = "0D point-source"
            arm.fiber_description = 'bias frame'
            arm.fib_obstype = ['BIAS', 'BIAS']
            arm.fib_src = ['BIAS', 'BIAS']
            arm.add_image(np.zeros(arm.CCD_DIMS, dtype=np.uint16))

        else:
            print "Initializing Fiber %s" % arm.fib_char[i]
            wavemap = False

            if simmode is 'C':
                arm.fiber_description = 'Laser Comb Spectrum'
                fn = "MODELS/COMB/comb.npy"
                wavelengths, intensities = np.load(fn)
                inds = (wavelengths >= arm.wmin) & (wavelengths <= arm.wmax)
                wavelengths = wavelengths[inds]
                intensities = intensities[inds]
                if arm.SAMPLING == "grid":
                    wf.wavegrid(
                        arm,
                        wavelengths=wavelengths,
                        intensities=intensities,
                        assign=True,
                        telluric=False)
                elif arm.SAMPLING == "mc":
                    arm.wavelengths = wavelengths
                    arm.intensities = intensities
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])
                arm.infiles[i] = fn
                arm.fib_obstype[i] = 'WAVE'
                arm.fib_src[i] = 'COMB'

            elif simmode is 'F':
                arm.fiber_description = 'flatfield spectrum simulation'
                arm.wavelengths = wf.calculate_wavelengths(arm, mode='CCD', nwaves=arm.nw)
                arm.intensities = np.ones(arm.wavelengths.size)
                arm.fib_obstype[i] = 'FLAT'
                arm.fib_src[i] = 'HAL'

            elif simmode is 'L':
                if len(arm.infiles[i]) is not 1:
                    raise ValueError
                print 'loading %s' % (arm.infiles[i][0])
                arm.fiber_description = 'emission line list spectrum simulation'
                wavelengths, intensities = np.loadtxt(arm.infiles[i][0], unpack=True)
                #wavelengths, intensities = np.load(arm.infiles[i][0])
                wavelengths *= 1.0e-7    # convert Ang to mm @MZ
                inds = np.where((wavelengths >= arm.wmin) & (wavelengths <= arm.wmax))
                arm.wavelengths = wavelengths[inds]
                arm.intensities = intensities[inds]
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])

            elif simmode is 'P':
                arm.fiber_description = """PHOENIX model spectrum simulation: T_eff=3000 K, [Fe/H]=0.0, log(g)=5.0"""
                wavelengths = np.load('phx_wavelengths.npy')
                intensities = np.load('phx_intensities.npy')
                inds = (wavelengths >= arm.wmin) & (wavelengths <= arm.wmax)
                wavelengths = wavelengths[inds]
                intensities = intensities[inds]
                if arm.SAMPLING == 'grid':
                    wf.wavegrid(
                        arm,
                        wavelengths=wavelengths,
                        intensities=intensities,
                        assign=True,
                        telluric=arm.tell)
                elif arm.SAMPLING == 'mc':
                    arm.wavelengths = wavelengths
                    arm.intensities = intensities
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])
                if arm.tell:
                    arm.intensities = wf.convolve_telluric_lines(
                        arm,
                        arm.wavelengths,
                        arm.intensities)
                arm.infiles[i] = ['phx_wavelengths.npy', 'phx_intensities.npy']
                arm.catg = 'SCI'
                arm.fib_obstype[i] = 'STAR'
                arm.fib_src[i] = 'OBJ'
                arm.sciobject = 'PHOE_lte03000-5.00-0.0'

            elif simmode is 'S':
                if len(arm.infiles[i]) == 1:
                    try:
                        wavelengths, intensities = pyfits.getdata(arm.infiles[i][0])
                    except IOError:
                        wavelengths, intensities = np.loadtxt(arm.infiles[i][0], unpack=True)
                elif len(arm.infiles[i]) == 2:
                    try:
                        wavelengths = pyfits.getdata(arm.infiles[i][0])
                        intensities = pyfits.getdata(arm.infiles[i][1])
                    except IOError:
                        wavelengths = np.loadtxt(arm.infiles[i][0])
                        intensities = np.loadtxt(arm.infiles[i][1])
                elif len(arm.infiles[i]) not in [1, 2]:
                    raise AttributeError("Please specify an input wavelength file and input flux/counts file.")
                arm.fiber_description = 'object spectrum simulation'
                wavelengths *= 1.0e-7    # convert Ang to mm
                inds = (wavelengths >= arm.wmin) & (wavelengths <= arm.wmax)
                wavelengths = wavelengths[inds]
                intensities = intensities[inds]
                # wavelength sampling, add telluric lines
                if arm.SAMPLING == 'grid':
                    wf.wavegrid(
                        arm,
                        wavelengths=wavelengths,
                        intensities=intensities,
                        assign=True,
                        telluric=arm.tell)
                elif arm.SAMPLING == 'mc':
                    arm.wavelengths = wavelengths
                    arm.intensities = intensities
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])
                if arm.tell:
                    arm.intensities = wf.convolve_telluric_lines(
                        arm,
                        arm.wavelengths,
                        arm.intensities)
                arm.fib_obstype[i] = 'STAR'
                arm.fib_src[i] = 'OBJ'

            elif simmode is 'T':
                from importtharlines import import_thar_lines
                arm.fiber_description = 'Thorium Argon line list'
                wavelengths, intensities = import_thar_lines()
                inds = (wavelengths >= arm.wmin) & (wavelengths <= arm.wmax)
                arm.wavelengths = wavelengths[inds]
                arm.intensities = intensities[inds]
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])
                arm.infiles[i] = "Lovis,2007; Kerber,2008"
                arm.fib_obstype[i] = 'WAVE'
                arm.fib_src[i] = 'ThAr'

            elif simmode is 'U':
                arm.fiber_description = 'Uranium Neon calibration lamp'
                fn = "MODELS/UNE/une.npy"
                arm.infiles[i] = fn
                wavelengths, intensities = np.load(fn)
                inds = (wavelengths >= arm.wmin) & (wavelengths <= arm.wmax)
                wavelengths = wavelengths[inds]
                intensities = intensities[inds]
                if arm.fib_rv[i]:
                    wavelengths = redshift(wavelengths, arm.fib_rv[i])
                if arm.SAMPLING == "grid":
                    wf.wavegrid(
                        arm,
                        wavelengths=wavelengths,
                        intensities=intensities,
                        assign=True,
                        telluric=False)
                elif arm.SAMPLING == "mc":
                    arm.wavelengths = wavelengths
                    arm.intensities = intensities
                arm.fib_obstype[i] = 'WAVE'
                arm.fib_src[i] = 'UNe'

            elif simmode is 'W':
                arm.fiber_description = 'wavelength mapping'
                arm.wavelengths = wf.calculate_wavelengths(arm, mode='CCD', nwaves=arm.nw)
                arm.intensities = arm.wavelengths
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])
                arm.wavemap = True
                arm.fib_obstype[i] = 'WAVE'
                arm.fib_src[i] = 'WAVE'

            elif simmode is 'X':
                arm.fiber_description = 'Fabry-Perot spectrum simulation'
                # arm.wavelengths = np.load('phx_wavelengths.npy')
                # print arm.wavelengths
                arm.wavelengths = wf.calculate_wavelengths(arm, mode='CCD', nwaves=arm.nw)
                arm.intensities = fabry_perot(arm.wavelengths)
                if arm.fib_rv[i]:
                    arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i])
                arm.fib_obstype[i] = 'WAVE'
                arm.fib_src[i] = 'FP'

            ## APPLY RV SHIFT TO WAVELENGTHS
            #if arm.fib_rv[i]:
                #arm.wavelengths = redshift(arm.wavelengths, arm.fib_rv[i]) # @MZ arm.fib_rv[fiber] not fiber

            # RUN SIMULATION
            sim_general(arm, i)

            print "Fiber %s simulation finished\n" % arm.fib_char[i]

    # PROCESS ARRAYS INTO A SINGLE IMAGE / ADD NOISE
    arm.add_full_sim_time(sim_start)          # start time
    arm.add_full_sim_time(time.time() - t1)   # finish time
    arm.image = arm.images[0] + arm.images[1] # final image array
     # ADD NOISE
    # if arm.noise:
    add_noise(arm)

    # prevents value overflow wraparound
    arm.image[arm.image > np.iinfo(np.uint16).max] = np.iinfo(np.uint16).max
    arm.image = np.array(arm.image, dtype=np.uint16) # noise converts array to int64 at some step - let's bring it back to uint16
    print "Converting image array to %s" % arm.image.dtype

    # WRITE FITS FILE
    write_to_fits(arm)

    # WAVELENGTH SOLUTIONS; CREATE WAVETRACE REGION FILES
    for i in xrange(len(arm.wt)):
        if arm.wt[i]:
            wlt.wavelength_trace(arm, i)

    # SPAWN FITS OUTPUT IMAGE IN SAO-DS9
    if args.ds9 or args.ds9_mosaic:
        if args.ds9_mosaic:
            ds9_param = '-mosaicimage iraf'
        else:
            ds9_param = '-mecube'
        if arm.WT_FLAG:
            if args.no_compress:
                call = 'ds9 %s.fits %s -region %s.reg' % (ds9_param, rm.outfile, arm.outfile)
            else:
                call = 'ds9 %s.fits.gz %s -region %s.reg' % (ds9_param, arm.outfile, arm.outfile)
        else:
            if args.no_compress:
                call = 'ds9 %s %s.fits' % (ds9_param, arm.outfile)
            else:
                call = 'ds9 %s %s.fits.gz' % (ds9_param, arm.outfile)
        print "Executing '%s'" % call
        try:
            subprocess.check_call(call.split())
        except OSError, e:
            print e
            print "Shell call failed. You can just run the following line in the shell:"
            print call

if __name__ == "__main__":
    main()

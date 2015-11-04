"""
wavelengthtrace.py

Performs wavelength solution of a point source, and finds the average
for each pixel.
Writes solutions to disk, and writes a DS9 region file according to user input.

23 Aug, 2013:
- by default, now always computes wavelengths on the fly, which seems just
  as fast as loading from disk and finding_nearest point
  NOTE possibly _solve_wavelengths() and _write_region_read() can be deprecated in favor of
    _write_region_solve()
- order overlap bug fixed

08 Mar, 2013:
- now does not write to *.gz by default (to keep the txt files readable)
- write_region() now loops over every order to trace repeated wavelengths
- <>_wavelength_solutions_fiber_<>.txt order changed, now:
  wavelength[ang], order, x[pix], y[pix]

30 Jan, 2013:
- self changed to arm for consistency with the rest of the program
"""
import numpy as np
import time
import os
from physics import find_nearest, n_sell
from printer import Printer
from printtiming import print_timing
from simgeneral_c import sim_general
from slitfunctions import point_source
import cfunctions
import wavefuncs


def _solve_wavelengths(arm, fiber, wave_trace=True, write=True):
    """
    Finds wavelength solutions using a point source slit function.
    """
    solution_file = 'DATA/%s_wavelength_solutions_fiber_%s.txt' % (arm.ARM.lower(), arm.fib_char[fiber].lower())

    if arm.fib_char[fiber] in ['A', 'B']: #@MZ
        offset = arm.fib_OFFSET[fiber]

    # CREATE POINT SOURCE SLIT FUNCTION
    slitx, slity = point_source(arm, offset=offset)
    nslit = slitx.size
    wavelengths = wavefuncs.calculate_wavelengths(arm, mode="wavetrace", nwaves=arm.dwt) # reminder: wavelengths is 2d array in this case
    intensities = np.ones(wavelengths.shape) # SET ALL INTENSITIES TO 1
    orders = arm.OSET
    norders = arm.NOSET

    # INITIALIZE IMAGE ARRAYS
    image = np.zeros(arm.CCD_DIMS)
    counts = np.zeros(arm.CCD_DIMS, dtype=np.uint)
    m_list = np.zeros(arm.CCD_DIMS, dtype=np.uint) # records order of each pixel

    # OTHER INITIALIZATIONS
    returnx = np.empty(0)
    returny = np.empty(0)
    LOC_FLAG = 0
    BLAZE_FLAG = 0

    for i,m in enumerate(orders):
        waves, _weights = wavefuncs.feed_wavelengths(arm, m, wavelengths=wavelengths, intensities=intensities)
        nwaves = waves.size
        n_g_sell = n_sell(arm.ARM, waves)

        cfunctions.compute(arm.ARM_FLAG,
                           BLAZE_FLAG,
                           LOC_FLAG,
                           nwaves,
                           nslit,
                           m,
                           arm.XD_0,
                           arm.YD_0,
                           np.ascontiguousarray(n_g_sell, dtype=np.float64),
                           np.ascontiguousarray(slitx, dtype=np.float64),
                           np.ascontiguousarray(slity, dtype=np.float64),
                           np.ascontiguousarray(waves, dtype=np.float64),
                           np.ascontiguousarray(waves, dtype=np.float64),
                           image,
                           np.ascontiguousarray(counts),
                           np.ascontiguousarray(m_list),
                           returnx,
                           returny)
        output = '  * Wavetrace order %i (%i of %i), %.2f%% complete.' % (m, i+1, norders, ((i+1.) * 100.0 / norders))
        Printer(output)
    print '\n'
    inds = np.where(counts != 0)
    image[inds] *= 1.0e7 / counts[inds] # average each pixel, convert mm to Angstrom

    # saves solutions to disk
    if write:
        newimage = image[inds]
        m_list = m_list[inds]
        y = inds[0]
        x = inds[1]
        n = x.size
        print "  * Saving wavelength solutions to '%s'" % solution_file
        filename = open(solution_file, 'w')
        for i,lamb in enumerate(newimage):
            filename.write( '%s %s %s %s\n' % (m_list[i], lamb, x[i], y[i]) )
            output1 = '  * %0.2f%% Complete. point %s of %s.' % \
                      (100.0 * float(i) / n, i+1, n)
            Printer(output1)
        print "\n"
        filename.close()
        return m_list, newimage, x, y
    return 0


def _write_region_read(arm, fiber, m_solutions=None, lamb_solutions=None, x=None,
    y=None):
    """
    Writes .reg file, using regions from a pre-written wavelength solution file.

    lamb_solutions are in Ang
    """
    solution_file = 'DATA/%s_wavelength_solutions_fiber_%s.txt' % \
        (arm.ARM.lower(), arm.fib_char[fiber].lower())
    print "  * Writing region file to '%s.reg'" % arm.outfile
    dlamb = arm.dwt         # Ang
    wtlist = arm.wtlist     # Ang
    if arm.WT_FLAG:
        outfile = open('%s.reg' % arm.outfile, 'a')
    elif not arm.WT_FLAG:
        outfile = open('%s.reg' % arm.outfile, 'w')
        outfile.write('# Region file format: DS9 version 4.0\n')
        outfile.write('# Filename: %s.fits\n' % arm.outfile.replace('FITS/', ''))
        outfile.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
        outfile.write('physical\n')
    # import wavelist, or create a wavelist using the ccd limits, ang -> mm
    if wtlist:
        input_waves = np.loadtxt(wtlist, unpack=True) * 1.e-7
    else:
        input_waves = np.arange(arm.wmin*1.e7, arm.wmax*1.e7, dlamb) * 1.e-7# step size of dlamb in Angstrom
    _weights = np.zeros(input_waves.size) # dummy variable for feed_wavelengths()
    orders = arm.OSET
    norders = arm.NOSET

    # find wavelengths for each order
    for i,m in enumerate(orders):
        waves, _intensities = wavefuncs.feed_wavelengths(arm, m,
            wavelengths=input_waves, intensities=_weights)
        inds = np.where(m_solutions == m) # use solutions for a given order
        for wavelength in waves:
            xi = x[inds]
            yi = y[inds]
            j, lamb = find_nearest(lamb_solutions[inds], wavelength*1.e7) # mm to ang
            # FITS indices start at 1, not 0
            outfile.write('point(%s,%s) # point=cross text={%.4f}\n' % \
                (xi[j]+1, yi[j]+1, lamb))
        output = '  * Region order %s (%s of %s), %.2f%% complete.' % \
            (m, i+1, norders, ((i+1.) * 100.0 / norders))
        Printer(output)
    print "\n"
    outfile.close()
    arm.set_wt_flag(True) # set flag to enable appending for region file
    return 0

def _write_region_solve(arm, fiber, m_solutions=None, lamb_solutions=None, x=None,
    y=None):
    """
    Writes .reg file using regions from a pre-written wavelength solution file.

    lamb_solutions are in Ang
    NOTE possibly other functions in this module can be deprecated in favor of
    this one.
    """
    solution_file = 'DATA/%s_wavelength_solutions_fiber_%s.txt' % \
        (arm.ARM.lower(), arm.fib_char[fiber].lower())
    print "  * Writing region file to '%s.reg'" % arm.outfile
    dlamb = arm.dwt         # Ang
    wtlist = arm.wtlist     # Ang
    if arm.WT_FLAG:
        outfile = open('%s.reg' % arm.outfile, 'a')
    elif not arm.WT_FLAG:
        outfile = open('%s.reg' % arm.outfile, 'w')
        outfile.write('# Region file format: DS9 version 4.0\n')
        outfile.write('# Filename: %s.fits\n' % arm.outfile.replace('FITS/', ''))
        outfile.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
        outfile.write('physical\n')

    # INITIALIZATIONS
    if arm.fib_char[fiber] in ['A', 'B']: #@MZ
        offset = arm.fib_OFFSET[fiber]
    slitx, slity = point_source(arm, offset=offset)
    nslit = slitx.size
    # import wavelist, or create a wavelist using the ccd limits, ang -> mm
    if wtlist:
        input_waves = np.loadtxt(wtlist, unpack=True) * 1.e-7
    else:
        input_waves = np.arange(arm.wmin*1.e7, arm.wmax*1.e7, dlamb) * 1.e-7# step size of dlamb in Angstrom
    _weights = np.zeros(input_waves.size) # dummy variable for feed_wavelengths()
    orders = arm.OSET
    norders = arm.NOSET
    # INITIALIZE IMAGE ARRAYS
    image = np.zeros(arm.CCD_DIMS)
    counts = np.zeros(arm.CCD_DIMS, dtype=np.uint)
    m_list = np.zeros(arm.CCD_DIMS, dtype=np.uint) # records order of each pixel

    # OTHER INITIALIZATIONS
    LOC_FLAG = 2
    BLAZE_FLAG = 0

    # find wavelengths for each order
    for i,m in enumerate(orders):
        waves, _intensities = wavefuncs.feed_wavelengths(arm, m,
            wavelengths=input_waves, intensities=_weights)
        nwaves = waves.size
        n_g_sell = n_sell(arm.ARM, waves)
        returnx = np.empty(nwaves)
        returny = np.empty(nwaves)
        cfunctions.compute(arm.ARM_FLAG,
                           BLAZE_FLAG,
                           LOC_FLAG,
                           nwaves,
                           nslit,
                           m,
                           arm.XD_0,
                           arm.YD_0,
                           np.ascontiguousarray(n_g_sell, dtype=np.float64),
                           np.ascontiguousarray(slitx, dtype=np.float64),
                           np.ascontiguousarray(slity, dtype=np.float64),
                           np.ascontiguousarray(waves, dtype=np.float64),
                           np.ascontiguousarray(waves, dtype=np.float64),
                           image,
                           np.ascontiguousarray(counts),
                           np.ascontiguousarray(m_list),
                           returnx,
                           returny)
        x = returnx
        y = returny
        print x
        print y
        # FITS indices start at 1, not 0
        outfile.write("".join(['point(%i,%i) # point=cross text={%.0f}\n' % \
            tup for tup in zip(x+1, y+1, waves*1.e7)]))  # ds9 has 1-based indexing ? @CM yes, I believe so
        output = '  * Region order %s (%s of %s), %.2f%% complete.' % \
            (m, i+1, norders, ((i+1.) * 100.0 / norders))
        Printer(output)
    print "\n"
    outfile.close()
    arm.set_wt_flag(True) # set flag to enable appending for region file
    return 0
    
# :::::::::::::::::::::::::::: MAIN FUNCTION :::::::::::::::::::::::::::::::::
    
@print_timing
def wavelength_trace(arm, fiber, infile=None):
    print "Fiber %s - wavelength trace" % (arm.fib_char[fiber])
    #solution_file = 'DATA/%s_wavelength_solutions_fiber_%s.txt' % \
        #(arm.ARM.lower(), arm.fib_char[fiber].lower())

    ## check if solutions already exist
    #if arm.solve or not os.path.isfile(solution_file):
        #m, lamb_solutions, x, y = _solve_wavelengths(arm, fiber, write=True)
    #else:
        ## loading is potential bottleneck
        #print "  * Loading wavelength solutions for %s arm..." % (arm.ARM)
        #m, lamb_solutions, x, y = np.loadtxt(solution_file, unpack=True)  # @CM

    ## write region file
    #_write_region_read(arm, fiber, m_solutions=m, lamb_solutions=lamb_solutions,
        #x=x, y=y)

    # write region file
    _write_region_solve(arm, fiber)
    return 0
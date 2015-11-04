"""
wavefuncs.py

Module containing wavelength functions.
"""
#import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
try:
    import pyfits
except ImportError:
    import astropy.io.fits as pyfits
import os
from physics import *
import cfunctions
from printer import Printer
import defaults

def _compute(self, slitx, slity, m, waves, returnx, returny):
    """Wrapper for cfunctions.compute"""
    nw = waves.size
    ns = slitx.size
    ng = n_sell(self.ARM, waves)
    cfunctions.compute(
        self.ARM_FLAG,
        0, # blaze flag
        1, # location flag
        self.GAP_FLAG,
        nw,
        ns,
        m,
        self.gap,
        self.gapoffset,
        self.XD_0,
        self.YD_0,
        np.ascontiguousarray(ng, dtype=np.float64),
        np.ascontiguousarray(slitx, dtype=np.float64),
        np.ascontiguousarray(slity, dtype=np.float64),
        np.ascontiguousarray(waves, dtype=np.float64),
        np.empty(0), # weights
        np.empty((0,0)), # ccd
        np.empty((0,0), dtype=np.uint), # counts
        np.empty((0,0), dtype=np.uint), # order list
        np.ascontiguousarray(returnx, dtype=np.float64),
        np.ascontiguousarray(returny, dtype=np.float64))
    return returnx, returny

# INITIALIZATION METHODS
def find_blaze_offset(self):
    """
    Finds the echelle blaze offset of the detector.
    Sets values to self.XD_0
    """
    print "  * Solving offset on detector for %s arm..." % (self.ARM)

    orders = self.corner_orders if self.ARM == "VIS" else self.zemax_orders
    norders = orders.size
    wavelengths = self.corner_wavelengths if self.ARM == "VIS" else self.zemax_wavelengths
    nw = wavelengths[0].size
    s = 0.0
    slitx = np.array(s)
    slity = np.array(s)
    xd = np.empty(wavelengths.shape)
    yd = np.empty(wavelengths.shape)
    ng = np.empty(wavelengths.shape)    # ADDED
    returnx = np.empty(wavelengths[0].size) # ADDED
    returny = np.empty(wavelengths[0].size) # ADDED

    for i in xrange(norders):
        lamb = wavelengths[i]
        returnx, returny = _compute(self, slitx, slity, orders[i], lamb, returnx, returny)
        xd[i] = returnx # mm
        yd[i] = returny # mm
    xdmin = xd.min()
    xdmax = xd.max()
    xdrange = np.absolute(xdmax - xdmin)
    ydmin = yd.min()
    ydmax = yd.max()
    ydrange = np.absolute(ydmax - ydmin)
    # SET X OFFSET CORRECTION TO CENTER THE FSR ON DETECTOR [mm]
    #self.XD_0 = (xdmin + xdmax) / 2.0 # mm
    #self.XD_0 = (np.abs(xdmax) - np.abs(xdmin)) / 2.0 # mm
#    if self.ARM == "NIR":
#        self.XD_0 = -(xd[0].min() + xd[0].max()) / 2.0
    XD_0 = -(xdmin + xdmax) / 2.0 # mm
    return XD_0


def find_ccd_wavelength_limits(self, dlamb=1.0e-6):
    """
    Find the CCD wavelength range of each order.
    Sets values to self.xleftedge, self.xrightedge, self.ccd_waves_min,
    self.ccd_waves_max, and sets self.YD_0
    """

    self.xleftedge = -self.NXPIX/2.0 * self.DPIX    # left edge of ccd [mm]
    self.xrightedge = self.NXPIX/2.0 * self.DPIX    # right edge of ccd [mm]

    orders = self.ORDERS
    norders = orders.size

    # FOR LEFT EDGE
    xle = np.empty(self.NORDERS)
    yle = np.empty(self.NORDERS)
    lle = np.empty(self.NORDERS)
    print 'Solving wavelengths at the left edge of the detector'\
        ' for %s arm...' % self.ARM

    slitx = np.array((0.0))
    slity = np.array((0.0))
    ng = np.array((0.0))
    returnx = np.array((0.0))
    returny = np.array((0.0))

    for i in xrange(norders):
        lamb = self.bl_ech[i]
        xd = 0.0
        yd = 0.0
        while np.abs(xd) < np.abs(self.xleftedge):
            returnx, returny = _compute(self, slitx, slity, orders[i], lamb, returnx, returny)
            xd = returnx
            yd = returny
            lamb -= dlamb
        xle[i] = xd
        yle[i] = yd
        lle[i] = lamb
        output = '  %0.2f%% Complete.' % ((i+1.0)/self.NORDERS * 100.0)
        Printer(output)
    print ''
    self.ccd_waves_min = lle

    # FOR RIGHT EDGE
    xre = np.empty(self.NORDERS)
    yre = np.empty(self.NORDERS)
    lre = np.empty(self.NORDERS)
    print 'Solving wavelengths at the right edge of the detector'\
        ' for %s arm...' % self.ARM
    for i in xrange(norders):
        lamb = self.bl_ech[i]
        xd = 0.0
        yd = 0.0
        while np.abs(xd) < np.abs(self.xrightedge):
            returnx, returny = _compute(self, slitx, slity, orders[i], lamb, returnx, returny)
            xd = returnx
            yd = returny
            lamb += dlamb
        xre[i] = xd
        yre[i] = yd
        lre[i] = lamb
        output = '  %0.2f%% Complete.' % ((i+1.0) / self.NORDERS * 100.0)
        PP = Printer(output)
    print ''
    self.ccd_waves_max = lre

    # find y-offset [mm], YD_0 is center of y-coordinates
    ys = np.concatenate((yle, yre))
    #self.YD_0 = -(ys.min() + ys.max()) / 2.0 # mm
    YD_0 = -(ys.max() + ys.min()) / 2.0 # mm
    return YD_0


def calculate_wavelengths(self, orders=None, nwaves=None, mode=None, buff=1.e-7):
    """
    Calculates the wavelengths to be used for the simulation.

    Flatfield finds wavelengths at the center of each pixel.
    """

    print "  * Calculating wavelengths for %s arm..." % self.ARM
    if orders is None:
        orders = self.ORDERS
    mi = orders - self.M0
    norders = orders.size

    if nwaves is None:
        nwaves = self.DEFAULT_NWF * self.NXPIX # default oversampling factor of 10
    else:
        nwaves *= self.NXPIX

    # DIVIDE WAVELENGTH RANGE OF EACH ORDER INTO EVENLY SPACED INTERVALS
    # WAVELENGTHS is a 2D ARRAY
    if mode is 'wavetrace':
        wmin = self.ccd_waves_min
        wmax = self.ccd_waves_max
        wavelengths = np.empty((norders, nwaves))
        for i in xrange(orders.size):
            wavelengths[i, :] = np.linspace(wmin[i], wmax[i], nwaves)
            output = '  * %0.2f%% Complete.' % ((i+1.0)/norders * 100.0)
            PP = Printer(output)
        print ''

    elif mode is 'FSR':
        wmin = self.fsr_waves_min
        wmax = self.fsr_waves_max
        wavelengths = np.empty((norders, nwaves))
        for i in xrange(orders.size):
            wavelengths[i, :] = np.linspace(wmin[i], wmax[i], nwaves)
            output = '  * %0.2f%% Complete.' % ((i+1.0)/norders * 100.0)
            PP = Printer(output)
        print ''
    elif mode is 'CCD':
        #wmin = self.ccd_waves_min.min()
        #wmax = self.ccd_waves_max.max()
        # DIVIDE ENTIRE WAVELENGTH RANGE INTO EVENLY SPACED INTERVALS
        # WAVELENGTHS is a 1D ARRAY
        #wavelengths = np.linspace(wmin, wmax, num=nwaves*norders)

        wmin = self.wmin - buff
        wmax = self.wmax + buff

        # determine number of samples based on desired resolution
        dw = self.nw * 1.e-7 # angstrom -> mm
        n = np.int(np.ceil( (wmax-wmin) / dw ))
        wavelengths = np.linspace(wmin, wmax, num=n)
        print "Evaluating wavelength grid to %i linearly spaced samples\n  dw = %4.e Angstrom" % (n, self.nw)
    else:
        raise ValueError
    return wavelengths


def feed_wavelengths(self, m, wavelengths=None, intensities=None):
    """
    calculates the wavelengths and intensities to be fed into
    the simulator for a given order
    """
    if wavelengths is None: wavelengths = self.wavelengths
    if intensities is None: intensities = self.intensities
    i = m - self.M0
    wmin = self.ccd_waves_min[i]
    wmax = self.ccd_waves_max[i]
    indices = (wavelengths >= wmin) & (wavelengths <= wmax)
    wavelengths = np.array(wavelengths[indices])#, order='C', dtype=np.double)
    intensities = np.array(intensities[indices])#, dtype=np.double)#order='C',
    return wavelengths, intensities


def wavegrid(self, wavelengths=None, intensities=None, factor=10.0,
    buff=1.0e-7, kind=3, assign=False, telluric=False):
    """
    Re-evaluates input spectra on a pre-determined wavelength grid.
    """
    if wavelengths is None:
        wavelengths = self.wavelengths
    if intensities is None:
        intensities = self.intensities
    wmin = self.wmin - buff
    wmax = self.wmax + buff

    # determine number of samples based on desired resolution
    dw = self.nw * 1.e-7 # angstrom -> mm
    n = np.int(np.ceil( (wmax-wmin) / dw ))

    #cond = np.logical_and(wavelengths >= wmin, wavelengths <= wmax)
    #n = wavelengths[cond].size * factor

    # new evenly spaced wavelength grid
    #wavelengths2 = np.linspace(wmin-buff, wmax+buff, num=n)
    # this still does not work quite right...resolution decreases towards
    # smaller wavelengths, which ENHANCES the gradient effect
    # if self.wlog:
    #     n = defaults.RESOLUTION * self.wlog
        #wavelengths2 = np.logspace(np.log(wmin),np.log(wmax), num=n, base=np.exp(1.0))
        #wavelengths2 = wavelengths2  # reverse
        #print wavelengths2
        #print "dw = %s" % (np.diff(wavelengths2)* 1.e7)

        # wavelengths2 = np.logspace(np.log(wmin),np.log(wmax),  num=n, base=np.exp(1.0))
        #print wavelengths2
        #print "dw = %s" % (np.diff(wavelengths2)* 1.e7)
        # print "Re-evaluating wavelength grid to %i logarithmically spaced samples" % (wavelengths2.size)
    # else:
    wavelengths2 = np.linspace(wmin, wmax, num=n)
    print "Re-evaluating wavelength grid to %i linearly spaced samples\n  dw = %4.e Angstrom" % (n, self.nw)
    #print "Re-evaluating wavelength grid to %i evenly spaced samples, dw = %4.e nm" % (n, np.diff(wavelengths2)[0]*1.e7)
    
    f = InterpolatedUnivariateSpline(wavelengths, intensities)
    intensities2 = f(wavelengths2)

    # telluric lines
    #if telluric:
        #print "Adding telluric species: %s" % self.tell
        #tellwaves, tellflux = load_telluric_lines(self)
        #ftell = InterpolatedUnivariateSpline(tellwaves, tellflux)
        #tellflux2 = ftell(wavelengths2)
    #else:
        #print "No telluric species added."
        #tellflux2 = np.ones(wavelengths2.shape)

    # preview
    #plt.plot(wavelengths2[::100], intensities2[::100], c="k", alpha=0.6)
    #plt.plot(wavelengths2[::100], intensities2[::100]*tellflux2[::100], label="w/ telluric lines", c="r", alpha=0.6)
    #plt.xlabel("Wavelength [Angstrom]")
    #plt.ylabel("Normalized Flux")
    #plt.legend()
    #plt.show()

    if assign:
        self.wavelengths = wavelengths2
        self.intensities = intensities2 #* tellflux2


def load_telluric_lines(self):
    """
    Adds telluric lines to spectra.
    """
    tell_dir = os.path.join("MODELS", "TELLURIC")
    wmin = self.wmin
    wmax = self.wmax

    # import all fluxes; wavelengths are same for all species
    fluxes = []
    species_fn = [os.path.join(tell_dir, "LBLRTM_%s_+0.0.npy" % (specie)) for specie in self.tell]
    for species in species_fn:
        data = np.load(species)
        fluxes.append(data[1])
    waves = data[0]
    all_spec = np.prod(fluxes, axis=0) # multiply all species together
    #plt.plot(waves, all_spec)
    #plt.xlabel("Angstrom")
    #plt.ylabel("Normalized flux")
    #plt.xlim([wmin, wmax])
    #plt.show()
    return waves, all_spec

    
def convolve_telluric_lines(self, wavelengths, intensities):
    """
    Multiplies spectra by telluric lines.

    Parameters
    ----------
    wavelengths : array_like
        Wavelength array in Angstroms.
    intensities : array_like
        Intensity array in arbitrary units.

    Returns
    -------
    intensities : array_like
        intensities convolved with telluric species.

    """
    print "Loading telluric species: %s" % self.tell
    tellwaves, tellflux = load_telluric_lines(self)
    ftell = InterpolatedUnivariateSpline(tellwaves, tellflux)
    tellflux = ftell(wavelengths)
    return intensities * tellflux
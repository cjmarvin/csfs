"""
debug.py

Debugging module for spectrograph simulation. This is to be used to check the
footprint diagram of the spectrographs, not for the optimization purposes
(but that could be implemented).
"""

import numpy as np
#import matplotlib.pyplot as plt
from printer import Printer
from physics import *
from slitfunctions import *
import cfunctions

def fmtcols(mylist, cols):
    lines = ("\t".join(mylist[i:i+cols]) for i in xrange(0,len(mylist),cols))
    return '\n'.join(lines)

def print_wave_table(arm):
    dvline1 = ":::::::::::::::::::::::::::::::::::::::::::::::::"
    header = "FSR Wavelengths"
    columns = "Order\tWmin [A]\tWblaze [A]\tWmax [A]"
    dvline2 = "-------------------------------------------------"
    print dvline1
    print header
    print columns
    print dvline2
    for m,wmin,wblaze,wmax in zip(arm.ORDERS, arm.fsr_waves_min, arm.bl_ech, arm.fsr_waves_max):
        print "%i:\t%.3f\t%.3f\t%.3f" % (m, wmin*1.e7, wblaze*1.e7, wmax*1.e7)

    dvline1 = ":::::::::::::::::::::::::::::::::::::::::::::::::"
    header = "CCD Wavelengths"
    columns = "Order\tWmin [A]\tWblaze [A]\tWmax [A]"
    dvline2 = "-------------------------------------------------"
    print dvline1
    print header
    print columns
    print dvline2
    for m,wmin,wblaze,wmax in zip(arm.ORDERS, arm.ccd_waves_min, arm.bl_ech, arm.ccd_waves_max):
        print "%i:\t%.3f\t%.3f\t%.3f" % (m, wmin*1.e7, wblaze*1.e7, wmax*1.e7)
    print ""


def compute(arm, orders, wavelengths, slit=0):
    """Compute wrapper for debug."""
    
    norders = orders.size
    nw = wavelengths.shape[0]
    
    if slit is 0:
        slitfunc = point_source
    else:
        slitfunc = half_moon_even
    slitx, slity = slitfunc(arm, offset=0.0)
    print "xd_0 = %.2f mm" % arm.XD_0
    print "yd_0 = %.2f mm" % arm.YD_0
    

    xd = np.empty(wavelengths.shape)
    yd = np.empty(wavelengths.shape)
    ng = np.empty(wavelengths.shape)
    returnx = np.empty(wavelengths[0].size)
    returny = np.empty(wavelengths[0].size)
    

    for i in xrange(norders):
        lamb = wavelengths[i]
        ng = n_sell(arm.ARM, lamb)
        print "n_sell = ", ng
        cfunctions.compute(arm.ARM_FLAG,
                           0, # blaze flag
                           1, # location flag
                           nw,
                           slitx.size,
                           orders[i],
                           arm.XD_0,
                           arm.YD_0,
                           np.ascontiguousarray(ng, dtype=np.float64),
                           np.ascontiguousarray(slitx, dtype=np.float64),
                           np.ascontiguousarray(slity, dtype=np.float64),
                           np.ascontiguousarray(lamb, dtype=np.float64),
                           np.empty(0), # weights
                           np.empty((0,0)), # ccd
                           np.empty((0,0), dtype=np.uint), # counts
                           np.empty((0,0), dtype=np.uint), # order list
                           np.ascontiguousarray(returnx, dtype=np.float64),
                           np.ascontiguousarray(returny, dtype=np.float64))
        xd[i] = returnx # mm
        yd[i] = returny # mm

    xd = xd.flatten()/ arm.DPIX
    yd = yd.flatten()/arm.DPIX
    xdmin = xd.min()
    xdmax = xd.max()
    xdrange = np.absolute(xdmax - xdmin)
    ydmin = yd.min()
    ydmax = yd.max()
    ydrange = np.absolute(ydmax - ydmin)
    return wavelengths.flatten(), xd, yd, xdrange, ydrange

def plot_corners(arm, slit=0):

    orders = arm.corner_orders
    wavelengths = arm.corner_wavelengths

    wavelengths, xd, yd, xdrange, ydrange = compute(arm, orders, wavelengths)

    print "::::: FSR WAVELENGTHS :::::"
    print "xdrange = %.2f pix" % xdrange
    print "ydrange = %.2f pix" % ydrange
    xdmin = xd.min()
    xdmax = xd.max()
    ydmin = yd.min()
    ydmax = yd.max()

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="auto")
    ax.set_title("%s post-xy offset correction" % arm.ARM)
    ax.scatter(xd, yd, label="FSR")
    #ax.grid()
    for i,wave in enumerate(wavelengths):
        ax.annotate("%.2f" % (wave*1.e7), xy=(xd[i], yd[i]))
    ax.vlines(-arm.NXPIX/2, ydmin, ydmax)
    ax.vlines(arm.NXPIX/2, ydmin, ydmax)
    ax.fill([-arm.NXPIX/2, arm.NXPIX/2, arm.NXPIX/2, -arm.NXPIX/2],
        [-arm.NYPIX/2, -arm.NYPIX/2, arm.NYPIX/2, arm.NYPIX/2],
        'r',
        alpha=0.1,
        edgecolor='g')
    plt.legend()
    plt.show()


def plot_zemax(arm):
    orders = arm.zemax_orders
    wavelengths = arm.zemax_wavelengths
    wavelengths, xd, yd, xdrange, ydrange = compute(arm, orders, wavelengths)

    print "::::: ZEMAX WAVELENGTHS :::::"
    print "xdrange = %.2f pix" % xdrange
    print "ydrange = %.2f pix" % ydrange
    print "xd_0 = %.2f mm" % arm.XD_0
    print "yd_0 = %.2f mm" % arm.YD_0
    
    xdmin = xd.min()
    xdmax = xd.max()
    ydmin = yd.min()
    ydmax = yd.max()
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="auto")
    ax.set_title("%s post-xy offset correction" % arm.ARM)
    ax.scatter(xd, yd, label="ZEMAX")
    #ax.grid()
    for i,wave in enumerate(wavelengths):
        ax.annotate("%.2f" % (wave*1.e7), xy=(xd[i], yd[i]))
    ax.vlines(-arm.NXPIX/2, ydmin, ydmax)
    ax.vlines(arm.NXPIX/2, ydmin, ydmax)
    ax.fill([-arm.NXPIX/2, arm.NXPIX/2, arm.NXPIX/2, -arm.NXPIX/2],
       [-arm.NYPIX/2, -arm.NYPIX/2, arm.NYPIX/2, arm.NYPIX/2],
       'r',
       alpha=0.1,
       edgecolor='g')
    plt.legend()
    plt.show()
    
    
# :::::::::::::::::::::::::: MAIN ::::::::::::::::::::::::::::::::::

def debug(arm):
    #print_wave_table(arm)
    #plot_corners(arm)
    plot_zemax(arm)
    #import_zemax_table(arm)
    #plot_slit()
    #plot_offset()
    #plot_fsr()

    return 0
"""plant.py

This function pickles the included FITS files into the NumPy NPY binary format,
allowing for much faster reading and loading.
"""
#import matplotlib.pylab as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import glob
import pprint
import sys


tell_dirname = "MODELS/TELLURIC/"
species = "CH4 CO2 H2O N2O O2 O3".split()
une_dirname = "MODELS/UNE/"
comb_dirname = "MODELS/COMB/"

#def _tellurictxt2fits():
    #"""
    #This helper function is only used in the construction of the program.
    #"""
    #species = "CH4 CO2 H2O N2O O2 O3".split()

    #species_fn = ["%sLBLRTM_%s_+0.0" % (tell_dirname, specie) for specie in species]
    #pprint.pprint(species_fn)
    #for fn in species_fn:
        #print "Loading %s" % fn
        #waves, flux = np.loadtxt(fn, unpack=True)
        #print "Converting ln(wavelengths) to wavelengths..."
        #waves = np.exp(waves)
        #print "Opening HDU..."
        #hdu = pyfits.PrimaryHDU([waves, flux])
        #print "Adding header entries..."
        #header = hdu.header
        #header["AXIS1U"] = ("Angstrom", "NAXIS1 units")
        #header["AXIS2U"] = ("Normalized Flux", "NAXIS2 units")
        #for h in hdu.header:
            #print h, hdu.header[h], hdu.header.comments[h]
        #print "Writing %s.fits to disk" % fn
        #hdu.writeto("%s.fits" % fn)

#def _unetxt2fits():
    #"""
    #This helper function is only used in the construction of the program.
    #"""
    
    #filename = "MODELS/UNE/une_200scan_8ma_001res.gz"
    #print "Loading '%s'..." % filename
    #wavenumbers, flux = np.loadtxt(filename, unpack=True)
    #print "Converting wavenumbers to wavelengths..."
    #wavelengths = 10000.0 / wavenumbers[::-1] * 1.e4 # sorted, -> angstrom
    #flux = flux[::-1]

    ## PREPARE FITS FILE
    #writefn = "MODELS/UNE/une.fits"
    #hdu = pyfits.PrimaryHDU([wavelengths, flux])
    #print "Adding header entries..."
    #header = hdu.header
    #header["AXIS1U"] = ("Angstrom", "NAXIS1 units")
    #header["AXIS2U"] = ("Arbitrary Flux", "NAXIS2 units")
    #print "Writing '%s' to disk..." % writefn
    #hdu.writeto(writefn)

#def _combtxt2fits():
    #"""
    #This helper function is only used in the construction of the program.
    #"""

    #wave_fn = "MODELS/COMB/comb_lamb.dat.gz"
    #flux_fn = "MODELS/COMB/comb_ints.dat.gz"
    #print "Loading '%s' and '%s'..." % (wave_fn, flux_fn)
    #waves = np.loadtxt(wave_fn)
    #flux = np.loadtxt(flux_fn)

    ## PREPARE FITS FILE
    #writefn = "MODELS/COMB/comb.fits"
    #hdu = pyfits.PrimaryHDU([waves, flux])
    #print "Adding header entries..."
    #header = hdu.header
    #header["AXIS1U"] = ("Angstrom", "NAXIS1 units")
    #header["AXIS2U"] = ("Relative Intensity", "NAXIS2 units")
    #print "Writing '%s' to disk..." % writefn
    #hdu.writeto(writefn)

        
def plant_telluric_spectra():
    """
    Loads fits.gz files and saves to an npy file.

    wavelengths : mm
    flux : normalized to unity

    NOTE Handling of writing files needed? Probably not since dependencies
    are handled in Makefile.
    """
    spec_fn = glob.glob(tell_dirname + "*fits.gz")
    for species in spec_fn:
        print "Loading '%s'..." % species
        data = pyfits.getdata(species)
        data[0] *= 1.0e-7    # convert Ang to mm
        plant_fn = species.replace("fits.gz", "npy")
        print "Writing '%s' to disk." % plant_fn
        np.save(plant_fn, data)

def plant_une_spectra():
    """
    Loads fits.gz files and saves to an npy file.

    wavelengths : mm
    flux : normalized to unity    
    """
    fn = une_dirname + "une.fits.gz"
    print "Loading '%s'..." % fn
    data = pyfits.getdata(fn)
    n = data.size * 3
    data[0] *= 1.0e-7    # convert Ang to mm
    newwaves = np.linspace(data[0].min(), data[0].max(), n)
    f = InterpolatedUnivariateSpline(data[0], data[1])
    out_arr = np.array((newwaves, f(newwaves)))
    plant_fn = fn.replace("fits.gz", "npy")
    print "Writing '%s' to disk." % plant_fn
    np.save(plant_fn, out_arr)


def plant_comb_spectra():
    """
    Loads fits.gz files and saves to an npy file.

    wavelengths : mm
    flux : normalized to unity
    """
    fn = comb_dirname + "comb.fits.gz"
    print "Loading '%s'..." % fn
    waves, flux = pyfits.getdata(fn)
    waves = waves.ravel()
    flux = flux.ravel()
    n = waves.size * 10
    waves *= 1.0e-7    # convert Ang to mm
    
    newwaves = np.linspace(waves.min(), waves.max(), n)
    f = InterpolatedUnivariateSpline(waves, flux)
    out_arr = np.array((newwaves, f(newwaves)))
    plant_fn = fn.replace("fits.gz", "npy")
    print "Writing '%s' to disk." % plant_fn
    np.save(plant_fn, out_arr)


# ============================== MAIN ========================================
    
def main():

    #_tellurictxt2fits
    #_unetxt2fits
    #_combtxt2fits

    process = {
        'tell' : plant_telluric_spectra,
        'une' : plant_une_spectra,
        'comb' : plant_comb_spectra}
        #plant model spectra

    for arg in sys.argv[1:]:
        if arg in process:
            process[arg]()

if __name__ == "__main__":
    main()
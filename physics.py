from numpy import sin, cos, argmin, abs, searchsorted, sqrt
__all__ = ["n_sell",
           "lam_blaze_ech",
           "lam_blaze_grism",
           "lam_range",
           "find_nearest",
           "redshift",
           "fabry_perot"]

c = 299792458.0     # [m/s]
h = 6.62606957e10-34 # m^2 kg / s

def n_sell(arm_id, lamb):
    """
    Returns wavelength-dependent refractive index according to the Sellmeier
    equation.

    Parameters
    ----------
    arm : {"NIR", "VIS"}
        Spectral arm. "VIS" arm initializes LF5 constants, while "NIR"
        initializes Infrasil 301 constants. See References.
    lamb : array_like (M)
        Input wavelength(s) [nm]

    Returns
    -------
    array_like (M)
        Wavelength-dependent refractive index.

    Notes
    ----------
    LF5 (VIS) dispersion constants taken from SCHOTT.
    ..[1] http://www.schott.com/advanced_optics/us/abbe_datasheets/schott_datasheet_all_us.pdf

    Infrasil 301 (NIR) dispersion constants taken from Heraeus.
    ..[2] http://heraeus-quarzglas.com/media/webmedia_local/downloads/broschren_mo/dataandproperties_optics_fusedsilica.pdf
    ..[3] CARMENES-TNO068
    """
    lamb = lamb * 1.0e3 # mm to microns
    if arm_id is "VIS":
        B1 = 1.28035628
        B2 = 0.163505973
        B3 = 0.893930112
        C1 = 0.00929854416
        C2 = 0.0449135769
        C3 = 110.493685
        
        return sqrt( (B1 * lamb**2)/(lamb**2 - C1) + (B2 * lamb**2)/(lamb**2 - C2) + (B3 * lamb**2)/(lamb**2 - C3) + 1.0 )
    elif arm_id is "NIR":
        B1 = 4.76523070e-1
        B2 = 6.27786368e-1
        B3 = 8.72274404e-1
        C1 = 2.84888095e-3
        C2 = 1.18369052e-2
        C3 = 9.56856012e1
        S1 = 1.10098155807840
        S2 = 0.00124893116320
        S3 = 0.78106788864000
        L1 =0.00789065877562
        L2 =0.07653985064959
        L3 =85.96849321128080
        #return sqrt( (B1 * lamb**2)/(lamb**2 - C1) + (B2 * lamb**2)/(lamb**2 - C2) + (B3 * lamb**2)/(lamb**2 - C3) + 1.0 )
        return sqrt((S1*lamb**2/(lamb**2 - L1)) + (S2*lamb**2/(lamb**2 - L2)) + (S3*lamb**2/(lamb**2 - L3)) + 1.0)
    else:
        raise ValueError

        
def _test_nsell():
    print "VIS %.4f = 1.5814 ?" % n_sell("VIS", 587.6e-6) #== 1.58143765521
    print "NIR %.4f = 1.4602 ?" % n_sell("NIR", 546.07e-6) #== 1.46018134651


# PHYSICS FUNCTIONS
def lam_blaze_ech(m, sigma, alpha):
	"Blaze Wavelength of Echelle Grating"
	return 2.0 * sigma * sin(alpha) / m

def lam_blaze_grism(n, sigma, alpha):
	"Blaze Wavelength of Grism Grating"
	return sigma * (n - 1.0) * sin(alpha)

def lam_range(lam_blaze, m):
	"The wavelength limits for each echelle order"
	return lam_blaze / m
#	return (m * lam_blaze) / (m**2 - 0.25)

def redshift(wavelengths, rv=0.0):
    return wavelengths + wavelengths * rv / c

def fabry_perot(wavelengths):
    """
    Fabry Perot etalon
    """
    l = 9.99          # FP length [mm]
    Fin = 8.0         # Finess
    nrefrac = 1.0     # refraction index
    dpha = 4.0 * 3.14159265 / wavelengths * nrefrac * l   # phase
    intensities =  1.0 / (1.0 + Fin*sin(dpha/2.0)**2)
    return intensities
    
def energy2counts(wavelengths, flux):
    """Convert energy per wavelength to counts per wavelength.
    
    Parameters
    ----------
    wavelengths : array_like
        Wavelengths in m.
    flux : array_like
        Energy per wavelength in some form of ergs/s.
        
    Returns
    -------
    array_like
        Counts per wavelength.
    """
    assert wavelengths.shape == flux.shape
    #waves = wavelengths * 1.e-10 # Angstrom -> m
    return flux * wavelengths / h / c
    

# MISCELLANEOUS FUNCTIONS
def find_nearest(array, target):
	"""Finds nearest value in an array. Returns the index and array value"""
	i = argmin(abs(array - target))
	return i, array[i]

if __name__ == "__main__":
    _test_nsell()
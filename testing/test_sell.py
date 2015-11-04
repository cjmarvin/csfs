from numpy import sqrt

def n_sell(arm, lamb):
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

    lamb must be in micrometers
    """
    if arm is "VIS":
        B1 = 1.28035628
        B2 = 0.163505973
        B3 = 0.893930112
        C1 = 0.00929854416
        C2 = 0.0449135769
        C3 = 110.493685
    elif arm is "NIR":
        B1 = 4.76523070e-1
        B2 = 6.27786368e-1
        B3 = 8.72274404e-1
        C1 = 2.84888095e-3
        C2 = 1.18369052e-2
        C3 = 9.56856012e1
    else:
        raise ValueError

    #lamb = lamb * 1.0e-3 # nm to microns
    lamb = lamb * 1.0e3 # mm to microns
    return sqrt( (B1 * lamb**2)/(lamb**2 - C1) + (B2 * lamb**2)/(lamb**2 - C2) + (B3 * lamb**2)/(lamb**2 - C3) + 1.0 )

def _test_nsell():
    print "VIS %.4f = 1.5814 ?" % n_sell("VIS", 587.6e-6) #== 1.58143765521
    print "NIR %.4f = 1.4602 ?" % n_sell("NIR", 546.07e-6) #== 1.46018134651

if __name__ == "__main__":
    _test_nsell()
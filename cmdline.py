import argparse
import os
import defaults
from version import __version__


class CustomFormatter(
    argparse.RawDescriptionHelpFormatter,
    argparse.ArgumentDefaultsHelpFormatter):
    """Keep whitespace in description and epilog, and add default values
    in argument help, respectively.
    """
    def _split_lines(self, text, width):
        """@Mathias: Include whitespace in help.
        """
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


class ArgParser(argparse.ArgumentParser):
    """This class customizes command-line argument behavior.

    Arguments can be loaded from files.
    """
    def __init__(self, *args, **kwargs):
        super(ArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, arg_line):
        """Treat anything after # as a comment."""
        for arg in arg_line.split():
            if not arg.strip():
                continue
            if arg[0] == "#":
                break
            yield arg


# multiline strings
description = """Forward simulates input spectra into a modelled CARMENES
spectrograph, outputting a 2D raw image."""

help_a = [
    "B : bias frame,",
    "C : laser comb calibration spectrum,",
    "F : ideal flatfield exposure (completely normalized continuum spectrum),",
    "L : line list,",
    "    note: input parameter also required, eg.",
    "    '-ainput </path/to/linelist>.txt'",
    "P : example PHOENIX model M dwarf spectrum,",
    "    T_eff=3000 K, [Fe/H]=0.0, log(g)=5.0 (Husser et al., 2013)",
    "S : input object spectrum,",
    "    note: input parameters also required, (see '-ainput') eg.",
    "    '-ainput </path/to/spectrum.fits>' or",
    "    '-ainput </path/to/wavelengths.fits> </path/to/flux.fits>'",
    "T : Thorium Argon calibration lamp line list (Lovis 2007; Kerber 2008),",
    "U : Uranium Neon calibration lamp spectrum (Sarmiento & Reiners 2013),",
    "W : wavemap (pixel values are wavelength instead of counts),",
    "X : Fabry-Perot calibration spectrum,",
    "0 : none\n"]

help_slit = [
    "0 : 0-D point source at absolute slit image center,",
    "1 : 2-D half-moon sampled from an even grid,",
    "2 : 2-D half-moon sampled from random uniform grid,",
    "3 : 2-D half-moon sampled from a Gaussian distribution",
    "    Gaussian depends on the following number of arguments [a,b,c,d] of type <float>:",
    "      0 args -> centered at (0,0), with standard deviation (width) of (1,1);",
    "      1 args -> centered at (0,0), width is argument (a, a);",
    "      2 args -> centered at (a,a), width is (b,b);",
    "      3 args -> centered at (a,b), width is (c,c);",
    "      4 args -> centered at (a,b), width is (c,d);",
    "    (see examples).\n"]

help_ainput = [
    "For fiber argument 'L': the expected parameter is of the form",
    "    '</path/to/linelist>.txt'.",
    "The file must be an ASCII-text file, with no header or footer, and 2 columns:",
    "    wavelengths [angstroms] and flux [arbitrary units].",
    "",
    "For fiber argument 'S': the expected parameters must be one of the following:",
    "    1) '</path/to/spectrum>.fits' or",
    "    2) '</path/to/wavelengths>.fits </path/to/flux>.fits'.",
    "The wavelengths fits-file must be in units of angstroms,",
    "and the flux fits-file can be of any units (arbitrary).\n"]

help_sampling = [
    "cdf  : Full Monte Carlo sampling. Cumulative Distribution Function (CDF)",
    "       for wavelength sampling and Rejection Sampling for slit sampling",
    "rej  : Full Monte Carlo sampling. Rejection Sampling for both",
    "       wavelengths and slit.\n"]

help_example = '' # dummy

help_examples = """
Example Usage
=============

Running CSFS from the command line must follow these basic guidelines:
>>> python csfs.py <spectral-arm> <slit-function> --<sampling-method> -a <a-fiber-type> -b <b-fiber-type>

One sampling method must be specified. For extraction purposes,
we recommend using:
    --monte-carlo
(-sm cdf is optional as it is actually the default setting for
--monte-carlo sampling)

Examples
--------

>>> python csfs.py N 1 -a S -b C -ainput spectrum.fits --grid
Half-moon 2d image, uniformly sampled. Fiber A is a given input spectrum
comprised of "spectrum.fits". Fiber B is a laser comb spectrum. Wavelengths
and fiber samples are on a grid.


>>> python csfs.py V 2 -a S -ainput wavelength.fits flux.fits -ds9 -mc -sm rej
Half-moon 2d image, uniform-randomly sampled. The spectrum is comprised of
"wavelength.fits" and "flux.fits". The FITS image will be opened in SAO-DS9
upon completion. Wavelengths are sampled using the spectrum as a Probability
Distribution Function, using Rejection Sampling to pick random samples.


>>> python csfs.py V 2 -a P -b 0 -noise -tell -blaze -mc -sm cdf
Include noise telluric lines. Also includes the echelle blaze function in the
simulation. Wavelengths are sampled using the spectrum as a Probability
Distribution Function, using its Cumulative Distribution Function (CDF) to pick
random samples.


Gaussian Fiber Examples
-----------------------
>>> python csfs.py V 3 -g -a P
This creates a slit cross section which has a randomly sampled Gaussian
distribution centered at (0,0) and a width of 1 in both the x and y direction.


>>>  python csfs.py V 3 0.5 -g -a P -g
The "0.5" argument given after the 3 indicates the width (standard deviation)
of the Gaussian in both the x and y direction. It should be a positive, real
number.


>>>  python csfs.py V 3 -0.5 0.7 -g -a P -g
The first argument here is the offset location for (x,y) in terms of the unit
circle. This means that the center (mean) of the Gaussian will be at
(-0.5, -0.5). The second argument 0.7 is the width in both x and y directions.


>>>  python csfs.py V 3 -0.5 0.2 0.7 -g -a P -g
Here, the center of the distribution is the first two arguments, so that the
central location is (x,y)=(-0.5, 0.2). The width in both x- and y- directions
is 0.7.


>>>  python csfs.py V 3 -0.5 0.2 0.7 0.3 -g -a P -g
The last of the possibilities is with 4 following arguments. The first two
define the central location (x,y)=(-0.5, 0.2) The third, 0.7, defines the
width in the x-direction. The fourth, 0.3, defines the width in the y-direction.


>>> python csfs.py V 3 -g -a P --preview
Instead of fully running the simulation, Images of the slit cross section will
pop up to show the slit image, and the program will exit.
The '--preview' option has been included to help visualize and test different
arguments.
"""

help_examples_grid = """Example Usage:"""

help_examples_mc = """Example Usages:"""


def parse_command_line():
    """
    Command line control of CSFS.

    """
    parser = ArgParser(
        prog="CARMENES Spectrum Forward Simulator",
        description=description,
        formatter_class=CustomFormatter,
        fromfile_prefix_chars='@',
        add_help=True,
        epilog=help_examples)

    parser.add_argument('-v', '--version',
        action="version",
        version='%(prog)s v{version}'.format(version=__version__))

    # ROOT PROGRAM ARGUMENTS
    parser.add_argument('arm',
        action='store',
        choices="N V".split(),
        help="""Spectral arm: {'N': NIR, 'V': VIS}.""")

    parser.add_argument('slit',
        action='store',
        nargs='+',
        # default=False,
        # metavar='', # this causes an AssertionError! an argparse bug...
        help="R|Slit function/PSF/cross-section:\n%s" % ('\n'.join('%s' % line for line in help_slit)))

    # =========================================================================

    # REQUIRE SA SAMPLING METHOD
    megroup_sampling = parser.add_mutually_exclusive_group(required=True)

    # SUBPROGRAM GRID
    # -------------------------------------------------------------------------
    megroup_sampling.add_argument("-g", "--grid",
        dest="grid",
        action="store_true",
        default=False,
        help="""Grid sampling of ray wavelengths and slit positions.""")

    group_grid = parser.add_argument_group(
        description='Grid sampling arguments.',
        title="Options work ONLY if '-g, --grid' sampling is selected.")

    group_grid.add_argument("-p", "--perturb",
        action="store_true",
        default=False,
        dest="perturb",
        help="""Perturb slit locations for a quasi-random simulation.
             This helps to avoid aliasing artifacts in the extraction of the image.
             Only works with '-g, --grid' option.""")

    group_grid.add_argument('--maxsnr',
        action='store',
        default=defaults.DMAX_SNR,
        type=float,
        metavar='[DN]',
        help="""Maximum Signal-to-Noise Ratio of simulated image
        (before noise). Renormalizes image to this value.
        Only works with '-g, --grid' option.""")

    group_grid.add_argument('-ns',
        action='store',
        dest='ns',
        default=defaults.DEFAULT_NS,
        type=int,
        help="""R|Number of samples to create the slit function grid before
the cross-section filter.
For slit=1, ns is defined along the fiber diameter.
The effective number of slit samples will be ~0.8*ns.
ns should be >= 60.

For slit=2, ns is defined along the fiber radius.
The number of effective slit samples will be ~Pi*ns^2. ns should be >= 10,000.
Only works with '-g, --grid' option.\n""")


    # -------------------------------------------------------------------------

    # SUBPROGRAM MONTE CARLO
    # -------------------------------------------------------------------------
    megroup_sampling.add_argument("-mc", "--monte-carlo",
        action="store_true",
        dest="mc",
        default=False,
        help="""Monte Carlo (MC) sampling of ray wavelengths and slit positions.""")

    group_mc = parser.add_argument_group(
        description="Options work ONLY if '-mc, --monte-carlo' sampling is selected",
        title='Monte Carlo sampling arguments.')

    group_mc.add_argument('-sm', '--mc-method',
        action="store",
        default='cdf',
        dest="sm",
        choices="cdf rej".split(),
        help=""""R|Method for MC wavelength sampling:\n%s.
Only works with '-mc, --monte-carlo' option.\n""" % ('\n'.join('%s' % line for line in help_sampling)))

    group_mc.add_argument("-nr", "--nrays",
        dest="nr",
        type=float,
        default=int(1e8),
        metavar='',
        nargs='+',
        help="""R|Approximate number of rays to simulate.
If 2 numbers are specified, then the first input goes to Fiber A,
and the second goes to Fiber B.
Due to order overlapping of wavelengths, this parameter is
only a lower-limit approximation.
Only works with '-mc, --monte-carlo' option.\n""")

    group_mc.add_argument('--interp-flag',
        action='store',
        default=False,
        type=int,
        help="""R|Interpolation flag for CDF function.\n0 : do not interpolate\n3 : cubic interpolation.""")

    # -------------------------------------------------------------------------

    # COMMON PARSER ARGUMENTS

    # FIBER OPTIONS
    group_fiber = parser.add_argument_group(title='Fiber arguments')

    group_fiber.add_argument('-a',
        action='store',
        default='0',
        #nargs='*',
        #required=True,
        choices='B C F L P S T U W X 0'.split(),
        metavar='',
        help="R|Fiber A input:\n%s" % ('\n'.join('%s' % line for line in help_a)))

    group_fiber.add_argument('-b',
        action='store',
        default='0',
        #nargs='*',
        #required=True,
        choices='B C F L P S T U W X 0'.split(),
        metavar='',
        help="""R|Fiber B input. Similar to Fiber A. See '-a'.\n""")

    group_fiber.add_argument('-A', '--ainput',
        dest='afile',
        action='store',
        nargs='+',
        metavar='',
        help="R|Input file(s) for Fiber A.\n%s" % ('\n'.join('%s' % line for line in help_ainput)))

    group_fiber.add_argument('-B', '--binput',
        dest='bfile',
        action='store',
        nargs='+',
        metavar='',
        help="""R|Input file(s) for Fiber B. See '-ainput'.\n""")

    group_fiber.add_argument('-rva',
        action='store',
        type=float,
        default=0.0,
        metavar='',
        help="""R|Radial velocity shift of Fiber A in [m/s].\n""")

    group_fiber.add_argument('-rvb',
        action='store',
        type=float,
        default=0.0,
        metavar='',
        help="""R|Radial velocity shift of Fiber B in [m/s].\n""")

    # deprecated
    #parser.add_argument('-wlog',
        #action='store',
        #default=False,
        #type=float,
        #help="""Sampling of wavelength grid in logarithmic space.""")

    # WAVELENGTH TRACE OPTIONS
    group_wtrace = parser.add_argument_group(title='Wavelength trace arguments',
        description="WARNING: Not currently supported!")

    group_wtrace.add_argument('-wta',
        action='store_true',
        default=False,
        # metavar='',
        help="""Create a wavelength region file (.reg) for Fiber A. The region
        file is a text file which includes the detector-positions of
        wavelengths. When opened in a FITS viewer with its associated FITS
        image, certain wavelengths are labeled at their respective
        positions.""")

    group_wtrace.add_argument('-wtb',
        action='store_true',
        default=False,
        # metavar='',
        help="""Create a wavelength region file (.reg) for Fiber B.
        See '-wta'.""")

    group_wtrace.add_argument('--wtlist',
        action='store',
        metavar='',
        help="""Input file of wavelengths for wavelength trace region file.
        Must be of the form '</path/to/wavelengthlist>.txt'.
        The file is an ASCII txt file, with 1 column
        (wavelengths [Angstrom]).""")

    group_wtrace.add_argument('-dwt',
        action='store',
        type=float,
        default=defaults.DEFAULT_DWT,
        metavar='',
        help="""Wavelength sampling [Angstrom] for wavelength trace region
        file.""")

    # SPECTRAL OPTIONS
    group_spectral = parser.add_argument_group(title='Spectral arguments')

    group_grid.add_argument('-dw',
        action='store',
        dest='nw',
        default=defaults.DEFAULT_DW,
        type=float,
        help="""Sampling of wavelength grid [Angstrom]. This defines the
        sampling resolution of functional fiber inputs {C, F U, W, X}.
        With the '-g, --grid' option, this also regrids the input wavelengths
        to the specified resolution.""")

    group_spectral.add_argument('--oset',
        action='store',
        nargs='+',
        type=int,
        metavar='',
        help="""Specify a list of orders to use in the simulation.
        The list can consist of either all integers, or can
        be 3 integers as follows: (min, max-1, increment).
        The latter is similar to the python/numpy range/arange function.
        Eg., '-oset 59 60 61 62 63' will simulate orders 59-63, whereas
        '-oset 59 80 ' will only use orders 59, 62, 65, 68, 71, 74, 77 (80 is
        not included because of the Python style).""")

    #If no argument is
    #specified, the gap will be the default value in 'armparams.h'
    #centered at 0.
    #If 1 argument is given in [mm], this will be the offset with respect
    #to the center at 0.
    group_spectral.add_argument('--gap',
        action="store_true",
        default=False,
        help="""Include gap between the two NIR detectors [2.54 mm].
        If spectral arm is VIS, this will have no effect and default to False.
        """)

    group_spectral.add_argument('--blaze',
        action='store_true',
        default=False,
        help="""R|Include echelle blaze function in computation.\n""")

    group_spectral.add_argument('--noise',
        const=True,
        default=False,
        nargs='?',
        metavar='DN',
        help="""Include noise.
        Currently shot noise, readout noise, and a bias of %s [DN] are included.
        A value for the readout noise can also be given as an argument.""" % defaults.sig_rn)

    group_spectral.add_argument('--tell',
        #action='store_const',
        choices=defaults.TELLURIC_LINES.split(),
        default=False,
        nargs='*',
        metavar=defaults.TELLURIC_LINES,
        help="""Include telluric lines computed using LBLRTM (Line-By-Line
        Radiative Transfer Model, based on Clough et al. 2005)
        Independent species can also be given as an argument.""")

    group_spectral.add_argument('--crop',
        action="store",
        default=0,
        type=int,
        help="""Number of pixels to crop on border of CCD. These cropped pixels
        will only be affected by readout noise if noise is selected, otherwise
        they will be set to 0. Must be an integer value. Only applies to NIR arm.""")

    # OUTPUT OPTIONS
    group_output = parser.add_argument_group(title='Output arguments')

    group_output.add_argument('--no-compress',
        dest="no_compress",
        action='store_true',
        default=False,
        help="""Do not compress output FITS file to .gz format.""")

    group_output.add_argument('--subwindows-ext',
        dest="subwindows_ext",
        action='store_true',
        default=False,
        help="""Include creation of SUBWINDOWS extension. Does not apply to VIS arm.""")

    group_output.add_argument('--single-ext',
        action='store_true',
        default=False,
        help="""Suppress creation of TOM and JERRY extensions, outputting to a single Primary HDU.
        Only applies to NIR arm.""")

    group_output.add_argument('--out',
        action='store',
        metavar='output/filename/filename.fits',
        help="""Specify a name for FITS file (and region file if it applies).
        The default format is '<date><time><a_type><b_type>.fits.gz'.
        Output is saved to a folder of this name in 'output/' unless another
        path is specifed. If --no-full-output is passed, then only a fits image
        will save to 'output/', with no folder created.""")

    group_output.add_argument('--full-output',
        action='store_true',
        default=False,
        dest="full_output",
        help="""Outputs auxiliary files (wavelengths.dat & slit_cross_section.pdf""")

    group_output.add_argument('--path',
        action='store',
        # default=os.path.join(os.getcwd(), "output"),
        default="output/",
        metavar='/path/to/output/',
        help="""Specify the output path.""")

    group_output.add_argument('--ds9',
        action='store_true',
        default=False,
        help="""Open FITS output file (and associated region file) in SAO-DS9
        after simulation finishes.""")

    group_output.add_argument('--ds9-mosaic',
        action='store_true',
        default=False,
        help="""Open FITS output file (and associated region file) in SAO-DS9
        after simulation finishes with --mosaicimage iraf option.
        This will merge all extensions into one.""")

    # OTHER OPTIONS
    group_other = parser.add_argument_group(title='Other arguments')

    group_other.add_argument('--reset',
        action='store_true',
        default=False,
        help="""Delete and re-solve all existing offset data, wavelength data,
        and wavelength solutions (equivalent to "rm DATA/*").
        NOTE: Usage was mostly for debugging purposes, and should be
        deprecated in favor of 'make clean' in the shell.""")

    group_other.add_argument('--solve',
        dest='solve',
        action='store_true',
        default=False,
        help="""For the current run, re-solve all existing offset data,
        wavelength data, and wavelength solutions instead of importing data.
        Mostly used for debugging purposes.""")

    group_other.add_argument('--debug',
        action='store_true',
        default=False,
        help="""Debug mode. Used for checking the optical model - not for
        debugging the program. CURRENTLY UNDER CONSTRUCTION.""")

    group_other.add_argument('--preview',
        action='store_true',
        default=False,
        help="""Preview fiber illumination and exit. Use only with '-g/--grid'
        option.""")

    args = parser.parse_args()
    return args

# for testing purposes
if __name__ == "__main__":
    args = parse_command_line()
    print args
    for k, v in args.__dict__.iteritems():
        print(k, v)

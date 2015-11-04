# __CARMENES Spectrum Forward Simulator__#
__Version 5.6.0__

![VIS Channel](http://carmenes.caha.es/ext/gallery/instrument/CARMENES-20131002-csfs-rainbowfibreanoaxis.png)

#Prerequisites/Requirements
- GCC
- GNU Scientific Library [(GSL)](http://www.gnu.org/software/gsl)
- GNU Make
- SAO-DS9 (for viewing output __FITS__)
- Python 2.6 or 2.7
    + NumPy 1.7 or later
    + Scipy 0.12 or later
    + Matplotlib 1.1 or later
    + Astropy 0.3 (or Pyfits 2.0) or later

In order to facilitate the installation of required Python packages,
an installation of the [Anaconda Python Distribution](https://store.continuum.io/cshop/anaconda/ "Anaconda Python Distribution")
is highly recommended, as it contains all required packages.
The full version is completely free, and installation does not require root priveledges.

# __Download__
You can get the bleeding edge version using Git:

    git clone https://cmarvin@bitbucket.org/cmarvin/csfs.git <name-of-destination>


Or go to the download page:
[https://bitbucket.org/cmarvin/csfs/downloads](https://bitbucket.org/cmarvin/csfs/downloads).

## Update
If using Git, you can easily update by

    git pull

# __Setup__

## Compiling
Before running, the C source code must be compiled.
This can be done simply by typing in the terminal:

	make

Any changes made in __armparams.h__ or __functions.c__ code requires
re-compilation by running `make` again.

## Cleaning the CSFS directory

	make clean

This command will remove all hidden **Mac OS X** files/directories,
and files with the following extensions:

+ __so__
+ __pyc__
+ __~__

Files in the __DATA__ and __FITS__ directories will be left untouched.

Once compiled, __csfs.py__ should run normally.

# __Examples__

## get help / list options

`python csfs.py -h` or `python csfs.py --help`


## debug plot _*currently under construction*_

The debug plot is a way to cross-check the simulation using certain
wavelengths and slit positions.
This avoids going through the process of producing a FITS file.
It can be used to create plots of the slit function, as well as footprints of
the spectral format.

	python csfs.py V 1 -debug -solve

# sampling methods

At least one sampling method must be specified. For extraction purpose,
it is recommended to use *monte-carlo* sampling.
The *grid* option has been kept for archival reasons.

## grid

This is the original sampling implementation.
It samples wavelengths on an even grid.
Samples from the fiber PSF are either chosen by using the `slit` argument as follows:

- `0` = point source
- `1` = even grid
- `2` = uniform random distribution
- `3` = gaussian random distribution

Select it by entering either `-g` or `--grid`.
Since the fiber image in this case is static, the locations can be perturbed
for each ray by providing the `--perturb` argument. This semi-randomizes the
simulation, or makes it quasi-Monte Carlo.

## Monte Carlo (MC)

Monte Carlo sampling involves picking samples at random.
Wavelength samples are chosen by using the spectrum as a Probability
Distribution Function (PDF).
Fiber PSF samples are chosen at random as the `slit` argument, according to a distribution:

- `0`, `1`, `2` = uniform random distribution
- `3` = gaussian random distribution

Select it by entering either `-mc` or `--monte-carlo`.
The different MC sampling methods are chosen using `-sm` or `--mc-method`:
- `cdf` = Cumulative Distribution Function (CDF)
- `rej` = Rejection Sampling

## input spectra

Two input FITS files are needed for the *input spectra* `S` option.
They are specified with the `-A/--ainput` command (or `-B/--binput` for fiber B):

	python csfs.py V 1 -monte--carlo -a S --ainput /path/to/wavelength.fits /path/to/flux.fits --ds9

The `--ds9` option opens the output FITS file in SAO-DS9.
*__TIP__*
To keep the command-line short, you can specify your input spectra files as
variables in the Unix Shell:

	IN1="/path/to/wavelength.fits"
	IN2="/path/to/flux.fits"

Then they will save space on the command-line:
	
	python csfs.py V 1 -mc -a S --ainput $IN1 $IN2

## telluric lines

Telluric lines can be added via the *input spectra* `S` option.
The option `--tell` can be given which will include all species in the telluric
absorption.

	python csfs.py V 1 -mc -a S -A $IN1 $IN2 --tell

Also, isolated species can be specified as optional arguments.
In this example, only H2O and O3 lines are included in the simulation.

	python csfs.py V 1 -mc -a S -A $IN1 $IN2 --tell H2O O3

## noise

Shot noise and readout noise can also be added to the simulated image.

	python csfs.py V 1 -g --grid -a S -A $IN1 $IN2 --noise

The default rms readout noise value can be found/changed in `defaults.py`.
Also, an rms readout noise value can be specified as an optional argument.

	python csfs.py V 1 -g --monte-carlo -a S -A $IN1 $IN2 --noise 30

Note that shot noise is automatically included if the `-mc/--monte-carlo`
is specified.


## command line files

Sometimes the command line input can be long and tiring to type.
Therefore, you can create a text file with line-separated arguments, and
specify `@`+<file-name> as an argument itself.
For example, if the text below is saved as _ex.args_, then you would simply
call

	python csfs.py @ex.args
to save a lot of command line typing.

	V
	3
	-a P
	#-b T
	--monte-carlo
	--ds9
	--noise
	--tell
	--blaze
	--no-full-output
	-nr 1e9
	
You can also use `#` to comment out lines.

## gaussian fiber cross-section

	python csfs.py V 3 -g -a P
This creates a slit cross section which has a randomly sampled Gaussian distribution centered at (0,0) and a width of 1 in both the x and y direction.

	python csfs.py V 3 0.5 -g -a P
The `0.5` argument given after the `3` indicates the width (standard deviation) of the Gaussian in both the x and y direction. It should be a positive, real number.

	python csfs.py V 3 -0.5 0.7 -g -a P
The first argument here is the offset location for (x,y) in terms of the unit circle. This means that the center (mean) of the Gaussian will be at (-0.5, -0.5). The second argument `0.7` is the width in both x and y directions.

	python csfs.py V 3 -0.5 0.2 0.7 -g -a P
Here, the center of the distribution is the first two arguments, so that the central location is (x,y)=(-0.5, 0.2). The width in both x- and y- directions is 0.7.

	python csfs.py V 3 -0.5 0.2 0.7 0.3 -g -a P -g
The last of the possibilities is with 4 following arguments. The first two define the central location (x,y)=(-0.5, 0.2) The third, `0.7`, defines the width in the x-direction. The fourth, `0.3`, defines the width in the y-direction.

	python csfs.py V 3 -g -a P --preview

Instead of fully running the simulation, Images of the slit cross section will pop up to show the slit image, and the program will exit.
The `--preview` option has been included to help visualize and test different arguments.



To-Do List
----------

- ~~finish telluric line implementation~~
- ~~_SNR_ option~~
- ~~export of data and plots~~ more or less
- ~~change output to folders of data products~~
- ~~cubic spline interpolation of wavelength grid~~
- cleaning of code; structure similar to CriForSS
- ~~_UNe_ calibration spectrum~~
- ~~add shot noise~~
- ~~add readout noise~~
- add dark current
- ~~add intialize or setup feature to unzip included telluric line spectra and convert to binary _npy_ files implemented in Makefile~~
- NIR arm option for 2 CCDs (with _gap_ option)


Changes from v_3.1.2
--------------------
- NIR arm implemented with Sellmeier equation
- compute_c changed to a general function, encompassing wavetrace, offset, etc arguments changed as well (order to be more clear, ARM flag added)
- rv_shift function corrected
- debug option added
- SpectralArm object now passed as a parameter in functions
- cmdline.py module added, separating command line parsing as a separate module
  1) *resetall* clears wavelimits, and offset (erases DATA folder)
  2) *debug* plots, print
  3) *solve* re-solves all data without erasing (similar to reset without erasing DATA)

- Makefile added
- armparams.h replaces vis.h
    all NIR and VIS specific parameters are given here
- defaults.py added
    all other simulation parameters reside here
- wavefuncs.py added
    wavelength methods from spectralarm.py are moved here
- writefits.py updated
    change in pyfits header convention
    truncated card fixed
- noise.py created but not implemented
- wavelengthtrace.py order overlap bug fixed

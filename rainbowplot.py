"""rainbowplot.py
This is a self-running script to produce a colorized echellogram, hence
the name 'rainbow' plot.
In order to create this, 2 images are needed - the spectrum you want
to colorize, and a wavemap which provides the color information.
You can either produce these images on-the-fly or specify an image to use.

Usage:
  rainbowplot.py [--spectrum=SPECTRUM] [--wavemap=WAVEMAP] [-hvpwo OUTFILE]
                 [--cmin=CMIN --cmax=CMAX] [--scale=SCALE] [--lowpass=LOWPASS]
                 [--highpass=HIGHPASS]
  rainbowplot.py --version

Options:
  -h --help                    show this help message and exit
  -v --version                 show version and exit
  --spectrum=SPECTRUM          input spectrum or CSFS script
                               [default: rainbow/default_rainbow_spectrum_pho-fab.fits.gz]
  --wavemap=WAVEMAP            input wavemap or CSFS wavemap script
                               [default: rainbow/default_viswavemap.fits.gz]
  -o OUTFILE --output=OUTFILE  output filename
                               [default: rainbow/rainbowplot.pdf]
  -p --plot                    open plot
  -w --write                   save data array to disk
  -s SCALE --scale=SCALE       scaling function of image intensity. can be
                               linear, log, or p<N> (for power)
                               [default: linear]
  --cmin=CMIN                  color of longest wavelength,
                               can be from 0-360 (hue scale)
                               or a color {red, yellow, green, cyan, blue,
                               magenta} [default: red]
  --cmax=CMAX                  color of shortest wavelength,
                               can be from 0-360 (hue scale)
                               or a color {red, yellow, green, cyan, blue,
                               magenta} [default: blue]
  --highpass=HIGHPASS          minimum percentage value to scale [default: 0.0]
  --lowpass=LOWPASS            maximum percentage value to scale [default: 1.0]

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import os
import shutil
import subprocess
from docopt import docopt
import urllib

# ddir = "rainbow"
ddir = "RAINBOW"
_url = "http://www.astro.physik.uni-goettingen.de/~cmarvin/csfs/library/rainbow"


def rgba_plot():
    colors_scaled = colors_ma / colors_ma.max()
    kwargs = {'interpolation': 'none', 'vmin': colors_ma.min()/colors_ma.max(), 'vmax': 1.0, 'origin': 'lower'}
    RGBA = plt.cm.jet(colors)
    RGBA[..., 3] = -intensity + 1.0

    plt.imshow(RGBA, origin='lower')
    print colors_scaled.min()
    print np.where(RGBA == colors_scaled.min())
    plt.colorbar()
    plt.show()


def hsv_plot(args):
    fac = 1.0
    pltkwargs = {
        'interpolation' : 'none',
        'vmin' : 0.0,
        'vmax' : 1.0,
        'origin' : 'lower',
        'aspect' : 'equal'
        }

    hue = {
        'red' :     0.0,
        'yellow' :  60.0,
        'green' :   120.0,
        'cyan' :    180.0,
        'blue' :    240.0,
        'magenta' : 300.0
        }

    if 'python' in args['--spectrum']:
        print 'executing %s' % args['--spectrum']
        subprocess.check_call(('%s --path rainbow --out tmpspectrum' % args['--spectrum']).split())
        spectrum = 'rainbow/tmpspectrum/tmpspectrum.fits.gz'
    else:
        spectrum = args['--spectrum']

    if 'python' in args['--wavemap']:
        print 'executing %s' % args['--wavemap']
        subprocess.check_call(('%s --path rainbow --out tmpwavemap' % args['--wavemap']).split())
        wavemap = 'rainbow/tmpwavemap/tmpwavemap.fits.gz'
    else:
        wavemap = args['--wavemap']

    try:
        cmin = float([args['--cmin']])
    except TypeError:
        cmin = hue[args['--cmin']]
    try:
        cmax = float([args['--cmax']])
    except TypeError:
        cmax = hue[args['--cmax']]

    if not os.path.exists('rainbow'):
        os.makedirs('rainbow')

    if not os.path.isfile(spectrum):
        s_url = os.path.join(_url, os.path.basename(spectrum))
        print "{} not found.".format(spectrum)
        print "Downloading {}".format(s_url)
        urllib.urlretrieve(s_url, spectrum)

    if not os.path.isfile(wavemap):
        s_url = os.path.join(_url, os.path.basename(wavemap))
        print "{} not found.".format(wavemap)
        print "Downloading {}".format(s_url)
        urllib.urlretrieve(s_url, wavemap)

    print "Loading spectrum file: %s" % spectrum
    intensity = pyfits.getdata(spectrum)
    print "Loading wavemap file: %s" % wavemap
    colors = pyfits.getdata(wavemap)

    print 'scaling intensity'
    _scale_func = {
        'log'  : np.log10,
        'sqrt' : np.sqrt,
        }
    vmin = intensity[intensity > 0].min() * float(args['--highpass'])
    vmax = intensity.max() * float(args['--lowpass'])
    intensity[intensity <= vmin] = 0
    intensity[intensity >= vmax] = intensity.max()
    if args['--scale'] != 'linear':
        if args['--scale'].startswith('p'):
            p = float(args['--scale'][1:])
            intensity = intensity**p
        else:
            intensity = _scale_func[args['--scale']](intensity)
    intensity /= intensity.max() # normalize to 1
    mask = ~(colors > 0)
    colors_ma = np.ma.array(colors, mask=mask)

    print 'scaling colors'
    colors_scaled = (colors_ma - colors_ma.max() ) / (colors_ma.min() - colors_ma.max())
    colors_scaled *= (cmax - cmin) / 360.
    HSV = np.ones(intensity.shape + (3,))
    # print colors_scaled.min(), colors_scaled.max()
    HSV[..., 0] = colors_scaled
    HSV[..., 2] = intensity

    print 'creating image'
    HSV = matplotlib.colors.hsv_to_rgb(HSV)
    rgb_text = 'R G B A'.split()
    hsv_text = 'H S V'.split()
    # write the array to disk
    if args['--write']:
        i = 0
        with file('rainbow/rainbow_array.txt', 'w') as outfile:
            outfile.write('# Array shape: {0}\n'.format(HSV.shape))
            # Iterating through a ndimensional array produces slices along
            # the last axis. This is equivalent to data[i,:,:] in this case
            for data_slice in HSV.reshape(3,4096,4096):
                # Writing out a break to indicate different slices...
                print 'writing %s array to disk...' % hsvstr[i]
                outfile.write('#%s\n' % rgb_text[i])
                # The formatting string indicates that I'm writing out
                # the values in left-justified columns 7 characters in width
                # with 2 decimal places.
                np.savetxt(outfile, data_slice)#, fmt='%-7.2f')
                i += 1
    # set up plot
    print "Setting up plot"
    fig = plt.figure(frameon=False)
    fig.set_size_inches(8, 8)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    print "Opening plot"
    ax.imshow(HSV, **pltkwargs)#, **kwargs)#, origin='lower', cmap=plt.cm.jet)
    if args['--plot']:
        plt.show()
    if args['--output']:
        print "Saving plot as %s" % args['--output']
        fig.savefig(args['--output'], dpi=300)#, bbox_inches='tight', pad_inches=0)

    # REMOVE TMP
    print "Cleaning rainbow/tmp directories"
    shutil.rmtree("rainbow/tmpspectrum", ignore_errors=True)
    shutil.rmtree("rainbow/tmpwavemap", ignore_errors=True)

if __name__ == '__main__':
    args = docopt(__doc__)
    # print args
    hsv_plot(args)

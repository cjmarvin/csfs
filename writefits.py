"""
writefits.py

09-04-2014
----------
Image headers according to
    Document: FDR-11A
    Issue number: 1.6
    Issue date: 01 Feb 2013

14-08-2013
----------
 * header.update() method changed to header['keyword'] = (value, comment)
 * Line 50: 'Reference of data associated w/ header' changed to 'Reference where the data associated with the header are published'
 * Line 60: 'Organization responsible for creating FITS file' changed to 'Identifies the organization responsible for creating the FITS file'

"""

import numpy as np
# import matplotlib.pyplot as plt
import os
import subprocess
import time
try:
    import pyfits
except ImportError:
    import astropy.io.fits as pyfits
import astrodate
import noise
from version import __version__


def recurs_makedirs(out_dir):
    """
    Makes a directory.

    If out_dir exists, creates out_dir+1 so that the original out_dir is not
    overwritten.
    """
    try:
        os.makedirs(out_dir)
    except OSError:
        if out_dir[-1].isdigit():
            out_dir = out_dir[:-1] + str(int(out_dir[-1]) + 1)
        else:
            out_dir += "1"
        return recurs_makedirs(out_dir)
    else:
        return out_dir


def setup_output(self):
    """
    Sets up output folder.

    Create necessary folder. Write slit function array.
    """

    # CREATE OUTPUT FOLDER; RENAME IF NECESSARY
    if self.outfile:
        outdirname = self.outfile
        outfile = self.outfile
    else:
        _out_type = "%.3s%.3s" % (
            self.simmodetypes[0].replace(" ", "").lower(),
            self.simmodetypes[1].replace(" ", "").lower())
        if 'BIAS' in self.simmodetypes:
            _out_type = 'bias'
        outdirname = 'csfs_%s_%.3s_%.6s' % (
            time.strftime("%Y_%m_%dT%H_%M_%S"),
            self.ARM.lower(),
            _out_type)
        outfile = 'csfs.%s-%s-%.6s' % (
            time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()),
            self.ARM.lower(),
            _out_type)
    outdir = os.path.join(self.outpath, outdirname) # add path

    if self.full_output:
        self.outdir = recurs_makedirs(outdir) # makedir and rename (+1) if necessary
        # self.outfile = os.path.join(outdir, outfile)
        self.outfile = os.path.join(self.outpath, outfile)
        # SAVE DATA
        # write slit
        if self.SAMPLING == "grid":
            headerline = "x[mm] y[mm]".split()
            headerfmt = "#{:<18} {:<18}\n"
            datafmt = "{:<19} {:<18}\n"
            out_slit_fn = os.path.join(self.outdir, "slit_cross_section.dat")
            data = self.slit_cross_section
            with open(out_slit_fn, "w") as f:
                    f.write(headerfmt.format(*headerline))
                    for line in data:
                        f.write(datafmt.format(*line))
            # write slit image
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect="equal")
            ax.scatter(self.slit_cross_section[:,0], self.slit_cross_section[:,1], alpha=0.5, edgecolors="none")
            ax.set_xlabel(r"$x$ [mm]")
            ax.set_ylabel(r"$y$ [mm]")
            ax.set_title("Fiber Cross-section")
            plt.savefig(os.path.join(self.outdir, "slit_cross_section.pdf"))

        # write wavelengths
        headerline = "ORDER CCD_MIN[A] FSR_MIN[A] BLAZE[A] FSR_MAX[A] CCD_MAX[A]".split()
        headerfmt = "#{:<7} {:<15} {:<18} {:<18} {:<18} {:<18}\n"
        datafmt = "{:<8} {:<15} {:<18} {:<18} {:<18} {:<18}\n"
        out_waves_fn = os.path.join(self.outdir, "wavelengths.dat")
        data = zip(
            self.ORDERS,
            self.ccd_waves_min*1.e7,
            self.fsr_waves_min*1.e7,
            self.bl_ech*1.e7,
            self.fsr_waves_max*1.e7,
            self.ccd_waves_max*1.e7)
        with open(out_waves_fn, "w") as f:
            f.write(headerfmt.format(*headerline))
            for line in data:
                f.write(datafmt.format(*line))
    else:
        self.outfile = os.path.join(self.outpath, outfile)
    # print self.outfile
    return 0


def write_to_fits(self, gzip=True):
    """
    Write simulation results to fits file.

    Write header and data to FITS HDU. Compress if necessary.
    """

    # create folders and files
    setup_output(self)

    # FOR CONVENIENCE AND CLARITY, ASSIGN HEADER VARIABLES UP HERE
    DATE = time.strftime("%Y-%m-%d", time.localtime())
    FILENAME = os.path.basename(self.outfile)
    OBJECT = self.sciobject
    if OBJECT is None:
        OBJECT = ("%s,%s") % tuple(self.fib_src)
    CATG = self.catg
    TECH = self.tech
    TYPE = ("%s,%s") % tuple(self.fib_obstype)
    SOURCE = ("%s,%s") % tuple(self.fib_src)
    DATE_OBS = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime())
    MJD_OBS = astrodate.AstroDate().mjd
    EFFTIME = DATE_OBS
    CCDPSIZ = self.DPIX * 1.e3 # [microns]
    CCDMEAN = self.bias # [DN]
    CCDRON = self.sig_rn # [e-]
    ORIGSECX = self.NXPIX
    ORIGSECY = self.NYPIX

    # SIMU keywords
    SMETH = self.SAMPLING
    PERTURB = self.perturb
    MCMETH = self.sm

    FIBA = self.fib_simmode[0]
    FIBB = self.fib_simmode[1]
    FIBA_SRC = self.fib_src[0]
    FIBB_SRC = self.fib_src[1]
    INFILEA1 = None
    INFILEA2 = None
    INFILEB1 = None
    INFILEB2 = None
    if self.infiles[0]:
        if len(self.infiles[0]) == 2:
            INFILEA1 = os.path.basename(self.infiles[0][0])
            INFILEA2 = os.path.basename(self.infiles[0][1])
        else:
            INFILEA1 = os.path.basename(self.infiles[0][0])
    if self.infiles[1]:
        if len(self.infiles[1]) == 2:
            INFILEB1 = os.path.basename(self.infiles[1][0])
            INFILEB2 = os.path.basename(self.infiles[1][1])
        else:
            INFILEB1 = os.path.basename(self.infiles[1][0])
    RVA = self.fib_rv[0]
    RVB = self.fib_rv[1]
    SLITTYP = self.slittyp
    if self.SAMPLING == "grid":
        NS_EFF = self.fib_nslit[0]
        NS_IN = self.ns
    if self.SAMPLING == "mc":
        NRAYS_IN = self.nr
        NRAYS_INMA = self.fib_nrays[0]
        NRAYS_INMB = self.fib_nrays[1]
        MEAN_NRAYS_A = self.fib_mean_rays_per_pixel[0]
        MEAN_NRAYS_B = self.fib_mean_rays_per_pixel[1]
        MIN_NRAYS_A = self.fib_min_rays_per_pixel[0]
        MIN_NRAYS_B = self.fib_min_rays_per_pixel[1]
        MAX_NRAYS_A = self.fib_max_rays_per_pixel[0]
        MAX_NRAYS_B = self.fib_max_rays_per_pixel[1]
        TOT_NRAYS_A = self.fib_nrays_tot[0]
        TOT_NRAYS_B = self.fib_nrays_tot[1]
    DW = self.nw
    OMIN = self.ORDERS.min()
    OMAX = self.ORDERS.max()
    NIRGAP = self.gap
    BLAZE = self.blaze
    PN = self.noise if self.SAMPLING == "grid" else True
    BIAS = self.bias
    RON = self.sig_rn
    if self.tell:
        CH4 = 'CH4' in self.tell
        CO2 = 'CO2' in self.tell
        H2O = 'H2O' in self.tell
        N2O = 'N2O' in self.tell
        O2 = 'O2' in self.tell
        O3 = 'O3' in self.tell
    else:
        CH4 = False
        CO2 = False
        H2O = False
        N2O = False
        O2 = False
        O3 = False
    TIME_A = np.float("%.2f" % self.fib_sim_time[0])
    TIME_B = np.float("%.2f" % self.fib_sim_time[1])

    # ================== FITS FUNCTIONS ==================================

    # create PrimaryHDU object to encapsulate data
    if self.ARM == 'VIS' or self.single_ext:
        hdu = pyfits.PrimaryHDU(data=self.image)
    elif self.ARM == 'NIR':
        hdu = pyfits.PrimaryHDU()
    header = hdu.header

    # MANDATORY KEYWORDS
    header['PCCOUNT'] = ('0', 'Number of pixels following the data in the heap')
    header['GCOUNT'] = (1, 'Number of groups')

    # GENERAL KEYWORDS
    #header['CHANNEL'] = (self.ARM, 'spectral arm')
    header['DATE'] = (DATE, 'Date on which the HDU was created')
    header['FILENAME'] = (FILENAME, 'Original name of this file')
    header['ICSVER'] = ('1.0', 'ICS version')
    header['PHASE'] = ('SIMULATION', 'CARMENES phase')
    header['AUTHOR'] = ('Quirrenbach et al.', 'Reference author')
    header['REFERENC'] = ('2012SPIE.8446E..0RQ', 'Reference of data associated w/ header')#'Reference where the data associated with the header are published'
    header['OBSERVAT'] = ('CAHA', 'Centro Astronomico Hispano Aleman A.I.E.')
    header['TELESCOP'] = ('CA-3.5', '3.5 m Calar Alto telescope')
    header['INSTRUME'] = ('CARMENES', 'instrument used')
    header['TIMETYPE'] = ('Guaranteed', 'Guaranteed or open time')
    header['PROG-NUM'] = ('CARMENES', 'Calar Alto programme identification')
    header['PROG-PI'] = ('CARMENES', 'Programme principal investigator identification')
    header['OBSERVER'] = ('Service', 'Service or visiting observer name')
    header['OPERATOR'] = ('', 'Telescope operator')
    header['OBJECT'] = (OBJECT, 'Name for the object observed')
    header['ORIGIN'] = ('CARMENES', 'Organization responsible for creating FITS file')#'Identifies the organization responsible for creating the FITS file'
    header['OBS-MODE'] = ('NIRmaster+VIS', 'Observation mode')
    #header['RA'] = (5.022730, '00:20:05.4 RA (J2000) pointing')
    #header['DEC'] = (-64.87395, '-64:52:26.2 DEC (J2000) pointing')
    header['RA'] = (0.0, '00:00:00.00 Requested RA (J2000) telescope coordinate')
    header['DEC'] = (0.0, '+00:00:00.00 Requested DEC (J2000) telescope coordinate')
    header['EQUINOX'] = (2000., 'standard FK5 (years)')
    header['RADECSYS'] = ('FK5', 'coordinate reference frame')
    header['SPECTRNG'] = ('520-1700 nm (VIS+NIR)', 'Spectral range')
    header['SPECTRES'] = ('82000', 'Spectral resolution Delta lambda / lambda')
    header['APERTURE'] = ('1.5', 'Aperture of fibres on sky in arcsec')

    # EXTENSION HEADERS
    header['CATG'] = (CATG, 'Observation category')
    header['TECH'] = (TECH, 'Technical type')
    header['TYPE'] = (TYPE, 'Observation type in fibA and fibB')
    header['SOURCE'] = (SOURCE, 'Observation source in fibA and fibB')
    header['EXPTIME'] = (0.0, 'Exposition time')
    if self.ARM == 'NIR':
        header['DETSIZE'] = ('[1:4096,1:2048]', '[px] x-range, yrange of full frame')
    header['RO_TIME'] = (0.0, 'Read-out time')
    header['SEQNUM'] = (0, 'Exposition number')
    header['EXPNUM'] = (0, 'Number of exposures')
    header['ARCVER'] = ('', 'ARC version') # only for VIS
    header['DATE-OBS'] = (DATE_OBS, 'date (yyyy-mm-dd) of observation') # NOTE check this
    header['MJD-OBS'] = (MJD_OBS, 'Modified Julian Date at start of observation')
    header['OBSEPOCH'] = ('', 'Observation epoch')
    header['EFFTIME'] = (EFFTIME, 'Effective time')
    header['AIRMASS'] = (0.0, 'Airmass')
    header['AMSTART'] = (0.0, 'Airmass at start')
    header['AMEND'] = (0.0, 'Airmass end')

    # METEO
    header['HIERARCH CAHA GEN AMBI WIND SPEED'] = (0.0, '[m/s] Wind speed')
    header['HIERARCH CAHA GEN AMBI WIND DIR'] = (0.0, '[deg] Wind direction')
    header['HIERARCH CAHA GEN AMBI RHUM'] = (0.0, '[%] Relative humidity')
    header['HIERARCH CAHA GEN AMBI PRES'] = (0.0, '[1e2 Pa] Air pressure')
    header['HIERARCH CAHA GEN AMBI TEMP'] = (0.0, '[degC] Air temperature')

    # TELESCOPE
    header['HIERARCH CAHA TEL GEOELEV'] = (2168, 'Elevation above sea [m]')
    header['HIERARCH CAHA TEL GEOLAT'] = (37.2236, 'Geographic latitude [deg]')
    header['HIERARCH CAHA TEL GEOLON'] = (-2.5463, 'Geographic longitude [deg]')
    header['HIERARCH CAHA TEL DOME.AZ'] = (147.00, 'Dome azumut [deg]')
    header['HIERARCH CAHA TEL DOME.EL_LOW_E'] = (70.00, 'Dome door low limit')
    header['HIERARCH CAHA TEL DOME.EL_UPP_E'] = (93.00, 'Dome door upper limit')
    header['HIERARCH CAHA TEL FOCU ID'] = ('CASS', 'Instrument focal position')
    header['HIERARCH CAHA TEL MIRR S1 COLLAREA'] = (2.94, 'Effective collecting area of primary mirror [m**2]')
    header['HIERARCH CAHA TEL FOCU F_RATIO'] = (8.00, 'Telescope F-ratio')
    header['HIERARCH CAHA TEL FOCU LEN'] = (17.610, 'Telescope focal length [m]')
    header['HIERARCH CAHA TEL FOCU SCALE'] = (0.0826, 'Telescope focus scale at instrument focal position [mm/arcsec]')
    header['HIERARCH CAHA TEL FOCU VALUE'] = (23.060, 'Focus relative setting')
    header['HIERARCH CAHA TEL POS SET RA'] = (388065.0, '07:11:11:0 RA preset [sec]')
    header['HIERARCH CAHA TEL POS SET DEC'] = (43.499167, '+43:29:57.0 DEC preset [sec]')
    header['HIERARCH CAHA TEL POS SET EQUINOX'] = (2000.0, 'Equinox at present')
    header['HIERARCH CAHA TEL POS AZ_START'] = (0.0, 'Telescope azimuth at start of observation')
    header['HIERARCH CAHA TEL POS EL_START'] = (0.0, 'Telescope elevation at start of observation')
    header['HIERARCH CAHA TEL POS HA_START'] = (6.662358, '00:26:39.0 HA [deg]')
    header['HIERARCH CAHA TEL SLATEL'] = ('', 'Telescope name known to SLALIB')

    # FRONTEND
    header['HIERARCH CAHA INS FRONTEND PICKMIRR'] = ('', 'Pick-up mirror position')
    header['HIERARCH CAHA INS FRONTEND ADCANG1'] = (0.0, 'Rotating angle of first ADC prism [deg]')
    header['HIERARCH CAHA INS FRONTEND ADCANG2'] = (0.0, 'Rotating angle of second ADC prism [deg]')
    header['HIERARCH CAHA INS FRONTEND ACGFOCUS'] = (0.0, 'Acquisition & guiding camera focus [mm]')
    header['HIERARCH CAHA INS FRONTEND NIRMIRROR'] = ('A+B', 'Position of calibration pickup mirror in front of NIR fibre input unit')
    header['HIERARCH CAHA INS FRONTEND VISMIRROR'] = ('', 'Position of calibration pickup mirror in front of VIS fibre input unit')
    header['HIERARCH CAHA INS FRONTEND REFLECT'] = ('Beamsplitter', 'Beamsplitter or Mirror in optical path')

    # A&G
    header['HIERARCH CAHA INS ACQGUIDE GUIDING'] = ('', 'Acquisition and guide computer status')
    header['HIERARCH CAHA INS ACQGUIDE NIMAGES'] = (1, 'Number of images taken by A&G system during exposure')

    # CRYOSTATS
    header['HIERARCH CAHA CRYO PRESS1'] = (0.0, 'Pressure in sensor 1 [mbar]')
    header['HIERARCH CAHA CRYO PRESS2'] = (0.0, 'Pressure in sensor 2 [mbar]')
    header['HIERARCH CAHA CRYO TEMPMON1'] = (0.0, 'Filter box [K]')
    header['HIERARCH CAHA CRYO TEMPMON2'] = (0.0, 'Motor [K]')
    header['HIERARCH CAHA CRYO TEMPMON3'] = (0.0, 'Optics [K]')
    header['HIERARCH CAHA CRYO TEMPMON4'] = (0.0, 'Fan out [K]')
    header['HIERARCH CAHA CRYO TEMPMON5'] = (0.0, 'Detector plate [K]')
    header['HIERARCH CAHA CRYO TEMPMON6'] = (0.0, 'Cold plate [K]')
    header['HIERARCH CAHA CRYO TEMPMON7'] = (0.0, 'Inner shield [K]')
    header['HIERARCH CAHA CRYO TEMPMON8'] = (0.0, 'Outer shield [K]')

    # DETECTORS
    header['HIERARCH CAHA DET BIASEC'] = ('[2049:2070,1:2048]', 'Bias section. Only for VIS') # only for VIS # WHAT IS THIS?
    header['HIERARCH CAHA DET TRIMSEC'] = ('4:2045,4:2045', 'Effective section for cutting the detector. Only for VIS.') # only for VIS
    header['HIERARCH CAHA DET CAMERA'] = ('', 'Camera name')
    header['HIERARCH CAHA DET PIXSCALE'] = (0.0, 'arcsec/pixel')
    header['HIERARCH CAHA DET ENOISE'] = (0.0, 'Electrons/read')
    header['HIERARCH CAHA DET ELECTRON'] = ('', '') # only for NIR
    header['HIERARCH CAHA DET FILTER'] = ('', 'Filter macro name of filter combinations') # only for NIR
    header['HIERARCH CAHA DET TPLNAME'] = ('', 'Macro/template name') # only for NIR
    header['HIERARCH CAHA DET TIMER0'] = ('', 'milliseconds') # only for NIR
    header['HIERARCH CAHA DET TIMER1'] = (0, 'milliseconds') # only for NIR
    header['HIERARCH CAHA DET TIMER2'] = (0, 'microseconds') # only for NIR
    header['HIERARCH CAHA DET PTIME'] = (0, 'Pixel-time (units)') # only for NIR
    header['HIERARCH CAHA DET READMODE'] = ('', 'Read cycle-type')
    header['HIERARCH CAHA DET IDLEMODE'] = ('', 'Idle to read transition')
    header['HIERARCH CAHA DET SAVEMODE'] = (0, 'Cave cycle-type')
    header['HIERARCH CAHA DET CPAR1'] = ('', 'Cycle type parameter')
    header['HIERARCH CAHA DET ITIME'] = (0.0, '(on chip) integration time [s]') # only for NIR
    header['HIERARCH CAHA DET HCOADDS'] = (0, 'Number of hardware coadds') # only for NIR
    header['HIERARCH CAHA DET PCOADDS'] = (0, 'Number of coadded plateus/periods') # only for NIR
    header['HIERARCH CAHA DET SCOADDS'] = (0, '# of software coadds') # only for NIR
    header['HIERARCH CAHA DET NCOADDS'] = (0, 'Effective coadds (total)') # only for NIR
    header['HIERARCH CAHA DET EXPTIME'] = (0.0, 'Total integration time [s]')
    header['HIERARCH CAHA DET FRAMENUM'] = (0, '') # only for NIR
    header['HIERARCH CAHA DET SKYFRAME'] = ('', '') # only for NIR
    header['HIERARCH CAHA DET SAVEAREA'] = ('', '') # only for NIR
    header['HIERARCH CAHA DET SOFTWARE'] = ('', 'Software version') # only for NIR
    header['HIERARCH CAHA DET CCDNAME'] = ('', 'Name of CCD detector') # only for VIS
    header['HIERARCH CAHA DET CCDPSIZ'] = (CCDPSIZ, 'Pixel size [mum]') # only for VIS
    header['HIERARCH CAHA DET CCDORI'] = (0.0, 'CCD orientation') # only for VIS
    header['HIERARCH CAHA DET CCDSPEED'] = (0.0, 'CCD readout speed') # only for VIS
    header['HIERARCH CAHA DET CCDBINX'] = (1.000, 'Binning factor along X axis') # only for VIS
    header['HIERARCH CAHA DET CCDBINY'] = (1.000, 'Binning factor along Y axis') # only for VIS
    header['HIERARCH CAHA DET CCDGAIN'] = (1.0, 'CCD gain e-/ADU') # only for VIS
    header['HIERARCH CAHA DET CCDMEAN'] = (CCDMEAN, 'Bias of CCD at selected Gain (DN)') # only for VIS
    header['HIERARCH CAHA DET CCDSAT'] = (0, 'Saturation of CCD (DN)') # only for VIS
    header['HIERARCH CAHA DET CCDRON'] = (CCDRON, 'Readout noise of CCD at selected gain (e-)') # only for VIS
    header['HIERARCH CAHA DET CCDTEMP'] = (0.0, '[K] CCD temperature') # only for VIS
    header['HIERARCH CAHA DET DATASEC'] = ('', 'Image portion of frame') # only for VIS
    header['HIERARCH CAHA DET CCDSEC'] = ('', 'Orientation of full frame') # only for VIS
    header['HIERARCH CAHA DET ORIGSECX'] = (ORIGSECX, 'Original size full frame nx') # only for VIS
    header['HIERARCH CAHA DET ORIGSECY'] = (ORIGSECY, 'Original size full frame ny') # only for VIS
    #header['HIERARCH CAHA DET CCDPSIZ'] = (self.DPIX*1.e3, 'Pixel size [mum]')
    ##header.update('HIERARCH CAHA DET CCD1 PSZY', self.DPIX*1.e3, 'pixel size in microns')
    #header['HIERARCH CAHA DET CCDRON '] = ( 0.000, 'Readout noise of CCD at selected gain (e-)')
    #header['HIERARCH CAHA DET CCDGAIN'] = (1.000, 'CCD gain e-/ADU')
    #header['HIERARCH CAHA DET CCDBINX'] = (1.000, 'Binning factor along X axis')
    #header['HIERARCH CAHA DET CCDBINY'] = (1.000, 'Binning factor along Y axis')

    # CALIBRATION UNITS
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM1'] = ('', 'Daily hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM2'] = ('', 'First master hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM3'] = ('', 'Second master hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM4'] = ('', 'Third master hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM5'] = ('', 'First super-master hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM6'] = ('', 'Second super-master hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET NUM7'] = ('', 'Third super-master hollow-cathode lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET HALOGEN'] = ('', 'Flat-field halogen lamp')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT1'] = (0.00, 'Daily hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT2'] = (0.00, 'First master hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT3'] = (0.00, 'Second master hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT4'] = (0.00, 'Third master hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT5'] = (0.00, 'First super-master hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT6'] = (0.00, 'Second super-master hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET CURRENT7'] = (0.00, 'Third super-master hollow-cathode lamp current in 1e 3 A')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE1'] = (0, 'Daily hollow-cathode lamp age [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE2'] = (0, 'First master hollow-cathode lamp age [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE3'] = (0, 'Second master hollow-cathode lamp age [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE4'] = (0, 'Third master hollow-cathode lampage [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE5'] = (0, 'First super-master hollow-cathode lamp age [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE6'] = (0, 'Second super-master hollow-cathode lamp age [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET AGE7'] = (0, 'Third super-master hollow-cathode lamp age [mA h]')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM1'] = ('', 'Daily hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM2'] = ('', 'First master hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM3'] = ('', 'Second master hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM4'] = ('', 'Third master hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM5'] = ('', 'First super-master hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM6'] = ('', 'Second super-master hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM7'] = ('', 'Third super-master hollow-cathode lamp code')
    header['HIERARCH CAHA INS CALUNIT SOCKET LAMPNUM8'] = ('', 'Flat-field halogen lamp code')
    header['HIERARCH CAHA INS CALUNIT OCTAGON'] = (1, 'Socket number selected by the mirror at centre of octagon')
    header['HIERARCH CAHA INS CALUNIT WHEEL-A'] = ('', 'Position of filter wheel of A fibre')
    header['HIERARCH CAHA INS CALUNIT WHEEL-B'] = ('', 'Position of filter wheel of B fibre')
    header['HIERARCH CAHA INS CALUNIT NOTCHFIL'] = ('', 'Notch filter')
    header['HIERARCH CAHA INS CALUNIT 3RD-SHUT'] = ('', 'Third-input shutter')
    header['HIERARCH CAHA INS CALUNIT LGHTSNSR'] = ('', 'Light sensor output')

    # CHAMBER
    header['HIERARCH CAHA INS CHAMBER TEMP1'] = (0.0, '')
    header['HIERARCH CAHA INS CHAMBER TEMP2'] = (0.0, '')
    header['HIERARCH CAHA INS CHAMBER POXYGEN1'] = (0.0, 'Partial pressure of O2 (Master) [units TBD]')
    header['HIERARCH CAHA INS CHAMBER POXYGEN2'] = (0.0, 'Partial pressure of O2 (Slave) [units TBD]')

    # EXPOSURE METERS
    header['HIERARCH CAHA INS EXPMETER DETECTOR'] = ('ON', '')
    header['HIERARCH CAHA INS EXPMETER STATUS'] = ('Online', 'Connected to ICS')
    header['HIERARCH CAHA INS EXPMETER WHEEL'] = ('Open', 'Position of the expmeter filter wheel')
    header['HIERARCH CAHA INS EXPMETER FILE'] = ('', 'ExpMeter ascii filename')

    # FIBRES
    header['HIERARCH CAHA INS FIBRE SHAKER'] = ('ON', 'Optical fibre shaker')
    header['HIERARCH CAHA INS FIBRE PDU'] = (0, '')

    # TANK
    header['HIERARCH CAHA INS TANK TEMP1'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP2'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP3'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP4'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP5'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP6'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP7'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP8'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP9'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP10'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP11'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP12'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP13'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP14'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP15'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK TEMP16'] = (0.0, '[K]')
    header['HIERARCH CAHA INS TANK PRESS1'] = (0.0, '[Pa]')
    header['HIERARCH CAHA INS TANK PRESS2'] = (0.0, '[Pa]')
    header['HIERARCH CAHA INS TANK SORPPUMP STATUS'] = ('', '')
    header['HIERARCH CAHA INS TANK SORPPUMP TEMP1'] = (0, '')
    header['HIERARCH CAHA INS TANK SORPPUMP TEMP2'] = (0, '')

    # TARGET
    # NOTE: ADDED, NOT IN IMAGE HEADERS DOCUMENT FDR-11A
    header['HIERARCH CAHA TARG NAME'] = ('', 'Name of target')
    header['HIERARCH CAHA TARG RA'] = (0.0, '00:00:00.000 (J2000)')
    header['HIERARCH CAHA TARG DEC'] = (0.0, '00:00:00.00 (J2000)')
    header['HIERARCH CAHA TARG MURA'] = (0.0, '[mas/a]')
    header['HIERARCH CAHA TARG MUDEC'] = (0.0, '[mas/a]')
    header['HIERARCH CAHA TARG RV'] = (0, '[km/s] Radial velocity')
    header['HIERARCH CAHA TARG SPEC TYP'] = ('MV', 'Spectral Type')
    header['HIERARCH CAHA TARG MASK'] = ('', 'Spectral Mask')

    # SIMULATION

    # STUPID HIERARCH BUG IN ASTROPY
    a = len('HIERARCH CAHA SIMU IN A1')
    c = len('FibB infile1')
    limit = 73
    if INFILEA1:
        b = len(INFILEA1)
        tot = a+b+c
        if a+b+c > limit:
            newlength = limit-tot
            INFILEA1 = INFILEA1[:newlength]
    if INFILEA2:
        b = len(INFILEA2)
        tot = a+b+c
        if a+b+c > limit:
            newlength = limit-tot
            INFILEA2 = INFILEA2[:newlength]
    if INFILEB1:
        b = len(INFILEB1)
        tot = a+b+c
        if a+b+c > limit:
            newlength = limit-tot
            INFILEB1 = INFILEB1[:newlength]
    if INFILEB2:
        b = len(INFILEB2)
        tot = a+b+c
        if a+b+c > limit:
            newlength = limit-tot
            INFILEB2 = INFILEB2[:newlength]

    header['HIERARCH CAHA SIMU VERSION'] = (__version__, 'Simulation version')
    header['HIERARCH CAHA SIMU SMP METH'] = (SMETH, 'Sampling method')
    if self.SAMPLING == "grid":
        header['HIERARCH CAHA SIMU PRTB WIDTH'] = (PERTURB, 'Perturbation width (unit-circle units)')
    if self.SAMPLING == "mc":
        header['HIERARCH CAHA SIMU MC METH'] = (MCMETH, 'Monte-Carlo sampling method')
    header['HIERARCH CAHA SIMU FIBA'] = (FIBA, 'Simulation type Fiber A')
    header['HIERARCH CAHA SIMU FIBA SRC'] = (FIBA_SRC, 'Fiber A source')
    header['HIERARCH CAHA SIMU INA1'] = (INFILEA1, 'FibA infile1')
    header['HIERARCH CAHA SIMU INA2'] = (INFILEA2, 'FibA infile2')
    header['HIERARCH CAHA SIMU RVA'] = (RVA, '[m/s] Input radial velocity Fiber A')
    header['HIERARCH CAHA SIMU FIBB'] = (FIBB, 'Simulation type Fiber B')
    header['HIERARCH CAHA SIMU FIBB SRC'] = (FIBB_SRC, 'Fiber B source')
    header['HIERARCH CAHA SIMU INB1'] = (INFILEB1, 'FibB infile1')
    header['HIERARCH CAHA SIMU INB2'] = (INFILEB2, 'FibB infile2')
    header['HIERARCH CAHA SIMU RVB'] = (RVB, '[m/s] Input radial velocity Fiber B')
    header['HIERARCH CAHA SIMU SLITTYP'] = (SLITTYP, 'Slit function cross-section type')
    if self.SAMPLING == "grid":
        header['HIERARCH CAHA SIMU NS IN'] = (NS_IN, 'Input number of slit cross-section samples')
        header['HIERARCH CAHA SIMU NS EFF'] = (NS_EFF, 'Effective number of slit cross-section samples')
    elif self.SAMPLING == "mc":
        header['HIERARCH CAHA SIMU NRA IN'] = (NRAYS_IN[0], 'Input number of rays fib A(at cmdline)')
        header['HIERARCH CAHA SIMU NRB IN'] = (NRAYS_IN[1], 'Input number of rays fib B(at cmdline)')
        for i,N in enumerate(NRAYS_INMA):
            header['HIERARCH CAHA SIMU NRA%s' % (i+OMIN)] = (N, 'Input number of rays order %s fib A' % (i+OMIN))
        for i,N in enumerate(NRAYS_INMB):
            header['HIERARCH CAHA SIMU NRB%s' % (i+OMIN)] = (N, 'Input number of rays order %s fib B' % (i+OMIN))
        header['HIERARCH CAHA SIMU MEANNRA'] = (MEAN_NRAYS_A, 'Mean rays per raw pixel fiber A')
        header['HIERARCH CAHA SIMU MEANNRB'] = (MEAN_NRAYS_B, 'Mean rays per raw pixel fiber B')
        header['HIERARCH CAHA SIMU MINNRA'] = (MIN_NRAYS_A, 'Min rays per raw pixel fiber A')
        header['HIERARCH CAHA SIMU MINNRB'] = (MIN_NRAYS_B, 'Min rays per raw pixel fiber B')
        header['HIERARCH CAHA SIMU MAXNRA'] = (MAX_NRAYS_A, 'Max rays per raw pixel fiber A')
        header['HIERARCH CAHA SIMU MAXNRB'] = (MAX_NRAYS_B, 'Max rays per raw pixel fiber B')
        header['HIERARCH CAHA SIMU TOTNRA'] = (TOT_NRAYS_A, 'Total rays fiber A')
        header['HIERARCH CAHA SIMU TOTNRB'] = (TOT_NRAYS_B, 'Total rays fiber B')
    header['HIERARCH CAHA SIMU DW'] = (DW, '[nm] Input wavelength grid sampling')
    header['HIERARCH CAHA SIMU OMIN'] = (OMIN, 'Minimum echelle order')
    header['HIERARCH CAHA SIMU OMAX'] = (OMAX, 'Maximum echelle order')
    #header['HIERARCH CAHA SIMU NWPIX'] = ('', 'Average wavelengths per pixel')
    header['HIERARCH CAHA SIMU NIRGAP'] = (NIRGAP, '[mm] Gap length of NIR detectors')
    header['HIERARCH CAHA SIMU BLAZE'] = (BLAZE, 'Blaze function simulated')
    header['HIERARCH CAHA SIMU PN'] = (PN, 'Shot noise simulated')
    header['HIERARCH CAHA SIMU BIAS'] = (BIAS, '[DN] Input bias')
    header['HIERARCH CAHA SIMU RON'] = (RON, '[DN] Input readout noise')
    header['HIERARCH CAHA SIMU TELL CH4'] = (CH4, 'Telluric species included')
    header['HIERARCH CAHA SIMU TELL CO2'] = (CO2, 'Telluric species included')
    header['HIERARCH CAHA SIMU TELL H2O'] = (H2O, 'Telluric species included')
    header['HIERARCH CAHA SIMU TELL N2O'] = (N2O, 'Telluric species included')
    header['HIERARCH CAHA SIMU TELL O2'] = (O2, 'Telluric species included')
    header['HIERARCH CAHA SIMU TELL O3'] = (O3, 'Telluric species included')
    header['HIERARCH CAHA SIMU TIME FIBA'] = (TIME_A, '[s] Fiber A simulation time')
    header['HIERARCH CAHA SIMU TIME FIBB'] = (TIME_B, '[s] Fiber B simulation time')


    if self.ARM == 'NIR' and not self.single_ext:
        hdu_tom = pyfits.ImageHDU(np.asarray(self.image[:self.NYPIX, :self.NXPIX/2], dtype=np.uint16), name='TOM')
        tom_header = hdu_tom.header
        tom_header['DETSEC'] = ('[1:2048,1:2048]', '[px] section of DETSIZE ')
        hdu_jerry = pyfits.ImageHDU(np.asarray(self.image[:self.NYPIX, self.NXPIX/2:self.NXPIX], dtype=np.uint16), name='JERRY')
        jerry_header = hdu_jerry.header
        jerry_header['DETSEC'] = ('[2049:4096,1:2048]', '[px] section of DETSIZE ')
        _hdu_list = [hdu, hdu_tom, hdu_jerry]
        if self.subwindows_ext:
            hdu_subwin = pyfits.ImageHDU(np.zeros((32, 32), dtype=np.uint16), name='SUBWINDOWS')
            _hdu_list.append(hdu_subwin)
    else:
        _hdu_list = [hdu]   # temp HDU list
    hdulist = pyfits.HDUList(_hdu_list)
    hdulist.writeto('%s.fits' % self.outfile, clobber=True)

    if gzip:
        call_gzip = "gzip %s.fits" % self.outfile
        subprocess.check_call(call_gzip.split())
        print "Data written to '%s.fits.gz'." % self.outfile
    else:
        print "Data written to '%s.fits'." % self.outfile
    print hdulist.info()
    return 0

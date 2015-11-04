"""
spectralarm.py

The SpectralArm class creates an object instance that encompasses most of the
self-describing attributes needed to describe a spectral channel (VIS or NIR)
and each instance.

In the program, any simulator information should be stored as a SpectralArm
attribute in order to ensure easy transmission of global variables through
different functions of CSFS.
"""

import numpy as np
#import matplotlib.pyplot as plt
import ctypes
import os
import sys
from physics import *
from printer import Printer
import cfunctions
from importparams import parse_params
import debug
import wavefuncs
import defaults

# descriptions
arm_dict = {
    'N' : 'NIR',
    'V' : 'VIS'
    }

sim_dict = {
    'B' : 'BIAS',
    'C' : 'LASER_COMB',
    'F' : 'FLATFIELD',
    'L' : 'LINE_LIST',
    'P' : 'PHOENIX_MODEL',
    'S' : 'SPECTRUM',
    'T' : 'THAR',
    'U' : 'UNE',
    'W' : 'WAVEMAP',
    'X' : 'FABRY-PEROT',
    '0' : 'NONE'
    }

meth_dict = {
    'cdf'  : 'Cumulative_Distribution_Function',
    'rs'   : 'Rejection_Sampling',
    'grid' : 'Grid_Sampling',
    'ps'   : 'Perturbation_Slit',
    'psw'  : 'Perturbation_Slit_Random_Wavelengths'
    }


class SpectralArm(object):
    """
    Class for a spectral arm. Initialization is either for NIR or VIS arm.

    """
    def __init__(self, args):

        params = parse_params() # get spectral arm parameters

        # INITIALIZE DESCRIPTIVE VALUES
        self.fib_char = ['A', 'B']
        self.catg = 'CAL'         # SCI|CAL|TST
        self.tech = 'ECHELLE'     # ECHELLE|IMAGE
        self.fib_obstype = ['OFF', 'OFF']     # fibA,fibB|BIAS|DARK
        # where fibA/B = STAR|SKY|OFF|FLAT|WAVE
        self.fib_src = ['OFF', 'OFF'] # fibA,fibB
        # where fibA/B = OBJ|STD|SKY|OFF|HAL|DOME|ThNe1|ThNe2|UNe1|UNe2|...|FP
        self.sciobject = None

        self.ARM = arm_dict[args.arm]
        print 'Initializing %s arm...' % (self.ARM)


        # GET COMMAND LINE ARGUMENTS
        # ---------------------------------------------------------------------
        self.blaze = args.blaze         # blaze flag
        self.noise = args.noise
        if self.noise:
            if self.noise == True:
                self.sig_rn = defaults.sig_rn
            else:
                self.sig_rn = self.noise
            self.bias = defaults.bias
        else:
            self.bias = 0.0
            self.sig_rn = 0.0
        self.SIM_LIST = defaults.SIM_LIST
        self.fib_simmode = [args.a, args.b] # @MZ
        self.fib_rv = [args.rva, args.rvb] # rv shifts
        self.simmodetypes = [sim_dict[mode] for mode in self.fib_simmode] # fiber sim description
        self.infiles = [args.afile, args.bfile] # @MZ
        self.slit = args.slit # slit function
        self.solve = args.solve
        self.preview = args.preview
        self.outfile = args.out
        self.outpath = args.path
        self.full_output = args.full_output
        self.subwindows_ext = args.subwindows_ext
        self.crop = args.crop
        self.single_ext = args.single_ext
        # telluric lines
        if args.tell:
            self.tell = sorted(args.tell)
        else:
            self.tell = args.tell
        # wavelength trace
        self.WT_FLAG = False            # wavelength trace flag
        self.wtlist = args.wtlist   # input wavelength file
        self.dwt = args.dwt
        self.wt = [args.wta, args.wtb]
        # initialize simulation times
        self.fib_sim_time = [0, 0]   # simulation times for each fiber
        self.full_sim_times = []    # start and finish simulation times
        self.wavemap = False
        self.nw = args.nw # user input nw
        self.interp_flag = args.interp_flag

        # sampling method dependent
        if args.grid:
            self.SAMPLING = 'grid'
            self.MAX_SNR = args.maxsnr    # ccd normalization value
            #self.NSLIT = args.ns            # slit samples
            self.ns = args.ns   # user entered ns
            self.nr = False
            self.fib_nslit = [0, 0]     # actual nslit used in simulation for each fiber; just initialized
            self.fib_nrays = False
            self.perturb = args.perturb
            self.sm = False
        elif args.mc:
            self.SAMPLING = 'mc'
            self.MAX_SNR = False    # ccd normalization value
            #self.NSLIT = False            # slit samples
            # self.nw = False # user input nw @CM remove for functional wavelength inputs
            self.ns = False   # user entered ns
            try:
                if len(args.nr) == 1:
                    self.nr = [int(args.nr[0]), int(args.nr[0])]
                elif len(args.nr) == 2:
                    self.nr = [int(args.nr[0]), int(args.nr[1])]
                else:
                    raise ValueError("Too many -nr/--nrays inputs.")
            except TypeError:
                self.nr = [int(args.nr), int(args.nr)]
            self.fib_nslit = False
            self.fib_nrays = [[],[]]
            self.fib_mean_rays_per_pixel = [0, 0]
            self.fib_min_rays_per_pixel = [0, 0]
            self.fib_max_rays_per_pixel = [0, 0]
            self.fib_nrays_tot = [0, 0]
            self.perturb = False
            self.sm = args.sm
        # ---------------------------------------------------------------------


        # GET SPECTROGRAPH PARAMETERS
        # ---------------------------------------------------------------------
        self.TAU_0 = np.radians(params["TAU_0"])
        self.DFIBER = params["DFIBER"]
        self.HEIGHT_RATIO = params["HEIGHT_RATIO"]
        self.fib_OFFSET = [params["OFFSET_A"], params["OFFSET_B"]]

        # arm dependent parameters
        self.ARM_FLAG = 0 if self.ARM == "NIR" else 1
        self.FN = params["NIR_FN"] if self.ARM == "NIR" else params["VIS_FN"]
        self.DPIX = params["NIR_DPIX"] if self.ARM == "NIR" else params["VIS_DPIX"]
        self.NXPIX = params["NIR_NXPIX"] if self.ARM == "NIR" else params["VIS_NXPIX"]
        self.NYPIX = params["NIR_NYPIX"] if self.ARM == "NIR" else params["VIS_NYPIX"]
        self.F_COL = params["NIR_F_COL"] if self.ARM == "NIR" else params["VIS_F_COL"]
        self.M0 = params["NIR_M_LOWER"] if self.ARM == "NIR" else params["VIS_M_LOWER"]
        self.MN = params["NIR_M_UPPER"] if self.ARM == "NIR" else params["VIS_M_UPPER"]
        self.ALPHA_E_DEG = params["NIR_ALPHA_E_DEG"] if self.ARM == "NIR" else params["VIS_ALPHA_E_DEG"]
        self.GAMMA_E_DEG = params["NIR_GAMMA_E_DEG"] if self.ARM == "NIR" else params["VIS_GAMMA_E_DEG"]
        self.SIGMA_E = params["NIR_SIGMA_E"] if self.ARM == "NIR" else params["VIS_SIGMA_E"]
        self.ALPHA_G_DEG = params["NIR_ALPHA_G_DEG"] if self.ARM == "NIR" else params["VIS_ALPHA_G_DEG"]
        self.SIGMA_G = params["NIR_SIGMA_G"] if self.ARM == "NIR" else params["VIS_SIGMA_G"]
        self.N_G = params["NIR_N_G"] if self.ARM == "NIR" else params["VIS_N_G"]
        self.F_CAM = params["NIR_F_CAM"] if self.ARM == "NIR" else params["VIS_F_CAM"]
        self.ALPHA_E = np.radians(self.ALPHA_E_DEG)
        self.GAMMA_E = np.radians(self.GAMMA_E_DEG)
        self.ALPHA_G = np.radians(self.ALPHA_G_DEG)
        self.ORDERS = np.arange(self.M0, self.MN+1, 1)
        self.NORDERS = self.ORDERS.size
        self.CCD_DIMS = (self.NYPIX, self.NXPIX)    # CCD dimensions
        # initialize detector offset
        self.XD_0 = 0.0
        self.YD_0 = 0.0

        # NIR detector gap
        if args.gap and self.ARM == "NIR":
            self.GAP_FLAG = 1
            try:
                if len(args.gap) == 1:
                    self.gap = args.gap
                    self.gapoffset = 0.0
                elif len(args.gap) == 2:
                    self.gap = args.gap[0]
                    self.gapoffset = args.gap[1]
            except TypeError:
                self.gap = params["NIR_DGAP"]
                self.gapoffset = 0.0
        else:
            self.GAP_FLAG = 0
            self.gap = False
            self.gapoffset = False

        # user input spectral orders
        if args.oset:
            self.OSET = args.oset
            if len(self.OSET) == 3 and self.OSET[2] < self.OSET[0]:
                self.OSET = np.arange(self.OSET[0], self.OSET[1], self.OSET[2])
            else:
                self.OSET.sort()
                self.OSET = np.array(self.OSET, dtype=np.uint16)
            self.NOSET = self.OSET.size
            print "  * %i orders will be simulated: " % self.NOSET, self.OSET
        else:
            self.OSET = self.ORDERS
            self.NOSET = self.NORDERS
            print "  * All %i orders will be simulated." % self.NOSET

        # ---------------------------------------------------------------------

        # CALCULATE ARM-DEPENDENT SPECTROGRAPH PROPERTIES
        # ---------------------------------------------------------------------
        print '  * Finding blaze wavelengths for %s arm...' % self.ARM
        self.bl_ech = lam_blaze_ech(self.ORDERS, self.SIGMA_E, self.ALPHA_E)
        self.bl_cd = lam_blaze_grism(self.N_G, self.SIGMA_G, self.ALPHA_G)

        print '  * Finding central order...'
        self.m_central = find_nearest(self.bl_ech, self.bl_cd)[0] + self.M0
        self.lamb_central = self.bl_ech[self.m_central - self.M0]

        print '  * Finding FSR wavelength limits...'
        self.fsr_waves_min = self.bl_ech - self.bl_ech / (2 * self.ORDERS)
        self.fsr_waves_max = self.bl_ech + self.bl_ech / (2 * self.ORDERS)

        # set corner orders and wavelengths
        self.corner_orders = self.ORDERS[[0, self.m_central - self.M0, -1]]

        self.corner_wavelengths = np.array((
            (self.fsr_waves_min[0], self.bl_ech[0], self.fsr_waves_max[0]),
            (self.fsr_waves_min[self.m_central - self.M0],
                self.bl_ech[self.m_central - self.M0],
                self.fsr_waves_max[self.m_central - self.M0]),
            (self.fsr_waves_min[-1], self.bl_ech[-1], self.fsr_waves_max[-1])))

        # zemax data from FDR
        self.zemax_orders, self.zemax_wavelengths = self.import_zemax_table()

        # ---------------------------------------------------------------------

        # HANDLE SOLVED/SOLVING DATA
        # ---------------------------------------------------------------------
        # clear offset data, wavelength data, and wavelength solutions
        reset_message = "Are you sure you want to reset the detector offset "\
            "and the known wavelength limits of the detector for the %s "\
            "arm? This cannot be undone." % self.ARM
        possible_answers = ("Y", "YES")
        data_list = [
            "DATA/ccd_waves_min_%s.txt" % self.ARM.lower() ,
            "DATA/ccd_waves_max_%s.txt" % self.ARM.lower() ,
            "DATA/xy_offset_%s.txt" % self.ARM.lower() ,
            "DATA/fsr_waves_min_%s.txt" % self.ARM.lower() ,
            "DATA/fsr_waves_max_%s.txt" % self.ARM.lower() ,
            "DATA/blaze_waves_%s.txt" % self.ARM.lower()]
        data_bool = [not os.path.isfile(datafile) for datafile in data_list]

        if args.reset:
            answer = raw_input("%s\n" % reset_message)
            if answer.upper() in possible_answers:
                print "Pre-existing data has been deleted and will be "\
                    "re-calculated during this simulation.\n"
                os.system("rm DATA/*")
            else:
                print "Pre-existing data will not be deleted.\n"

        if args.solve:
            print '  * Solving CCD wavelength limits...'
            print "  * xd_0 = %.2f mm, yd_0 = %.2f mm" % (self.XD_0, self.YD_0)
            self.find_det_offset()
            print "  * xd_0 = %.2f mm, yd_0 = %.2f mm" % (self.XD_0, self.YD_0)
        elif any(data_bool):
            print '  * Solving CCD wavelength limits...'
            self.find_det_offset()
            data_data = (
                self.ccd_waves_min,
                self.ccd_waves_max,
                np.array((self.XD_0, self.YD_0)),
                self.fsr_waves_min,
                self.fsr_waves_max,
                self.bl_ech)
            for datafile,dataitem in zip(data_list,data_data):
                np.savetxt(datafile, dataitem)
        else:
            print '  * Loading CCD wavelength limits...'
            self.ccd_waves_min = np.loadtxt('DATA/ccd_waves_min_%s.txt' % self.ARM.lower())
            self.ccd_waves_max = np.loadtxt('DATA/ccd_waves_max_%s.txt' % self.ARM.lower())
            self.XD_0, self.YD_0 = np.loadtxt('DATA/xy_offset_%s.txt' % self.ARM.lower())

        # ---------------------------------------------------------------------

        # wavelength coverage of spectral arm
        self.wmin = self.ccd_waves_min.min()
        self.wmax = self.ccd_waves_max.max()

        # pixel bins
        self.ccd_xbins = np.arange((self.NXPIX+1)) - self.NXPIX/2
        self.ccd_ybins = np.arange((self.NYPIX+1)) - self.NYPIX/2

        # CCD center locations
        self.xcenters = (np.arange(1, self.NXPIX+1, 1) - (self.NXPIX+1)/2.0) * self.DPIX / 2.0

        # debug - WHERE SHOULD THIS GO? STILL NEEDED?
        self.debug_mode = args.debug
        if self.debug_mode:
            sys.exit(debug.debug(self))

        # image arrays
        fiber_a = np.zeros(self.CCD_DIMS)#, dtype=np.uint16)
        fiber_b = np.zeros(self.CCD_DIMS)#, dtype=np.uint16)
        self.zero_images = [fiber_a, fiber_b]
        self.images = []

        print "%s arm initialization complete.\n" % self.ARM

# =============================================================================
# :::::::::::::::::::::::::::: METHODS ::::::::::::::::::::::::::::::
# =============================================================================

    def find_det_offset(self):
        """Sets detector offset values.
        """
        self.XD_0 = wavefuncs.find_blaze_offset(self)
        self.YD_0 = wavefuncs.find_ccd_wavelength_limits(self)

    def import_zemax_table(self):
        """
        Imports ZEMAX data given in the FDR internal documents.
        """
        infile = "testing/zemax_%s_wavelengths.txt" % self.ARM.lower()
        m, li, lc, lf = np.loadtxt(infile, unpack=True, skiprows=2)
        m = np.int_(m)  # float to int
        li *= 1.e-3
        lc *= 1.e-3
        lf *= 1.e-3
        wavelengths = np.array(zip(li, lc, lf))
        return m, wavelengths

    # ATTRIBUTE SETTING METHODS
    # NOTE these methods are not necessary, and can be deprecated in the future
    def add_image(self, image):
        self.images.append(image)

    def add_full_sim_time(self, time):
        self.full_sim_times.append(time)

    def set_wt_flag(self, flag):
        self.WT_FLAG = flag

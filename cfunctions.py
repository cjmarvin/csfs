"""
cfunctions.py

Wraps c-code into python. Each block is for a different function.
TODO change parameters for possible dummy variables to keyword arguments.
"""

import numpy as np
import ctypes

# load c library
_lib = np.ctypeslib.load_library('functions', '.')
_lib2 = np.ctypeslib.load_library('sampling_cum_sim', '.')

# define numpy array types
array_1d_double = np.ctypeslib.ndpointer(dtype=np.double,
                                         ndim=1,
                                         flags='C_CONTIGUOUS')
array_1d_int = np.ctypeslib.ndpointer(dtype=np.int,
                                      ndim=1,
                                      flags='C_CONTIGUOUS')
array_2d_double = np.ctypeslib.ndpointer(dtype=np.double,
                                         ndim=2,
                                         flags='C_CONTIGUOUS')
array_2d_int = np.ctypeslib.ndpointer(dtype=np.int,
                                      ndim=2,
                                      flags='C_CONTIGUOUS')
array_2d_uint = np.ctypeslib.ndpointer(dtype=np.uint,
                                      ndim=2,
                                      flags='C_CONTIGUOUS')
# array_2d_ulong = np.ctypeslib.ndpointer(dtype=np.u,
#                                         ndim=2,
#                                         flags='C_CONTIGUOUS')

# ORIGINAL COMPUTE
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute
_lib.compute_c_general.argtypes = [
    ctypes.c_int,    # ARM ID
    ctypes.c_int,    # BLAZE_FLAG
    ctypes.c_int,    # LOC_FLAG
    ctypes.c_int,    # GAP_FLAG
    ctypes.c_ulong,  # nw
    ctypes.c_ulong,  # nslit
    ctypes.c_uint,   # m
    ctypes.c_double, # gap
    ctypes.c_double, # gapoffset
    ctypes.c_double, # xd_0
    ctypes.c_double, # yd_0
    array_1d_double, # n_sell
    array_1d_double, # slitx
    array_1d_double, # slity
    array_1d_double, # lamb
    array_1d_double, # weights
    array_2d_double, # ccd
    array_2d_uint,   # counts
    array_2d_uint,   # m_list
    array_1d_double, # returnx
    array_1d_double  # returny
    ]
_lib.compute_c_general.restype = None

def compute(
    arm_flag,
    blaze_flag,
    location_flag,
    gap_flag,
    nwaves,
    ns,
    m,
    gap,
    gapoffset,
    xd_0,
    yd_0,
    n_sell,
    slitx,
    slity,
    waves,
    weights,
    ccd,
    counts,
    m_list,
    returnx,
    returny):

    _lib.compute_c_general(
        arm_flag,
        blaze_flag,
        location_flag,
        gap_flag,
        nwaves,
        ns,
        m,
        gap,
        gapoffset,
        xd_0,
        yd_0,
        n_sell,
        slitx,
        slity,
        waves,
        weights,
        ccd,
        counts,
        m_list,
        returnx,
        returny)
    return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# COMPUTE GRID (compute c without general options)
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute
_lib.compute_c_grid.argtypes = [
    ctypes.c_int,    # ARM ID
    ctypes.c_int,    # BLAZE_FLAG
    ctypes.c_int,    # GAP_FLAG
    ctypes.c_ulong,  # nw
    ctypes.c_ulong,  # nslit
    ctypes.c_uint,   # m
    ctypes.c_double, # gap
    ctypes.c_double, # gapoffset
    ctypes.c_double, # xd_0
    ctypes.c_double, # yd_0
    array_1d_double, # n_sell
    array_1d_double, # slitx
    array_1d_double, # slity
    array_1d_double, # lamb
    array_1d_double, # weights
    array_2d_double, # ccd
    array_2d_uint   # counts
    ]
_lib.compute_c_grid.restype = None

def compute_grid(
    arm_flag,
    blaze_flag,
    gap_flag,
    nwaves,
    ns,
    m,
    gap,
    gapoffset,
    xd_0,
    yd_0,
    n_sell,
    slitx,
    slity,
    waves,
    weights,
    ccd,
    counts):

    _lib.compute_c_grid(
        arm_flag,
        blaze_flag,
        gap_flag,
        nwaves,
        ns,
        m,
        gap,
        gapoffset,
        xd_0,
        yd_0,
        n_sell,
        slitx,
        slity,
        waves,
        weights,
        ccd,
        counts)
    return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# COMPUTE GRID WITH SLIT PERTURBATIONS
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute
_lib.compute_c_grid_perturb.argtypes = [
    ctypes.c_int,    # ARM ID
    ctypes.c_int,    # BLAZE_FLAG
    ctypes.c_int,    # GAP_FLAG
    ctypes.c_ulong,  # nw
    ctypes.c_ulong,  # nslit
    ctypes.c_uint,   # m
    ctypes.c_double, # gap
    ctypes.c_double, # gapoffset
    ctypes.c_double, # xd_0
    ctypes.c_double, # yd_0
    ctypes.c_double, # perturb
    array_1d_double, # n_sell
    array_1d_double, # slitx
    array_1d_double, # slity
    array_1d_double, # lamb
    array_1d_double, # weights
    array_2d_double, # ccd
    array_2d_uint   # counts
    ]
_lib.compute_c_grid_perturb.restype = None

def compute_grid_perturb(
    arm_flag,
    blaze_flag,
    gap_flag,
    nwaves,
    ns,
    m,
    gap,
    gapoffset,
    xd_0,
    yd_0,
    perturb,
    n_sell,
    slitx,
    slity,
    waves,
    weights,
    ccd,
    counts):
    """c - wrapper"""

    _lib.compute_c_grid_perturb(
        arm_flag,
        blaze_flag,
        gap_flag,
        nwaves,
        ns,
        m,
        gap,
        gapoffset,
        xd_0,
        yd_0,
        perturb,
        n_sell,
        slitx,
        slity,
        waves,
        weights,
        ccd,
        counts)
    return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# COMPUTE FULL MONTE CARLO METHOD USING REJECTION SAMPLING FOR SED SAMPLING
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute
_lib.compute_c_mc_rejsamp.argtypes = [
    ctypes.c_int,    # ARM ID
    ctypes.c_int,    # BLAZE_FLAG
    ctypes.c_int,    # GAP_FLAG
    ctypes.c_int,    # SLIT_FLAG
    ctypes.c_ulong,  # nw
    ctypes.c_ulong,  # nphotons
    ctypes.c_uint,   # m
    ctypes.c_double, # gap
    ctypes.c_double, # gapoffset
    ctypes.c_double, # xd_0
    ctypes.c_double, # yd_0
    ctypes.c_double, # flux max
    ctypes.c_double, # offset
    ctypes.c_double, # slit_locx
    ctypes.c_double, # slit_locy
    ctypes.c_double, # slit_scalex
    ctypes.c_double, # slit_scaley
    array_1d_double, # lamb
    array_1d_double, # weights
    array_2d_uint,   # ccd
    array_2d_double  # wavemap values
    ]
_lib.compute_c_mc_rejsamp.restype = None

def compute_mc_rejsamp(
    arm_flag,
    blaze_flag,
    gap_flag,
    slit_flag,
    nwaves,
    nphotons,
    m,
    gap,
    gapoffset,
    xd_0,
    yd_0,
    fluxmax,
    offset,
    slit_locx,
    slit_locy,
    slit_scalex,
    slit_scaley,
    waves,
    weights,
    ccd,
    wavemap_values):
    """c - wrapper"""

    _lib.compute_c_mc_rejsamp(
        arm_flag,
        blaze_flag,
        gap_flag,
        slit_flag,
        nwaves,
        nphotons,
        m,
        gap,
        gapoffset,
        xd_0,
        yd_0,
        fluxmax,
        offset,
        slit_locx,
        slit_locy,
        slit_scalex,
        slit_scaley,
        waves,
        weights,
        ccd,
        wavemap_values)
    return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# COMPUTE FULL MONTE CARLO METHOD USING CDF FOR SED SAMPLING
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# sampling CDF
_lib.compute_c_mc_cdf.argtypes = [
    ctypes.c_int,    # ARM ID
    ctypes.c_int,    # BLAZE_FLAG
    ctypes.c_int,    # GAP_FLAG
    ctypes.c_int,    # SLIT_FLAG
    ctypes.c_ulong,  # nphotons
    ctypes.c_uint,   # m
    ctypes.c_double, # gap
    ctypes.c_double, # gapoffset
    ctypes.c_double, # xd_0
    ctypes.c_double, # yd_0
    ctypes.c_double, # offset
    ctypes.c_double, # slit_locx
    ctypes.c_double, # slit_locy
    ctypes.c_double, # slit_scalex
    ctypes.c_double, # slit_scaley
    array_1d_double, # lamb
    array_2d_uint,   # ccd
    array_2d_double  # wavelength values
    ]
_lib.compute_c_mc_cdf.restype = None

def compute_mc_cdf(
    arm_flag,
    blaze_flag,
    gap_flag,
    slit_flag,
    nphotons,
    m,
    gap,
    gapoffset,
    xd_0,
    yd_0,
    offset,
    slit_locx,
    slit_locy,
    slit_scalex,
    slit_scaley,
    waves,
    ccd,
    wavemap_values):
    """c - wrapper"""

    _lib.compute_c_mc_cdf(
        arm_flag,
        blaze_flag,
        gap_flag,
        slit_flag,
        nphotons,
        m,
        gap,
        gapoffset,
        xd_0,
        yd_0,
        offset,
        slit_locx,
        slit_locy,
        slit_scalex,
        slit_scaley,
        waves,
        ccd,
        wavemap_values)
    return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for wavelength_grid_test
_lib.wavelength_rejection_sampling.argtypes = [
    ctypes.c_ulong,  # nw
    ctypes.c_ulong,  # nwout
    array_1d_double, # waves
    array_1d_double, # weights
    array_1d_double, # env
    array_1d_double # outwaves
    ]
_lib.wavelength_rejection_sampling.restype = None
    
def wavelength_rejection_sampling(
    nw,
    nwout,
    waves,
    weights,
    env,
    outwaves
    ):
    """c - wrapper"""
    waves = np.ascontiguousarray(waves, dtype=np.float64)
    weights = np.ascontiguousarray(weights, dtype=np.float64)
    env = np.ascontiguousarray(env, dtype=np.float64)
    outwaves = np.ascontiguousarray(outwaves, dtype=np.float64)
    _lib.wavelength_rejection_sampling(
        nw,
        nwout,
        waves,
        weights,
        env,
        outwaves
        )
    return 0    
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# CDF SAMPLING FUNCTION
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# sampling CDF    
_lib2.random_wave.argtypes = [
    array_1d_double, # waves
    array_1d_double, # weights
    ctypes.c_int,    # nw
    ctypes.c_int,    # interpolation flag
    array_1d_double, # outwaves
    ctypes.c_int     # outwaves size
    ]
_lib2.random_wave.restype = None

def random_wave_cdf(
    waves,
    weights,
    nw,
    flag,
    outwaves,
    nwout
    ):
    """c - wrapper"""
    waves = np.ascontiguousarray(waves, dtype=np.float64)
    weights = np.ascontiguousarray(weights, dtype=np.float64)
    outwaves = np.ascontiguousarray(outwaves, dtype=np.float64)
    _lib2.random_wave(
        waves,
        weights,
        nw,
        flag,
        outwaves,
        nwout
        )
    return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



# COMPUTE WITH DETECTOR GAP (ONLY FOR NIR) 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute_gap
#_lib.compute_c_gap.argtypes = [ctypes.c_int,    # ARM ID
                           #ctypes.c_int,    # BLAZE_FLAG
                           #ctypes.c_int,    # LOC_FLAG
                           #ctypes.c_ulong,  # nw
                           #ctypes.c_ulong,  # nslit
                           #ctypes.c_uint,   # m
                           #ctypes.c_double, # gap
                           #ctypes.c_double, # gapoff
                           #ctypes.c_double, # xd_0
                           #ctypes.c_double, # yd_0
                           #array_1d_double, # n_sell
                           #array_1d_double, # slitx
                           #array_1d_double, # slity
                           #array_1d_double, # lamb
                           #array_1d_double, # weights
                           #array_2d_double, # ccd
                           #array_2d_uint,    # counts
                           #array_2d_uint,    # m_list
                           #array_1d_double, # returnx
                           #array_1d_double  # returny
                           #]
#_lib.compute_c_gap.restype = None

#def compute_gap(arm_flag,
            #blaze_flag,
            #location_flag,
            #nwaves,
            #ns,
            #m,
            #gap,
            #gapoff,
            #xd_0,
            #yd_0,
            #n_sell,
            #slitx,
            #slity,
            #waves,
            #weights,
            #ccd,
            #counts,
            #m_list,
            #returnx,
            #returny):
    #"""c - wrapper"""

    #_lib.compute_c_gap(arm_flag,
                   #blaze_flag,
                   #location_flag,
                   #nwaves,
                   #ns,
                   #m,
                   #gap,
                   #gapoff,
                   #xd_0,
                   #yd_0,
                   #n_sell,
                   #slitx,
                   #slity,
                   #waves,
                   #weights,
                   #ccd,
                   #counts,
                   #m_list,
                   #returnx,
                   #returny)
    #return 0
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# COMPUTE WITH PERTURBED SLIT LOCATIONS AND RANDOM WAVELENGTH (SED) SAMPLING
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# declare arguments and return types for compute
#_lib.compute_c_pdf.argtypes = [ctypes.c_int,    # ARM ID
                           #ctypes.c_int,    # BLAZE_FLAG
                           #ctypes.c_int,    # LOC_FLAG
                           #ctypes.c_ulong,  # nw
                           #ctypes.c_ulong,  # nslit
                           #ctypes.c_uint,   # m
                           #ctypes.c_double, # xd_0
                           #ctypes.c_double, # yd_0
                           #ctypes.c_double, # perturb
                           #ctypes.c_double, # flux max
                           #array_1d_double, # n_sell
                           #array_1d_double, # slitx
                           #array_1d_double, # slity
                           #array_1d_double, # lamb
                           #array_1d_double, # weights
                           #array_2d_double, # ccd
                           #array_2d_uint,    # counts
                           #array_2d_uint,    # m_list
                           #array_1d_double, # returnx
                           #array_1d_double  # returny
                           #]
#_lib.compute_c_pdf.restype = None

#def compute_pdf(arm_flag,
            #blaze_flag,
            #location_flag,
            #nwaves,
            #ns,
            #m,
            #xd_0,
            #yd_0,
            #perturb,
            #fluxmax,
            #n_sell,
            #slitx,
            #slity,
            #waves,
            #weights,
            #ccd,
            #counts,
            #m_list,
            #returnx,
            #returny):
    #"""c - wrapper"""

    #_lib.compute_c_pdf(arm_flag,
                   #blaze_flag,
                   #location_flag,
                   #nwaves,
                   #ns,
                   #m,
                   #xd_0,
                   #yd_0,
                   #perturb,
                   #fluxmax,
                   #n_sell,
                   #slitx,
                   #slity,
                   #waves,
                   #weights,
                   #ccd,
                   #counts,
                   #m_list,
                   #returnx,
                   #returny)
    #return 0
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

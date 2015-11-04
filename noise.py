"""
noise.py

Functions that simulate noise in spectral observations.
Default values are in defaults.py.
"""
import numpy as np
import defaults as df
#import matplotlib.pylab as plt

def det_bias(bias):
    print "Adding bias: %s [DN]" % bias
    return bias

def shot_noise(image):
    print "Adding shot noise"
    return np.random.poisson(lam=image)

def readout_noise(rn, shape):
    print "Adding readout noise: rms = %s [DN]" % rn
    return np.random.normal(scale=rn, size=shape)

def dark_current():
    pass

def add_noise(arm):
    mean = arm.image.mean()
    shape = arm.image.shape

    # shot noise
    if arm.SAMPLING == "grid" and arm.noise:
        arm.image = shot_noise(arm.image)

    # crop border (NIR arm)
    if arm.crop and arm.ARM == "NIR":
        c = arm.crop
        nx = arm.NXPIX
        ny = arm.NYPIX
        print 'Cropping %s pixels' % c
        image_tom = arm.image[:ny, :nx/2]
        image_jerry = arm.image[:ny, nx/2:nx]
        mask = np.logical_and.outer(
            [0]*c + [1]*(ny-c-c) + [0]*c,
            [0]*c + [1]*(ny-c-c) + [0]*c)
        image_tom[~mask] = 0
        image_jerry[~mask] = 0
        image_tom = shot_noise(image_tom)
        image_jerry = shot_noise(image_jerry)
        arm.image = np.concatenate((image_tom, image_jerry), axis=1)

    # bias
    if arm.noise:
        arm.image += det_bias(df.bias)

    # readout noise
    if arm.noise:
        if arm.noise == True:
            rn = df.sig_rn
        elif arm.noise:
            rn = arm.noise
        arm.image += readout_noise(rn, shape)

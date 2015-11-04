"""
slitfunctions.py

Generates slit functions according to parameters.

Changes:

08 Mar 2013:
- all returns are changed so that 1D ndarrays are returned. This ensures
  compatibility and flexibility with C wrapping in computing_functions.c
"""

__all__ = ["point_source", "half_moon_even", "half_moon_random"]

import numpy as np
# import matplotlib.pylab as plt
import sys

def preview(ax, x, y, circle=False):
    ax.scatter(x, y, alpha=0.5, edgecolors="none")
    if circle:
        an = np.linspace(0, 2*np.pi, 100)
        ax.plot(circle*np.cos(an), circle*np.sin(an), color="r", linewidth=3, alpha=0.6)
    #ax.set_xlabel(r"$x$ [mm]")
    #ax.set_ylabel(r"$y$ [mm]")
    #ax.set_title(title)
    return ax

def point_source(arm, offset=0.0):
    """
    initialized point source slit function.
    """
    x = 0.0
    y = 0.0
    y += offset

    # scale slit function
    x = x * arm.FN
    y = y * arm.FN

    # apply rotation
    xp = (x * np.cos(arm.TAU_0)) + (y * np.sin(arm.TAU_0))
    yp = (-x * np.sin(arm.TAU_0)) + (y * np.cos(arm.TAU_0))

    # convert to array for looping compatibility
    slitx = np.array((xp))
    slity = np.array((yp))
    #slit = np.array((xp, yp))

    # SET SLIT CROSS SECTION ARRAY
    arm.slit_cross_section = np.array(zip(x, y))

    # preview slit
#     fig = plt.figure()
#     ax = fig.add_subplot(111, aspect='equal')
#     ax.scatter(slit[:, 0], slit[:, 1], s=20, edgecolor=None)
#     plt.show()
    print "  * Point source initialized."
    arm.slittyp = "0D point-source"
    return slitx, slity

def half_moon_even(arm, n=30, offset=0.0):
    """
    creates an evenly sampled slit function.
    """
    # create circular grid
    X, Y = np.mgrid[-n-2:n+2, -n-2:n+2]
    # .5 shift ensures that the exact fiber edges are not sampled ie. at R
    X = X + 0.5
    Y = Y + 0.5
    #X += 0.5#X + 0.5
    #Y += 0.5#Y + 0.5
    message = "  * Created evenly sampled half-moon slit function with %s flux points."
    arm.slittyp = "2D half-moon, even grid"
    return slit_function(arm, X, Y, message, n, offset)


def half_moon_random(arm, n=30, offset=0.0):
    """
    creates a randomly sampled slit function.
    probably should use ns >= 10,000.
    """
    X = np.random.uniform(-n, n+1, size=n)
    Y = np.random.uniform(-n, n+1, size=n)
    message = "  * Created randomly sampled half-moon slit function created with %s flux points."
    arm.slittyp = "2D half-moon, random grid"
    return slit_function(arm, X, Y, message, n, offset)


def slit_function(arm, X, Y, message, n=30, offset=0.0):
    """
    Creates two-slicer slit function.
    """

    height_ratio = arm.HEIGHT_RATIO
    R = n # / 2 @ MZ
    ind = (X**2 + Y**2 <= R**2)
    X = X[ind]
    Y = Y[ind]

    # slice circle
    left_inds = np.where(X < 0.0)
    right_inds = np.where(X > 0.0)
    xl = X[left_inds]
    yl = Y[left_inds]
    xr = X[right_inds]
    yr = Y[right_inds]

    # shift slices into half moon positions
    X[left_inds] += R / 2.0
    Y[left_inds] += R * height_ratio / 2.0
    X[right_inds] -= R / 2.0
    Y[right_inds] -= R * height_ratio / 2.0

    x = X #np.concatenate((xl, xr)) @MZ
    y = Y #np.concatenate((yl, yr)) @MZ

    # scale slit function
    x *= arm.DFIBER / (2.0 * R)#n @ MZ
    y *= arm.DFIBER / (2.0 * R)#n @ MZ
    y += offset

    x = x * arm.FN
    y = y * arm.FN

    # SET DE-ROTATED SLIT IMAGE
    arm.slit_cross_section = np.array(zip(x, y))

    # apply rotation
    xp = (x * np.cos(arm.TAU_0)) + (y * np.sin(arm.TAU_0))
    yp = (-x * np.sin(arm.TAU_0)) + (y * np.cos(arm.TAU_0))

    # ??? HOW SHOULD WE RETURN THE SLIT ???
    #slit = np.array((xp, yp))
    #np.savetxt("DATA/lastslit.dat", zip(x,y)) # @MZ
    slit = np.array(zip(xp, yp))

    # preview slit
    #fig = plt.figure()
    #ax = fig.add_subplot(111, aspect='equal')
    #ax.scatter(slit[:, 0], slit[:, 1], s=20, edgecolor=None)
    #plt.show()

    print message % (slit.shape[0])
    return xp, yp

# ===================== NEW ALGORITHM FOR SLIT FUNCTIONS =====================

def populate(arm):
    """
    Sampling population of the fiber cross-section.
    """
    slit = arm.slit
    n = arm.ns


    # if slit[0] == "0":
        # arm.slittyp = "0D point-source"

    if slit[0] == "1":
        if n % 2 != 0:
            n += 1
        print n
        x, y = np.mgrid[-n/2-1:n/2+1, -n/2-1:n/2+1] / np.float(n/2) + 0.5/n
        arm.slittyp = "2D half-moon, even grid"
    elif slit[0] == "2":
        x = np.random.uniform(-1, 1, size=n)
        y = np.random.uniform(-1, 1, size=n)
        arm.slittyp = "2D half-moon, random grid"
    elif slit[0] == "3":
        if len(slit) == 1:
            locx = 0.0
            locy = 0.0
            scalex = 1.0
            scaley = 1.0
        elif len(slit) == 2:
            locx = 0.0
            locy = 0.0
            scalex = slit[1]
            scaley = slit[1]
        elif len(slit) == 3:
            locx = slit[1]
            locy = slit[1]
            scalex = slit[2]
            scaley = slit[2]
        elif len(slit) == 4:
            locx = slit[1]
            locy = slit[2]
            scalex = slit[3]
            scaley = slit[3]
        elif len(slit) == 5:
            locx = slit[1]
            locy = slit[2]
            scalex = slit[3]
            scaley = slit[4]
        else:
            raise ValueError("Incorrect arguments for slit:3")
        x = np.random.normal(loc=locx, scale=scalex, size=n)
        y = np.random.normal(loc=locy, scale=scaley, size=n)
        arm.slittyp = "2D half-moon, gaussian grid"
    else:
        raise ValueError("Incorrect argument for slit. Choose from {0,1,2,3}.")

    return x, y

def circle_filter(arm, x, y):
    """Only points that fall under the circle pass through.
    """
    #r = arm.DFIBER * arm.FN
    inds = (x**2 + y**2 <= 1.0)
    return x[inds], y[inds]

def slicer(arm, x, y, offset=0.0):
    """Returns a sliced 2d image.
    """
    #c = arm.DFIBER * arm.FN * 0.125
    li = np.where(x < 0.0)
    ri = np.where(x > 0.0)
    x[li] += 0.5
    y[li] += arm.HEIGHT_RATIO*0.5
    x[ri] -= 0.5
    y[ri] -= (arm.HEIGHT_RATIO*0.5)
    return x, y

def scale(arm, x, y, offset=0.0):
    """
    Rescales fiber image.
    """
    # scale fiber
    x *= arm.DFIBER * arm.FN * 0.5 # [mm]
    y *= arm.DFIBER * arm.FN * 0.5 # [mm]
    # fiber offset
    y += offset*arm.FN
    return x, y

# =============================================================================

def slit_image(arm, offset=0.0):
    """
    General slit cross-section generation algorithm.
    """

    #offset = 0.0 # DELETE
    # POPULATE FIBER SQUARE
    # - even, random, gaussian
    x, y = populate(arm)
    if arm.preview:
        import matplotlib.pyplot as plt
        fig = plt.figure("Fiber Cross-section")
        ax1 = fig.add_subplot(221, aspect="equal", title="Pre-Scaled Sample Population")
        ax2 = fig.add_subplot(222, aspect="equal", title="Pre-Scaled Fiber Sample Population")
        ax3 = fig.add_subplot(223, aspect="equal", title="Sliced image")
        ax4 = fig.add_subplot(224, aspect="equal", title="Rotated image", xlabel="[mm]", ylabel="[mm]")
        preview(ax1, x, y, circle=1.0)

    # FILTER FIBER
    x, y = circle_filter(arm, x, y)
    if arm.preview:
        preview(ax2, x, y, circle=1.0)

    # SLICE FIBER
    x, y = slicer(arm, x, y, offset=offset)
    if arm.preview:
        preview(ax3, x, y)

    # SCALE FIBER with OFFSET
    x, y = scale(arm, x, y, offset=offset)

    # SAVE CROSS SECTION
    arm.slit_cross_section = np.array(zip(x, y))
    arm.ns_eff = x.size
    # arm.slitperturb = np.abs(np.mean(np.diff(x)))
    arm.slitperturb = np.abs(np.mean(np.diff(np.concatenate((x,y)))))
    print "Perturbation width = %s" % arm.slitperturb

    # COMMAND LINE INFO
    print "Input number of fiber cross-section samples: %s" % arm.ns
    print "%s fiber cross-section created" % arm.slittyp
    print "Effective number of fiber cross-section samples: %s" % arm.ns_eff

    # ROTATE FIBER IMAGE
    x = (x * np.cos(arm.TAU_0)) + (y * np.sin(arm.TAU_0))
    y = (-x * np.sin(arm.TAU_0)) + (y * np.cos(arm.TAU_0))

    if arm.preview:
        preview(ax4, x, y)
        plt.show()
        sys.exit(0)

    return x, y

if __name__ == "__main__":
    slit_image()

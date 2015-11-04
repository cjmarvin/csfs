#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "armparams.h"

#if __STDC_VERSION__ < 199901L
#define restrict
#endif


// progress bar, inspired by Jinglin
inline void progress_bar(int k, int n) {			
	int res = 20;     // progress bar resolution
	int window = 50;  // progress bar width
	int counter;
	int progress;
	double ratio;
	
	if (n < res) return;
	if ( (k != n) && (k + 1) % (n/res) != 0) return;
  ratio = (double)(k + 1) / n;
	counter = ratio * window;
  printf("    [");
  for (progress=0; progress<window; ++progress) {
		if (progress < counter) {
    	printf("=");
		}
		else if (progress == counter) {
			printf(">");
		}
		else {
			printf(" ");
		}
	}
	printf("\b] %6.2f %%\r", (ratio*100.0) );
  fflush(stdout);
}


inline double cradians(double angle) {
    return angle * M_PI / 180.0;
}


inline void slit_image(
    double dfiber,
    double fn,
    double tau_0,
    double height_ratio,
    double offset,
    gsl_rng* r,
    double* outarr) {
  int i = 0;
  double x;
  double y;
  double tau = cradians(tau_0);

  // CIRCLE FILTER
  while (i != 1) {
    x = gsl_rng_uniform(r) * (2.0) - 1.0;
    y = gsl_rng_uniform(r) * (2.0) - 1.0;
    if (x*x + y*y < 1.0)
      i = 1;
  }
  // SLICE
  if (x < 0.0) {
    x += 0.5;
    y += (height_ratio * 0.5);
  }
  else if (x > 0.0) {
    x -= 0.5;
    y -= (height_ratio * 0.5);
  }
  // SCALE and OFFSET SHIFT
  x *= dfiber * fn * 0.5;
  y *= dfiber * fn * 0.5;
  y += offset * fn;
  // ROTATE
  outarr[0] = x * cos(tau) + y * sin(tau);
  outarr[1] = -x * sin(tau) + y * cos(tau);
}


inline void slit_image_gaussian(
    double dfiber,
    double fn,
    double tau_0,
    double height_ratio,
    double offset,
    double slit_locx,
		double slit_locy,
		double slit_scalex,
		double slit_scaley,
    gsl_rng* r,
    double* outarr) {
  int i = 0;
  double x;
  double y;
  double tau = cradians(tau_0);

  // CIRCLE FILTER
  while (i != 1) {
    x = gsl_ran_gaussian_ziggurat(r, slit_scalex) * (2.0) - 1.0 + slit_locx;
    y = gsl_ran_gaussian_ziggurat(r, slit_scaley) * (2.0) - 1.0 + slit_locy;
    if (x*x + y*y < 1.0)
      i = 1;
  }
  // SLICE
  if (x < 0.0) {
    x += 0.5;
    y += (height_ratio * 0.5);
  }
  else if (x > 0.0) {
    x -= 0.5;
    y -= (height_ratio * 0.5);
  }
  // SCALE and OFFSET SHIFT
  x *= dfiber * fn * 0.5;
  y *= dfiber * fn * 0.5;
  y += offset * fn;
  // ROTATE
  outarr[0] = x * cos(tau) + y * sin(tau);
  outarr[1] = -x * sin(tau) + y * cos(tau);
}


// n_sell()
// Returns wavelength-dependent refractive index according to the Sellmeier
// equation.
//
// Parameters
// ----------
// arm : {"NIR", "VIS"}
//     Spectral arm. "VIS" arm initializes LF5 constants, while "NIR"
//     initializes Infrasil 301 constants. See References.
// lamb : array_like (M)
//     Input wavelength(s) [nm]
//
// Returns
// -------
// array_like (M)
//     Wavelength-dependent refractive index.
//
// Notes
// ----------
// LF5 (VIS) dispersion constants taken from SCHOTT.
// ..[1] http://www.schott.com/advanced_optics/us/abbe_datasheets/schott_datasheet_all_us.pdf
//
// Infrasil 301 (NIR) dispersion constants taken from Heraeus.
// ..[2] http://heraeus-quarzglas.com/media/webmedia_local/downloads/broschren_mo/dataandproperties_optics_fusedsilica.pdf
// ..[3] CARMENES-TNO068
inline double nsell(int arm_id, double lamb) {
  double B1, B2, B3, C1, C2, C3, L1, L2, L3, S1, S2, S3;
  lamb = lamb * 1.0e3; // mm -> microns
  // VIS arm
  if (arm_id == 1) {
    B1 = 1.28035628;
    B2 = 0.163505973;
    B3 = 0.893930112;
    C1 = 0.00929854416;
    C2 = 0.0449135769;
    C3 = 110.493685;
    return sqrt( (B1 * lamb*lamb)/(lamb*lamb - C1) + (B2 * lamb*lamb)/(lamb*lamb - C2) + (B3 * lamb*lamb)/(lamb*lamb - C3) + 1.0 );
  }
  // NIR arm
  else if (arm_id == 0) {
    B1 = 4.76523070e-1;
    B2 = 6.27786368e-1;
    B3 = 8.72274404e-1;
    C1 = 2.84888095e-3;
    C2 = 1.18369052e-2;
    C3 = 9.56856012e1;
    S1 = 1.10098155807840;
    S2 = 0.00124893116320;
    S3 = 0.78106788864000;
    L1 =0.00789065877562;
    L2 =0.07653985064959;
    L3 =85.96849321128080;
    return sqrt((S1*lamb*lamb/(lamb*lamb - L1)) + (S2*lamb*lamb/(lamb*lamb - L2)) + (S3*lamb*lamb/(lamb*lamb - L3)) + 1.0);
  }
  return 0.0;
}


// Computes a population according to a Probability Distribution Function (PDF)
// using Rejection Sampling.
// Can be used universally, but named according to wavelength distribution.
// TODO generalize variable names and function name
void wavelength_rejection_sampling(
    unsigned long nw,
    unsigned long nwout,
    double fluxmax,
    double* restrict lamb,
    double* restrict weights,
    double* restrict out) {
  const double wmin = lamb[0];
  const double wmax = lamb[nw - 1];
  unsigned long i, k;   /* wavelength iterator */
  double u, fc, lc, alpha; // rejection sampling variables
  long seed = time(NULL);
  unsigned long pc = nwout / 100;

  // RANDOM NUMBER GENERATOR INITIALIZATION
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus); /* global rng generator */
  gsl_rng_set(r, seed); /* seed the rng */

  // CUBIC SPLINE OF SPECTRA
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_akima, nw);
  gsl_spline_init(spline, lamb, weights, nw);

  for (i=0; i<nwout; ++i) {
    // REJECTION SAMPLING ALGORITHM FOR WAVELENGTHS
    k = 0;
    while (k != 1) {
      lc = gsl_rng_uniform(r) * (wmax - wmin) + wmin;
      u = gsl_rng_uniform(r);
      fc = gsl_spline_eval(spline, lc, acc);
      alpha = fc / fluxmax;
      if (alpha >= u)
        k = 1;
    }
    out[i] = lc;
    if (i % pc == 0)
      printf("Wavelength sampling %.2f %% done.\n", (float)i / (float)nwout * 100.0 );
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free(acc);
  gsl_rng_free(r);
}


// compute_c_general
// General computation function with optional flags for different outputs.
// Memory can be saved by including dummy input arrays for unused variables.
void compute_c_general(
    int ARM,
    int BLAZE_FLAG,
    int LOC_FLAG,
    int GAP_FLAG,
    unsigned long nw,          /* size of w (and also x and y) */
    unsigned long nslit,         /* number of slit elements */
    unsigned short m,
    double gap,
    double gapoffset,
    double xd_0,
    double yd_0,
    double* restrict n_sell,
    double* restrict x,   /* slit location */
    double* restrict y,   /* slit location */
    double* restrict lamb,      /* wavelengths */
    double* restrict w,         /* intensities */
    double* restrict ccd,
    unsigned long* restrict counts,
    unsigned long* restrict m_list,
    double* restrict returnx,
    double* restrict returny) {
  
  /* assign constants according to spectral arm */
  const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
  const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
  const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
  const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
  const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
  const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
  const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
  const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
  const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
  const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;

  /* pre-calculate constants */
  const double ORDER = (double)m;
  const double F_COL_2 = F_COL * F_COL;
  const double MU_E0 = ALPHA_E - M_PI;
  const double NU_E0 = GAMMA_E;
  const double MU_E1 = ALPHA_E + M_PI;
  const double NU_E1 = -GAMMA_E;

  /* these grism surface angles are still unclear */
  const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
  const double NU_G1 = ALPHA_G;
  const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
  //const double NU_G0 = cradians(0.0);
  //const double NU_G1 = ALPHA_G;
  //const double NU_G2 = -ALPHA_G;
  const double COS_MU_E0 = cos(MU_E0);
  const double SIN_MU_E0 = sin(MU_E0);
  const double COS_MU_E1 = cos(MU_E1);
  const double SIN_MU_E1 = sin(MU_E1);
  const double COS_NU_E0 = cos(NU_E0);
  const double SIN_NU_E0 = sin(NU_E0);
  const double COS_NU_E1 = cos(NU_E1);
  const double SIN_NU_E1 = sin(NU_E1);
  const double COS_NU_G0 = cos(NU_G0);
  const double SIN_NU_G0 = sin(NU_G0);
  const double COS_NU_G1 = cos(NU_G1);
  const double SIN_NU_G1 = sin(NU_G1);
  const double COS_NU_G2 = cos(NU_G2);
  const double SIN_NU_G2 = sin(NU_G2);
  const double XD02 = (double)NXPIX/2.0;
  const double YD02 = (double)NYPIX/2.0;
  const double LGAP = (gapoffset - gap) / 2.0; /* [mm] */
  const double RGAP = (gapoffset + gap) / 2.0; /* [mm] */
  const double TOMLE = -XD02*DPIX + LGAP; /* [mm] */
  const double TOMRE = LGAP; /* [mm] */
  const double JERLE = RGAP; /* [mm] */
  const double JERRE = XD02*DPIX + RGAP; /* [mm] */
  /* initialize temporary variables */
  double veclen;
  double li;
  double wi;
  double n_g;
  double xi, xx, x1, x2, x21, x3, x4, x5, xd;
  double yi, yy, y1, y2, y21, y3, y4, y5, yd;
  double zz, z1, z2, z3, z21, z4, z5;
  double beta;
  double blaze_eff;
  // unsigned long xbin; /* x coordinate [pix] */
  // unsigned long ybin; /* y coordinate [pix] */
  long xbin; /* x coordinate [pix] */
  long ybin; /* y coordinate [pix] */
  unsigned long i;   /* wavelength iterator */
  unsigned long j;   /* slit iterator */

  for (j=0; j<nslit; ++j) {
    xi = x[j];       /* x coord */
    yi = y[j];       /* y coord */
    for (i=0; i<nw; ++i) {
      li = lamb[i];    /* wavelength */
      wi = w[i];       /* intensities */
      n_g = n_sell[i]; /* refractive indices */

      /* LAUNCH RAYS (create normalized vectors) */
      veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
      xx = xi / veclen;
      yy = yi / veclen;
      zz = F_COL / veclen;

      /* BLAZE FUNCTION */
      //xb = xx;
      /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
      /* INTO ECHELLE RF */
      x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
      y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
      z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
      /* NORMALIZATION AFTER ECHELLE RELATION */
      z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);

      /* OUT OF ECHELLE RF */
      x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
      y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
      z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

      /* BLAZE FUNCTION */
      beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
      blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);

      /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
      /* INTO PLANE RF */
      x21 = x2;
      y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
      z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
      /* PLANE RELATION */
      x3 = x21 / n_g;
      y3 = y21 / n_g;
      /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
      z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);

      /* INTO GRISM RF and GRISM RELATION */
      x4 = x3 * n_g;
      y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
      z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
      /* NORMALIZATION AFTER GRATING RELATION */
      z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);

      /* OUT OF GRISM RF */
      x5 = x4;
      y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
      z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

      /* PROJECTION ONTO DETECTOR */
      //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
      //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
      //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
      //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
      //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;

      xd = x5/z5 * F_CAM + xd_0;
      yd = (y5/z5 * F_CAM + yd_0);
      //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);

      switch (LOC_FLAG) {
        /* BIN PIXELS */
        case 0: {
          // NIR DETECTOR GAP
          if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
            xd -= LGAP;
	          xd /= DPIX;
	          yd /= DPIX;
	          xbin = (int)floor(xd + XD02);
	          if (xbin >= 0 && xbin < NXPIX) {
	            ybin = (int)floor(yd + YD02);
	            if (ybin >= 0 && ybin < NYPIX) {
	              if (BLAZE_FLAG == 1) {
	                ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	                counts[xbin+NXPIX*ybin] += 1;
	                m_list[xbin+NXPIX*ybin] = m;
	              }
	              else {
	                ccd[xbin+NXPIX*ybin] += wi;
	                counts[xbin+NXPIX*ybin] += 1;
	                m_list[xbin+NXPIX*ybin] = m;
	              }
	            }
	          }
          }
          else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
            xd -= RGAP;
	          xd /= DPIX;
	          yd /= DPIX;
	          xbin = (int)floor(xd + XD02);
	          if (xbin >= 0 && xbin < NXPIX) {
	            ybin = (int)floor(yd + YD02);
	            if (ybin >= 0 && ybin < NYPIX) {
	              if (BLAZE_FLAG == 1) {
	                ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	                counts[xbin+NXPIX*ybin] += 1;
	                m_list[xbin+NXPIX*ybin] = m;
	              }
	              else {
	                ccd[xbin+NXPIX*ybin] += wi;
	                counts[xbin+NXPIX*ybin] += 1;
	                m_list[xbin+NXPIX*ybin] = m;
	              }
	            }
	          }
          }
					else if (GAP_FLAG == 0) {
	          xd /= DPIX;
	          yd /= DPIX;
	          xbin = (int)floor(xd + XD02);
	          if (xbin >= 0 && xbin < NXPIX) {
	            ybin = (int)floor(yd + YD02);
	            if (ybin >= 0 && ybin < NYPIX) {
	              if (BLAZE_FLAG == 1) {
	                ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	                counts[xbin+NXPIX*ybin] += 1;
	                m_list[xbin+NXPIX*ybin] = m;
	              }
	              else {
	                ccd[xbin+NXPIX*ybin] += wi;
	                counts[xbin+NXPIX*ybin] += 1;
	                m_list[xbin+NXPIX*ybin] = m;
	              }
	            }
	          }
					}
          break;
        }
        /* RETURN EXACT LOCATIONS [mm] */
        case 1: {
          if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
            xd -= LGAP;
	          returnx[i+j*NXPIX] = xd;
	          returny[i+j*NXPIX] = yd;
          }
          else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
            xd -= RGAP;
	          returnx[i+j*NXPIX] = xd;
	          returny[i+j*NXPIX] = yd;
          }
					else if (GAP_FLAG == 0) {
	          returnx[i+j*NXPIX] = xd;
	          returny[i+j*NXPIX] = yd;
					}
          break;
        }
        /* RETURN LOCATIONS 0-CENTERED IN [pix] */
        case 2: {
          if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
            xd -= LGAP;
	          returnx[i+j*NXPIX] = xd/DPIX ;
	          returny[i+j*NXPIX] = yd/DPIX ;
          }
          else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
            xd -= RGAP;
	          returnx[i+j*NXPIX] = xd/DPIX ;
	          returny[i+j*NXPIX] = yd/DPIX ;
          }
					else if (GAP_FLAG == 0) {
	          returnx[i+j*NXPIX] = xd/DPIX ;
	          returny[i+j*NXPIX] = yd/DPIX ;
					}
          break;
        }
        /* RETURN LOCATIONS [pix] */
        case 3: {
          if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
            xd -= LGAP;
	          returnx[i+j*NXPIX] = xd/DPIX + (double)NXPIX/2.0;
	          returny[i+j*NXPIX] = yd/DPIX + (double)NYPIX/2.0;
          }
          else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
            xd -= RGAP;
	          returnx[i+j*NXPIX] = xd/DPIX + (double)NXPIX/2.0;
	          returny[i+j*NXPIX] = yd/DPIX + (double)NYPIX/2.0;
          }
					else if (GAP_FLAG == 0) {
	          returnx[i+j*NXPIX] = xd/DPIX + (double)NXPIX/2.0;
	          returny[i+j*NXPIX] = yd/DPIX + (double)NYPIX/2.0;
					}
          break;
        }
        default:
            break;
      }
			// progress_bar(i, nw);
    }
	}
	// printf("\n");
}


void compute_c_grid(
    int ARM,
    int BLAZE_FLAG,
    int GAP_FLAG,
    unsigned long nw,          /* size of w (and also x and y) */
    unsigned long nslit,         /* number of slit elements */
    unsigned short m,
    double gap,
    double gapoffset,
    double xd_0,
    double yd_0,
    double* restrict n_sell,
    double* restrict x,   /* slit location */
    double* restrict y,   /* slit location */
    double* restrict lamb,      /* wavelengths */
    double* restrict w,         /* intensities */
    double* restrict ccd,
		unsigned long* restrict counts) {

  /* assign constants according to spectral arm */
  const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
  const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
  const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
  const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
  const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
  const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
  const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
  const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
  const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
  const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;

  /* pre-calculate constants */
  const double ORDER = (double)m;
  const double F_COL_2 = F_COL * F_COL;
  const double MU_E0 = ALPHA_E - M_PI;
  const double NU_E0 = GAMMA_E;
  const double MU_E1 = ALPHA_E + M_PI;
  const double NU_E1 = -GAMMA_E;

  /* these grism surface angles are still unclear */
  const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
  const double NU_G1 = ALPHA_G;
  const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
  //const double NU_G0 = cradians(0.0);
  //const double NU_G1 = ALPHA_G;
  //const double NU_G2 = -ALPHA_G;
  const double COS_MU_E0 = cos(MU_E0);
  const double SIN_MU_E0 = sin(MU_E0);
  const double COS_MU_E1 = cos(MU_E1);
  const double SIN_MU_E1 = sin(MU_E1);
  const double COS_NU_E0 = cos(NU_E0);
  const double SIN_NU_E0 = sin(NU_E0);
  const double COS_NU_E1 = cos(NU_E1);
  const double SIN_NU_E1 = sin(NU_E1);
  const double COS_NU_G0 = cos(NU_G0);
  const double SIN_NU_G0 = sin(NU_G0);
  const double COS_NU_G1 = cos(NU_G1);
  const double SIN_NU_G1 = sin(NU_G1);
  const double COS_NU_G2 = cos(NU_G2);
  const double SIN_NU_G2 = sin(NU_G2);
  const double XD02 = (double)NXPIX/2.0;
  const double YD02 = (double)NYPIX/2.0;
  const double LGAP = (gapoffset - gap) / 2.0; /* [mm] */
  const double RGAP = (gapoffset + gap) / 2.0; /* [mm] */
  const double TOMLE = -XD02*DPIX + LGAP; /* [mm] */
  const double TOMRE = LGAP; /* [mm] */
  const double JERLE = RGAP; /* [mm] */
  const double JERRE = XD02*DPIX + RGAP; /* [mm] */
  /* initialize temporary variables */
  double veclen;
  double li;
  double wi;
  double n_g;
  double xi, xx, x1, x2, x21, x3, x4, x5, xd;
  double yi, yy, y1, y2, y21, y3, y4, y5, yd;
  double zz, z1, z2, z3, z21, z4, z5;
  double beta;
  double blaze_eff;
  // unsigned long xbin; /* x coordinate [pix] */
  // unsigned long ybin; /* y coordinate [pix] */
  long xbin; /* x coordinate [pix] */
  long ybin; /* y coordinate [pix] */
  unsigned long i;   /* wavelength iterator */
  unsigned long j;   /* slit iterator */

  for (j=0; j<nslit; ++j) {
    xi = x[j];       /* x coord */
    yi = y[j];       /* y coord */
    for (i=0; i<nw; ++i) {
      li = lamb[i];    /* wavelength */
      wi = w[i];       /* intensities */
      n_g = n_sell[i]; /* refractive indices */

      /* LAUNCH RAYS (create normalized vectors) */
      veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
      xx = xi / veclen;
      yy = yi / veclen;
      zz = F_COL / veclen;

      /* BLAZE FUNCTION */
      //xb = xx;
      /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
      /* INTO ECHELLE RF */
      x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
      y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
      z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
      /* NORMALIZATION AFTER ECHELLE RELATION */
      z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);

      /* OUT OF ECHELLE RF */
      x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
      y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
      z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

      /* BLAZE FUNCTION */
      beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
      blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);

      /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
      /* INTO PLANE RF */
      x21 = x2;
      y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
      z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
      /* PLANE RELATION */
      x3 = x21 / n_g;
      y3 = y21 / n_g;
      /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
      z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);

      /* INTO GRISM RF and GRISM RELATION */
      x4 = x3 * n_g;
      y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
      z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
      /* NORMALIZATION AFTER GRATING RELATION */
      z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);

      /* OUT OF GRISM RF */
      x5 = x4;
      y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
      z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

      /* PROJECTION ONTO DETECTOR */
      //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
      //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
      //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
      //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
      //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;

      xd = x5/z5 * F_CAM + xd_0;
      yd = (y5/z5 * F_CAM + yd_0);
      //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);

      // NIR DETECTOR GAP
      if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
        xd -= LGAP;
	      /* BIN PIXELS */
	      xd /= DPIX;
	      yd /= DPIX;
	      xbin = (int)floor(xd + XD02);
	      if (xbin >= 0 && xbin < NXPIX) {
	        ybin = (int)floor(yd + YD02);
	        if (ybin >= 0 && ybin < NYPIX) {
	          if (BLAZE_FLAG == 1) {
	            ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	          }
	          else {
	            ccd[xbin+NXPIX*ybin] += wi;
							counts[xbin+NXPIX*ybin] += 1;
	          }
	        }
	      }
      }
      else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
        xd -= RGAP;
	      /* BIN PIXELS */
	      xd /= DPIX;
	      yd /= DPIX;
	      xbin = (int)floor(xd + XD02);
	      if (xbin >= 0 && xbin < NXPIX) {
	        ybin = (int)floor(yd + YD02);
	        if (ybin >= 0 && ybin < NYPIX) {
	          if (BLAZE_FLAG == 1) {
	            ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	          }
	          else {
	            ccd[xbin+NXPIX*ybin] += wi;
							counts[xbin+NXPIX*ybin] += 1;
	          }
	        }
	      }
      }
			else if (GAP_FLAG == 0) {
	      /* BIN PIXELS */
	      xd /= DPIX;
	      yd /= DPIX;
	      xbin = (int)floor(xd + XD02);
	      if (xbin >= 0 && xbin < NXPIX) {
	        ybin = (int)floor(yd + YD02);
	        if (ybin >= 0 && ybin < NYPIX) {
	          if (BLAZE_FLAG == 1) {
	            ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	          }
	          else {
	            ccd[xbin+NXPIX*ybin] += wi;
							counts[xbin+NXPIX*ybin] += 1;
	          }
	        }
	      }
			}
    }
		progress_bar(j, nslit);
	}
	printf("\n");
}


void compute_c_grid_perturb(
    int ARM,
    int BLAZE_FLAG,
    int GAP_FLAG,
    unsigned long nw,          /* size of w (and also x and y) */
    unsigned long nslit,         /* number of slit elements */
    unsigned short m,
    double gap,
    double gapoffset,
    double xd_0,
    double yd_0,
    double perturb, /* perturbation width of slit location */
    double* restrict n_sell,
    double* restrict x,   /* slit location */
    double* restrict y,   /* slit location */
    double* restrict lamb,      /* wavelengths */
    double* restrict w,         /* intensities */
    double* restrict ccd,
		unsigned long* restrict counts) {

  /* assign constants according to spectral arm */
  const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
  const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
  const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
  const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
  const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
  const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
  const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
  const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
  const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
  const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;

  /* pre-calculate constants */
  const double ORDER = (double)m;
  const double F_COL_2 = F_COL * F_COL;
  const double MU_E0 = ALPHA_E - M_PI;
  const double NU_E0 = GAMMA_E;
  const double MU_E1 = ALPHA_E + M_PI;
  const double NU_E1 = -GAMMA_E;
  /* these grism surface angles are still unclear */
  const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
  const double NU_G1 = ALPHA_G;
  const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
  //const double NU_G0 = cradians(0.0);
  //const double NU_G1 = ALPHA_G;
  //const double NU_G2 = -ALPHA_G;
  const double COS_MU_E0 = cos(MU_E0);
  const double SIN_MU_E0 = sin(MU_E0);
  const double COS_MU_E1 = cos(MU_E1);
  const double SIN_MU_E1 = sin(MU_E1);
  const double COS_NU_E0 = cos(NU_E0);
  const double SIN_NU_E0 = sin(NU_E0);
  const double COS_NU_E1 = cos(NU_E1);
  const double SIN_NU_E1 = sin(NU_E1);
  const double COS_NU_G0 = cos(NU_G0);
  const double SIN_NU_G0 = sin(NU_G0);
  const double COS_NU_G1 = cos(NU_G1);
  const double SIN_NU_G1 = sin(NU_G1);
  const double COS_NU_G2 = cos(NU_G2);
  const double SIN_NU_G2 = sin(NU_G2);
  const double XD02 = (double)NXPIX/2.0;
  const double YD02 = (double)NYPIX/2.0;
  const double LGAP = (gapoffset - gap) / 2.0; /* [mm] */
  const double RGAP = (gapoffset + gap) / 2.0; /* [mm] */
  const double TOMLE = -XD02*DPIX + LGAP; /* [mm] */
  const double TOMRE = LGAP; /* [mm] */
  const double JERLE = RGAP; /* [mm] */
  const double JERRE = XD02*DPIX + RGAP; /* [mm] */
  /* initialize temporary variables */
  double veclen;
  double li;
  double wi;
  double n_g;
  double xi, xx, x1, x2, x21, x3, x4, x5, xd;
  double yi, yy, y1, y2, y21, y3, y4, y5, yd;
  double zz, z1, z2, z3, z21, z4, z5;
  double beta;
  double blaze_eff;
  long xbin; /* x coordinate [pix] */
  long ybin; /* y coordinate [pix] */
  /* index variables */
  unsigned long i;   /* wavelength iterator */
  unsigned long j;   /* slit iterator */

  // RANDOM NUMBER GENERATOR INITIALIZATION
  long seed = time(NULL);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus); /* global generator */
  gsl_rng_set(r, seed); /* seed the rng */

  for (i=0; i<nw; ++i) {
    li = lamb[i];    /* wavelength */
    wi = w[i];       /* intensities */
    n_g = n_sell[i]; /* refractive indices */

    for (j=0; j<nslit; ++j) {
      xi = x[j] + gsl_ran_gaussian_ziggurat(r, perturb); /* randomize x coord */
      yi = y[j] + gsl_ran_gaussian_ziggurat(r, perturb); /* randomize y coord */

      /* LAUNCH RAYS (create normalized vectors) */
      veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
      xx = xi / veclen;
      yy = yi / veclen;
      zz = F_COL / veclen;

      /* BLAZE FUNCTION */
      //xb = xx;

      /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
      /* INTO ECHELLE RF */
      x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
      y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
      z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
      /* NORMALIZATION AFTER ECHELLE RELATION */
      z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);

      /* OUT OF ECHELLE RF */
      x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
      y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
      z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

      /* BLAZE FUNCTION */
      beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
      blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);

      /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
      /* INTO PLANE RF */
      x21 = x2;
      y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
      z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
      /* PLANE RELATION */
      x3 = x21 / n_g;
      y3 = y21 / n_g;
      /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
      z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);

      /* INTO GRISM RF and GRISM RELATION */
      x4 = x3 * n_g;
      y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
      z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
      /* NORMALIZATION AFTER GRATING RELATION */
      z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);

      /* OUT OF GRISM RF */
      x5 = x4;
      y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
      z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
      /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

      /* PROJECTION ONTO DETECTOR */
      //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
      //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
      //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
      //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
      //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;

      xd = x5/z5 * F_CAM + xd_0;
      yd = (y5/z5 * F_CAM + yd_0);
      //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);

      // NIR DETECTOR GAP
      if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
        xd -= LGAP;
	      xd /= DPIX;
	      yd /= DPIX;
	      xbin = (int)floor(xd + XD02);
	      if (xbin >= 0 && xbin < NXPIX) {
	        ybin = (int)floor(yd + YD02);
	        if (ybin >= 0 && ybin < NYPIX) {
	          if (BLAZE_FLAG == 1) {
	            ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	          }
	          else {
	            ccd[xbin+NXPIX*ybin] += wi;
							counts[xbin+NXPIX*ybin] += 1;
	          }
	        }
	      }
      }
      else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
        xd -= RGAP;
	      xd /= DPIX;
	      yd /= DPIX;
	      xbin = (int)floor(xd + XD02);
	      if (xbin >= 0 && xbin < NXPIX) {
	        ybin = (int)floor(yd + YD02);
	        if (ybin >= 0 && ybin < NYPIX) {
	          if (BLAZE_FLAG == 1) {
	            ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	          }
	          else {
	            ccd[xbin+NXPIX*ybin] += wi;
							counts[xbin+NXPIX*ybin] += 1;
	          }
	        }
	      }
      }
			else if (GAP_FLAG == 0) {
	      xd /= DPIX;
	      yd /= DPIX;
	      xbin = (int)floor(xd + XD02);
	      if (xbin >= 0 && xbin < NXPIX) {
	        ybin = (int)floor(yd + YD02);
	        if (ybin >= 0 && ybin < NYPIX) {
	          if (BLAZE_FLAG == 1) {
	            ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
	          }
	          else {
	            ccd[xbin+NXPIX*ybin] += wi;
							counts[xbin+NXPIX*ybin] += 1;
	          }
	        }
	      }
			}
			progress_bar(i, nw);
    }
	}
	printf("\n");
  gsl_rng_free(r);
}


/*
  Spectrum treated as Probability Distribution Function (PDF) using
  Rejection Sampling techniques
*/
void compute_c_mc_rejsamp(
    int ARM,
    int BLAZE_FLAG,
    int GAP_FLAG,
    int SLIT_FLAG,
    unsigned long nw,          /* size of w (and also x and y) */
    unsigned long nphotons,    /* number of photons to simulate */
    unsigned short m,
    double gap,
    double gapoffset,    
    double xd_0,
    double yd_0,
    double fluxmax,
    double offset,
    double slit_locx,					// x- Gaussian location
		double slit_locy,					// y- Gaussian location
		double slit_scalex,				// x- Gaussian scale
		double slit_scaley,				// y- Gaussian scale
    double* restrict lamb,      /* wavelengths */
    double* restrict w,         /* intensities */
    unsigned long* restrict ccd,
		double* restrict wavemap) {

  /* assign constants according to spectral arm */
  const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
  const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
  const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
  const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
  const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
  const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
  const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
  const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
  const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
  const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;
  const double FN = (ARM == 0) ? NIR_FN : VIS_FN;

  /* pre-calculate constants */
  const double ORDER = (double)m;
  const double F_COL_2 = F_COL * F_COL;
  const double MU_E0 = ALPHA_E - M_PI;
  const double NU_E0 = GAMMA_E;
  const double MU_E1 = ALPHA_E + M_PI;
  const double NU_E1 = -GAMMA_E;
  /* these grism surface angles are still unclear */
  const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
  const double NU_G1 = ALPHA_G;
  const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
  //const double NU_G0 = cradians(0.0);
  //const double NU_G1 = ALPHA_G;
  //const double NU_G2 = -ALPHA_G;
  const double COS_MU_E0 = cos(MU_E0);
  const double SIN_MU_E0 = sin(MU_E0);
  const double COS_MU_E1 = cos(MU_E1);
  const double SIN_MU_E1 = sin(MU_E1);
  const double COS_NU_E0 = cos(NU_E0);
  const double SIN_NU_E0 = sin(NU_E0);
  const double COS_NU_E1 = cos(NU_E1);
  const double SIN_NU_E1 = sin(NU_E1);
  const double COS_NU_G0 = cos(NU_G0);
  const double SIN_NU_G0 = sin(NU_G0);
  const double COS_NU_G1 = cos(NU_G1);
  const double SIN_NU_G1 = sin(NU_G1);
  const double COS_NU_G2 = cos(NU_G2);
  const double SIN_NU_G2 = sin(NU_G2);
  const double XD02 = (double)NXPIX/2.0;
  const double YD02 = (double)NYPIX/2.0;
  const double LGAP = (gapoffset - gap) / 2.0; /* [mm] */
  const double RGAP = (gapoffset + gap) / 2.0; /* [mm] */
  const double TOMLE = -XD02*DPIX + LGAP; /* [mm] */
  const double TOMRE = LGAP; /* [mm] */
  const double JERLE = RGAP; /* [mm] */
  const double JERRE = XD02*DPIX + RGAP; /* [mm] */  
  const double wmin = lamb[0];
  const double wmax = lamb[nw - 1];

  /* initialize temporary variables */
  double veclen;
  double li;
  double n_g;
  double xi, xx, x1, x2, x21, x3, x4, x5, xd;
  double yi, yy, y1, y2, y21, y3, y4, y5, yd;
  double zz, z1, z2, z3, z21, z4, z5;
  double beta;
  double blaze_eff;
  double u, fc, lc, alpha; // rejection sampling variables
  long xbin; /* x coordinate [pix] */
  long ybin; /* y coordinate [pix] */

  /* index variables */
  unsigned long i;   /* wavelength iterator */
  unsigned int k;   /* rej samp flag */

  // RANDOM NUMBER GENERATOR INITIALIZATION
  long seed = time(NULL);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus); /* global rng generator */
  gsl_rng_set(r, seed); /* seed the rng */

  // CUBIC SPLINE OF SPECTRA
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_akima, nw);
  gsl_spline_init(spline, lamb, w, nw);

  for (i=0; i<nphotons; ++i) {
    // REJECTION SAMPLING ALGORITHM FOR WAVELENGTHS
    k = 0;
    while (k != 1) {
      lc = gsl_rng_uniform(r) * (wmax - wmin) + wmin;
      u = gsl_rng_uniform(r);
      fc = gsl_spline_eval(spline, lc, acc);
      alpha = fc / fluxmax;
      if (alpha >= u)
        k = 1;
    }
    li = lc;    /* wavelength */
    n_g = nsell(ARM, lc);
    double outarr[2] = {0.0, 0.0};
    if (SLIT_FLAG == 0) {
      slit_image(DFIBER, FN, TAU_0, HEIGHT_RATIO, offset, r, outarr);
    }
    else if (SLIT_FLAG == 1) {
      slit_image_gaussian(DFIBER, FN, TAU_0, HEIGHT_RATIO, offset, slit_locx, slit_locy, slit_scalex, slit_scaley, r, outarr);
    }
    xi = outarr[0];
    yi = outarr[1];

    /* LAUNCH RAYS (create normalized vectors) */
    veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
    xx = xi / veclen;
    yy = yi / veclen;
    zz = F_COL / veclen;

    /* BLAZE FUNCTION */
    //xb = xx;
    /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
    /* INTO ECHELLE RF */
    x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
    y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
    z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
    /* NORMALIZATION AFTER ECHELLE RELATION */
    z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);

    /* OUT OF ECHELLE RF */
    x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
    y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
    z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
    /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    /* BLAZE FUNCTION */
    beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
    blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);

    /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
    /* INTO PLANE RF */
    x21 = x2;
    y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
    z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
    /* PLANE RELATION */
    x3 = x21 / n_g;
    y3 = y21 / n_g;
    /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
    z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);

    /* INTO GRISM RF and GRISM RELATION */
    x4 = x3 * n_g;
    y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
    z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
    /* NORMALIZATION AFTER GRATING RELATION */
    z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);

    /* OUT OF GRISM RF */
    x5 = x4;
    y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
    z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
    /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    /* PROJECTION ONTO DETECTOR */
    //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
    //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
    //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
    //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
    //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;
    xd = x5/z5 * F_CAM + xd_0;
    yd = (y5/z5 * F_CAM + yd_0);
    //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);
    
    // NIR DETECTOR GAP
    if (GAP_FLAG == 1 && xd > TOMLE && xd < TOMRE && ARM == 0) {
      xd -= LGAP;
	    xd /= DPIX;
	    yd /= DPIX;
	    xbin = (int)floor(xd + XD02);
	    if (xbin >= 0 && xbin < NXPIX) {
	      ybin = (int)floor(yd + YD02);
	      if (ybin >= 0 && ybin < NYPIX) {
	        if (BLAZE_FLAG == 1) {
	          u = gsl_rng_uniform(r);
	          if (u <= blaze_eff) {
	            ccd[xbin+NXPIX*ybin] += 1;
	          }
	        }
	        else {
	          ccd[xbin+NXPIX*ybin] += 1;
	        }
	      }
	    }
    }
    else if (GAP_FLAG == 1 && xd > JERLE && xd < JERRE && ARM == 0) {
      xd -= RGAP;
	    xd /= DPIX;
	    yd /= DPIX;
	    xbin = (int)floor(xd + XD02);
	    if (xbin >= 0 && xbin < NXPIX) {
	      ybin = (int)floor(yd + YD02);
	      if (ybin >= 0 && ybin < NYPIX) {
	        if (BLAZE_FLAG == 1) {
	          u = gsl_rng_uniform(r);
	          if (u <= blaze_eff) {
	            ccd[xbin+NXPIX*ybin] += 1;
	          }
	        }
	        else {
	          ccd[xbin+NXPIX*ybin] += 1;
						wavemap[xbin+NXPIX*ybin] += li;
	        }
	      }
	    }
    }
    else if (GAP_FLAG == 0) {
	    xd /= DPIX;
	    yd /= DPIX;
	    xbin = (int)floor(xd + XD02);
	    if (xbin >= 0 && xbin < NXPIX) {
	      ybin = (int)floor(yd + YD02);
	      if (ybin >= 0 && ybin < NYPIX) {
	        if (BLAZE_FLAG == 1) {
	          u = gsl_rng_uniform(r);
	          if (u <= blaze_eff) {
	            ccd[xbin+NXPIX*ybin] += 1;
	          }
	        }
	        else {
	          ccd[xbin+NXPIX*ybin] += 1;
						wavemap[xbin+NXPIX*ybin] += li;
	        }
	      }
	    }
		}
		progress_bar(i, nphotons);
	}
	printf("\n");
  gsl_spline_free (spline);
  gsl_interp_accel_free(acc);
  gsl_rng_free(r);
}


/*
  Spectrum treated as Probability Distribution Function (PDF) using
  Cumulative Distribution Function (CDF) techniques.
*/
void compute_c_mc_cdf(
    int ARM,
    int BLAZE_FLAG,
    int GAP_FLAG,
    int SLIT_FLAG,
    unsigned long nphotons,    /* number of photons to simulate */
    unsigned short m,
    double gap,
    double gapoffset,
    double xd_0,
    double yd_0,
    double offset,
    double slit_locx,					// x- Gaussian location
		double slit_locy,					// y- Gaussian location
		double slit_scalex,				// x- Gaussian scale
		double slit_scaley,				// y- Gaussian scale
    double* restrict lamb,      /* wavelengths */
    unsigned long* restrict ccd,
		double* restrict wavemap) {
  
  /* assign constants according to spectral arm */
  const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
  const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
  const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
  const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
  const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
  const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
  const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
  const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
  const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
  const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;
  const double FN = (ARM == 0) ? NIR_FN : VIS_FN;
  /* pre-calculate constants */
  const double ORDER = (double)m;
  const double F_COL_2 = F_COL * F_COL;
  const double MU_E0 = ALPHA_E - M_PI;
  const double NU_E0 = GAMMA_E;
  const double MU_E1 = ALPHA_E + M_PI;
  const double NU_E1 = -GAMMA_E;
  /* these grism surface angles are still unclear */
  const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
  const double NU_G1 = ALPHA_G;
  const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
  //const double NU_G0 = cradians(0.0);
  //const double NU_G1 = ALPHA_G;
  //const double NU_G2 = -ALPHA_G;
  const double COS_MU_E0 = cos(MU_E0);
  const double SIN_MU_E0 = sin(MU_E0);
  const double COS_MU_E1 = cos(MU_E1);
  const double SIN_MU_E1 = sin(MU_E1);
  const double COS_NU_E0 = cos(NU_E0);
  const double SIN_NU_E0 = sin(NU_E0);
  const double COS_NU_E1 = cos(NU_E1);
  const double SIN_NU_E1 = sin(NU_E1);
  const double COS_NU_G0 = cos(NU_G0);
  const double SIN_NU_G0 = sin(NU_G0);
  const double COS_NU_G1 = cos(NU_G1);
  const double SIN_NU_G1 = sin(NU_G1);
  const double COS_NU_G2 = cos(NU_G2);
  const double SIN_NU_G2 = sin(NU_G2);
  const double XD02 = (double)NXPIX/2.0; // pix
  const double YD02 = (double)NYPIX/2.0; // pix
  const double LGAP = gapoffset - (gap / 2.0); /* [mm] */
  const double RGAP = gapoffset + (gap / 2.0); /* [mm] */
  const double TOMLE = -XD02*DPIX + LGAP; /* [mm] */
  const double TOMRE = LGAP; /* [mm] */
  const double JERLE = RGAP; /* [mm] */
  const double JERRE = XD02*DPIX + RGAP; /* [mm] */
  /* initialize temporary variables */
  double veclen;
  double li;
  double n_g;
  double xi, xx, x1, x2, x21, x3, x4, x5, xd;
  double yi, yy, y1, y2, y21, y3, y4, y5, yd;
  double zz, z1, z2, z3, z21, z4, z5;
  double beta;
  double blaze_eff;
  double u; // random number
  long xbin; /* x coordinate [pix] */
  long ybin; /* y coordinate [pix] */
  unsigned long i;   /* wavelength iterator */
  // RANDOM NUMBER GENERATOR INITIALIZATION
  long seed = time(NULL);
  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus); /* global rng generator */
  gsl_rng_set(r, seed); /* seed the rng */

  for (i=0; i<nphotons; ++i) {
    li = lamb[i];    /* wavelength */
    n_g = nsell(ARM, li);
    double outarr[2] = {0.0, 0.0};
    if (SLIT_FLAG == 0) {
      slit_image(DFIBER, FN, TAU_0, HEIGHT_RATIO, offset, r, outarr);
    }
    else if (SLIT_FLAG == 1) {
      slit_image_gaussian(DFIBER, FN, TAU_0, HEIGHT_RATIO, offset, slit_locx, slit_locy, slit_scalex, slit_scaley, r, outarr);
    }
    xi = outarr[0];
    yi = outarr[1];

    /* LAUNCH RAYS (create normalized vectors) */
    veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
    xx = xi / veclen;
    yy = yi / veclen;
    zz = F_COL / veclen;
    /* BLAZE FUNCTION */
    //xb = xx;
    /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
    /* INTO ECHELLE RF */
    x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
    y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
    z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
    /* NORMALIZATION AFTER ECHELLE RELATION */
    z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);

    /* OUT OF ECHELLE RF */
    x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
    y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
    z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
    /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    /* BLAZE FUNCTION */
    beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
    blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);

    /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
    /* INTO PLANE RF */
    x21 = x2;
    y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
    z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
    /* PLANE RELATION */
    x3 = x21 / n_g;
    y3 = y21 / n_g;
    /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
    z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);

    /* INTO GRISM RF and GRISM RELATION */
    x4 = x3 * n_g;
    y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
    z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
    /* NORMALIZATION AFTER GRATING RELATION */
    z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);

    /* OUT OF GRISM RF */
    x5 = x4;
    y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
    z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
    /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    /* PROJECTION ONTO DETECTOR */
    //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
    //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
    //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
    //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
    //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;

    xd = x5/z5 * F_CAM + xd_0;
    yd = (y5/z5 * F_CAM + yd_0);
    //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);

    // NIR DETECTOR GAP
    if ( (GAP_FLAG == 1) && (xd > TOMLE && xd < TOMRE) && (ARM == 0)) {
      xd -= LGAP;
	    xd /= DPIX;
	    yd /= DPIX;
	    xbin = (int)floor(xd + XD02);
	    if (xbin >= 0 && xbin < NXPIX) {
	      ybin = (int)floor(yd + YD02);
	      if (ybin >= 0 && ybin < NYPIX) {
	        if (BLAZE_FLAG == 1) {
	          u = gsl_rng_uniform(r);
	          if (u <= blaze_eff) {
	            ccd[xbin+NXPIX*ybin] += 1;
	          }
	        }
	        else {
	          ccd[xbin+NXPIX*ybin] += 1;
	        }
	      }
	    }
    }
    else if ((GAP_FLAG == 1) && (xd > JERLE && xd < JERRE) && (ARM == 0)) {
      xd -= RGAP;
	    xd /= DPIX;
	    yd /= DPIX;
	    xbin = (int)floor(xd + XD02);
	    if (xbin >= 0 && xbin < NXPIX) {
	      ybin = (int)floor(yd + YD02);
	      if (ybin >= 0 && ybin < NYPIX) {
	        if (BLAZE_FLAG == 1) {
	          u = gsl_rng_uniform(r);
	          if (u <= blaze_eff) {
	            ccd[xbin+NXPIX*ybin] += 1;
	          }
	        }
	        else {
	          ccd[xbin+NXPIX*ybin] += 1;
						wavemap[xbin+NXPIX*ybin] += li;
	        }
	      }
	    }
    }
		else if (GAP_FLAG == 0) {
	    xd /= DPIX;
	    yd /= DPIX;
	    xbin = (int)floor(xd + XD02);
	    if (xbin >= 0 && xbin < NXPIX) {
	      ybin = (int)floor(yd + YD02);
	      if (ybin >= 0 && ybin < NYPIX) {
	        if (BLAZE_FLAG == 1) {
	          u = gsl_rng_uniform(r);
	          if (u <= blaze_eff) {
	            ccd[xbin+NXPIX*ybin] += 1;
	          }
	        }
	        else {
	          ccd[xbin+NXPIX*ybin] += 1;
						wavemap[xbin+NXPIX*ybin] += li;
	        }
	      }
	    }
    }
		progress_bar(i, nphotons);
	}
	printf("\n");
  gsl_rng_free(r);
}


// void compute_c_gap(int ARM,
//                int BLAZE_FLAG,
//                int LOC_FLAG,
//                unsigned long nw,          /* size of w (and also x and y) */
//                unsigned long nslit,         /* number of slit elements */
//                unsigned short m,
//                double gap,
//                double gapoffset,
//                double xd_0,
//                double yd_0,
//                double *n_sell,
//                double *x,   /* slit location */
//                double *y,   /* slit location */
//                double *lamb,      /* wavelengths */
//                double *w,         /* intensities */
//                double *ccd,
//                unsigned long *counts,
//                unsigned long *m_list,
//                double *returnx,
//                double *returny) {
// 
//     /* assign constants according to spectral arm */
//     const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
//     const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
//     const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
//     const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
//     const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
//     const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
//     const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
//     const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
//     const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
//     const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;
// 
//     /* pre-calculate constants */
//     const double ORDER = (double)m;
//     const double F_COL_2 = F_COL * F_COL;
// 
//     const double MU_E0 = ALPHA_E - M_PI;
//     const double NU_E0 = GAMMA_E;
// 
//     const double MU_E1 = ALPHA_E + M_PI;
//     const double NU_E1 = -GAMMA_E;
// 
//     /* these grism surface angles are still unclear */
//     const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
//     const double NU_G1 = ALPHA_G;
//     const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
//     //const double NU_G0 = cradians(0.0);
//     //const double NU_G1 = ALPHA_G;
//     //const double NU_G2 = -ALPHA_G;
// 
//     const double COS_MU_E0 = cos(MU_E0);
//     const double SIN_MU_E0 = sin(MU_E0);
// 
//     const double COS_MU_E1 = cos(MU_E1);
//     const double SIN_MU_E1 = sin(MU_E1);
// 
//     const double COS_NU_E0 = cos(NU_E0);
//     const double SIN_NU_E0 = sin(NU_E0);
// 
//     const double COS_NU_E1 = cos(NU_E1);
//     const double SIN_NU_E1 = sin(NU_E1);
// 
//     const double COS_NU_G0 = cos(NU_G0);
//     const double SIN_NU_G0 = sin(NU_G0);
// 
//     const double COS_NU_G1 = cos(NU_G1);
//     const double SIN_NU_G1 = sin(NU_G1);
// 
//     const double COS_NU_G2 = cos(NU_G2);
//     const double SIN_NU_G2 = sin(NU_G2);
//     const double XD02 = (double)NXPIX/2.0;
//     const double YD02 = (double)NYPIX/2.0;
//     const double LGAP = (gapoff - gap) / 2.0; /* [mm] */
//     const double RGAP = (gapoff + gap) / 2.0; /* [mm] */
//     const double TOMLE = -XD02*DPIX + LGAP; /* [mm] */
//     const double TOMRE = LGAP; /* [mm] */
//     const double JERLE = RGAP; /* [mm] */
//     const double JERRE = XD02*DPIX + RGAP; /* [mm] */
// 
//     /* initialize temporary variables */
//     double veclen;
//     double li;
//     double wi;
//     double n_g;
// 
//     double xi, xx, x1, x2, x21, x3, x4, x5, xd;
//     double yi, yy, y1, y2, y21, y3, y4, y5, yd;
//     double zz, z1, z2, z3, z21, z4, z5;
// 
//     double beta;
//     double blaze_eff;
// 
//     unsigned long xbin; /* x coordinate [pix] */
//     unsigned long ybin; /* y coordinate [pix] */
// 
//     /* index variables */
// //     unsigned long nind = nw * nslit; /* 1D index */
// //     unsigned long ind;
//     unsigned long i;   /* wavelength iterator */
//     unsigned long j;   /* slit iterator */
// 
//     //printf("\nTOM LEFT = %.2f\n", TOMLE);
//     //printf("TOM RIGHT = %.2f\n", TOMRE);
//     //printf("JERRY LEFT = %.2f\n", JERLE);
//     //printf("JERRY RIGHT = %.2f\n", JERRE);
// 
//     for (j=0; j<nslit; ++j)
//     {
//       xi = x[j];       /* x coord */
//       yi = y[j];       /* y coord */
//       for (i=0; i<nw; ++i)
//       {
// //         i = ind / nslit; /* slit iterator */
// //         j = ind % nslit; /* wave iterator */
// 
//         li = lamb[i];    /* wavelength */
//         wi = w[i];       /* intensities */
//         n_g = n_sell[i]; /* refractive indices */
// 
//         /* LAUNCH RAYS (create normalized vectors) */
//         veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
//         xx = xi / veclen;
//         yy = yi / veclen;
//         zz = F_COL / veclen;
// 
//         /* BLAZE FUNCTION */
//         //xb = xx;
// 
//         /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
//         /* INTO ECHELLE RF */
//         x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
//         y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
//         z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
//         /* NORMALIZATION AFTER ECHELLE RELATION */
//         z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);
// 
//         /* OUT OF ECHELLE RF */
//         x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
//         y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
//         z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
//         /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// 
// 
//         /* BLAZE FUNCTION */
//         beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
//         blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);
// 
// 
//         /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
//         /* INTO PLANE RF */
//         x21 = x2;
//         y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
//         z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
//         /* PLANE RELATION */
//         x3 = x21 / n_g;
//         y3 = y21 / n_g;
//         /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
//         z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);
// 
//         /* INTO GRISM RF and GRISM RELATION */
//         x4 = x3 * n_g;
//         y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
//         z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
//         /* NORMALIZATION AFTER GRATING RELATION */
//         z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);
// 
//         /* OUT OF GRISM RF */
//         x5 = x4;
//         y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
//         z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
//         /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// 
// 
//         /* PROJECTION ONTO DETECTOR */
//         //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
//         //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
//         //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
//         //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
//         //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;
// 
//         xd = x5/z5 * F_CAM + xd_0;
//         yd = (y5/z5 * F_CAM + yd_0);
//         //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);
// 
//         /* BIN PIXELS */
//         //xd /= DPIX;
//         //yd /= DPIX;
//         //xbin = (int)floor(xd + XD02);
// 
//         if (xd > TOMLE && xd < TOMRE)
//         {
//             xd -= LGAP;
//             xd /= DPIX;
//             yd /= DPIX;
//             xbin = (int)floor(xd + XD02);
//             ybin = (int)floor(yd + YD02);
//             if (ybin >= 0 && ybin < NYPIX)
//             {
//                 if (BLAZE_FLAG == 1)
//                 {
//                     ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
//                     counts[xbin+NXPIX*ybin] += 1;
//                     m_list[xbin+NXPIX*ybin] = m;
//                 }
//                 else
//                 {
//                     ccd[xbin+NXPIX*ybin] += wi;
//                     counts[xbin+NXPIX*ybin] += 1;
//                     m_list[xbin+NXPIX*ybin] = m;
//                 }
//             }
//         }
// 
//         else if (xd > JERLE && xd < JERRE)
//         {
//             xd -= RGAP;
//             xd /= DPIX;
//             yd /= DPIX;
//             xbin = (int)floor(xd + XD02);
//             ybin = (int)floor(yd + YD02);
//             if (ybin >= 0 && ybin < NYPIX)
//             {
//                 if (BLAZE_FLAG == 1)
//                 {
//                     ccd[xbin+NXPIX*ybin] += wi * blaze_eff;
//                     counts[xbin+NXPIX*ybin] += 1;
//                     m_list[xbin+NXPIX*ybin] = m;
//                 }
//                 else
//                 {
//                     ccd[xbin+NXPIX*ybin] += wi;
//                     counts[xbin+NXPIX*ybin] += 1;
//                     m_list[xbin+NXPIX*ybin] = m;
//                 }
//             }
//         }
// 
// 
//     }}
// }


// void
// compute_c_pdf(int ARM,
//               int BLAZE_FLAG,
//               int LOC_FLAG,
//               unsigned long nw,          /* size of w (and also x and y) */
//               unsigned long nslit,         /* number of slit elements */
//               unsigned short m,
//               double xd_0,
//               double yd_0,
//               double perturb, /* perturbation width of slit location */
//               double fluxmax,
//               double* restrict n_sell,
//               double* restrict x,   /* slit location */
//               double* restrict y,   /* slit location */
//               double* restrict lamb,      /* wavelengths */
//               double* restrict w,         /* intensities */
//               double* restrict ccd,
//               unsigned long* restrict counts,
//               unsigned long* restrict m_list,
//               double* restrict returnx,
//               double* restrict returny)
//         /* Spectrum treated as Probability Distribution Function (PDF) using
//              Rejection Sampling techniques */
//         {
//     /* assign constants according to spectral arm */
//     const double F_COL = (ARM == 0) ? NIR_F_COL : VIS_F_COL;
//     const double ALPHA_E = (ARM == 0) ? cradians(NIR_ALPHA_E_DEG) : cradians(VIS_ALPHA_E_DEG);
//     const double GAMMA_E = (ARM == 0) ? cradians(NIR_GAMMA_E_DEG) : cradians(VIS_GAMMA_E_DEG);
//     const double SIGMA_E = (ARM == 0) ? NIR_SIGMA_E : VIS_SIGMA_E;
//     const double ALPHA_G = (ARM == 0) ? cradians(NIR_ALPHA_G_DEG) : cradians(VIS_ALPHA_G_DEG);
//     const double SIGMA_G = (ARM == 0) ? NIR_SIGMA_G : VIS_SIGMA_G;
//     const double F_CAM = (ARM == 0) ? NIR_F_CAM : VIS_F_CAM;
//     const double DPIX = (ARM == 0) ? NIR_DPIX : VIS_DPIX;
//     const int NXPIX = (ARM == 0) ? NIR_NXPIX : VIS_NXPIX;
//     const int NYPIX = (ARM == 0) ? NIR_NYPIX : VIS_NXPIX;
//
//     /* pre-calculate constants */
//     const double ORDER = (double)m;
//     const double F_COL_2 = F_COL * F_COL;
//     const double MU_E0 = ALPHA_E - M_PI;
//     const double NU_E0 = GAMMA_E;
//     const double MU_E1 = ALPHA_E + M_PI;
//     const double NU_E1 = -GAMMA_E;
//     /* these grism surface angles are still unclear */
//     const double NU_G0 = (ARM == 0) ? cradians(-8.4) : cradians(1.6);
//     const double NU_G1 = ALPHA_G;
//     const double NU_G2 = (ARM == 0) ? cradians(-11.8) : -ALPHA_G;
//     //const double NU_G0 = cradians(0.0);
//     //const double NU_G1 = ALPHA_G;
//     //const double NU_G2 = -ALPHA_G;
//     const double COS_MU_E0 = cos(MU_E0);
//     const double SIN_MU_E0 = sin(MU_E0);
//     const double COS_MU_E1 = cos(MU_E1);
//     const double SIN_MU_E1 = sin(MU_E1);
//     const double COS_NU_E0 = cos(NU_E0);
//     const double SIN_NU_E0 = sin(NU_E0);
//     const double COS_NU_E1 = cos(NU_E1);
//     const double SIN_NU_E1 = sin(NU_E1);
//     const double COS_NU_G0 = cos(NU_G0);
//     const double SIN_NU_G0 = sin(NU_G0);
//     const double COS_NU_G1 = cos(NU_G1);
//     const double SIN_NU_G1 = sin(NU_G1);
//     const double COS_NU_G2 = cos(NU_G2);
//     const double SIN_NU_G2 = sin(NU_G2);
//     const double XD02 = (double)NXPIX/2.0;
//     const double YD02 = (double)NYPIX/2.0;
//         const double wmin = lamb[0];
//         const double wmax = lamb[nw - 1];
//
//     /* initialize temporary variables */
//     double veclen;
//     double li;
//     double wi;
//     double n_g;
//     double xi, xx, x1, x2, x21, x3, x4, x5, xd;
//     double yi, yy, y1, y2, y21, y3, y4, y5, yd;
//     double zz, z1, z2, z3, z21, z4, z5;
//     double beta;
//     double blaze_eff;
//         double u, fc, lc, alpha; // rejection sampling variables
//     long xbin; /* x coordinate [pix] */
//     long ybin; /* y coordinate [pix] */
//
//     /* index variables */
//     unsigned long i;   /* wavelength iterator */
//     unsigned long j;   /* slit iterator */
//         unsigned int k;   /* rej samp iterator */
//
//         // RANDOM NUMBER GENERATOR INITIALIZATION
//         long seed = time(NULL);
//         gsl_rng* r = gsl_rng_alloc(gsl_rng_taus); /* global rng generator */
//         gsl_rng_set(r, seed); /* seed the rng */
//
//         // CUBIC SPLINE OF SPECTRA
//         gsl_interp_accel* acc = gsl_interp_accel_alloc();
//         gsl_spline* spline = gsl_spline_alloc(gsl_interp_akima, nw);
//         gsl_spline_init(spline, lamb, w, nw);
//
//         for (i=0; i<nw; ++i)
//     {
//       for (j=0; j<nslit; ++j)
//       {
//                 // REJECTION SAMPLING ALGORITHM FOR WAVELENGTHS
//                 k = 0;
//               while (k != 1)
//                 {
//                     lc = gsl_rng_uniform(r) * (wmax - wmin) + wmin;
//                     u = gsl_rng_uniform(r);
//                     fc = gsl_spline_eval(spline, lc, acc);
//                     alpha = fc / fluxmax;
//                     if (alpha >= u)
//                         k = 1;
//                 }
//           li = lc;    /* wavelength */
//                 n_g = nsell(ARM, lc);
//
//                 xi = x[j] + gsl_ran_gaussian_ziggurat(r, perturb); /* randomize x coord */
//                 yi = y[j] + gsl_ran_gaussian_ziggurat(r, perturb); /* randomize y coord */
//
//         /* LAUNCH RAYS (create normalized vectors) */
//         veclen = sqrt( xi*xi + yi*yi + F_COL_2 );
//         xx = xi / veclen;
//         yy = yi / veclen;
//         zz = F_COL / veclen;
//
//         /* BLAZE FUNCTION */
//         //xb = xx;
//
//         /* :::::::::::::::::::: ECHELLE :::::::::::::::::::::::::::::::::::::*/
//         /* INTO ECHELLE RF */
//         x1 = (ORDER * li / SIGMA_E) - (COS_MU_E0 * xx) + (SIN_MU_E0*SIN_NU_E0 * yy) + (SIN_MU_E0*COS_NU_E0 * zz);
//         y1 = -(COS_NU_E0 * yy) + (SIN_NU_E0 * zz);
//         z1 = -(SIN_MU_E0 * xx) - (COS_MU_E0*SIN_NU_E0 * yy) - (COS_MU_E0*COS_NU_E0 * zz);
//         /* NORMALIZATION AFTER ECHELLE RELATION */
//         z1 = (z1 / fabs(z1)) * sqrt(1.0 - y1*y1 - x1*x1);
//
//         /* OUT OF ECHELLE RF */
//         x2 = (COS_MU_E1 * x1) - (SIN_MU_E1 * z1);
//         y2 = -(SIN_MU_E1*SIN_NU_E1) * x1 + (COS_NU_E1 * y1 - COS_MU_E1*SIN_NU_E1 * z1);
//         z2 = (SIN_MU_E1*COS_NU_E1) * x1 + (SIN_NU_E1 * y1 + COS_MU_E1*COS_NU_E1 * z1);
//         /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
//
//         /* BLAZE FUNCTION */
//         beta = M_PI * cos(ALPHA_E) * SIGMA_E * (xx + x2) / li;
//         blaze_eff = (sin(beta)*sin(beta)) / (beta*beta);
//
//         /* :::::::::::::::::::: GRISM :::::::::::::::::::::::::::::::::::::::*/
//         /* INTO PLANE RF */
//         x21 = x2;
//         y21 = y2 * COS_NU_G0 - z2 * SIN_NU_G0;
//         z21 = y2 * SIN_NU_G0 + z2 * COS_NU_G0;
//         /* PLANE RELATION */
//         x3 = x21 / n_g;
//         y3 = y21 / n_g;
//         /* NORMALIZATION AFTER PLANE RELATION (REFRACTION) */
//         z3 = (z21 / fabs(z21)) * sqrt(1.0 - x3*x3 - y3*y3);
//
//         /* INTO GRISM RF and GRISM RELATION */
//         x4 = x3 * n_g;
//         y4 = li / (-SIGMA_G) + (n_g * COS_NU_G1 * y3) - (n_g * SIN_NU_G1 * z3);
//         z4 = (SIN_NU_G1 * y3) + (COS_NU_G1 * z3);
//         /* NORMALIZATION AFTER GRATING RELATION */
//         z4 = (z4 / fabs(z4)) * sqrt(1.0 - x4*x4 - y4*y4);
//
//         /* OUT OF GRISM RF */
//         x5 = x4;
//         y5 = (COS_NU_G2 * y4) - (SIN_NU_G2 * z4);
//         z5 = (SIN_NU_G2 * y4) + (COS_NU_G2 * z4);
//         /* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
//
//         /* PROJECTION ONTO DETECTOR */
//         //xd = (xd_0/F_CAM  + x5/z5) * (F_CAM / DPIX);
//         //yd = -yd_0/DPIX + y5 * (F_CAM/DPIX) / z5;
//         //xd = (xd_0 + F_CAM * x5 / z5) / DPIX;
//         //yd = -(yd_0 + F_CAM * y5 / z5) / DPIX;
//         //yd = -(yd_0 + F_CAM * y5 / sqrt(x5*x5 + z5*z5)) / DPIX;
//
//         xd = x5/z5 * F_CAM + xd_0;
//         yd = (y5/z5 * F_CAM + yd_0);
//         //yd = -(y5/sqrt(x5*x5 + z5*z5) * F_CAM + yd_0);
//
//         xd /= DPIX;
//         yd /= DPIX;
//         xbin = (int)floor(xd + XD02);
//         if (xbin >= 0 && xbin < NXPIX)
//         {
//             ybin = (int)floor(yd + YD02);
//             if (ybin >= 0 && ybin < NYPIX)
//             {
//                 if (BLAZE_FLAG == 1)
//                 {
//                                       u = gsl_rng_uniform(r);
//                                         if (u <= blaze_eff)
//                                         {
//                                             ccd[xbin+NXPIX*ybin] += 1;
//                     //counts[xbin+NXPIX*ybin] += 1;
//                     //m_list[xbin+NXPIX*ybin] = m;
//                                         }
//                 }
//                 else
//                 {
//                     ccd[xbin+NXPIX*ybin] += 1;
//                     //counts[xbin+NXPIX*ybin] += 1;
//                     //m_list[xbin+NXPIX*ybin] = m;
//                 }
//             }
//         }
//     }}
//         gsl_spline_free (spline);
//         gsl_interp_accel_free(acc);
//         gsl_rng_free(r);
// }
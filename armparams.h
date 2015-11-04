/* ::::::::::::: COMMON ARM PARAMETERS :::::::::::::::::::::::::::: */

const double HEIGHT_RATIO = 2.206;
const double OFFSET_A = -0.175; /* [mm] */
const double OFFSET_B =  0.175; /* [mm] */
const double DFIBER = 0.1;
const double TAU_0 = 9.01;               /* [deg] */
const double NORM_CCD_VAL = 40000.0;

/*::::::::::::::::::: NIR ARM :::::::::::::::::::::::::::::::::*/

/* COLLIMATOR */
const double NIR_FN = (10.218 / 3.5);
const double NIR_F_COL = 1565.0; /* collimator focal length [mm] */

/* ECHELLE */
const int NIR_M_LOWER = 36;  /* min ech order */
const int NIR_M_UPPER = 64;  /* max ech order */
const double NIR_ALPHA_E_DEG = 75.2;            /* [deg] */
const double NIR_GAMMA_E_DEG = 1.2;             /* [deg] */
const double NIR_SIGMA_E = (1.0 / 31.6);
//const int NIR_ORDER_LOWER = 36;
//const int NIR_ORDER_UPPER = 64;

/* GRISM */
const double NIR_ALPHA_G_DEG = 13.85;            /* [deg] */
const double NIR_NU_0 = 6.4; /* prism entrance angle [deg] 4.0+2.4 or 4.0-2.4 */
const double NIR_NU_1 = 4.0; /* grating entrance angle [deg] */
const double NIR_NU_2 = 4.0; /* grism exit angle [deg] */
const double NIR_SIGMA_G = (1.0 / 81.0);
const double NIR_N_G = 1.45780315;
const double NIR_NB1 = 1.28035628;
const double NIR_NB2 = 0.163505973;
const double NIR_NB3 = 0.893930112;
const double NIR_NC1 = 0.00929854416;
const double NIR_NC2 = 0.0449135769;
const double NIR_NC3 = 110.493685;

/* DETECTOR */
const double NIR_F_CAM = 548.276;

/* CCD */
const double NIR_DPIX = 18.0e-3; /* [mm] */
const double NIR_DGAP = 2.514; /* [mm] */
const int NIR_NXPIX = 4096;
const int NIR_NYPIX = 2048;

/* ::::::::::::::::::::: VIS ARM ::::::::::::::::::::::::::::::::::::: */

/* COLLIMATOR */
const double VIS_FN = (10.274/3.5);
const double VIS_F_COL = (0.5 * 3179.987); /* collimator focal length [mm] */

/* ECHELLE */
const int VIS_M_LOWER = 59;          /* min ech order */
const int VIS_M_UPPER = 111;         /* max ech order */
const double VIS_ALPHA_E_DEG = 75.2;            /* [deg] */
const double VIS_GAMMA_E_DEG = 1.2;             /* [deg] */
const double VIS_SIGMA_E = (1.0 / 31.6);
//const int VIS_ORDER = 59;
//const int VIS_ORDER = 111;

/* GRISM */
const double VIS_ALPHA_G_DEG = 17.8;            /* [deg] */
const double VIS_NU_0 = -8.4; /* prism entrance angle [deg] 4.0+2.4 or 4.0-2.4 */
const double VIS_NU_1 = 17.8; /* grating entrance angle [deg] */
const double VIS_NU_2 = -11.8; /* prism exit angle [deg] */
const double VIS_SIGMA_G = (1.0 / 223.0);
const double VIS_N_G = 1.58144;
const double VIS_NB1 = 1.28035628;
const double VIS_NB2 = 0.163505973;
const double VIS_NB3 = 0.893930112;
const double VIS_NC1 = 0.00929854416;
const double VIS_NC2 = 0.0449135769;
const double VIS_NC3 = 110.493685;

/* DETECTOR */
const double VIS_F_CAM = 455.0;

/* CCD */
const double VIS_DPIX = 15.0e-3;
const int VIS_NXPIX = 4096;
const int VIS_NYPIX = 4096;
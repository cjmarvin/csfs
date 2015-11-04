//based on spline_running_schnell
//dist array problem solved
//include RSS

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


/*
 *	Natural Splines
 *	---------------
 *	Here the end-conditions are determined by setting the second
 *	derivative of the spline at the end-points to equal to zero.
 *
 *	There are n-2 unknowns (y[i]'' at x[2], ..., x[n-1]) and n-2
 *	equations to determine them.  Either Choleski or Gaussian
 *	elimination could be used.
 */

/*-------------------------------------------------------------------------*/

static void natural_spline(
  int     n, 
  double *x, 
  double *y, 
  double *b, 
  double *c, 
  double *d
){
    int nm1, i;
    double t;

    x--; y--; b--; c--; d--;

    if(n < 2) {
	errno = EDOM;
	return;
    }

    if(n < 3) {
	t = (y[2] - y[1]);
	b[1] = t / (x[2] - x[1]);
	b[2] = b[1];
	c[1] = c[2] = d[1] = d[2] = 0.0;
	return;
    }

    nm1 = n - 1;

    /* Set up the tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */

    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1]) / d[1];
    for( i = 2 ; i < n ; i++) {
	d[i] 	= x[i+1] - x[i];
	b[i] 	= 2.0 * (d[i-1] + d[i]);
	c[i+1] 	= (y[i+1] - y[i]) / d[i];	
	c[i] 	= c[i+1] - c[i];
    }

    /* Gaussian elimination */

    for(i = 3; i < n; i++) {
	t = d[i-1] / b[i-1];
	b[i] = b[i] - t * d[i-1];
	c[i] = c[i] - t * c[i-1];
    }

    /* Backward substitution */

    c[nm1] = c[nm1] / b[nm1];
    for(i = n-2; i > 1; i--)
	c[i] = (c[i] - d[i] * c[i+1]) / b[i];

    /* End conditions */

    c[1] = c[n] = 0.0;

    /* Get cubic coefficients */

    b[1] = (y[2] - y[1]) / d[1] - d[i] * c[2];
    c[1] = 0.0;
    d[1] = c[2] / d[1];
    b[n] = (y[n] - y[nm1]) / d[nm1] + d[nm1] * c[nm1];
    for(i = 2; i < n; i++) {
	b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2.0 * c[i]);
	d[i] = (c[i+1] - c[i]) / d[i];
	c[i] = 3.0 * c[i];
    }
    c[n] = 0.0;
    d[n] = 0.0;

    return;
}

/*-------------------------------------------------------------------------*/

void spline_eval(
  int    *nu,	//evaluation: number of points to be evaluated
  double *u,	//evaluation: input array
  double *v,	//evaluation: output array
  int    *n,	//interporation: number of interporating points
  double *x,	//interporation: input array
  double *y,	//interporation: output array
  double *b,	//interporation coefficient
  double *c,	//interporation coefficients
  double *d 	//interporation coefficient
){
  
/* Evaluate  v[l] := spline(u[l], ...),       l = 1,..,nu, i.e. 0:(nu-1)
 * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(*n-1)
 */
    const int n_1 = *n - 1;
    int i, j, k, l;
    double ul, dx, tmp;

    for(l = 0; l < *nu; l++)
      v[l] = u[l];

    i = 0;
    for(l = 0; l < *nu; l++) {
       ul = v[l];
       if(ul < x[i] || (i < n_1 && x[i+1] < ul)) {
           /* reset i  such that  x[i] <= ul <= x[i+1] : */
           i = 0;
           j = *n;
           do {
              k = (i+j)/2;
              if(ul < x[k]) j = k;
              else i = k;
           }
           while(j > i+1);
       }
       dx = ul - x[i];
       v[l] = y[i] + dx*(b[i] + dx*(c[i] + dx*d[i]));
    }
}

/*-------------------------------------------------------------------------*/

void random_wave(double *wav_org, double *dist, int N, int flag, double *y, int NS){

      
    /*---CDF; Data extension---*/
    
    int n, N_plus = N + 1;			//number of extended data
    double *wav  = malloc(sizeof(double) * N_plus);
    double *CDF  = malloc(sizeof(double) * N_plus);
    double delta = (wav_org[1] - wav_org[0]) / 2;

    CDF[0] = 0;
    for(n = 0; n < N; n++) {
       CDF[n+1] = CDF[n] + dist[n];
       wav[n] = wav_org[n] - delta;
    }
    wav[N] = wav_org[N-1] + delta;
   
    double total = CDF[N];
    
    /* wav range */

    double lambda_min	= wav[0];
    double lambda_max	= wav[N];
    double range	= lambda_max - lambda_min;    
   
    printf("Reading data completed. Sampling started...\n");
    
    int ns;
    double *b, *c, *d;
    b = malloc(sizeof(double) * (N+1));		//to be investigated
    c = malloc(sizeof(double) * (N+1));
    d = malloc(sizeof(double) * (N+1));
    natural_spline(N_plus, CDF, wav, b, c, d);	//wav as a function of CDF
    double *x = malloc(sizeof(double) * NS);
    
    srand48((long)time(NULL));
    for (ns = 0; ns < NS; ns ++)
      x[ns] = total * drand48();
    spline_eval(&NS, x, y, &N_plus, CDF, wav, b, c, d);	//returns y to be the sampling wavelength
}

/*-------------------------------------------------------------------------*/

int main(int argc, char *argv[]){

    clock_t start = clock();
    
    int intpol_deg, NS = pow(10, 7);			//number of sampling points
    
    /*---interpolation mode---*/    
    
    FILE *file;
    if (argc == 2){
      file = fopen("./wav_int_sim.dat", "r");
      intpol_deg = (int) atoi(argv[1]);
    }
    else if (argc == 3){
      file = fopen(argv[1], "r");
      intpol_deg = (int) atoi(argv[2]);
    }
  
    /*---read in data---*/

    int n, N = 500001; 				//number of data
    double *y = malloc(sizeof(double) * NS);
    double *wav_org	= malloc(sizeof(double) * N);
    double *dist_org	= malloc(sizeof(double) * N);
    
    for(n = 0; n < N; n++)
       fscanf(file, "%lf%lf", &wav_org[n], &dist_org[n]);  

    /*---interporation---*/
    
    if (intpol_deg == 3)
      random_wave(wav_org, dist_org, N, intpol_deg, y, NS);
    
    clock_t ends = clock(); 
    printf("Running time = %7.3f s for %.1le samplings\n",(double) (ends - start) / CLOCKS_PER_SEC, (double)NS); 
    
    
    return 0;
}

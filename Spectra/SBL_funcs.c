/* (C) 2019 Jeffrey M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "SBL_funcs.h"

#include "spectra_aux.h"

/***************************************************************************
 *                        H E L P E R   F U N C T I O N S                  *
 ***************************************************************************/

/* Minimum of two numbers */
static inline double min(double x, double y)
{
  if (x < y)
    return x;
  else
    return y;
}

/* Square of real number */
static inline double square(double x)
{
  return x*x;
}

double theta_x;
double delta_m; /* Delta m_{41}^2 */

/**********************************
*  INPUT FOR PROBABILITY ENGINES  *
***********************************/

/* 
  Some of the input is obnoxiously long; it lives in "spectra_aux.h"
*/

static const int DANSS_up_long = 301;
static const int DANSS_down_long = 301;

static const int NEOS_long = 301;

static const int bugey_15_long = 401;
static const int bugey_40_long = 401;

/*********************************
*  DEFINING PROBABILITY ENGINES  *
**********************************/

static const double Dm21 = 7.39e-5; /* Recall: units are eV^2 */

static const double th12 = 0.59027;

  /**************************
  *  DANSS - UPPER POSITION *
  ***************************/

int DANSS_up_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

/* 
  We define the average probability for a given value of 
			q = 1.267 * (Dm241 [eV^2])/(E[GeV])
  via linear interpolation – we have numerically integrated over the DANSS core and
  detector (with some approximation) in Mathematica to get the probability, defined
  in the array DANSS_up_prob.
*/	

  int i, j;
  double l = length[0];
  double logq = log10(1.267*delta_m/E); /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (DANSS_up_prob[1][0]-DANSS_up_prob[0][0]);
  int pos = 0;
  double numerator, denominator;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  if (logq <= DANSS_up_prob[0][0]){
  /* For logq too small, quadratic oscillations */
    factor = DANSS_up_prob[0][1]*pow(10.0, 2.0*(logq - DANSS_up_prob[0][0])); 
  }
  else if (logq >= DANSS_up_prob[DANSS_up_long-1][0]){
  /* For logq too large, oscillations average out */
    factor = 0.5;
  }
  else{
    pos = floor( (logq-DANSS_up_prob[0][0])/step );
    numerator = (DANSS_up_prob[pos+1][1]-DANSS_up_prob[pos][1]);
    denominator = (DANSS_up_prob[pos+1][0]-DANSS_up_prob[pos][0]);;
    factor = DANSS_up_prob[pos][1] + (logq-DANSS_up_prob[pos][0])*numerator/denominator;
  }

  P[0][0] = 1.0-square(sin(2.0*theta_x))*factor ;
  return 0;
}

  /**************************
  *  DANSS - LOWER POSITION *
  ***************************/

int DANSS_down_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

/* 
  This is analogous to DANSS_up_probability_matrix, but for the lower position. See
  description above for more details.
*/	

  int i, j;
  double l = length[0];
  double logq = log10(1.267*delta_m/E); /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (DANSS_down_prob[1][0]-DANSS_down_prob[0][0]);
  int pos = 0;
  double numerator, denominator;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  if (logq <= DANSS_down_prob[0][0]){
  /* For logq too small, quadratic oscillations */
    factor = DANSS_down_prob[0][1]*pow(10.0, 2.0*(logq - DANSS_down_prob[0][0])); 
  }
  else if (logq >= DANSS_down_prob[DANSS_down_long-1][0]){
  /* For logq too large, oscillations average out */
    factor = 0.5;
  }
  else{
    pos = floor( (logq-DANSS_down_prob[0][0])/step );
    numerator = (DANSS_down_prob[pos+1][1]-DANSS_down_prob[pos][1]);
    denominator = (DANSS_down_prob[pos+1][0]-DANSS_down_prob[pos][0]);;
    factor = DANSS_down_prob[pos][1] + (logq-DANSS_down_prob[pos][0])*numerator/denominator;
  }

  P[0][0] = 1.0-square(sin(2.0*theta_x))*factor ;
  return 0;
}

  /*********
  *  NEOS  *
  **********/

int NEOS_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

/* 
  We define the average probability for a given value of 
			q = 1.267 * (Dm241 [eV^2])/(E[GeV])
  via linear interpolation – we have numerically integrated over the NEOS core and
  detector (with some approximation) in Mathematica to get the probability, defined
  in the array NEOS_prob.
*/	

  int i, j;
  double l = length[0];
  double logq = log10(1.267*delta_m/E); /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (NEOS_prob[1][0]-NEOS_prob[0][0]);
  int pos = 0;
  double numerator, denominator;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  if (logq <= NEOS_prob[0][0]){
   /* For logq too small, quadratic oscillations */
    factor = NEOS_prob[0][1]*pow(10.0, 2.0*(logq - NEOS_prob[0][0]));
  }
  else if (logq >= NEOS_prob[NEOS_long-1][0]){
  /* For logq too large, oscillations average out */
    factor = 0.5;
  }
  else{
    pos = floor( (logq-NEOS_prob[0][0])/step );
    numerator = (NEOS_prob[pos+1][1]-NEOS_prob[pos][1]);
    factor = NEOS_prob[pos][1] + (logq-NEOS_prob[pos][0])*numerator/step;
  }

  P[0][0] = 1.0-square(sin(2.0*theta_x))*factor ;
  return 0;
}

  /**********
  *  BUGEY  *
  ***********/

    /******************
    *  15 m position  *
    *******************/

int bugey_15_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

/* 
  We define the average probability for a given value of 
			q = 1.267 * (Dm241 [eV^2])/(E[GeV])
  via linear interpolation – we have numerically integrated over the core and
  detector (with some approximation) in Mathematica to get the probability.
*/	

  int i, j;
  double l = length[0];
  double logq = log10(1.267*delta_m/E); /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (bugey_15_prob[1][0]-bugey_15_prob[0][0]);
  int pos = 0;
  double numerator, denominator;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  if (logq <= bugey_15_prob[0][0]){
   /* For logq too small, quadratic oscillations */
    factor = bugey_15_prob[0][1]*pow(10.0, 2.0*(logq - bugey_15_prob[0][0]));
  }
  else if (logq >= bugey_15_prob[bugey_15_long-1][0]){
  /* For logq too large, oscillations average out */
    factor = 0.5;
  }
  else{
    pos = floor( (logq-bugey_15_prob[0][0])/step );
    numerator = (bugey_15_prob[pos+1][1]-bugey_15_prob[pos][1]);
    factor = bugey_15_prob[pos][1] + (logq-bugey_15_prob[pos][0])*numerator/step;
  }

  P[0][0] = 1.0-square(sin(2.0*theta_x))*factor ;
  return 0;
}

    /******************
    *  40 m position  *
    *******************/

int bugey_40_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

/* 
  We define the average probability for a given value of 
			q = 1.267 * (Dm241 [eV^2])/(E[GeV])
  via linear interpolation – we have numerically integrated over the core and
  detector (with some approximation) in Mathematica to get the probability.
*/	

  int i, j;
  double l = length[0];
  double logq = log10(1.267*delta_m/E); /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (bugey_40_prob[1][0]-bugey_40_prob[0][0]);
  int pos = 0;
  double numerator, denominator;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  if (logq <= bugey_40_prob[0][0]){
   /* For logq too small, quadratic oscillations */
    factor = bugey_40_prob[0][1]*pow(10.0, 2.0*(logq - bugey_40_prob[0][0]));
  }
  else if (logq >= bugey_40_prob[bugey_40_long-1][0]){
  /* For logq too large, oscillations average out */
    factor = 0.5;
  }
  else{
    pos = floor( (logq-bugey_40_prob[0][0])/step );
    numerator = (bugey_40_prob[pos+1][1]-bugey_40_prob[pos][1]);
    factor = bugey_40_prob[pos][1] + (logq-bugey_40_prob[pos][0])*numerator/step;
  }

  P[0][0] = 1.0-square(sin(2.0*theta_x))*factor ;
  return 0;
}

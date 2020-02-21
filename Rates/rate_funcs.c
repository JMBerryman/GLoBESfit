/* (C) 2019 Jeffrey M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "rate_funcs.h"
#include "rate_combo_aux.h"

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

/**********************************
*  INPUT FOR PROBABILITY ENGINES  *
***********************************/

static double theta_x;
static double delta_m; /* Delta m_{41}^2 */
static double delta_atm; /* Delta m_{31}^2 */
static double theta_13;

static const double Dm21 = 7.39e-5;
static const double Dm31 = 2.525e-3;

static const double th12 = 0.59027;

static const int DB_EH1_AD1_long = 2001;
static const int DB_EH1_AD2_long = 2001;
static const int DB_EH2_AD3_long = 2001;
static const int DB_EH2_AD8_long = 2001;

static const int RENO_ND_long = 1751;
static const int RENO_FD_long = 1751;

/**************************************************************************
*   This function is the result of                                        *
*       \int_La^Lb dL \sin^2(q*L)/L^2 / \int_La^Lb dL 1/L^2		  *
***************************************************************************/

static inline double lsin(double q, double La, double Lb)
{ 
  double result;

  result = (Lb*square(sin(La*q))-La*square(sin(Lb*q))+La*Lb*q*(gsl_sf_Si(2.0*Lb*q)-gsl_sf_Si(2.0*La*q)))/(Lb-La) ;

  return result;
}

/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/

int combo_set_oscillation_parameters(glb_params p, void *user_data)
{
  theta_13 = glbGetOscParams(p, MY_THETA_13);
  theta_x = glbGetOscParams(p, MY_THETA_X);
  delta_m =  glbGetOscParams(p, MY_DELTA_M);
  delta_atm =  glbGetOscParams(p, MY_DELTA_ATM);
  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/

int combo_get_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p,theta_13,MY_THETA_13); 
  glbSetOscParams(p,theta_x,MY_THETA_X); 
  glbSetOscParams(p,delta_m,MY_DELTA_M);
  glbSetOscParams(p,delta_atm,MY_DELTA_ATM);
  return 0;
}

/*********************************
*  DEFINING PROBABILITY ENGINES  *
**********************************/

  /****************************
  *  TWO-FLAVOR PROBABILITIES *
  *****************************/

int combo_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{

/* 
  "user_data" contains the width of each experiment's detector
*/	

  int i, j;
  double La,Lb;
  double l=length[0];
  double wide = ((double *)user_data)[0];
  
  La = l-wide;
  Lb = l+wide;

  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  P[0][0] = 1.0-square(sin(2.0*theta_x))*lsin(1.267*delta_m/E,La,Lb) ;
  return 0;
}

  /**************************
  *  STANDARD PROBABILITIES *
  ***************************/

int standard_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  int i,j;
  double La,Lb;
  double l=length[0];
  double wide = ((double *)user_data)[0];
  
  La = l-wide;
  Lb = l+wide;

  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  double Dm41 = delta_m;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21*l/E;
  double DD31 = 1.267*Dm31*l/E;
  double DD41 = 1.267*Dm41*l/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  P[0][0] = 1.0;

  for (i=0; i<3; i++){
    for(j=i+1; j<4; j++){
      P[0][0] += -4.0*square(u[i]*u[j])*lsin(lambda[j]-lambda[i],La,Lb);
    }
  }
  return 0;
}

  /**************************
  *  DAYA BAY PROBABILITIES *
  ***************************/

    /************
    *  EH1 AD1  *
    *************/

int DB_EH1_AD1_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  double numerator = 0.0;

  double Dm41 = delta_m;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21/E; /* Recall: E in GeV! */
  double DD31 = 1.267*Dm31/E;
  double DD41 = 1.267*Dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (DB_EH1_AD1_prob[1][0]-DB_EH1_AD1_prob[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  P[0][0] = 1.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = -1.0*log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= DB_EH1_AD1_prob[0][0]){
 	/* For logq too small, oscillations average out */
        factor = 0.5;
      }
      else if (logq >= DB_EH1_AD1_prob[DB_EH1_AD1_long-1][0]){
	/* For logq too large, quadratic oscillations */
        factor = DB_EH1_AD1_prob[DB_EH1_AD1_long-1][1]*pow(10.0, -2.0*(logq-DB_EH1_AD1_prob[DB_EH1_AD1_long-1][0]));
      }
      else{
        pos = floor( (logq-DB_EH1_AD1_prob[0][0])/step );
        numerator = (DB_EH1_AD1_prob[pos+1][1]-DB_EH1_AD1_prob[pos][1]);
        factor = DB_EH1_AD1_prob[pos][1] + (logq-DB_EH1_AD1_prob[pos][0])*numerator/step;
      }
      P[0][0] += -4.0*u[k]*u[k]*u[j]*u[j] * factor;
    }
  }

  return 0;
}  
    /************
    *  EH1 AD2  *
    *************/

int DB_EH1_AD2_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  double inner = 0.0;
  double numerator = 0.0;
  double denominator = 0.0;
  double Dm41 = delta_m;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21/E;
  double DD31 = 1.267*Dm31/E;
  double DD41 = 1.267*Dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (DB_EH1_AD2_prob[1][0]-DB_EH1_AD2_prob[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = -1.0*log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= DB_EH1_AD2_prob[0][0]){
 	/* For logq too small, oscillations average out */
        factor = 0.5;
      }
      else if (logq >= DB_EH1_AD2_prob[DB_EH1_AD2_long-1][0]){
	/* For logq too large, quadratic oscillations */
        factor = DB_EH1_AD2_prob[DB_EH1_AD2_long-1][1]*pow(10.0, -2.0*(logq-DB_EH1_AD2_prob[DB_EH1_AD2_long-1][0]));
      }
      else{
        pos = floor( (logq-DB_EH1_AD2_prob[0][0])/step );
        numerator = (DB_EH1_AD2_prob[pos+1][1]-DB_EH1_AD2_prob[pos][1]);
        denominator = (DB_EH1_AD2_prob[pos+1][0]-DB_EH1_AD2_prob[pos][0]);
        factor = DB_EH1_AD2_prob[pos][1] + (logq-DB_EH1_AD2_prob[pos][0])*numerator/denominator;
      }
      inner += 4.0*square(u[k]*u[j]) * factor;
    }
  }

  P[0][0] = 1.0 - inner ;
  return 0;
}  

    /************
    *  EH2 AD3  *
    *************/
   
int DB_EH2_AD3_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  double inner = 0.0;
  double numerator = 0.0;
  double denominator = 0.0;

  double Dm41 = delta_m;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21/E;
  double DD31 = 1.267*Dm31/E;
  double DD41 = 1.267*Dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (DB_EH2_AD3_prob[1][0]-DB_EH2_AD3_prob[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = -1.0*log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= DB_EH2_AD3_prob[0][0]){
 	/* For logq too small, oscillations average out */
        factor = 0.5;
      }
      else if (logq >= DB_EH2_AD3_prob[DB_EH2_AD3_long-1][0]){
	/* For logq too large, quadratic oscillations */
        factor = DB_EH2_AD3_prob[DB_EH2_AD3_long-1][1]*pow(10.0, -2.0*(logq-DB_EH2_AD3_prob[DB_EH2_AD3_long-1][0]));
      }
      else{
        pos = floor( (logq-DB_EH2_AD3_prob[0][0])/step );
        numerator = (DB_EH2_AD3_prob[pos+1][1]-DB_EH2_AD3_prob[pos][1]);
        denominator = (DB_EH2_AD3_prob[pos+1][0]-DB_EH2_AD3_prob[pos][0]);
        factor = DB_EH2_AD3_prob[pos][1] + (logq-DB_EH2_AD3_prob[pos][0])*numerator/denominator;
      }
      inner += 4.0*square(u[k]*u[j]) * factor;
    }
  }

  P[0][0] = 1.0 - inner ;
  return 0;
}  

    /************
    *  EH2 AD8  *
    *************/

int DB_EH2_AD8_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  double inner = 0.0;
  double numerator = 0.0;
  double denominator = 0.0;

  double Dm41 = delta_m;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21/E;
  double DD31 = 1.267*Dm31/E;
  double DD41 = 1.267*Dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (DB_EH2_AD8_prob[1][0]-DB_EH2_AD8_prob[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = -1.0*log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= DB_EH2_AD8_prob[0][0]){
 	/* For logq too small, oscillations average out */
        factor = 0.5;
      }
      else if (logq >= DB_EH2_AD8_prob[DB_EH2_AD8_long-1][0]){
	/* For logq too large, quadratic oscillations */
        factor = DB_EH2_AD8_prob[DB_EH2_AD8_long-1][1]*pow(10.0, -2.0*(logq-DB_EH2_AD8_prob[DB_EH2_AD8_long-1][0]));
      }
      else{
        pos = floor( (logq-DB_EH2_AD8_prob[0][0])/step );
        numerator = (DB_EH2_AD8_prob[pos+1][1]-DB_EH2_AD8_prob[pos][1]);
        denominator = (DB_EH2_AD8_prob[pos+1][0]-DB_EH2_AD8_prob[pos][0]);
        factor = DB_EH2_AD8_prob[pos][1] + (logq-DB_EH2_AD8_prob[pos][0])*numerator/denominator;
      }
      inner += 4.0*square(u[k]*u[j]) * factor;
    }
  }

  P[0][0] = 1.0 - inner ;
  return 0;
}

  /*********
  *  RENO  *
  **********/

    /******************
    *  NEAR DETECTOR  *
    *******************/
   
int RENO_ND_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  double numerator = 0.0;

  double Dm41 = delta_m;
  double Dm31 = delta_atm;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21/E;
  double DD31 = 1.267*Dm31/E;
  double DD41 = 1.267*Dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (RENO_ND_prob[1][0]-RENO_ND_prob[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  P[0][0] = 1.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = -1.0*log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= RENO_ND_prob[0][0]){
 	/* For logq too small, oscillations average out */
        factor = 0.5;
      }
      else if (logq >= RENO_ND_prob[RENO_ND_long-1][0]){
	/* For logq too large, quadratic oscillations */
        factor = RENO_ND_prob[RENO_ND_long-1][1]*pow(10.0, -2.0*(logq - RENO_ND_prob[RENO_ND_long-1][0]));
      }
      else{
        pos = floor( (logq-RENO_ND_prob[0][0])/step );
        numerator = (RENO_ND_prob[pos+1][1]-RENO_ND_prob[pos][1]);
        factor = RENO_ND_prob[pos][1] + (logq-RENO_ND_prob[pos][0])*numerator/step;
      }
      P[0][0] += - 4.0*square(u[k]*u[j]) * factor;
    }
  }
  return 0;
}  

    /*****************
    *  FAR DETECTOR  *
    ******************/

int RENO_FD_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  double numerator = 0.0;

  double Dm41 = delta_m;
  double Dm31 = delta_atm;
  double th14 = theta_x;
  double th13 = theta_13;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = 1.267*Dm21/E;
  double DD31 = 1.267*Dm31/E;
  double DD41 = 1.267*Dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (RENO_FD_prob[1][0]-RENO_FD_prob[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  P[0][0] = 1.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = -1.0*log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= RENO_FD_prob[0][0]){
 	/* For logq too small, oscillations average out */
        factor = 0.5;
      }
      else if (logq >= RENO_FD_prob[RENO_FD_long-1][0]){
	/* For logq too large, quadratic oscillations */
        factor = RENO_FD_prob[RENO_FD_long-1][1]*pow(10.0, -2.0*(logq - RENO_FD_prob[RENO_FD_long-1][0]));
      }

      else{
        pos = floor( (logq - RENO_FD_prob[0][0])/step );
        numerator = (RENO_FD_prob[pos+1][1] - RENO_FD_prob[pos][1]);
        factor = RENO_FD_prob[pos][1] + (logq - RENO_FD_prob[pos][0])*numerator/step;
      }

      P[0][0] += -4.0*u[k]*u[k]*u[j]*u[j] * factor;
    }
  }

  return 0;
}

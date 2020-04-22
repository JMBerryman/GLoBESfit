/* GLoBESfit -- GLoBES fitting tools
*  (C) 2019-2020 The GLoBESfit Team
*
* GLoBESfit is mainly intended for academic purposes. Proper
* credit must be given if you use GLoBESfit or parts of it. Please
* read the section 'Credit' in the README file.
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* (C) 2019, 2020 Patrick Huber, J. M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "glf_types.h"
#include "glf_probability.h"

#define GLF_LE_KM_GEV 1.0/(GLB_EV_TO_KM_FACTOR*4.0*1E9)


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

static double theta_14;
static double delta_m; /* Delta m_{41}^2 */
static double delta_atm; /* Delta m_{31}^2 */
static double theta_13;
static double theta_12;
static double delta_solar; /* Delta m_{41}^2 */


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

int glf_get_oscillation_parameters(glb_params p, void *user_data)
{
  theta_13 = glbGetOscParams(p, GLF_THETA_13);
  theta_14 = glbGetOscParams(p, GLF_THETA_14);
  delta_m =  glbGetOscParams(p, GLF_DELTA_M);
  delta_atm =  glbGetOscParams(p, GLF_DELTA_ATM);
  theta_12 =  glbGetOscParams(p, GLF_THETA_12);
  delta_solar =  glbGetOscParams(p, GLF_DELTA_SOLAR);
  return 0;
}

/***************************************************************************
 * Write oscillation parameters from internal data structures into p.      *
 ***************************************************************************/

int glf_set_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p,theta_13,GLF_THETA_13); 
  glbSetOscParams(p,theta_14,GLF_THETA_14); 
  glbSetOscParams(p,delta_m,GLF_DELTA_M);
  glbSetOscParams(p,delta_atm,GLF_DELTA_ATM);
  glbSetOscParams(p,theta_12,GLF_THETA_12);
  glbSetOscParams(p,delta_solar,GLF_DELTA_SOLAR);
  return 0;
}

/*********************************
*  DEFINING PROBABILITY ENGINES  *
**********************************/

  /****************************
  *  TWO-FLAVOR PROBABILITIES *
  *****************************/

int glf_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
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

  P[0][0] = 1.0-square(sin(2.0*theta_14))*lsin(GLF_LE_KM_GEV*delta_m/E,La,Lb) ;
  return 0;
}

  /**************************
  *  STANDARD PROBABILITIES *
  ***************************/

int glf_standard_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
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

  double th12 = theta_12;
  double th13 = theta_13;
  double th14 = theta_14;

  double dm21 = delta_solar;
  double dm31 = delta_atm;
  double dm41 = delta_m;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = GLF_LE_KM_GEV*dm21*l/E;
  double DD31 = GLF_LE_KM_GEV*dm31*l/E;
  double DD41 = GLF_LE_KM_GEV*dm41*l/E;

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

/*******************************************************************
 * Unified treatment of oscillation probabilities                  *
 *******************************************************************/

int glf_four_state_probability_matrix(
double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  glf_distance_data *distances = (glf_distance_data *) user_data;

  size_t dlength= distances->length;
  double (*data)[2]= distances->data;

  double numerator = 0.0;

  double th12 = theta_12;
  double th13 = theta_13;
  double th14 = theta_14;

  double dm21 = delta_solar;
  double dm31 = delta_atm;
  double dm41 = delta_m;

  double c12 = cos(th12);
  double s12 = sin(th12);

  double c13=cos(th13);
  double s13=sin(th13);

  double c14=cos(th14);
  double s14=sin(th14);

  double DD21 = GLF_LE_KM_GEV*dm21/E; /* Recall: E in GeV! */
  double DD31 = GLF_LE_KM_GEV*dm31/E;
  double DD41 = GLF_LE_KM_GEV*dm41/E;

  double lambda[4] = {0.0, DD21, DD31, DD41};

  double u[4]; /* The elements of the "e" row of the PMNS matrix */
  u[0] = c12*c13*c14;
  u[1] = c13*s12*c14;
  u[2] = s13*c14;
  u[3] = s14;

  int i, j, k;
  double logq; /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (data[1][0]-data[0][0]);
  int pos = 0;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  P[0][0] = 1.0;

  for (j=0; j<3; j++){
    for (k=j+1; k<4; k++){
      logq = log10(fabs(lambda[k]-lambda[j]));
  
      if (logq <= data[0][0]){
	/* For logq too small, quadratic oscillations */
        factor = data[0][1]*pow(10.0, 2.0*(logq-data[0][0]));
      }
      else if (logq >= data[dlength-1][0]){
 	/* For logq too large, oscillations average out */
        factor = 0.5;
      }
      else{
        pos = floor( (logq-data[0][0])/step );
        numerator = (data[pos+1][1]-data[pos][1]);
        factor = data[pos][1] + (logq-data[pos][0])*numerator/step;
      }
      P[0][0] += -4.0*u[k]*u[k]*u[j]*u[j] * factor;
    }
  }

  
  return 0;
}


int glf_two_state_probability_matrix(
double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
                          double filter_sigma, void *user_data)
{
  glf_distance_data *distances = (glf_distance_data *) user_data;

  size_t dlength= distances->length;
  double (*data)[2]= distances->data;

  int i, j;
  double l = length[0];
  double logq = log10(1.267*delta_m/E); /* NB: working with log10 of q! */
  double factor = 0.0;

  double step = (data[1][0]-data[0][0]);
  int pos = 0;
  double numerator, denominator;
  
  /* Set all probabilities to zero initially */
  for (i=0; i < 3; i++)
    for (j=0; j < 3; j++)
      P[i][j] = 0.0;

  if (logq <= data[0][0]){
   /* For logq too small, quadratic oscillations */
    factor = data[0][1]*pow(10.0, 2.0*(logq - data[0][0]));
  }
  else if (logq >= data[dlength-1][0]){
  /* For logq too large, oscillations average out */
    factor = 0.5;
  }
  else{
    pos = floor( (logq-data[0][0])/step );
    numerator = (data[pos+1][1]-data[pos][1]);
    factor = data[pos][1] + (logq-data[pos][0])*numerator/step;
  }

  P[0][0] = 1.0-square(sin(2.0*theta_14))*factor ;

  
  return 0;
}


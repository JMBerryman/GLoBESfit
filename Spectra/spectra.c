/* (C) 2019 Jeffrey M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "spectra.h"
#include "spectra_data.h"
#include "spectra_aux.h"
#include "spectra_aux2.h"

extern double systematic[10];

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
double delta_atm; /* Delta m_{31}^2 */
double theta_13;

/***************************************************************************
 * Store oscillation parameters in internal data structures.               *
 * For more sophisticated probability engines, this would be the right     *
 * place to pre-compute the mixing matrix and parts of the Hamiltonian in  *
 * order to speed up the calls to the actual probability matrix function.  *
 ***************************************************************************/

int spectra_get_oscillation_parameters(glb_params p, void *user_data)
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

int spectra_set_oscillation_parameters(glb_params p, void *user_data)
{
  glbSetOscParams(p, theta_13, MY_THETA_13); 
  glbSetOscParams(p, theta_x, MY_THETA_X); 
  glbSetOscParams(p, delta_m, MY_DELTA_M);
  glbSetOscParams(p, delta_atm, MY_DELTA_ATM);
  return 0;
}

  /**************************
  *  OTHER RELEVANT INPUTS  *
  ***************************/

/*
  The order of the experiments in the corresponding .glb file is:

  0 - DANSS (upper position)
  1 – DANSS (lower position)

  2 - Daya Bay EH1 AD1
  3 - Daya Bay EH1 AD2
  4 - Daya Bay EH2 AD3
  5 - Daya Bay EH2 AD8

  6 - Daya Bay EH3 AD4
  7 - Daya Bay EH3 AD5
  8 - Daya Bay EH3 AD6
  9 - Daya Bay EH3 AD7

  10 - NEOS

  11 - Double Chooz ND
  12 - Double Chooz FD

  13 - Bugey, 15 m
  14 - Bugey, 40 m

  15 - RENO ND
  16 - RENO FD
*/

/**************************************************
*                                                 *
*        THE CALCULATION OF THE CHI-SQUARED       *
*                                                 *
***************************************************/

double spectra_chi(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following

  There will be a "delta" for each experiment; assume no correlated uncertainties

  "user_data" contains flags indicated which experiments are to be included in the fit

*/

  double deltaDANSS[24] = {0.0};
  double deltaDB[52] = {0.0};
  double deltaNEOS[60] = {0.0};

  double deltaDC[26] = {0.0};
  double deltaBugey[25] = {0.0};
  double deltaRENO[25] = {0.0};

  int i, j, k;
  double numerator, numerator2, denominator;
  int YesNo[10];
  for (k=0; k<10; k++){
    YesNo[k] = ((int *)user_data)[k];
  }

  int Nnuisance=4; /* The number of nuisance parameters */

  int CheckZeros;
  /* A quantity that we use to look for zero entries after shifting energy scales */

  /* Various important quantities used below */
  double emin, emax; 
  /* Get emin/emax for given experiment; used for energy-scale variation */

  /**********
  *  DANSS  *
  ***********/

  /* Get emin and emax for DANSS */
  glbGetEminEmax(exp, &emin, &emax);

/*
  Note that we've included two bins on either side of the analysis window
  to help with the energy scale calculation; these are not included in our chi-squared
  analysis.
*/

  double DANSS_up_U235[32];
  double DANSS_up_U238[32];
  double DANSS_up_Pu239[32];
  double DANSS_up_Pu241[32];

  double DANSS_down_U235[32];
  double DANSS_down_U238[32];
  double DANSS_down_Pu239[32];
  double DANSS_down_Pu241[32];

  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp, 0),
                      DANSS_up_U235, 32, emin, emax);
  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp, 1),
                      DANSS_up_U238, 32, emin, emax);
  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp, 2),
                      DANSS_up_Pu239, 32, emin, emax);
  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp, 3),
                      DANSS_up_Pu241, 32, emin, emax);

  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp+1, 0),
                      DANSS_down_U235, 32, emin, emax);
  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp+1, 1),
                      DANSS_down_U238, 32, emin, emax);
  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp+1, 2),
                      DANSS_down_Pu239, 32, emin, emax);
  glbShiftEnergyScale(0.01*x[0], glbGetSignalFitRatePtr(exp+1, 3),
                      DANSS_down_Pu241, 32, emin, emax);

CheckZeros = 0;

if (YesNo[0] == 1){
  for (i=4; i<28; i++){

/*
  If these elements are zero, then something has gone with the minimization;
  take "deltaEXP" to be 1000. to force minimizer to turn around!
*/

    if(DANSS_up_U235[i]==0.0||DANSS_down_U235[i]==0.0){
      CheckZeros++;
    }
    if(DANSS_up_U238[i]==0.0||DANSS_down_U238[i]==0.0){
      CheckZeros++;
    }
    if(DANSS_up_Pu239[i]==0.0||DANSS_down_Pu239[i]==0.0){
      CheckZeros++;
    }
    if(DANSS_up_Pu241[i]==0.0||DANSS_down_Pu241[i]==0.0){
      CheckZeros++;
    }

    if (CheckZeros==0){
      numerator = DANSS_down_U235[i]*DANSSF[0];
      numerator += DANSS_down_U238[i]*DANSSF[1];
      numerator += DANSS_down_Pu239[i]*DANSSF[2];
      numerator += DANSS_down_Pu241[i]*DANSSF[3];

      denominator = DANSS_up_U235[i]*DANSSF[0];
      denominator += DANSS_up_U238[i]*DANSSF[1];
      denominator += DANSS_up_Pu239[i]*DANSSF[2];
      denominator += DANSS_up_Pu241[i]*DANSSF[3];

      deltaDANSS[i-4] = DANSSratio[i-4] - numerator/denominator;
    }
    else{
      deltaDANSS[i-4] = 1000.;
    }
  }
}
  /*************
  *  DAYA BAY  *
  **************/

  /*
    We introduce a couple of fudge factors for Daya Bay in order to reproduce their 
    determinations of theta_13 and Dm31(Dmee)
  */
  double DB12_fudge = 0.9933;
  double DB13_fudge = 0.9973;

/*
  WATCH CAREFULLY HERE! The four near detectors at Daya Bay are also used to determine
  the normalization for NEOS, so they are formatted differently from the four far
  detectors! Namely, while the far detectors are binned in 0.2 MeV bins from 2.08-8.68
  MeV, the near detectors are binned in 0.1 MeV bins from 1.78-8.68 MeV bins. The finer
  resolution is driven by the binning used at NEOS; just combine the appropriate bins to
  get the numbers of events useful at Daya Bay!

  I've also included one additional bin at low energies for the far detectors, which I'm
  not including in the calculation, because it gives a different result for the lowest bin
  that I am including. I'm not sure why this is the case...

  I also need the three-neutrino numbers for the near detectors; these are to be used with
  the NEOS analysis.
*/

    /************
    *  EH1 AD1  *
    *************/

  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+2, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+2, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+2, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+2, 3);
   
    /************
    *  EH1 AD2  *
    *************/

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+3, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+3, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+3, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+3, 3);

    /************
    *  EH2 AD3  *
    *************/

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
    /************
    *  EH2 AD8  *
    *************/

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);

    /************
    *  EH3 AD4  *
    *************/

  double *EH3_AD4_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *EH3_AD4_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *EH3_AD4_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *EH3_AD4_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);
   
    /************
    *  EH3 AD5  *
    *************/

  double *EH3_AD5_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *EH3_AD5_U238 = glbGetSignalFitRatePtr(exp+7, 1);
  double *EH3_AD5_Pu239 = glbGetSignalFitRatePtr(exp+7, 2);
  double *EH3_AD5_Pu241 = glbGetSignalFitRatePtr(exp+7, 3);

    /************
    *  EH3 AD6  *
    *************/

  double *EH3_AD6_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *EH3_AD6_U238 = glbGetSignalFitRatePtr(exp+8, 1);
  double *EH3_AD6_Pu239 = glbGetSignalFitRatePtr(exp+8, 2);
  double *EH3_AD6_Pu241 = glbGetSignalFitRatePtr(exp+8, 3);
   
    /************
    *  EH3 AD7  *
    *************/

  double *EH3_AD7_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *EH3_AD7_U238 = glbGetSignalFitRatePtr(exp+9, 1);
  double *EH3_AD7_Pu239 = glbGetSignalFitRatePtr(exp+9, 2);
  double *EH3_AD7_Pu241 = glbGetSignalFitRatePtr(exp+9, 3);

  int index;

/*
  REMINDER: We need to take the appropriate combinations of the Daya Bay bins
  because of the fine binning with which we have calculated their spectra!

  As a reminder, the Daya Bay spectra are broken into 227 bins (!!), each of width
  0.05 MeV. We restructure them to match the binning at Daya Bay, as follows:

    Bins 0-11: Daya Bay Bin 0 (0.7 - 1.3 MeV prompt)
    Bins 12 + 5*(i-1) + (0-4): Daya Bay Bin i ( = 1-24 )
                                (1.3 + 0.25*(i-1) + (0-0.25) MeV prompt)
    Bins 132-225: Daya Bay Bin 25 (7.3 - 12.0 MeV prompt)

  Below, we calculate the total number of events at each AD -- we initialize these here
*/

  double TotalAD1[26] = {0.0};
  double TotalAD2[26] = {0.0};
  double TotalAD3[26] = {0.0};
  double TotalAD8[26] = {0.0};

  double TotalAD4[26] = {0.0};
  double TotalAD5[26] = {0.0};
  double TotalAD6[26] = {0.0};
  double TotalAD7[26] = {0.0};

if (YesNo[1] == 1){
  for (i=0; i<12; i++){
    TotalAD1[0] += EH1_AD1_U235[i]*DayaBayF1[0];
    TotalAD1[0] += EH1_AD1_U238[i]*DayaBayF1[1];
    TotalAD1[0] += EH1_AD1_Pu239[i]*DayaBayF1[2];
    TotalAD1[0] += EH1_AD1_Pu241[i]*DayaBayF1[3];

    TotalAD2[0] += EH1_AD2_U235[i]*DayaBayF2[0];
    TotalAD2[0] += EH1_AD2_U238[i]*DayaBayF2[1];
    TotalAD2[0] += EH1_AD2_Pu239[i]*DayaBayF2[2];
    TotalAD2[0] += EH1_AD2_Pu241[i]*DayaBayF2[3];

    TotalAD3[0] += EH2_AD3_U235[i]*DayaBayF3[0];
    TotalAD3[0] += EH2_AD3_U238[i]*DayaBayF3[1];
    TotalAD3[0] += EH2_AD3_Pu239[i]*DayaBayF3[2];
    TotalAD3[0] += EH2_AD3_Pu241[i]*DayaBayF3[3];

    TotalAD8[0] += EH2_AD8_U235[i]*DayaBayF8[0];
    TotalAD8[0] += EH2_AD8_U238[i]*DayaBayF8[1];
    TotalAD8[0] += EH2_AD8_Pu239[i]*DayaBayF8[2];
    TotalAD8[0] += EH2_AD8_Pu241[i]*DayaBayF8[3];

    TotalAD4[0] += EH3_AD4_U235[i]*DayaBayF4[0];
    TotalAD4[0] += EH3_AD4_U238[i]*DayaBayF4[1];
    TotalAD4[0] += EH3_AD4_Pu239[i]*DayaBayF4[2];
    TotalAD4[0] += EH3_AD4_Pu241[i]*DayaBayF4[3];

    TotalAD5[0] += EH3_AD5_U235[i]*DayaBayF5[0];
    TotalAD5[0] += EH3_AD5_U238[i]*DayaBayF5[1];
    TotalAD5[0] += EH3_AD5_Pu239[i]*DayaBayF5[2];
    TotalAD5[0] += EH3_AD5_Pu241[i]*DayaBayF5[3];

    TotalAD6[0] += EH3_AD6_U235[i]*DayaBayF6[0];
    TotalAD6[0] += EH3_AD6_U238[i]*DayaBayF6[1];
    TotalAD6[0] += EH3_AD6_Pu239[i]*DayaBayF6[2];
    TotalAD6[0] += EH3_AD6_Pu241[i]*DayaBayF6[3];

    TotalAD7[0] += EH3_AD7_U235[i]*DayaBayF7[0];
    TotalAD7[0] += EH3_AD7_U238[i]*DayaBayF7[1];
    TotalAD7[0] += EH3_AD7_Pu239[i]*DayaBayF7[2];
    TotalAD7[0] += EH3_AD7_Pu241[i]*DayaBayF7[3];
  }

  for (j=1; j<25; j++){
    for (i=0; i<5; i++){
      TotalAD1[j] += EH1_AD1_U235[12+5*(j-1)+i]*DayaBayF1[0];
      TotalAD1[j] += EH1_AD1_U238[12+5*(j-1)+i]*DayaBayF1[1];
      TotalAD1[j] += EH1_AD1_Pu239[12+5*(j-1)+i]*DayaBayF1[2];
      TotalAD1[j] += EH1_AD1_Pu241[12+5*(j-1)+i]*DayaBayF1[3];

      TotalAD2[j] += EH1_AD2_U235[12+5*(j-1)+i]*DayaBayF2[0];
      TotalAD2[j] += EH1_AD2_U238[12+5*(j-1)+i]*DayaBayF2[1];
      TotalAD2[j] += EH1_AD2_Pu239[12+5*(j-1)+i]*DayaBayF2[2];
      TotalAD2[j] += EH1_AD2_Pu241[12+5*(j-1)+i]*DayaBayF2[3];

      TotalAD3[j] += EH2_AD3_U235[12+5*(j-1)+i]*DayaBayF3[0];
      TotalAD3[j] += EH2_AD3_U238[12+5*(j-1)+i]*DayaBayF3[1];
      TotalAD3[j] += EH2_AD3_Pu239[12+5*(j-1)+i]*DayaBayF3[2];
      TotalAD3[j] += EH2_AD3_Pu241[12+5*(j-1)+i]*DayaBayF3[3];

      TotalAD8[j] += EH2_AD8_U235[12+5*(j-1)+i]*DayaBayF8[0];
      TotalAD8[j] += EH2_AD8_U238[12+5*(j-1)+i]*DayaBayF8[1];
      TotalAD8[j] += EH2_AD8_Pu239[12+5*(j-1)+i]*DayaBayF8[2];
      TotalAD8[j] += EH2_AD8_Pu241[12+5*(j-1)+i]*DayaBayF8[3];

      TotalAD4[j] += EH3_AD4_U235[12+5*(j-1)+i]*DayaBayF4[0];
      TotalAD4[j] += EH3_AD4_U238[12+5*(j-1)+i]*DayaBayF4[1];
      TotalAD4[j] += EH3_AD4_Pu239[12+5*(j-1)+i]*DayaBayF4[2];
      TotalAD4[j] += EH3_AD4_Pu241[12+5*(j-1)+i]*DayaBayF4[3];

      TotalAD5[j] += EH3_AD5_U235[12+5*(j-1)+i]*DayaBayF5[0];
      TotalAD5[j] += EH3_AD5_U238[12+5*(j-1)+i]*DayaBayF5[1];
      TotalAD5[j] += EH3_AD5_Pu239[12+5*(j-1)+i]*DayaBayF5[2];
      TotalAD5[j] += EH3_AD5_Pu241[12+5*(j-1)+i]*DayaBayF5[3];

      TotalAD6[j] += EH3_AD6_U235[12+5*(j-1)+i]*DayaBayF6[0];
      TotalAD6[j] += EH3_AD6_U238[12+5*(j-1)+i]*DayaBayF6[1];
      TotalAD6[j] += EH3_AD6_Pu239[12+5*(j-1)+i]*DayaBayF6[2];
      TotalAD6[j] += EH3_AD6_Pu241[12+5*(j-1)+i]*DayaBayF6[3];

      TotalAD7[j] += EH3_AD7_U235[12+5*(j-1)+i]*DayaBayF7[0];
      TotalAD7[j] += EH3_AD7_U238[12+5*(j-1)+i]*DayaBayF7[1];
      TotalAD7[j] += EH3_AD7_Pu239[12+5*(j-1)+i]*DayaBayF7[2];
      TotalAD7[j] += EH3_AD7_Pu241[12+5*(j-1)+i]*DayaBayF7[3];
    }
  }

  for (i=132; i<226; i++){
    TotalAD1[25] += EH1_AD1_U235[i]*DayaBayF1[0];
    TotalAD1[25] += EH1_AD1_U238[i]*DayaBayF1[1];
    TotalAD1[25] += EH1_AD1_Pu239[i]*DayaBayF1[2];
    TotalAD1[25] += EH1_AD1_Pu241[i]*DayaBayF1[3];

    TotalAD2[25] += EH1_AD2_U235[i]*DayaBayF2[0];
    TotalAD2[25] += EH1_AD2_U238[i]*DayaBayF2[1];
    TotalAD2[25] += EH1_AD2_Pu239[i]*DayaBayF2[2];
    TotalAD2[25] += EH1_AD2_Pu241[i]*DayaBayF2[3];

    TotalAD3[25] += EH2_AD3_U235[i]*DayaBayF3[0];
    TotalAD3[25] += EH2_AD3_U238[i]*DayaBayF3[1];
    TotalAD3[25] += EH2_AD3_Pu239[i]*DayaBayF3[2];
    TotalAD3[25] += EH2_AD3_Pu241[i]*DayaBayF3[3];

    TotalAD8[25] += EH2_AD8_U235[i]*DayaBayF8[0];
    TotalAD8[25] += EH2_AD8_U238[i]*DayaBayF8[1];
    TotalAD8[25] += EH2_AD8_Pu239[i]*DayaBayF8[2];
    TotalAD8[25] += EH2_AD8_Pu241[i]*DayaBayF8[3];

    TotalAD4[25] += EH3_AD4_U235[i]*DayaBayF4[0];
    TotalAD4[25] += EH3_AD4_U238[i]*DayaBayF4[1];
    TotalAD4[25] += EH3_AD4_Pu239[i]*DayaBayF4[2];
    TotalAD4[25] += EH3_AD4_Pu241[i]*DayaBayF4[3];

    TotalAD5[25] += EH3_AD5_U235[i]*DayaBayF5[0];
    TotalAD5[25] += EH3_AD5_U238[i]*DayaBayF5[1];
    TotalAD5[25] += EH3_AD5_Pu239[i]*DayaBayF5[2];
    TotalAD5[25] += EH3_AD5_Pu241[i]*DayaBayF5[3];

    TotalAD6[25] += EH3_AD6_U235[i]*DayaBayF6[0];
    TotalAD6[25] += EH3_AD6_U238[i]*DayaBayF6[1];
    TotalAD6[25] += EH3_AD6_Pu239[i]*DayaBayF6[2];
    TotalAD6[25] += EH3_AD6_Pu241[i]*DayaBayF6[3];

    TotalAD7[25] += EH3_AD7_U235[i]*DayaBayF7[0];
    TotalAD7[25] += EH3_AD7_U238[i]*DayaBayF7[1];
    TotalAD7[25] += EH3_AD7_Pu239[i]*DayaBayF7[2];
    TotalAD7[25] += EH3_AD7_Pu241[i]*DayaBayF7[3];
  }

/*
  Now return to business as usual: the ratio between EH2/EH1 and EH3/EH1
*/

  for (i=0; i<26; i++){
    denominator = TotalAD1[i] + TotalAD2[i];
    numerator = TotalAD3[i] + TotalAD8[i];
    numerator2 = TotalAD4[i] + TotalAD5[i] + TotalAD6[i] + TotalAD7[i];

    deltaDB[i] = DayaBayRatio21[i] - DB12_fudge*numerator/denominator ;
    deltaDB[i+26] = DayaBayRatio31[i] - DB13_fudge*numerator2/denominator ;
  }
}
  /*********
  *  NEOS  *
  **********/

  /* Get emin and emax for NEOS */
  glbGetEminEmax(exp+10, &emin, &emax);

  double NEOS_U235[70];
  double NEOS_U238[70];
  double NEOS_Pu239[70];
  double NEOS_Pu241[70];

  double NEOS_U235_0[70];
  double NEOS_U238_0[70];
  double NEOS_Pu239_0[70];
  double NEOS_Pu241_0[70];

  glbShiftEnergyScale(0.01*x[1], glbGetSignalFitRatePtr(exp+10, 0),
                      NEOS_U235, 70, emin, emax);
  glbShiftEnergyScale(0.01*x[1], glbGetSignalFitRatePtr(exp+10, 1),
                      NEOS_U238, 70, emin, emax);
  glbShiftEnergyScale(0.01*x[1], glbGetSignalFitRatePtr(exp+10, 2),
                      NEOS_Pu239, 70, emin, emax);
  glbShiftEnergyScale(0.01*x[1], glbGetSignalFitRatePtr(exp+10, 3),
                      NEOS_Pu241, 70, emin, emax);

  glbShiftEnergyScale(0.01*x[1], glbGetRuleRatePtr(exp+10, 0),
                      NEOS_U235_0, 70, emin, emax);
  glbShiftEnergyScale(0.01*x[1], glbGetRuleRatePtr(exp+10, 1),
                      NEOS_U238_0, 70, emin, emax);
  glbShiftEnergyScale(0.01*x[1], glbGetRuleRatePtr(exp+10, 2),
                      NEOS_Pu239_0, 70, emin, emax);
  glbShiftEnergyScale(0.01*x[1], glbGetRuleRatePtr(exp+10, 3),
                      NEOS_Pu241_0, 70, emin, emax);

CheckZeros = 0;

  double NEOSnum = 0.0;
  double NEOSdenom = 0.0;
  double DBnum = 0.0;
  double DBdenom = 0.0;

if (YesNo[2] == 1){
    for (i=0; i<60; i++){

/*
  If these elements are zero, then something has gone with the minimization;
  take "deltaEXP" to be 1000. to force minimizer to turn around!
*/

      if(NEOS_U235[i]==0.0||NEOS_U235_0[i]==0.0){
        CheckZeros++;
      }
      if(NEOS_U238[i]==0.0||NEOS_U238_0[i]==0.0){
        CheckZeros++;
      }
      if(NEOS_Pu239[i]==0.0||NEOS_Pu239_0[i]==0.0){
        CheckZeros++;
      }
      if(NEOS_Pu241[i]==0.0||NEOS_Pu241_0[i]==0.0){
        CheckZeros++;
      }

      if(CheckZeros == 0){
        /* We take the double ratio (NEOS/NEOS_0)/(DB/DB_0) */
        NEOSnum = NEOS_U235[i]*NEOSF[0];
        NEOSnum += NEOS_U238[i]*NEOSF[1];
        NEOSnum += NEOS_Pu239[i]*NEOSF[2];
        NEOSnum += NEOS_Pu241[i]*NEOSF[3];

        NEOSdenom = NEOS_U235_0[i]*NEOSF[0];
        NEOSdenom += NEOS_U238_0[i]*NEOSF[1];
        NEOSdenom += NEOS_Pu239_0[i]*NEOSF[2];
        NEOSdenom += NEOS_Pu241_0[i]*NEOSF[3];

        DBnum = 0.0;
        DBdenom = 0.0;

        for (j=6+2*i; j<8+2*i; j++){
          DBnum += EH1_AD1_U235[j]* NEOSF[0];
          DBnum += EH1_AD1_U238[j]* NEOSF[1];
          DBnum += EH1_AD1_Pu239[j]* NEOSF[2];
          DBnum += EH1_AD1_Pu241[j]* NEOSF[3];

          DBnum += EH1_AD2_U235[j]* NEOSF[0];
          DBnum += EH1_AD2_U238[j]* NEOSF[1];
          DBnum += EH1_AD2_Pu239[j]* NEOSF[2];
          DBnum += EH1_AD2_Pu241[j]* NEOSF[3];
  
          DBdenom += EH1_AD1_U235_0[j]* NEOSF[0];
          DBdenom += EH1_AD1_U238_0[j]* NEOSF[1];
          DBdenom += EH1_AD1_Pu239_0[j]* NEOSF[2];
          DBdenom += EH1_AD1_Pu241_0[j]* NEOSF[3];

          DBdenom += EH1_AD2_U235_0[j]* NEOSF[0];
          DBdenom += EH1_AD2_U238_0[j]* NEOSF[1];
          DBdenom += EH1_AD2_Pu239_0[j]* NEOSF[2];
          DBdenom += EH1_AD2_Pu241_0[j]* NEOSF[3];
        }
        deltaNEOS[i] = NEOS_Ratio[i] - NEOSnum*DBdenom/(NEOSdenom*DBnum);
      }
      else{
        deltaNEOS[i] = 1000.;
      }
    }
  }
  /*****************
  *  DOUBLE CHOOZ  *
  ******************/

  /*
    We introduce a fudge factor for Double Chooz in order to reproduce their 
    determinations of theta_13 and Dm31/Dmee
  */
  double DC_fudge = 1.0026;

    /******************
    *  NEAR DETECTOR  *
    *******************/

  double *DC_ND_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *DC_ND_U238 = glbGetSignalFitRatePtr(exp+11, 1);
  double *DC_ND_Pu239 = glbGetSignalFitRatePtr(exp+11, 2);
  double *DC_ND_Pu241 = glbGetSignalFitRatePtr(exp+11, 3);
   
    /*****************
    *  FAR DETECTOR  *
    ******************/

  double *DC_FD_U235 = glbGetSignalFitRatePtr(exp+12, 0);
  double *DC_FD_U238 = glbGetSignalFitRatePtr(exp+12, 1);
  double *DC_FD_Pu239 = glbGetSignalFitRatePtr(exp+12, 2);
  double *DC_FD_Pu241 = glbGetSignalFitRatePtr(exp+12, 3);

  double TotND, TotFD;

if (YesNo[3] == 1){
  for (i=0; i<26; i++){
    TotND = DC_ND_U235[i]*DCF[0];
    TotND += DC_ND_U238[i]*DCF[1];
    TotND += DC_ND_Pu239[i]*DCF[2];
    TotND += DC_ND_Pu241[i]*DCF[3];

    TotFD = DC_FD_U235[i]*DCF[0];
    TotFD += DC_FD_U238[i]*DCF[1];
    TotFD += DC_FD_Pu239[i]*DCF[2];
    TotFD += DC_FD_Pu241[i]*DCF[3];

    deltaDC[i] = DC_Ratio[i] - DC_fudge*TotFD/TotND;
  }
}

  /**********
  *  BUGEY  *
  ***********/

  /* Get emin and emax for Bugey */
  glbGetEminEmax(exp+13, &emin, &emax);

  double Bugey_15_U235[35]; 
  double Bugey_15_U238[35];
  double Bugey_15_Pu239[35];
  double Bugey_15_Pu241[35];

  double Bugey_40_U235[35];
  double Bugey_40_U238[35];
  double Bugey_40_Pu239[35];
  double Bugey_40_Pu241[35];

  /******************
  *  15 m position  *
  *******************/

  glbShiftEnergyScale(0.01*x[2], glbGetSignalFitRatePtr(exp+13, 0),
                      Bugey_15_U235, 35, emin, emax);
  glbShiftEnergyScale(0.01*x[2], glbGetSignalFitRatePtr(exp+13, 1),
                      Bugey_15_U238, 35, emin, emax);
  glbShiftEnergyScale(0.01*x[2], glbGetSignalFitRatePtr(exp+13, 2),
                      Bugey_15_Pu239, 35, emin, emax);
  glbShiftEnergyScale(0.01*x[2], glbGetSignalFitRatePtr(exp+13, 3),
                      Bugey_15_Pu241, 35, emin, emax);

  /******************
  *  40 m position  *
  *******************/

  glbShiftEnergyScale(0.01*x[3], glbGetSignalFitRatePtr(exp+14, 0),
                      Bugey_40_U235, 35, emin, emax);
  glbShiftEnergyScale(0.01*x[3], glbGetSignalFitRatePtr(exp+14, 1),
                      Bugey_40_U238, 35, emin, emax);
  glbShiftEnergyScale(0.01*x[3], glbGetSignalFitRatePtr(exp+14, 2),
                      Bugey_40_Pu239, 35, emin, emax);
  glbShiftEnergyScale(0.01*x[3], glbGetSignalFitRatePtr(exp+14, 3),
                      Bugey_40_Pu241, 35, emin, emax);

CheckZeros = 0;

if (YesNo[4] == 1){
  for (i=5; i<30; i++){

/*
  If these elements are zero, then something has gone with the minimization;
  take "deltaEXP" to be 1000. to force minimizer to turn around!
*/

    if(Bugey_15_U235[i]==0.0||Bugey_40_U235[i]==0.0){
      CheckZeros++;
    }
    if(Bugey_15_U238[i]==0.0||Bugey_40_U238[i]==0.0){
      CheckZeros++;
    }
    if(Bugey_15_Pu239[i]==0.0||Bugey_40_Pu239[i]==0.0){
      CheckZeros++;
    }
    if(Bugey_15_Pu241[i]==0.0||Bugey_40_Pu241[i]==0.0){
      CheckZeros++;
    }

    if (CheckZeros == 0){
      numerator = Bugey_40_U235[i]*BugeyF[0];
      numerator += Bugey_40_U238[i]*BugeyF[1];
      numerator += Bugey_40_Pu239[i]*BugeyF[2];
      numerator += Bugey_40_Pu241[i]*BugeyF[3];

      denominator = Bugey_15_U235[i]*BugeyF[0];
      denominator += Bugey_15_U238[i]*BugeyF[1];
      denominator += Bugey_15_Pu239[i]*BugeyF[2];
      denominator += Bugey_15_Pu241[i]*BugeyF[3];

      deltaBugey[i-5] = Bugey_Ratio[i-5] - numerator/denominator;
    }
    else{
      deltaBugey[i-5] = 1000.;
    }
  }
}

  /*********
  *  RENO  *
  **********/

  /*
    We introduce a fudge factor for RENO in order to reproduce their determinations
    of theta_13 and Dm31/Dmee
  */
  double RENO_fudge = 1.00884;

    /******************
    *  NEAR DETECTOR  *
    *******************/

  double *RENO_ND_U235 = glbGetSignalFitRatePtr(exp+15, 0);
  double *RENO_ND_U238 = glbGetSignalFitRatePtr(exp+15, 1);
  double *RENO_ND_Pu239 = glbGetSignalFitRatePtr(exp+15, 2);
  double *RENO_ND_Pu241 = glbGetSignalFitRatePtr(exp+15, 3);
   
    /*****************
    *  FAR DETECTOR  *
    ******************/

  double *RENO_FD_U235 = glbGetSignalFitRatePtr(exp+16, 0);
  double *RENO_FD_U238 = glbGetSignalFitRatePtr(exp+16, 1);
  double *RENO_FD_Pu239 = glbGetSignalFitRatePtr(exp+16, 2);
  double *RENO_FD_Pu241 = glbGetSignalFitRatePtr(exp+16, 3);

if (YesNo[5] == 1){
  for (i=0; i<22; i++){ /* The first 22 bins are totally normal...  */
    TotND = RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD = RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];

    deltaRENO[i] = RENO_Ratio[i] - RENO_fudge*TotFD/TotND;
  }

  for (i=0; i<2; i++){ /* The next 2 bins are actually two bins combined...  */
    TotND = RENO_ND_U235[22+2*i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[22+2*i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[22+2*i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[22+2*i]* RENO_ND_F[3];

    TotND += RENO_ND_U235[22+2*i+1]* RENO_ND_F[0];
    TotND += RENO_ND_U238[22+2*i+1]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[22+2*i+1]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[22+2*i+1]* RENO_ND_F[3];

    TotFD = RENO_FD_U235[22+2*i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[22+2*i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[22+2*i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[22+2*i]* RENO_FD_F[3];

    TotFD += RENO_FD_U235[22+2*i+1]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[22+2*i+1]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[22+2*i+1]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[22+2*i+1]* RENO_FD_F[3];

    deltaRENO[i+22] = RENO_Ratio[i+22] - RENO_fudge*TotFD/TotND;
  }

  TotND = 0.0;
  TotFD = 0.0;

  for (i=26; i<29; i++){ /* The last bin is the sum of three bins  */
    TotND += RENO_ND_U235[i]* RENO_ND_F[0];
    TotND += RENO_ND_U238[i]* RENO_ND_F[1];
    TotND += RENO_ND_Pu239[i]* RENO_ND_F[2];
    TotND += RENO_ND_Pu241[i]* RENO_ND_F[3];

    TotFD += RENO_FD_U235[i]* RENO_FD_F[0];
    TotFD += RENO_FD_U238[i]* RENO_FD_F[1];
    TotFD += RENO_FD_Pu239[i]* RENO_FD_F[2];
    TotFD += RENO_FD_Pu241[i]* RENO_FD_F[3];
  }
  deltaRENO[24] = RENO_Ratio[24] - RENO_fudge*TotFD/TotND;
}

  /********************************
  *  Putting the pieces together  *
  *********************************/

  /* Adding the contributions from each experiment */

  double chi2 = 0.0;

  /* Starting with DANSS  */

if (YesNo[0] == 1){
  for (i=0; i<24; i++){
    for (j=0; j<24; j++){
      chi2 += deltaDANSS[i]*deltaDANSS[j]*DANSSCovariance[i][j];
    }
  }
}


  /* Next is Daya Bay piece (EH2/EH1 and EH3/EH1) */

if (YesNo[1] == 1){
  for (i=0; i<52; i++){
    for (j=0; j<52; j++){
      chi2 += deltaDB[i]*deltaDB[j]*DayaBayCovariance[i][j];
    }
  }
}

  /* Following that is NEOS/Daya Bay */

if (YesNo[2] == 1){
  for (i=0; i<60; i++){
    for (j=0; j<60; j++){
      chi2 += deltaNEOS[i]* deltaNEOS[j]* NEOS_Covariance[i][j];
    }
  }
}

  /* Following that is Double Chooz */

if (YesNo[3] == 1){
  for (i=0; i<26; i++){
    for (j=0; j<26; j++){
      chi2 += deltaDC[i]* deltaDC[j]*DC_Covariance[i][j];
    }
  }
}

  /* Following that is Bugey */

if (YesNo[4] == 1){
  for (i=0; i<25; i++){
    for (j=0; j<25; j++){
      chi2 += deltaBugey[i]*deltaBugey[j]*Bugey_Covariance[i][j];
    }
  }
}

  /* Following that is RENO */

if (YesNo[5] == 1){
  for (i=0; i<25; i++){
    for (j=0; j<25; j++){
      chi2 += deltaRENO[i]*deltaRENO[j]*RENO_Covariance[i][j];
    }
  }
}

/*
  Adding nuisance parameters; helps insure convergence of minimization to include
  nuisance parameters for experiments that are not included in fit...
*/

  for (i=0; i<Nnuisance; i++){
    if(fabs(x[i]) > errors[i])
      chi2 += square( (fabs(x[i]) ) / errors[i]);
  }

  /* Save the systematics parameters as starting values for the next step */

  for (i=0; i < Nnuisance; i++)
    systematic[i] = x[i];

  return chi2;
}  
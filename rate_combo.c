/* (C) 2019 Jeffrey M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "rate_combo.h"
#include "rate_funcs.h"

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

/************************************
*  INPUTS FOR THE RATE MEASUREMENTS *
*************************************/

/*
  This is the array containing the experimentally measured ratios;
  the order is given as follows:

  {Bugey4, Rovno91, Bugey3(15), Bugey3(40), Bugey(95), Gosgen(38), Gosgen(46),
	Gosgen(65), ILL, Krasnoyarsk87(33), Krasnoyarsk87(92), Krasnoyarsk94,
	Krasnoyarsk99, SRP-18, SRP-24, Rovno88-1I, Rovno-2I, Rovno-1S, Rovno-2S(25), 
	Rovno-2S(18), Nucifer, Palo Verde, Double Chooz, Chooz}
*/

/*
  These are the ratios of observed IBD events relative to the predicted number of events.
  It is critical to note that "predicted" means "having absolutely no oscillations" --
  not even known three-neutrino oscillations. Therefore, medium- and long-baseline
  experiments (PV, DC, Chooz, DB, RENO) won't prefer R = 1.0; they would prefer some
  value that depends on theta13 and Dm2ee. But because the oscillation engines include
  three-neutrino oscillations, this is the appropriate definition of Rexp.
*/

static const double Rexp[24] = {0.941, 0.939, 0.941, 0.947, 0.872, 0.972, 0.984, 0.940, 
		0.824, 0.936, 0.951, 0.945, 0.964, 0.917, 0.978, 0.905, 0.969, 0.960, 
		1.018, 0.930, 1.046, 0.971, 0.934, 0.976};

/*
  We have obtained the summation method calculation of the reactor flux; here we show the 
  same ratios as above with respect to these new predictions (see 1904.09358)
*/

static const double RexpSM[24] = {0.979, 0.983, 0.979, 0.985, 0.907, 1.013, 1.024, 0.975, 
		0.881, 1.001, 1.018, 1.011, 1.032, 0.981, 1.047, 0.946, 1.012, 1.003, 
		1.060, 0.971, 1.115, 1.014, 0.969, 1.013};

static const double RexpHKSS[24] = {0.9327, 0.9309, 0.9331, 0.9387, 0.8641, 0.9624, 
		0.9745, 0.9307, 0.8156, 0.9271, 0.9421, 0.9361, 0.9555, 0.9082, 0.9689, 
		0.8967, 0.9594, 0.9507, 1.0079, 0.9208, 1.0363, 0.9626, 0.9259, 0.9680};

/*
  These are arrays containing the fuel composition of each of the reactors, ordered:
  { U235, U238, Pu239, Pu241}
*/

static const double BugeyF[4] = {0.538, 0.078, 0.328, 0.056};
static const double Rovno_91F[4] = {0.614, 0.074, 0.274, 0.038};

/* Giunti uses different values than Vogel & Rovno91 paper - why? */

static const double Rovno88_1IF[4] = {0.607, 0.074, 0.277, 0.042};
static const double Rovno88_2IF[4] = {0.603, 0.076, 0.276, 0.045};
static const double Rovno88_1SF[4] = {0.606, 0.074, 0.277, 0.043};
static const double Rovno88_2S_25F[4] = {0.557, 0.076, 0.313, 0.054};
static const double Rovno88_2S_18F[4] = {0.606, 0.074, 0.274, 0.046};

static const double Gosgen38F[4] = {0.619, 0.067, 0.272, 0.042};
static const double Gosgen46F[4] = {0.584, 0.068, 0.298, 0.050};
static const double Gosgen65F[4] = {0.543, 0.070, 0.329, 0.058};

static const double NuciferF[4] = {0.926, 0.008, 0.061, 0.005};

/* These values are from the Nucifer paper; Giunti, et al., flipped 239 and 238! */

static const double PVF[4] = {0.600, 0.070, 0.270, 0.060};
static const double DCF[4] = {0.520, 0.087, 0.333, 0.060};
static const double ChoozF[4] = {0.496, 0.087, 0.351, 0.066};

/* 
  The (inverses of the) covariance matrix for the data (cov_matrix) and the
  correlation matrix of the normalization of the HM fluxes (V_inv). These have been
  calculated externally using Mathematica, taking into account the correlations described 
  in 1703.00860
*/

static const double cov_matrix[24][24]={{6802.72, -1700.68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-1700.68, 1700.68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 2650.2, -2270.55, -26.2913, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -2270.55, 2520.6, -17.3163, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -26.2913, -17.3163, 46.3025, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 377.17, -20.2866, -12.4825, -60.0549, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -20.2866, 377.17, -12.4825, -60.0549, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -12.4825, -12.4825, 236.878, -36.9523, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -60.0549, -60.0549, -36.9523, 148.146, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 411.167, -16.6083, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -16.6083, 24.7001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 566.893, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1111.11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1275.51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1189.06, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 262.458, -56.5211, -13.4043, -13.4043, -15.9841, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -56.5211, 262.458, -13.4043, -13.4043, -15.9841, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -13.4043, -13.4043, 201.095, -27.8423, -33.2009, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -13.4043, -13.4043, -27.8423, 201.095, -33.2009, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -15.9841, -15.9841, -33.2009, -33.2009, 233.409, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 87.3439, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 342.936, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10628.1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 976.562}};

/************************
*  INPUTS FOR DAYA BAY  *
*************************/

/*
  These are arrays containing the fuel fractions for each DB data point, ordered:
  { U235, U238, Pu239, Pu241}
*/

static const double DayaBay0[4] = {0.5113, 0.0767, 0.3445, 0.0675};
static const double DayaBay1[4] = {0.5279, 0.0766, 0.3326, 0.0629};
static const double DayaBay2[4] = {0.5418, 0.0764, 0.3219, 0.0599};
static const double DayaBay3[4] = {0.5553, 0.0762, 0.3113, 0.0572};
static const double DayaBay4[4] = {0.5699, 0.0760, 0.2992, 0.0549};
static const double DayaBay5[4] = {0.5849, 0.0758, 0.2878, 0.0515};
static const double DayaBay6[4] = {0.6033, 0.0757, 0.2744, 0.0466};
static const double DayaBay7[4] = {0.6304, 0.0754, 0.2525, 0.0417};

static const double *DayaBayFs[8] = {DayaBay0, DayaBay1, DayaBay2, DayaBay3,
					DayaBay4, DayaBay5, DayaBay6, DayaBay7};

/* These are the ratios of data/prediction for each of the bins */

static const double DBdata[8] = {0.9599, 0.9579, 0.9566, 0.9555, 0.9572, 0.9539,
		0.9528, 0.9558};

/*
  The ratio of the DB data with respect to the new flux calculation
*/

static const double DBdataSM[8] = {0.9970, 0.9960, 0.9955, 0.9952, 0.9978, 0.9953,
		0.9954, 1.0001};

static const double DBdataHKSS[8] = {0.9522, 0.9503, 0.9490, 0.9479, 0.9495, 0.9462, 
		0.9451, 0.9481};

/* 
  The (inverse of the) covariance matrix for Daya Bay (cov_matrix_DB), taken from the
  supplementary material to arXiv:1704.01082, including both systematic and statistical
  correlations. The ordering follows the ordering of the flux fractions laid out above.
*/

static const double cov_matrix_DB[8][8]={{209834., -46543.7, -42974., -62448.2, -16613.7, -17214.6, -30518.2, 6929.51}, {-46543.7, 223976., -72201.8, 10043.2, -28424.1, -29452.3, -56218., -603.994}, {-42974., -72201.8, 275079., -78655.1, -34661.9, -35915.7, 43006.7, -53588.2}, {-62448.2, 10043.2, -78655.1, 302389., -32260., -33427., -73774.7, -31704.2}, {-16613.7, -28424.1, -34661.9, -32260., 277938., -64534.4, -84319.8, -16783.6}, {-17214.6, -29452.3, -35915.7, -33427., -64534.4, 285657., -87369.9, -17390.7}, {-30518.2, -56218., 43006.7, -73774.7, -84319.8, -87369.9, 286871., 2671.42}, {6929.51, -603.994, -53588.2, -31704.2, -16783.6, -17390.7, 2671.42, 110349.}};

/********************
*  INPUTS FOR RENO  *
*********************/

/*
  These are arrays containing the fuel fractions for each RENO data point, ordered:
  { U235, U238, Pu239, Pu241}
*/

static const double RENO0[4] = {0.527, 0.074, 0.333, 0.066};
static const double RENO1[4] = {0.546, 0.073, 0.320, 0.061};
static const double RENO2[4] = {0.557, 0.075, 0.310, 0.058};
static const double RENO3[4] = {0.568, 0.074, 0.302, 0.056};
static const double RENO4[4] = {0.579, 0.074, 0.295, 0.052};
static const double RENO5[4] = {0.588, 0.075, 0.286, 0.051};
static const double RENO6[4] = {0.599, 0.073, 0.278, 0.050};
static const double RENO7[4] = {0.620, 0.072, 0.262, 0.046};

static const double *RENOFs[8] = {RENO0, RENO1, RENO2, RENO3, RENO4, RENO5, 
		RENO6, RENO7};

/* These are the ratios of data/prediction for each of the bins */

static const double RENOdata[8] = {0.9477, 0.9531, 0.9491, 0.9455, 0.9504, 0.9480,
		0.9493, 0.9447};

/*
  The ratio of the RENO data with respect to the new flux calculation
*/

static const double RENOdataSM[8] = {0.9858, 0.9923, 0.9901, 0.9857, 0.9900, 0.9881, 
		0.9871, 0.9814};

static const double RENOdataHKSS[8] = {0.9416, 0.9468, 0.9438, 0.9388, 0.9420, 0.9393,
		0.9373, 0.9303};

/* 
  The (inverse of the) covariance matrix for RENO (cov_matrix_RENO), taken from the
  supplementary material to arXiv:1806.00574, including both systematic and statistical
  correlations. The ordering follows the ordering of the flux fractions laid out above.
*/

static const double cov_matrix_RENO[8][8]={{112371., -16135.6, -16302., -16387.4, -14796.4, -16616.9, -16681.5, -15143.6}, {-16135.6, 102807., -14681.8, -14758.8, -13325.9, -14965.5, -15023.7, -13638.6}, {-16302., -14681.8, 103716., -14911., -13463.3, -15119.8, -15178.6, -13779.2}, {-16387.4, -14758.8, -14911., 104182., -13533.9, -15199.1, -15258.2, -13851.5}, {-14796.4, -13325.9, -13463.3, -13533.9, 95380.4, -13723.4, -13776.8, -12506.6}, {-16616.9, -14965.5, -15119.8, -15199.1, -13723.4, 105427., -15471.8, -14045.4}, {-16681.5, -15023.7, -15178.6, -15258.2, -13776.8, -15471.8, 105777., -14100.1}, {-15143.6, -13638.6, -13779.2, -13851.5, -12506.6, -14045.4, -14100.1, 97325.3}};

/**************************
*  OTHER RELEVANT INPUTS  *
***************************/

/* 
  The (inverse of the) correlation matrix of the normalization of the HM fluxes (V_inv).
  This has been calculated externally using Mathematica using Table III in
  arXiv:1703.00860

  I'm assuming that this also applies to Daya Bay -- I can't think of why it wouldn't...

  I'm assuming these theoretical errors still apply to the summation method calculation;
  this is almost certainly not correct, but let's roll with it!
*/

static const double V_inv_HM[4][4] = {{19.0983, 0., -6.04002, -12.8977}, 
	{0., 1., 0., 0.}, 
	{-6.04002, 0., 8.16908, -1.65811}, 
	{-12.8977, 0., -1.65811, 14.969}};

static const double V_inv_SM[4][4] = {{22.1693, 0., -6.79549, -15.2312}, 
	{0., 1., 0., 0.}, 
	{-6.79549, 0., 8.86624, -1.59453}, 
	{-15.2312, 0., -1.59453, 17.2477}};

static const double V_inv_HKSS[4][4] = {{24.1592, -0.0109924, -7.75505, -16.2475}, 
	{-0.0109924, 1.00357, 0.0300266, -0.0764211}, 
	{-7.75505, 0.0300266, 10.3186, -2.0983}, 
	{-16.2475, -0.0764211, -2.0983, 18.7623}};

/*
  As of present, this code contains the ability to toggle off blocks of experiments
  that are mutually correlated. The reason for not *quite* doing this at the level
  of individual experiments across the board is that throwing out different experiments
  *before* computing the inverse of the covariance matrix isn't the same as throwing out
  those parts of the inverse of the covariance matrix if that experiment is correlated
  with other ones.

  Until such a time as I can come up with a better scheme, this is how we're doing it!

  The blocks of experiments are:
	0: Bugey-4 + Rovno91
	1: Bugey-3 (15 + 40 + 95)
	2: Gosgen (38 + 46 + 65) + ILL
	3: Krasnoyarsk87 (33 + 92)
	4: Krasnoyarsk94
	5: Krasnoyarsk99
	6: Savannah River (18)
	7: Savannah River (24)
	8: Rovno88 (1I + 2I + 1S + 2S (25) + 3S (18))
	9: Nucifer
	10: Palo Verde
	11: Double Chooz
	12: Chooz
	13: Daya Bay -- fuel evolution
	14: RENO -- fuel evolution

  Note that Palo Verde, Double Chooz and Chooz are corrected for 3-nu oscillations
  in their common oscillation engine; th12, th13, Dm21 and Dm31 are hard-coded in there, 
  for the time being!

*/

/*****************************************************************
*                                                                *
*        THE CALCULATION OF A CHI-SQUARED WITH SYSTEMATICS       *
*                                                                *
******************************************************************/

double combo_rate_chi(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following, for instance, 1703.00860

  "delta" is given by (Rexp - Rpred)

  "user_data" contains the list of experiments being considered in a particular
  invocation of the code. 
*/

  int k;

  double delta[40] = {0.0};
  int YesNo[26];
  for (k=0; k<26; k++){
    YesNo[k] = ((int *)user_data)[k];
  }

/************
*  Bugey-4  *
*************/

  double *Bugey_4_U235 = glbGetSignalFitRatePtr(exp, 0);
  double *Bugey_4_U238 = glbGetSignalFitRatePtr(exp, 1);
  double *Bugey_4_Pu239 = glbGetSignalFitRatePtr(exp, 2);
  double *Bugey_4_Pu241 = glbGetSignalFitRatePtr(exp, 3);

  double *Bugey_4_U235_0 = glbGetRuleRatePtr(exp, 0);
  double *Bugey_4_U238_0 = glbGetRuleRatePtr(exp, 1);
  double *Bugey_4_Pu239_0 = glbGetRuleRatePtr(exp, 2);
  double *Bugey_4_Pu241_0 = glbGetRuleRatePtr(exp, 3);
   
  double Bugey_4_RatioN = 0.0;
  double Bugey_4_RatioD = 0.0;

if(YesNo[0]==1||YesNo[2]==1){
  Bugey_4_RatioN = x[0]*BugeyF[0]*Bugey_4_U235[0];
  Bugey_4_RatioN += x[1]*BugeyF[1]*Bugey_4_U238[0];
  Bugey_4_RatioN += x[2]*BugeyF[2]*Bugey_4_Pu239[0];
  Bugey_4_RatioN += x[3]*BugeyF[3]*Bugey_4_Pu241[0];

  Bugey_4_RatioD = BugeyF[0]*Bugey_4_U235_0[0];
  Bugey_4_RatioD += BugeyF[1]*Bugey_4_U238_0[0];
  Bugey_4_RatioD += BugeyF[2]*Bugey_4_Pu239_0[0];
  Bugey_4_RatioD += BugeyF[3]*Bugey_4_Pu241_0[0];
}
if(YesNo[0]==1){
  delta[0] = Rexp[0] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/************
*  Rovno91  *
*************/

  double *Rovno91_U235 = glbGetSignalFitRatePtr(exp+1, 0);
  double *Rovno91_U238 = glbGetSignalFitRatePtr(exp+1, 1);
  double *Rovno91_Pu239 = glbGetSignalFitRatePtr(exp+1, 2);
  double *Rovno91_Pu241 = glbGetSignalFitRatePtr(exp+1, 3);

  double *Rovno91_U235_0 = glbGetRuleRatePtr(exp+1, 0);
  double *Rovno91_U238_0 = glbGetRuleRatePtr(exp+1, 1);
  double *Rovno91_Pu239_0 = glbGetRuleRatePtr(exp+1, 2);
  double *Rovno91_Pu241_0 = glbGetRuleRatePtr(exp+1, 3);

  double Rovno91_RatioN = 0.0;
  double Rovno91_RatioD = 0.0;

if(YesNo[1]==1){
  Rovno91_RatioN = x[0]*Rovno_91F[0]*Rovno91_U235[0];
  Rovno91_RatioN += x[1]*Rovno_91F[1]*Rovno91_U238[0];
  Rovno91_RatioN += x[2]*Rovno_91F[2]*Rovno91_Pu239[0];
  Rovno91_RatioN += x[3]*Rovno_91F[3]*Rovno91_Pu241[0];

  Rovno91_RatioD = Rovno_91F[0]*Rovno91_U235_0[0];
  Rovno91_RatioD += Rovno_91F[1]*Rovno91_U238_0[0];
  Rovno91_RatioD += Rovno_91F[2]*Rovno91_Pu239_0[0];
  Rovno91_RatioD += Rovno_91F[3]*Rovno91_Pu241_0[0];

  delta[1] = Rexp[1] - Rovno91_RatioN/Rovno91_RatioD;
}

/*******************
*  Bugey-3 (15 m)  *
********************/

/*
  No need to recalculate rates -- reuse the output for Bugey-4!
*/

if(YesNo[2]==1){
  delta[2] = Rexp[2] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/*******************
*  Bugey-3 (40 m)  *
********************/

  double *Bugey_3_40_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *Bugey_3_40_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);
   
  double *Bugey_3_40_U235_0 = glbGetRuleRatePtr(exp+2, 0);
  double *Bugey_3_40_U238_0 = glbGetRuleRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239_0 = glbGetRuleRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241_0 = glbGetRuleRatePtr(exp+2, 3);
   
  double Bugey_3_40_RatioN = 0.0;
  double Bugey_3_40_RatioD = 0.0;

if(YesNo[3]==1){   
  Bugey_3_40_RatioN = x[0]*BugeyF[0]*Bugey_3_40_U235[0];
  Bugey_3_40_RatioN += x[1]*BugeyF[1]*Bugey_3_40_U238[0];
  Bugey_3_40_RatioN += x[2]*BugeyF[2]*Bugey_3_40_Pu239[0];
  Bugey_3_40_RatioN += x[3]*BugeyF[3]*Bugey_3_40_Pu241[0];

  Bugey_3_40_RatioD = BugeyF[0]*Bugey_3_40_U235_0[0];
  Bugey_3_40_RatioD += BugeyF[1]*Bugey_3_40_U238_0[0];
  Bugey_3_40_RatioD += BugeyF[2]*Bugey_3_40_Pu239_0[0];
  Bugey_3_40_RatioD += BugeyF[3]*Bugey_3_40_Pu241_0[0];

  delta[3] = Rexp[3] - Bugey_3_40_RatioN/Bugey_3_40_RatioD;
}

/*******************
*  Bugey-3 (95 m)  *
********************/

  double *Bugey_3_95_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *Bugey_3_95_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);
   
  double *Bugey_3_95_U235_0 = glbGetRuleRatePtr(exp+3, 0);
  double *Bugey_3_95_U238_0 = glbGetRuleRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239_0 = glbGetRuleRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241_0 = glbGetRuleRatePtr(exp+3, 3);
   
  double Bugey_3_95_RatioN = 0.0;
  double Bugey_3_95_RatioD = 0.0;

if(YesNo[4]==1){   
  Bugey_3_95_RatioN = x[0]*BugeyF[0]*Bugey_3_95_U235[0];
  Bugey_3_95_RatioN += x[1]*BugeyF[1]*Bugey_3_95_U238[0];
  Bugey_3_95_RatioN += x[2]*BugeyF[2]*Bugey_3_95_Pu239[0];
  Bugey_3_95_RatioN += x[3]*BugeyF[3]*Bugey_3_95_Pu241[0];

  Bugey_3_95_RatioD = BugeyF[0]*Bugey_3_95_U235_0[0];
  Bugey_3_95_RatioD += BugeyF[1]*Bugey_3_95_U238_0[0];
  Bugey_3_95_RatioD += BugeyF[2]*Bugey_3_95_Pu239_0[0];
  Bugey_3_95_RatioD += BugeyF[3]*Bugey_3_95_Pu241_0[0];

  delta[4] = Rexp[4] - Bugey_3_95_RatioN/Bugey_3_95_RatioD;
}

/********************
*  Gosgen (37.9 m)  *
*********************/

  double *Gosgen_38_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *Gosgen_38_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
  double *Gosgen_38_U235_0 = glbGetRuleRatePtr(exp+4, 0);
  double *Gosgen_38_U238_0 = glbGetRuleRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239_0 = glbGetRuleRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241_0 = glbGetRuleRatePtr(exp+4, 3);
   
  double Gosgen_38_RatioN = 0.0;
  double Gosgen_38_RatioD = 0.0;

if(YesNo[5]==1){
   
  Gosgen_38_RatioN = x[0]*Gosgen38F[0]*Gosgen_38_U235[0];
  Gosgen_38_RatioN += x[1]*Gosgen38F[1]*Gosgen_38_U238[0];
  Gosgen_38_RatioN += x[2]*Gosgen38F[2]*Gosgen_38_Pu239[0];
  Gosgen_38_RatioN += x[3]*Gosgen38F[3]*Gosgen_38_Pu241[0];

  Gosgen_38_RatioD = Gosgen38F[0]*Gosgen_38_U235_0[0];
  Gosgen_38_RatioD += Gosgen38F[1]*Gosgen_38_U238_0[0];
  Gosgen_38_RatioD += Gosgen38F[2]*Gosgen_38_Pu239_0[0];
  Gosgen_38_RatioD += Gosgen38F[3]*Gosgen_38_Pu241_0[0];

  delta[5] = Rexp[5] - Gosgen_38_RatioN/Gosgen_38_RatioD;
}

/********************
*  Gosgen (45.9 m)  *
*********************/

  double *Gosgen_46_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *Gosgen_46_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);
   
  double *Gosgen_46_U235_0 = glbGetRuleRatePtr(exp+5, 0);
  double *Gosgen_46_U238_0 = glbGetRuleRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239_0 = glbGetRuleRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241_0 = glbGetRuleRatePtr(exp+5, 3);
   
  double Gosgen_46_RatioN = 0.0;
  double Gosgen_46_RatioD = 0.0;

if(YesNo[6]==1){
   
  Gosgen_46_RatioN = x[0]*Gosgen46F[0]*Gosgen_46_U235[0];
  Gosgen_46_RatioN += x[1]*Gosgen46F[1]*Gosgen_46_U238[0];
  Gosgen_46_RatioN += x[2]*Gosgen46F[2]*Gosgen_46_Pu239[0];
  Gosgen_46_RatioN += x[3]*Gosgen46F[3]*Gosgen_46_Pu241[0];

  Gosgen_46_RatioD = Gosgen46F[0]*Gosgen_46_U235_0[0];
  Gosgen_46_RatioD += Gosgen46F[1]*Gosgen_46_U238_0[0];
  Gosgen_46_RatioD += Gosgen46F[2]*Gosgen_46_Pu239_0[0];
  Gosgen_46_RatioD += Gosgen46F[3]*Gosgen_46_Pu241_0[0];

  delta[6] = Rexp[6] - Gosgen_46_RatioN/Gosgen_46_RatioD;
}

/********************
*  Gosgen (64.7 m)  *
*********************/

  double *Gosgen_65_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *Gosgen_65_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);

  double *Gosgen_65_U235_0 = glbGetRuleRatePtr(exp+6, 0);
  double *Gosgen_65_U238_0 = glbGetRuleRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239_0 = glbGetRuleRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241_0 = glbGetRuleRatePtr(exp+6, 3);
   
  double Gosgen_65_RatioN = 0.0;
  double Gosgen_65_RatioD = 0.0;

if(YesNo[7]==1){
   
  Gosgen_65_RatioN = x[0]*Gosgen65F[0]*Gosgen_65_U235[0];
  Gosgen_65_RatioN += x[1]*Gosgen65F[1]*Gosgen_65_U238[0];
  Gosgen_65_RatioN += x[2]*Gosgen65F[2]*Gosgen_65_Pu239[0];
  Gosgen_65_RatioN += x[3]*Gosgen65F[3]*Gosgen_65_Pu241[0];

  Gosgen_65_RatioD = Gosgen65F[0]*Gosgen_65_U235_0[0];
  Gosgen_65_RatioD += Gosgen65F[1]*Gosgen_65_U238_0[0];
  Gosgen_65_RatioD += Gosgen65F[2]*Gosgen_65_Pu239_0[0];
  Gosgen_65_RatioD += Gosgen65F[3]*Gosgen_65_Pu241_0[0];

  delta[7] = Rexp[7] - Gosgen_65_RatioN/Gosgen_65_RatioD;
}

/*****************
*  ILL (8.76 m)  *
******************/

  double *ILL_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *ILL_U235_0 = glbGetRuleRatePtr(exp+7, 0);

if(YesNo[8]==1){
  delta[8] = Rexp[8] - x[0]*ILL_U235[0]/ILL_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (32.8 m)  *
****************************/

  double *Krasnoyarsk_33_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *Krasnoyarsk_33_U235_0 = glbGetRuleRatePtr(exp+8, 0);

if(YesNo[9]==1){
  delta[9] = Rexp[9] - x[0]* Krasnoyarsk_33_U235[0]/Krasnoyarsk_33_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (92.3 m)  *
****************************/

  double *Krasnoyarsk_92_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *Krasnoyarsk_92_U235_0 = glbGetRuleRatePtr(exp+9, 0);

if(YesNo[10]==1){
  delta[10] = Rexp[10] - x[0]* Krasnoyarsk_92_U235[0]/Krasnoyarsk_92_U235_0[0];
}

/***************************
*  Krasnoyarsk94 (57.0 m)  *
****************************/

  double *Krasnoyarsk_57_U235 = glbGetSignalFitRatePtr(exp+10, 0);
  double *Krasnoyarsk_57_U235_0 = glbGetRuleRatePtr(exp+10, 0);

if(YesNo[11]==1){
  delta[11] = Rexp[11] - x[0]* Krasnoyarsk_57_U235[0]/Krasnoyarsk_57_U235_0[0];
}

/***************************
*  Krasnoyarsk99 (34.0 m)  *
****************************/

  double *Krasnoyarsk_34_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *Krasnoyarsk_34_U235_0 = glbGetRuleRatePtr(exp+11, 0);

if(YesNo[12]==1){
  delta[12] = Rexp[12] - x[0]* Krasnoyarsk_34_U235[0]/Krasnoyarsk_34_U235_0[0];
}

/********************
*  SRP-18 (18.2 m)  *
*********************/

  double *SRP_18_U235 = glbGetSignalFitRatePtr(exp+12, 0);
  double *SRP_18_U235_0 = glbGetRuleRatePtr(exp+12, 0);

if(YesNo[13]==1){
  delta[13] = Rexp[13] - x[0]* SRP_18_U235[0]/SRP_18_U235_0[0];
}

/********************
*  SRP-24 (23.8 m)  *
*********************/

  double *SRP_24_U235 = glbGetSignalFitRatePtr(exp+13, 0);
  double *SRP_24_U235_0 = glbGetRuleRatePtr(exp+13, 0);

if(YesNo[14]==1){
  delta[14] = Rexp[14] - x[0]*SRP_24_U235[0]/SRP_24_U235_0[0];
}

/***************
*  Rovno88-1I  *
****************/

  double *Rovno88_1I_U235 = glbGetSignalFitRatePtr(exp+14, 0);
  double *Rovno88_1I_U238 = glbGetSignalFitRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239 = glbGetSignalFitRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241 = glbGetSignalFitRatePtr(exp+14, 3);

  double *Rovno88_1I_U235_0 = glbGetRuleRatePtr(exp+14, 0);
  double *Rovno88_1I_U238_0 = glbGetRuleRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239_0 = glbGetRuleRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241_0 = glbGetRuleRatePtr(exp+14, 3);

  double Rovno88_1I_RatioN = 0.0;
  double Rovno88_1I_RatioD = 0.0;

if(YesNo[15]==1){
  Rovno88_1I_RatioN = x[0]* Rovno88_1IF[0]*Rovno88_1I_U235[0];
  Rovno88_1I_RatioN += x[1]* Rovno88_1IF[1]*Rovno88_1I_U238[0];
  Rovno88_1I_RatioN += x[2]* Rovno88_1IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1I_RatioN += x[3]* Rovno88_1IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1I_RatioD = Rovno88_1IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[3]*Rovno88_1I_Pu241_0[0];

  delta[15] = Rexp[15] - Rovno88_1I_RatioN/Rovno88_1I_RatioD;
}

/***************
*  Rovno88-2I  *
****************/

  double Rovno88_2I_RatioN = 0.0;
  double Rovno88_2I_RatioD = 0.0;

if(YesNo[16]==1){   
  Rovno88_2I_RatioN = x[0]* Rovno88_2IF[0]*Rovno88_1I_U235[0];
  Rovno88_2I_RatioN += x[1]* Rovno88_2IF[1]*Rovno88_1I_U238[0];
  Rovno88_2I_RatioN += x[2]* Rovno88_2IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_2I_RatioN += x[3]* Rovno88_2IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_2I_RatioD = Rovno88_2IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[3]*Rovno88_1I_Pu241_0[0];

  delta[16] = Rexp[16] - Rovno88_2I_RatioN/Rovno88_2I_RatioD;
}
/***************
*  Rovno88-1S  *
****************/

  double Rovno88_1S_RatioN = 0.0;
  double Rovno88_1S_RatioD = 0.0;

if(YesNo[17]==1){
  Rovno88_1S_RatioN = x[0]* Rovno88_1SF[0]*Rovno88_1I_U235[0];
  Rovno88_1S_RatioN += x[1]* Rovno88_1SF[1]*Rovno88_1I_U238[0];
  Rovno88_1S_RatioN += x[2]* Rovno88_1SF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1S_RatioN += x[3]* Rovno88_1SF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1S_RatioD = Rovno88_1SF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[3]*Rovno88_1I_Pu241_0[0];

  delta[17] = Rexp[17] - Rovno88_1S_RatioN/Rovno88_1S_RatioD;
}

/**********************
*  Rovno88-2S (25 m)  *
***********************/
   
  double Rovno88_2S_25_RatioN = 0.0;
  double Rovno88_2S_25_RatioD = 0.0;

  double *Rovno88_2S_25_U235 = glbGetSignalFitRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238 = glbGetSignalFitRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239 = glbGetSignalFitRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241 = glbGetSignalFitRatePtr(exp+15, 3);

  double *Rovno88_2S_25_U235_0 = glbGetRuleRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238_0 = glbGetRuleRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239_0 = glbGetRuleRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241_0 = glbGetRuleRatePtr(exp+15, 3);

if(YesNo[18]==1){
  Rovno88_2S_25_RatioN = x[0]*Rovno88_2S_25F[0]*Rovno88_2S_25_U235[0];
  Rovno88_2S_25_RatioN += x[1]*Rovno88_2S_25F[1]*Rovno88_2S_25_U238[0];
  Rovno88_2S_25_RatioN += x[2]*Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239[0];
  Rovno88_2S_25_RatioN += x[3]*Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241[0];

  Rovno88_2S_25_RatioD = Rovno88_2S_25F[0]*Rovno88_2S_25_U235_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[1]*Rovno88_2S_25_U238_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241_0[0];

  delta[18] = Rexp[18] - Rovno88_2S_25_RatioN/Rovno88_2S_25_RatioD;
}

/**********************
*  Rovno88-2S (18 m)  *
***********************/

  double Rovno88_2S_18_RatioN = 0.0;
  double Rovno88_2S_18_RatioD = 0.0;

if(YesNo[19]==1){
  Rovno88_2S_18_RatioN = x[0]*Rovno88_2S_18F[0]*Rovno88_1I_U235[0];
  Rovno88_2S_18_RatioN += x[1]*Rovno88_2S_18F[1]*Rovno88_1I_U238[0];
  Rovno88_2S_18_RatioN += x[2]*Rovno88_2S_18F[2]*Rovno88_1I_Pu239[0];
  Rovno88_2S_18_RatioN += x[3]*Rovno88_2S_18F[3]*Rovno88_1I_Pu241[0];

  Rovno88_2S_18_RatioD = Rovno88_2S_18F[0]*Rovno88_1I_U235_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[1]*Rovno88_1I_U238_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[3]*Rovno88_1I_Pu241_0[0];

  delta[19] = Rexp[19] - Rovno88_2S_18_RatioN/Rovno88_2S_18_RatioD;
}

/************
*  Nucifer  *
*************/

  double *Nucifer_U235 = glbGetSignalFitRatePtr(exp+16, 0);
  double *Nucifer_U238 = glbGetSignalFitRatePtr(exp+16, 1);
  double *Nucifer_Pu239 = glbGetSignalFitRatePtr(exp+16, 2);
  double *Nucifer_Pu241 = glbGetSignalFitRatePtr(exp+16, 3);

  double *Nucifer_U235_0 = glbGetRuleRatePtr(exp+16, 0);
  double *Nucifer_U238_0 = glbGetRuleRatePtr(exp+16, 1);
  double *Nucifer_Pu239_0 = glbGetRuleRatePtr(exp+16, 2);
  double *Nucifer_Pu241_0 = glbGetRuleRatePtr(exp+16, 3);

  double Nucifer_RatioN=0.0;
  double Nucifer_RatioD=0.0;

if(YesNo[20]==1){
  Nucifer_RatioN = x[0]*NuciferF[0]*Nucifer_U235[0];
  Nucifer_RatioN += x[1]*NuciferF[1]*Nucifer_U238[0];
  Nucifer_RatioN += x[2]*NuciferF[2]*Nucifer_Pu239[0];
  Nucifer_RatioN += x[3]*NuciferF[3]*Nucifer_Pu241[0];

  Nucifer_RatioD = NuciferF[0]*Nucifer_U235_0[0];
  Nucifer_RatioD += NuciferF[1]*Nucifer_U238_0[0];
  Nucifer_RatioD += NuciferF[2]*Nucifer_Pu239_0[0];
  Nucifer_RatioD += NuciferF[3]*Nucifer_Pu241_0[0];

  delta[20] = Rexp[20] - Nucifer_RatioN/Nucifer_RatioD;
}

/**********************************
*  Palo Verde -- 750 m and 890 m  *
***********************************/

  double PV_RatioN=0.0;
  double PV_RatioD=0.0;

  double *PV_750_U235 = glbGetSignalFitRatePtr(exp+17, 0);
  double *PV_750_U238 = glbGetSignalFitRatePtr(exp+17, 1);
  double *PV_750_Pu239 = glbGetSignalFitRatePtr(exp+17, 2);
  double *PV_750_Pu241 = glbGetSignalFitRatePtr(exp+17, 3);

  double *PV_750_U235_0 = glbGetRuleRatePtr(exp+17, 0);
  double *PV_750_U238_0 = glbGetRuleRatePtr(exp+17, 1);
  double *PV_750_Pu239_0 = glbGetRuleRatePtr(exp+17, 2);
  double *PV_750_Pu241_0 = glbGetRuleRatePtr(exp+17, 3);

  double *PV_890_U235 = glbGetSignalFitRatePtr(exp+18, 0);
  double *PV_890_U238 = glbGetSignalFitRatePtr(exp+18, 1);
  double *PV_890_Pu239 = glbGetSignalFitRatePtr(exp+18, 2);
  double *PV_890_Pu241 = glbGetSignalFitRatePtr(exp+18, 3);

  double *PV_890_U235_0 = glbGetRuleRatePtr(exp+18, 0);
  double *PV_890_U238_0 = glbGetRuleRatePtr(exp+18, 1);
  double *PV_890_Pu239_0 = glbGetRuleRatePtr(exp+18, 2);
  double *PV_890_Pu241_0 = glbGetRuleRatePtr(exp+18, 3);

if(YesNo[21]==1){
  PV_RatioN = x[0]*PVF[0]*(PV_750_U235[0] + PV_890_U235[0]);
  PV_RatioN += x[1]*PVF[1]*(PV_750_U238[0] + PV_890_U238[0]);
  PV_RatioN += x[2]*PVF[2]*(PV_750_Pu239[0] + PV_890_Pu239[0]);
  PV_RatioN += x[3]*PVF[3]*(PV_750_Pu241[0] + PV_890_Pu241[0]);

  PV_RatioD = PVF[0]*(PV_750_U235_0[0]+PV_890_U235_0[0]);
  PV_RatioD += PVF[1]*(PV_750_U238_0[0]+PV_890_U238_0[0]);
  PV_RatioD += PVF[2]*(PV_750_Pu239_0[0]+PV_890_Pu239_0[0]);
  PV_RatioD += PVF[3]*(PV_750_Pu241_0[0]+PV_890_Pu241_0[0]);

  delta[21] = Rexp[21] - PV_RatioN/PV_RatioD;
}

/************************************
*  Double Chooz -- 355 m and 469 m  *
*************************************/

  double DC_RatioN=0.0;
  double DC_RatioD=0.0;

  double *DC_355_U235 = glbGetSignalFitRatePtr(exp+19, 0);
  double *DC_355_U238 = glbGetSignalFitRatePtr(exp+19, 1);
  double *DC_355_Pu239 = glbGetSignalFitRatePtr(exp+19, 2);
  double *DC_355_Pu241 = glbGetSignalFitRatePtr(exp+19, 3);

  double *DC_355_U235_0 = glbGetRuleRatePtr(exp+19, 0);
  double *DC_355_U238_0 = glbGetRuleRatePtr(exp+19, 1);
  double *DC_355_Pu239_0 = glbGetRuleRatePtr(exp+19, 2);
  double *DC_355_Pu241_0 = glbGetRuleRatePtr(exp+19, 3);

  double *DC_469_U235 = glbGetSignalFitRatePtr(exp+20, 0);
  double *DC_469_U238 = glbGetSignalFitRatePtr(exp+20, 1);
  double *DC_469_Pu239 = glbGetSignalFitRatePtr(exp+20, 2);
  double *DC_469_Pu241 = glbGetSignalFitRatePtr(exp+20, 3);

  double *DC_469_U235_0 = glbGetRuleRatePtr(exp+20, 0);
  double *DC_469_U238_0 = glbGetRuleRatePtr(exp+20, 1);
  double *DC_469_Pu239_0 = glbGetRuleRatePtr(exp+20, 2);
  double *DC_469_Pu241_0 = glbGetRuleRatePtr(exp+20, 3);

if(YesNo[22]==1){
  DC_RatioN = x[0]*DCF[0]*(DC_355_U235[0] + DC_469_U235[0]);
  DC_RatioN += x[1]*DCF[1]*(DC_355_U238[0] + DC_469_U238[0]);
  DC_RatioN += x[2]*DCF[2]*(DC_355_Pu239[0] + DC_469_Pu239[0]);
  DC_RatioN += x[3]*DCF[3]*(DC_355_Pu241[0] + DC_469_Pu241[0]);

  DC_RatioD = DCF[0]*(DC_355_U235_0[0] + DC_469_U235_0[0]);
  DC_RatioD += DCF[1]*(DC_355_U238_0[0] + DC_469_U238_0[0]);
  DC_RatioD += DCF[2]*(DC_355_Pu239_0[0] + DC_469_Pu239_0[0]);
  DC_RatioD += DCF[3]*(DC_355_Pu241_0[0] + DC_469_Pu241_0[0]);

  delta[22] = Rexp[22] - DC_RatioN/DC_RatioD;
}

/******************************
*  Chooz -- 998 m and 1115 m  *
*******************************/

  double Chooz_RatioN=0.0;
  double Chooz_RatioD=0.0;

  double *Chooz_998_U235 = glbGetSignalFitRatePtr(exp+21, 0);
  double *Chooz_998_U238 = glbGetSignalFitRatePtr(exp+21, 1);
  double *Chooz_998_Pu239 = glbGetSignalFitRatePtr(exp+21, 2);
  double *Chooz_998_Pu241 = glbGetSignalFitRatePtr(exp+21, 3);

  double *Chooz_998_U235_0 = glbGetRuleRatePtr(exp+21, 0);
  double *Chooz_998_U238_0 = glbGetRuleRatePtr(exp+21, 1);
  double *Chooz_998_Pu239_0 = glbGetRuleRatePtr(exp+21, 2);
  double *Chooz_998_Pu241_0 = glbGetRuleRatePtr(exp+21, 3);

  double *Chooz_1115_U235 = glbGetSignalFitRatePtr(exp+22, 0);
  double *Chooz_1115_U238 = glbGetSignalFitRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239 = glbGetSignalFitRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241 = glbGetSignalFitRatePtr(exp+22, 3);

  double *Chooz_1115_U235_0 = glbGetRuleRatePtr(exp+22, 0);
  double *Chooz_1115_U238_0 = glbGetRuleRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239_0 = glbGetRuleRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241_0 = glbGetRuleRatePtr(exp+22, 3);

if(YesNo[23]==1){
  Chooz_RatioN = x[0]*ChoozF[0]*(Chooz_998_U235[0] + Chooz_1115_U235[0]);
  Chooz_RatioN += x[1]*ChoozF[1]*(Chooz_998_U238[0] + Chooz_1115_U238[0]);
  Chooz_RatioN += x[2]*ChoozF[2]*(Chooz_998_Pu239[0] + Chooz_1115_Pu239[0]);
  Chooz_RatioN += x[3]*ChoozF[3]*(Chooz_998_Pu241[0] + Chooz_1115_Pu241[0]);

  Chooz_RatioD = ChoozF[0]*(Chooz_998_U235_0[0] + Chooz_1115_U235_0[0]);
  Chooz_RatioD += ChoozF[1]*(Chooz_998_U238_0[0] + Chooz_1115_U238_0[0]);
  Chooz_RatioD += ChoozF[2]*(Chooz_998_Pu239_0[0] + Chooz_1115_Pu239_0[0]);
  Chooz_RatioD += ChoozF[3]*(Chooz_998_Pu241_0[0] + Chooz_1115_Pu241_0[0]);

  delta[23] = Rexp[23] - Chooz_RatioN/Chooz_RatioD;
}

/************
*  EH1 AD1  *
*************/

  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+23, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+23, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+23, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+23, 3);
   
  double EH1_AD1_Total = 0.0;
  double EH1_AD1_Total_0 = 0.0;

/************
*  EH1 AD2  *
*************/

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+24, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+24, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+24, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+24, 3);

  double EH1_AD2_Total = 0.0;
  double EH1_AD2_Total_0 = 0.0;

/************
*  EH2 AD3  *
*************/

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+25, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+25, 3);

  double *EH2_AD3_U235_0 = glbGetRuleRatePtr(exp+25, 0);
  double *EH2_AD3_U238_0 = glbGetRuleRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239_0 = glbGetRuleRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241_0 = glbGetRuleRatePtr(exp+25, 3);
   
  double EH2_AD3_Total = 0.0;
  double EH2_AD3_Total_0 = 0.0;

/************
*  EH2 AD8  *
*************/

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+26, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+26, 3);

  double *EH2_AD8_U235_0 = glbGetRuleRatePtr(exp+26, 0);
  double *EH2_AD8_U238_0 = glbGetRuleRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239_0 = glbGetRuleRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241_0 = glbGetRuleRatePtr(exp+26, 3);

  double EH2_AD8_Total = 0.0;
  double EH2_AD8_Total_0 = 0.0;

/*********
*  RENO  *
**********/

  double *RENO_U235 = glbGetSignalFitRatePtr(exp+27, 0);
  double *RENO_U238 = glbGetSignalFitRatePtr(exp+27, 1);
  double *RENO_Pu239 = glbGetSignalFitRatePtr(exp+27, 2);
  double *RENO_Pu241 = glbGetSignalFitRatePtr(exp+27, 3);

  double *RENO_U235_0 = glbGetRuleRatePtr(exp+27, 0);
  double *RENO_U238_0 = glbGetRuleRatePtr(exp+27, 1);
  double *RENO_Pu239_0 = glbGetRuleRatePtr(exp+27, 2);
  double *RENO_Pu241_0 = glbGetRuleRatePtr(exp+27, 3);

  double RENO_Total = 0.0;
  double RENO_Total_0 = 0.0;

/*********************************
*  Assembling Daya Bay and RENO  *
**********************************/

  double chi2 = 0.0;
  int i,j;

  double num, denom;

  for (i=0; i<8; i++){
  if (YesNo[24] == 1){
    EH1_AD1_Total = x[0]*DayaBayFs[i][0]*EH1_AD1_U235[0];
    EH1_AD1_Total += x[1]*DayaBayFs[i][1]*EH1_AD1_U238[0];
    EH1_AD1_Total += x[2]*DayaBayFs[i][2]*EH1_AD1_Pu239[0];
    EH1_AD1_Total += x[3]*DayaBayFs[i][3]*EH1_AD1_Pu241[0];

    EH1_AD2_Total = x[0]*DayaBayFs[i][0]*EH1_AD2_U235[0];
    EH1_AD2_Total += x[1]*DayaBayFs[i][1]*EH1_AD2_U238[0];
    EH1_AD2_Total += x[2]*DayaBayFs[i][2]*EH1_AD2_Pu239[0];
    EH1_AD2_Total += x[3]*DayaBayFs[i][3]*EH1_AD2_Pu241[0];

    EH2_AD3_Total = x[0]*DayaBayFs[i][0]*EH2_AD3_U235[0];
    EH2_AD3_Total += x[1]*DayaBayFs[i][1]*EH2_AD3_U238[0];
    EH2_AD3_Total += x[2]*DayaBayFs[i][2]*EH2_AD3_Pu239[0];
    EH2_AD3_Total += x[3]*DayaBayFs[i][3]*EH2_AD3_Pu241[0];

    EH2_AD8_Total = x[0]*DayaBayFs[i][0]*EH2_AD8_U235[0];
    EH2_AD8_Total += x[1]*DayaBayFs[i][1]*EH2_AD8_U238[0];
    EH2_AD8_Total += x[2]*DayaBayFs[i][2]*EH2_AD8_Pu239[0];
    EH2_AD8_Total += x[3]*DayaBayFs[i][3]*EH2_AD8_Pu241[0];

    EH1_AD1_Total_0 = DayaBayFs[i][0]*EH1_AD1_U235_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][1]*EH1_AD1_U238_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][2]*EH1_AD1_Pu239_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][3]*EH1_AD1_Pu241_0[0];

    EH1_AD2_Total_0 = DayaBayFs[i][0]*EH1_AD2_U235_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][1]*EH1_AD2_U238_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][2]*EH1_AD2_Pu239_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][3]*EH1_AD2_Pu241_0[0];

    EH2_AD3_Total_0 = DayaBayFs[i][0]*EH2_AD3_U235_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][1]*EH2_AD3_U238_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][2]*EH2_AD3_Pu239_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][3]*EH2_AD3_Pu241_0[0];

    EH2_AD8_Total_0 = DayaBayFs[i][0]*EH2_AD8_U235_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][1]*EH2_AD8_U238_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][2]*EH2_AD8_Pu239_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][3]*EH2_AD8_Pu241_0[0];

    num = EH1_AD1_Total+ EH1_AD2_Total+ EH2_AD3_Total+ EH2_AD8_Total;
    denom = EH1_AD1_Total_0+ EH1_AD2_Total_0+ EH2_AD3_Total_0+ EH2_AD8_Total_0;

    delta[i+24] = num/denom - DBdata[i];
  }

  if (YesNo[25] == 1){
    RENO_Total = x[0]*RENOFs[i][0]*RENO_U235[0];
    RENO_Total += x[1]*RENOFs[i][1]*RENO_U238[0];
    RENO_Total += x[2]*RENOFs[i][2]*RENO_Pu239[0];
    RENO_Total += x[3]*RENOFs[i][3]*RENO_Pu241[0];

    RENO_Total_0 = RENOFs[i][0]*RENO_U235_0[0];
    RENO_Total_0 += RENOFs[i][1]*RENO_U238_0[0];
    RENO_Total_0 += RENOFs[i][2]*RENO_Pu239_0[0];
    RENO_Total_0 += RENOFs[i][3]*RENO_Pu241_0[0];

    delta[i+32] = RENO_Total/RENO_Total_0 - RENOdata[i] ;
  }
  }

/********************************
*  Putting the pieces together  *
*********************************/

  /* Adding the contributions from each experiment */

  for (i=0; i<24; i++){
    for (j=0; j<24; j++){
      chi2 += delta[i]*cov_matrix[i][j]*delta[j];
    }
  }

  /* Now, adding Daya Bay... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+24]*cov_matrix_DB[i][j]*delta[j+24];
    }
  }

  /* Now, adding RENO... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+32]*cov_matrix_RENO[i][j]*delta[j+32];
    }
  }

  /* Adding in the systematic uncertainties on the ratios... */

  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      chi2 += (x[i]-1.0)*(x[j]-1.0)*V_inv_HM[i][j]/(errors[i]*errors[j]);
    }
  }

  for (i = 0; i<4; i++){ /* Export nuisance parameters to main file */
    systematic[i] = x[i];
  }

  return chi2;
}  

/********************************************************************
*                                                                   *
*        THE CALCULATION OF A CHI-SQUARED WITHOUT SYSTEMATICS       *
*                                                                   *
*********************************************************************/

double combo_rate_chi_nosys(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following, for instance, 1703.00860

  "delta" is given by (Rexp - Rpred)
*/

  int k;

  double delta[40] = {0.0};

  int YesNo[26];
  for (k=0; k<26; k++){
    YesNo[k] = ((int *)user_data)[k];
  }

/************
*  Bugey-4  *
*************/

  double *Bugey_4_U235 = glbGetSignalFitRatePtr(exp, 0);
  double *Bugey_4_U238 = glbGetSignalFitRatePtr(exp, 1);
  double *Bugey_4_Pu239 = glbGetSignalFitRatePtr(exp, 2);
  double *Bugey_4_Pu241 = glbGetSignalFitRatePtr(exp, 3);

  double *Bugey_4_U235_0 = glbGetRuleRatePtr(exp, 0);
  double *Bugey_4_U238_0 = glbGetRuleRatePtr(exp, 1);
  double *Bugey_4_Pu239_0 = glbGetRuleRatePtr(exp, 2);
  double *Bugey_4_Pu241_0 = glbGetRuleRatePtr(exp, 3);
   
  double Bugey_4_RatioN = 0.0;
  double Bugey_4_RatioD = 0.0;

if(YesNo[0]==1||YesNo[2]==1){
  Bugey_4_RatioN = BugeyF[0]*Bugey_4_U235[0];
  Bugey_4_RatioN += BugeyF[1]*Bugey_4_U238[0];
  Bugey_4_RatioN += BugeyF[2]*Bugey_4_Pu239[0];
  Bugey_4_RatioN += BugeyF[3]*Bugey_4_Pu241[0];

  Bugey_4_RatioD = BugeyF[0]*Bugey_4_U235_0[0];
  Bugey_4_RatioD += BugeyF[1]*Bugey_4_U238_0[0];
  Bugey_4_RatioD += BugeyF[2]*Bugey_4_Pu239_0[0];
  Bugey_4_RatioD += BugeyF[3]*Bugey_4_Pu241_0[0];
}
if(YesNo[0]==1){
  delta[0] = Rexp[0] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/************
*  Rovno91  *
*************/

  double *Rovno91_U235 = glbGetSignalFitRatePtr(exp+1, 0);
  double *Rovno91_U238 = glbGetSignalFitRatePtr(exp+1, 1);
  double *Rovno91_Pu239 = glbGetSignalFitRatePtr(exp+1, 2);
  double *Rovno91_Pu241 = glbGetSignalFitRatePtr(exp+1, 3);

  double *Rovno91_U235_0 = glbGetRuleRatePtr(exp+1, 0);
  double *Rovno91_U238_0 = glbGetRuleRatePtr(exp+1, 1);
  double *Rovno91_Pu239_0 = glbGetRuleRatePtr(exp+1, 2);
  double *Rovno91_Pu241_0 = glbGetRuleRatePtr(exp+1, 3);

  double Rovno91_RatioN = 0.0;
  double Rovno91_RatioD = 0.0;

if(YesNo[1]==1){
  Rovno91_RatioN = Rovno_91F[0]*Rovno91_U235[0];
  Rovno91_RatioN += Rovno_91F[1]*Rovno91_U238[0];
  Rovno91_RatioN += Rovno_91F[2]*Rovno91_Pu239[0];
  Rovno91_RatioN += Rovno_91F[3]*Rovno91_Pu241[0];

  Rovno91_RatioD = Rovno_91F[0]*Rovno91_U235_0[0];
  Rovno91_RatioD += Rovno_91F[1]*Rovno91_U238_0[0];
  Rovno91_RatioD += Rovno_91F[2]*Rovno91_Pu239_0[0];
  Rovno91_RatioD += Rovno_91F[3]*Rovno91_Pu241_0[0];

  delta[1] = Rexp[1] - Rovno91_RatioN/Rovno91_RatioD;
}

/*******************
*  Bugey-3 (15 m)  *
********************/

/*
  No need to recalculate rates -- reuse the output for Bugey-4!
*/

if(YesNo[2]==1){
  delta[2] = Rexp[2] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/*******************
*  Bugey-3 (40 m)  *
********************/

  double *Bugey_3_40_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *Bugey_3_40_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);
   
  double *Bugey_3_40_U235_0 = glbGetRuleRatePtr(exp+2, 0);
  double *Bugey_3_40_U238_0 = glbGetRuleRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239_0 = glbGetRuleRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241_0 = glbGetRuleRatePtr(exp+2, 3);
   
  double Bugey_3_40_RatioN = 0.0;
  double Bugey_3_40_RatioD = 0.0;

if(YesNo[3]==1){   
  Bugey_3_40_RatioN = BugeyF[0]*Bugey_3_40_U235[0];
  Bugey_3_40_RatioN += BugeyF[1]*Bugey_3_40_U238[0];
  Bugey_3_40_RatioN += BugeyF[2]*Bugey_3_40_Pu239[0];
  Bugey_3_40_RatioN += BugeyF[3]*Bugey_3_40_Pu241[0];

  Bugey_3_40_RatioD = BugeyF[0]*Bugey_3_40_U235_0[0];
  Bugey_3_40_RatioD += BugeyF[1]*Bugey_3_40_U238_0[0];
  Bugey_3_40_RatioD += BugeyF[2]*Bugey_3_40_Pu239_0[0];
  Bugey_3_40_RatioD += BugeyF[3]*Bugey_3_40_Pu241_0[0];

  delta[3] = Rexp[3] - Bugey_3_40_RatioN/Bugey_3_40_RatioD;
}

/*******************
*  Bugey-3 (95 m)  *
********************/

  double *Bugey_3_95_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *Bugey_3_95_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);
   
  double *Bugey_3_95_U235_0 = glbGetRuleRatePtr(exp+3, 0);
  double *Bugey_3_95_U238_0 = glbGetRuleRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239_0 = glbGetRuleRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241_0 = glbGetRuleRatePtr(exp+3, 3);
   
  double Bugey_3_95_RatioN = 0.0;
  double Bugey_3_95_RatioD = 0.0;

if(YesNo[4]==1){   
  Bugey_3_95_RatioN = BugeyF[0]*Bugey_3_95_U235[0];
  Bugey_3_95_RatioN += BugeyF[1]*Bugey_3_95_U238[0];
  Bugey_3_95_RatioN += BugeyF[2]*Bugey_3_95_Pu239[0];
  Bugey_3_95_RatioN += BugeyF[3]*Bugey_3_95_Pu241[0];

  Bugey_3_95_RatioD = BugeyF[0]*Bugey_3_95_U235_0[0];
  Bugey_3_95_RatioD += BugeyF[1]*Bugey_3_95_U238_0[0];
  Bugey_3_95_RatioD += BugeyF[2]*Bugey_3_95_Pu239_0[0];
  Bugey_3_95_RatioD += BugeyF[3]*Bugey_3_95_Pu241_0[0];

  delta[4] = Rexp[4] - Bugey_3_95_RatioN/Bugey_3_95_RatioD;
}

/********************
*  Gosgen (37.9 m)  *
*********************/

  double *Gosgen_38_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *Gosgen_38_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
  double *Gosgen_38_U235_0 = glbGetRuleRatePtr(exp+4, 0);
  double *Gosgen_38_U238_0 = glbGetRuleRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239_0 = glbGetRuleRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241_0 = glbGetRuleRatePtr(exp+4, 3);
   
  double Gosgen_38_RatioN = 0.0;
  double Gosgen_38_RatioD = 0.0;

if(YesNo[5]==1){
   
  Gosgen_38_RatioN = Gosgen38F[0]*Gosgen_38_U235[0];
  Gosgen_38_RatioN += Gosgen38F[1]*Gosgen_38_U238[0];
  Gosgen_38_RatioN += Gosgen38F[2]*Gosgen_38_Pu239[0];
  Gosgen_38_RatioN += Gosgen38F[3]*Gosgen_38_Pu241[0];

  Gosgen_38_RatioD = Gosgen38F[0]*Gosgen_38_U235_0[0];
  Gosgen_38_RatioD += Gosgen38F[1]*Gosgen_38_U238_0[0];
  Gosgen_38_RatioD += Gosgen38F[2]*Gosgen_38_Pu239_0[0];
  Gosgen_38_RatioD += Gosgen38F[3]*Gosgen_38_Pu241_0[0];

  delta[5] = Rexp[5] - Gosgen_38_RatioN/Gosgen_38_RatioD;
}

/********************
*  Gosgen (45.9 m)  *
*********************/

  double *Gosgen_46_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *Gosgen_46_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);
   
  double *Gosgen_46_U235_0 = glbGetRuleRatePtr(exp+5, 0);
  double *Gosgen_46_U238_0 = glbGetRuleRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239_0 = glbGetRuleRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241_0 = glbGetRuleRatePtr(exp+5, 3);
   
  double Gosgen_46_RatioN = 0.0;
  double Gosgen_46_RatioD = 0.0;

if(YesNo[6]==1){
   
  Gosgen_46_RatioN = Gosgen46F[0]*Gosgen_46_U235[0];
  Gosgen_46_RatioN += Gosgen46F[1]*Gosgen_46_U238[0];
  Gosgen_46_RatioN += Gosgen46F[2]*Gosgen_46_Pu239[0];
  Gosgen_46_RatioN += Gosgen46F[3]*Gosgen_46_Pu241[0];

  Gosgen_46_RatioD = Gosgen46F[0]*Gosgen_46_U235_0[0];
  Gosgen_46_RatioD += Gosgen46F[1]*Gosgen_46_U238_0[0];
  Gosgen_46_RatioD += Gosgen46F[2]*Gosgen_46_Pu239_0[0];
  Gosgen_46_RatioD += Gosgen46F[3]*Gosgen_46_Pu241_0[0];

  delta[6] = Rexp[6] - Gosgen_46_RatioN/Gosgen_46_RatioD;
}

/********************
*  Gosgen (64.7 m)  *
*********************/

  double *Gosgen_65_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *Gosgen_65_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);

  double *Gosgen_65_U235_0 = glbGetRuleRatePtr(exp+6, 0);
  double *Gosgen_65_U238_0 = glbGetRuleRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239_0 = glbGetRuleRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241_0 = glbGetRuleRatePtr(exp+6, 3);
   
  double Gosgen_65_RatioN = 0.0;
  double Gosgen_65_RatioD = 0.0;

if(YesNo[7]==1){
   
  Gosgen_65_RatioN = Gosgen65F[0]*Gosgen_65_U235[0];
  Gosgen_65_RatioN += Gosgen65F[1]*Gosgen_65_U238[0];
  Gosgen_65_RatioN += Gosgen65F[2]*Gosgen_65_Pu239[0];
  Gosgen_65_RatioN += Gosgen65F[3]*Gosgen_65_Pu241[0];

  Gosgen_65_RatioD = Gosgen65F[0]*Gosgen_65_U235_0[0];
  Gosgen_65_RatioD += Gosgen65F[1]*Gosgen_65_U238_0[0];
  Gosgen_65_RatioD += Gosgen65F[2]*Gosgen_65_Pu239_0[0];
  Gosgen_65_RatioD += Gosgen65F[3]*Gosgen_65_Pu241_0[0];

  delta[7] = Rexp[7] - Gosgen_65_RatioN/Gosgen_65_RatioD;
}

/*****************
*  ILL (8.76 m)  *
******************/

  double *ILL_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *ILL_U235_0 = glbGetRuleRatePtr(exp+7, 0);

if(YesNo[8]==1){
  delta[8] = Rexp[8] - ILL_U235[0]/ILL_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (32.8 m)  *
****************************/

  double *Krasnoyarsk_33_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *Krasnoyarsk_33_U235_0 = glbGetRuleRatePtr(exp+8, 0);

if(YesNo[9]==1){
  delta[9] = Rexp[9] - Krasnoyarsk_33_U235[0]/Krasnoyarsk_33_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (92.3 m)  *
****************************/

  double *Krasnoyarsk_92_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *Krasnoyarsk_92_U235_0 = glbGetRuleRatePtr(exp+9, 0);

if(YesNo[10]==1){
  delta[10] = Rexp[10] - Krasnoyarsk_92_U235[0]/Krasnoyarsk_92_U235_0[0];
}

/***************************
*  Krasnoyarsk94 (57.0 m)  *
****************************/

  double *Krasnoyarsk_57_U235 = glbGetSignalFitRatePtr(exp+10, 0);
  double *Krasnoyarsk_57_U235_0 = glbGetRuleRatePtr(exp+10, 0);

if(YesNo[11]==1){
  delta[11] = Rexp[11] - Krasnoyarsk_57_U235[0]/Krasnoyarsk_57_U235_0[0];
}

/***************************
*  Krasnoyarsk99 (34.0 m)  *
****************************/

  double *Krasnoyarsk_34_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *Krasnoyarsk_34_U235_0 = glbGetRuleRatePtr(exp+11, 0);

if(YesNo[12]==1){
  delta[12] = Rexp[12] - Krasnoyarsk_34_U235[0]/Krasnoyarsk_34_U235_0[0];
}

/********************
*  SRP-18 (18.2 m)  *
*********************/

  double *SRP_18_U235 = glbGetSignalFitRatePtr(exp+12, 0);
  double *SRP_18_U235_0 = glbGetRuleRatePtr(exp+12, 0);

if(YesNo[13]==1){
  delta[13] = Rexp[13] - SRP_18_U235[0]/SRP_18_U235_0[0];
}

/********************
*  SRP-24 (23.8 m)  *
*********************/

  double *SRP_24_U235 = glbGetSignalFitRatePtr(exp+13, 0);
  double *SRP_24_U235_0 = glbGetRuleRatePtr(exp+13, 0);

if(YesNo[14]==1){
  delta[14] = Rexp[14] - SRP_24_U235[0]/SRP_24_U235_0[0];
}

/***************
*  Rovno88-1I  *
****************/

  double *Rovno88_1I_U235 = glbGetSignalFitRatePtr(exp+14, 0);
  double *Rovno88_1I_U238 = glbGetSignalFitRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239 = glbGetSignalFitRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241 = glbGetSignalFitRatePtr(exp+14, 3);

  double *Rovno88_1I_U235_0 = glbGetRuleRatePtr(exp+14, 0);
  double *Rovno88_1I_U238_0 = glbGetRuleRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239_0 = glbGetRuleRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241_0 = glbGetRuleRatePtr(exp+14, 3);

  double Rovno88_1I_RatioN = 0.0;
  double Rovno88_1I_RatioD = 0.0;

if(YesNo[15]==1){
  Rovno88_1I_RatioN = Rovno88_1IF[0]*Rovno88_1I_U235[0];
  Rovno88_1I_RatioN += Rovno88_1IF[1]*Rovno88_1I_U238[0];
  Rovno88_1I_RatioN += Rovno88_1IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1I_RatioN += Rovno88_1IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1I_RatioD = Rovno88_1IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[3]*Rovno88_1I_Pu241_0[0];

  delta[15] = Rexp[15] - Rovno88_1I_RatioN/Rovno88_1I_RatioD;
}

/***************
*  Rovno88-2I  *
****************/

  double Rovno88_2I_RatioN = 0.0;
  double Rovno88_2I_RatioD = 0.0;

if(YesNo[16]==1){   
  Rovno88_2I_RatioN = Rovno88_2IF[0]*Rovno88_1I_U235[0];
  Rovno88_2I_RatioN += Rovno88_2IF[1]*Rovno88_1I_U238[0];
  Rovno88_2I_RatioN += Rovno88_2IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_2I_RatioN += Rovno88_2IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_2I_RatioD = Rovno88_2IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[3]*Rovno88_1I_Pu241_0[0];

  delta[16] = Rexp[16] - Rovno88_2I_RatioN/Rovno88_2I_RatioD;
}
/***************
*  Rovno88-1S  *
****************/

  double Rovno88_1S_RatioN = 0.0;
  double Rovno88_1S_RatioD = 0.0;

if(YesNo[17]==1){
  Rovno88_1S_RatioN = Rovno88_1SF[0]*Rovno88_1I_U235[0];
  Rovno88_1S_RatioN += Rovno88_1SF[1]*Rovno88_1I_U238[0];
  Rovno88_1S_RatioN += Rovno88_1SF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1S_RatioN += Rovno88_1SF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1S_RatioD = Rovno88_1SF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[3]*Rovno88_1I_Pu241_0[0];

  delta[17] = Rexp[17] - Rovno88_1S_RatioN/Rovno88_1S_RatioD;
}

/**********************
*  Rovno88-2S (25 m)  *
***********************/
   
  double Rovno88_2S_25_RatioN = 0.0;
  double Rovno88_2S_25_RatioD = 0.0;

  double *Rovno88_2S_25_U235 = glbGetSignalFitRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238 = glbGetSignalFitRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239 = glbGetSignalFitRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241 = glbGetSignalFitRatePtr(exp+15, 3);

  double *Rovno88_2S_25_U235_0 = glbGetRuleRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238_0 = glbGetRuleRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239_0 = glbGetRuleRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241_0 = glbGetRuleRatePtr(exp+15, 3);

if(YesNo[18]==1){
  Rovno88_2S_25_RatioN = Rovno88_2S_25F[0]*Rovno88_2S_25_U235[0];
  Rovno88_2S_25_RatioN += Rovno88_2S_25F[1]*Rovno88_2S_25_U238[0];
  Rovno88_2S_25_RatioN += Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239[0];
  Rovno88_2S_25_RatioN += Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241[0];

  Rovno88_2S_25_RatioD = Rovno88_2S_25F[0]*Rovno88_2S_25_U235_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[1]*Rovno88_2S_25_U238_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241_0[0];

  delta[18] = Rexp[18] - Rovno88_2S_25_RatioN/Rovno88_2S_25_RatioD;
}

/**********************
*  Rovno88-3S (18 m)  *
***********************/

  double Rovno88_2S_18_RatioN = 0.0;
  double Rovno88_2S_18_RatioD = 0.0;

if(YesNo[19]==1){
  Rovno88_2S_18_RatioN = Rovno88_2S_18F[0]*Rovno88_1I_U235[0];
  Rovno88_2S_18_RatioN += Rovno88_2S_18F[1]*Rovno88_1I_U238[0];
  Rovno88_2S_18_RatioN += Rovno88_2S_18F[2]*Rovno88_1I_Pu239[0];
  Rovno88_2S_18_RatioN += Rovno88_2S_18F[3]*Rovno88_1I_Pu241[0];

  Rovno88_2S_18_RatioD = Rovno88_2S_18F[0]*Rovno88_1I_U235_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[1]*Rovno88_1I_U238_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[3]*Rovno88_1I_Pu241_0[0];

  delta[19] = Rexp[19] - Rovno88_2S_18_RatioN/Rovno88_2S_18_RatioD;
}

/************
*  Nucifer  *
*************/

  double *Nucifer_U235 = glbGetSignalFitRatePtr(exp+16, 0);
  double *Nucifer_U238 = glbGetSignalFitRatePtr(exp+16, 1);
  double *Nucifer_Pu239 = glbGetSignalFitRatePtr(exp+16, 2);
  double *Nucifer_Pu241 = glbGetSignalFitRatePtr(exp+16, 3);

  double *Nucifer_U235_0 = glbGetRuleRatePtr(exp+16, 0);
  double *Nucifer_U238_0 = glbGetRuleRatePtr(exp+16, 1);
  double *Nucifer_Pu239_0 = glbGetRuleRatePtr(exp+16, 2);
  double *Nucifer_Pu241_0 = glbGetRuleRatePtr(exp+16, 3);

  double Nucifer_RatioN=0.0;
  double Nucifer_RatioD=0.0;

if(YesNo[20]==1){
  Nucifer_RatioN = NuciferF[0]*Nucifer_U235[0];
  Nucifer_RatioN += NuciferF[1]*Nucifer_U238[0];
  Nucifer_RatioN += NuciferF[2]*Nucifer_Pu239[0];
  Nucifer_RatioN += NuciferF[3]*Nucifer_Pu241[0];

  Nucifer_RatioD = NuciferF[0]*Nucifer_U235_0[0];
  Nucifer_RatioD += NuciferF[1]*Nucifer_U238_0[0];
  Nucifer_RatioD += NuciferF[2]*Nucifer_Pu239_0[0];
  Nucifer_RatioD += NuciferF[3]*Nucifer_Pu241_0[0];

  delta[20] = Rexp[20] - Nucifer_RatioN/Nucifer_RatioD;
}

/**********************************
*  Palo Verde -- 750 m and 890 m  *
***********************************/

  double PV_RatioN=0.0;
  double PV_RatioD=0.0;

  double *PV_750_U235 = glbGetSignalFitRatePtr(exp+17, 0);
  double *PV_750_U238 = glbGetSignalFitRatePtr(exp+17, 1);
  double *PV_750_Pu239 = glbGetSignalFitRatePtr(exp+17, 2);
  double *PV_750_Pu241 = glbGetSignalFitRatePtr(exp+17, 3);

  double *PV_750_U235_0 = glbGetRuleRatePtr(exp+17, 0);
  double *PV_750_U238_0 = glbGetRuleRatePtr(exp+17, 1);
  double *PV_750_Pu239_0 = glbGetRuleRatePtr(exp+17, 2);
  double *PV_750_Pu241_0 = glbGetRuleRatePtr(exp+17, 3);

  double *PV_890_U235 = glbGetSignalFitRatePtr(exp+18, 0);
  double *PV_890_U238 = glbGetSignalFitRatePtr(exp+18, 1);
  double *PV_890_Pu239 = glbGetSignalFitRatePtr(exp+18, 2);
  double *PV_890_Pu241 = glbGetSignalFitRatePtr(exp+18, 3);

  double *PV_890_U235_0 = glbGetRuleRatePtr(exp+18, 0);
  double *PV_890_U238_0 = glbGetRuleRatePtr(exp+18, 1);
  double *PV_890_Pu239_0 = glbGetRuleRatePtr(exp+18, 2);
  double *PV_890_Pu241_0 = glbGetRuleRatePtr(exp+18, 3);

if(YesNo[21]==1){
  PV_RatioN = PVF[0]*(PV_750_U235[0]+PV_890_U235[0]);
  PV_RatioN += PVF[1]*(PV_750_U238[0]+ PV_890_U238[0]);
  PV_RatioN += PVF[2]*(PV_750_Pu239[0]+ PV_890_Pu239[0]);
  PV_RatioN += PVF[3]*(PV_750_Pu241[0]+ PV_890_Pu241[0]);

  PV_RatioD = PVF[0]*(PV_750_U235_0[0]+PV_890_U235_0[0]);
  PV_RatioD += PVF[1]*(PV_750_U238_0[0]+ PV_890_U238_0[0]);
  PV_RatioD += PVF[2]*(PV_750_Pu239_0[0]+ PV_890_Pu239_0[0]);
  PV_RatioD += PVF[3]*(PV_750_Pu241_0[0]+ PV_890_Pu241_0[0]);

  delta[21] = Rexp[21] - PV_RatioN/PV_RatioD;
}

/************************************
*  Double Chooz -- 355 m and 469 m  *
*************************************/

  double DC_RatioN=0.0;
  double DC_RatioD=0.0;

  double *DC_355_U235 = glbGetSignalFitRatePtr(exp+19, 0);
  double *DC_355_U238 = glbGetSignalFitRatePtr(exp+19, 1);
  double *DC_355_Pu239 = glbGetSignalFitRatePtr(exp+19, 2);
  double *DC_355_Pu241 = glbGetSignalFitRatePtr(exp+19, 3);

  double *DC_355_U235_0 = glbGetRuleRatePtr(exp+19, 0);
  double *DC_355_U238_0 = glbGetRuleRatePtr(exp+19, 1);
  double *DC_355_Pu239_0 = glbGetRuleRatePtr(exp+19, 2);
  double *DC_355_Pu241_0 = glbGetRuleRatePtr(exp+19, 3);

  double *DC_469_U235 = glbGetSignalFitRatePtr(exp+20, 0);
  double *DC_469_U238 = glbGetSignalFitRatePtr(exp+20, 1);
  double *DC_469_Pu239 = glbGetSignalFitRatePtr(exp+20, 2);
  double *DC_469_Pu241 = glbGetSignalFitRatePtr(exp+20, 3);

  double *DC_469_U235_0 = glbGetRuleRatePtr(exp+20, 0);
  double *DC_469_U238_0 = glbGetRuleRatePtr(exp+20, 1);
  double *DC_469_Pu239_0 = glbGetRuleRatePtr(exp+20, 2);
  double *DC_469_Pu241_0 = glbGetRuleRatePtr(exp+20, 3);

if(YesNo[22]==1){
  DC_RatioN = DCF[0]*(DC_355_U235[0] + DC_469_U235[0]);
  DC_RatioN += DCF[1]*(DC_355_U238[0] + DC_469_U238[0]);
  DC_RatioN += DCF[2]*(DC_355_Pu239[0] + DC_469_Pu239[0]);
  DC_RatioN += DCF[3]*(DC_355_Pu241[0] + DC_469_Pu241[0]);

  DC_RatioD = DCF[0]*(DC_355_U235_0[0] + DC_469_U235_0[0]);
  DC_RatioD += DCF[1]*(DC_355_U238_0[0] + DC_469_U238_0[0]);
  DC_RatioD += DCF[2]*(DC_355_Pu239_0[0] + DC_469_Pu239_0[0]);
  DC_RatioD += DCF[3]*(DC_355_Pu241_0[0] + DC_469_Pu241_0[0]);

  delta[22] = Rexp[22] - DC_RatioN/DC_RatioD;
}

/******************************
*  Chooz -- 998 m and 1115 m  *
*******************************/

  double Chooz_RatioN=0.0;
  double Chooz_RatioD=0.0;

  double *Chooz_998_U235 = glbGetSignalFitRatePtr(exp+21, 0);
  double *Chooz_998_U238 = glbGetSignalFitRatePtr(exp+21, 1);
  double *Chooz_998_Pu239 = glbGetSignalFitRatePtr(exp+21, 2);
  double *Chooz_998_Pu241 = glbGetSignalFitRatePtr(exp+21, 3);

  double *Chooz_998_U235_0 = glbGetRuleRatePtr(exp+21, 0);
  double *Chooz_998_U238_0 = glbGetRuleRatePtr(exp+21, 1);
  double *Chooz_998_Pu239_0 = glbGetRuleRatePtr(exp+21, 2);
  double *Chooz_998_Pu241_0 = glbGetRuleRatePtr(exp+21, 3);

  double *Chooz_1115_U235 = glbGetSignalFitRatePtr(exp+22, 0);
  double *Chooz_1115_U238 = glbGetSignalFitRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239 = glbGetSignalFitRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241 = glbGetSignalFitRatePtr(exp+22, 3);

  double *Chooz_1115_U235_0 = glbGetRuleRatePtr(exp+22, 0);
  double *Chooz_1115_U238_0 = glbGetRuleRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239_0 = glbGetRuleRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241_0 = glbGetRuleRatePtr(exp+22, 3);

if(YesNo[23]==1){
  Chooz_RatioN = ChoozF[0]*(Chooz_998_U235[0] + Chooz_1115_U235[0]);
  Chooz_RatioN += ChoozF[1]*(Chooz_998_U238[0] + Chooz_1115_U238[0]);
  Chooz_RatioN += ChoozF[2]*(Chooz_998_Pu239[0] + Chooz_1115_Pu239[0]);
  Chooz_RatioN += ChoozF[3]*(Chooz_998_Pu241[0] + Chooz_1115_Pu241[0]);

  Chooz_RatioD = ChoozF[0]*(Chooz_998_U235_0[0] + Chooz_1115_U235_0[0]);
  Chooz_RatioD += ChoozF[1]*(Chooz_998_U238_0[0] + Chooz_1115_U238_0[0]);
  Chooz_RatioD += ChoozF[2]*(Chooz_998_Pu239_0[0] + Chooz_1115_Pu239_0[0]);
  Chooz_RatioD += ChoozF[3]*(Chooz_998_Pu241_0[0] + Chooz_1115_Pu241_0[0]);

  delta[23] = Rexp[23] - Chooz_RatioN/Chooz_RatioD;
}

/************
*  EH1 AD1  *
*************/

  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+23, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+23, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+23, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+23, 3);
   
  double EH1_AD1_Total = 0.0;
  double EH1_AD1_Total_0 = 0.0;

/************
*  EH1 AD2  *
*************/

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+24, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+24, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+24, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+24, 3);

  double EH1_AD2_Total = 0.0;
  double EH1_AD2_Total_0 = 0.0;

/************
*  EH2 AD3  *
*************/

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+25, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+25, 3);

  double *EH2_AD3_U235_0 = glbGetRuleRatePtr(exp+25, 0);
  double *EH2_AD3_U238_0 = glbGetRuleRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239_0 = glbGetRuleRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241_0 = glbGetRuleRatePtr(exp+25, 3);
   
  double EH2_AD3_Total = 0.0;
  double EH2_AD3_Total_0 = 0.0;

/************
*  EH2 AD8  *
*************/

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+26, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+26, 3);

  double *EH2_AD8_U235_0 = glbGetRuleRatePtr(exp+26, 0);
  double *EH2_AD8_U238_0 = glbGetRuleRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239_0 = glbGetRuleRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241_0 = glbGetRuleRatePtr(exp+26, 3);

  double EH2_AD8_Total = 0.0;
  double EH2_AD8_Total_0 = 0.0;

/*********
*  RENO  *
**********/

  double *RENO_U235 = glbGetSignalFitRatePtr(exp+27, 0);
  double *RENO_U238 = glbGetSignalFitRatePtr(exp+27, 1);
  double *RENO_Pu239 = glbGetSignalFitRatePtr(exp+27, 2);
  double *RENO_Pu241 = glbGetSignalFitRatePtr(exp+27, 3);

  double *RENO_U235_0 = glbGetRuleRatePtr(exp+27, 0);
  double *RENO_U238_0 = glbGetRuleRatePtr(exp+27, 1);
  double *RENO_Pu239_0 = glbGetRuleRatePtr(exp+27, 2);
  double *RENO_Pu241_0 = glbGetRuleRatePtr(exp+27, 3);

  double RENO_Total = 0.0;
  double RENO_Total_0 = 0.0;

/*********************************
*  Assembling Daya Bay and RENO  *
**********************************/

  double chi2 = 0.0;
  int i,j;

  double num, denom;

  for (i=0; i<8; i++){
  if (YesNo[24] == 1){
    EH1_AD1_Total = DayaBayFs[i][0]*EH1_AD1_U235[0];
    EH1_AD1_Total += DayaBayFs[i][1]*EH1_AD1_U238[0];
    EH1_AD1_Total += DayaBayFs[i][2]*EH1_AD1_Pu239[0];
    EH1_AD1_Total += DayaBayFs[i][3]*EH1_AD1_Pu241[0];

    EH1_AD2_Total = DayaBayFs[i][0]*EH1_AD2_U235[0];
    EH1_AD2_Total += DayaBayFs[i][1]*EH1_AD2_U238[0];
    EH1_AD2_Total += DayaBayFs[i][2]*EH1_AD2_Pu239[0];
    EH1_AD2_Total += DayaBayFs[i][3]*EH1_AD2_Pu241[0];

    EH2_AD3_Total = DayaBayFs[i][0]*EH2_AD3_U235[0];
    EH2_AD3_Total += DayaBayFs[i][1]*EH2_AD3_U238[0];
    EH2_AD3_Total += DayaBayFs[i][2]*EH2_AD3_Pu239[0];
    EH2_AD3_Total += DayaBayFs[i][3]*EH2_AD3_Pu241[0];

    EH2_AD8_Total = DayaBayFs[i][0]*EH2_AD8_U235[0];
    EH2_AD8_Total += DayaBayFs[i][1]*EH2_AD8_U238[0];
    EH2_AD8_Total += DayaBayFs[i][2]*EH2_AD8_Pu239[0];
    EH2_AD8_Total += DayaBayFs[i][3]*EH2_AD8_Pu241[0];

    EH1_AD1_Total_0 = DayaBayFs[i][0]*EH1_AD1_U235_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][1]*EH1_AD1_U238_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][2]*EH1_AD1_Pu239_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][3]*EH1_AD1_Pu241_0[0];

    EH1_AD2_Total_0 = DayaBayFs[i][0]*EH1_AD2_U235_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][1]*EH1_AD2_U238_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][2]*EH1_AD2_Pu239_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][3]*EH1_AD2_Pu241_0[0];

    EH2_AD3_Total_0 = DayaBayFs[i][0]*EH2_AD3_U235_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][1]*EH2_AD3_U238_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][2]*EH2_AD3_Pu239_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][3]*EH2_AD3_Pu241_0[0];

    EH2_AD8_Total_0 = DayaBayFs[i][0]*EH2_AD8_U235_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][1]*EH2_AD8_U238_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][2]*EH2_AD8_Pu239_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][3]*EH2_AD8_Pu241_0[0];

    num = EH1_AD1_Total+ EH1_AD2_Total+ EH2_AD3_Total+ EH2_AD8_Total;
    denom = EH1_AD1_Total_0+ EH1_AD2_Total_0+ EH2_AD3_Total_0+ EH2_AD8_Total_0;

    delta[i+24] = num/denom - DBdata[i];
  }

  if (YesNo[25] == 1){
    RENO_Total = RENOFs[i][0]*RENO_U235[0];
    RENO_Total += RENOFs[i][1]*RENO_U238[0];
    RENO_Total += RENOFs[i][2]*RENO_Pu239[0];
    RENO_Total += RENOFs[i][3]*RENO_Pu241[0];

    RENO_Total_0 = RENOFs[i][0]*RENO_U235_0[0];
    RENO_Total_0 += RENOFs[i][1]*RENO_U238_0[0];
    RENO_Total_0 += RENOFs[i][2]*RENO_Pu239_0[0];
    RENO_Total_0 += RENOFs[i][3]*RENO_Pu241_0[0];

    delta[i+32] = RENO_Total/RENO_Total_0 - RENOdata[i] ;
  }
  }

/********************************
*  Putting the pieces together  *
*********************************/

  /* Adding the contributions from each experiment */

  for (i=0; i<24; i++){
    for (j=0; j<24; j++){
      chi2 += delta[i]*cov_matrix[i][j]*delta[j];
    }
  }

  /* Now, adding Daya Bay... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+24]*cov_matrix_DB[i][j]*delta[j+24];
    }
  }

  /* Now, adding RENO... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+32]*cov_matrix_RENO[i][j]*delta[j+32];
    }
  }

  return chi2;
}

/********************************************************************
*                                                                   *
*    THE CALCULATION OF A CHI-SQUARED WITH HM SYSTEMATICS USING     *
*                SUMMATION METHOD ANTINEUTRINO FLUXES               *
*                                                                   *
*********************************************************************/

double combo_rate_chi_SM(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following, for instance, 1703.00860

  "delta" is given by (Rexp - Rpred)
*/

  int k;

  double delta[40] = {0.0};

  int YesNo[26];
  for (k=0; k<26; k++){
    YesNo[k] = ((int *)user_data)[k];
  }

/************
*  Bugey-4  *
*************/

  double *Bugey_4_U235 = glbGetSignalFitRatePtr(exp, 0);
  double *Bugey_4_U238 = glbGetSignalFitRatePtr(exp, 1);
  double *Bugey_4_Pu239 = glbGetSignalFitRatePtr(exp, 2);
  double *Bugey_4_Pu241 = glbGetSignalFitRatePtr(exp, 3);

  double *Bugey_4_U235_0 = glbGetRuleRatePtr(exp, 0);
  double *Bugey_4_U238_0 = glbGetRuleRatePtr(exp, 1);
  double *Bugey_4_Pu239_0 = glbGetRuleRatePtr(exp, 2);
  double *Bugey_4_Pu241_0 = glbGetRuleRatePtr(exp, 3);
   
  double Bugey_4_RatioN = 0.0;
  double Bugey_4_RatioD = 0.0;

if(YesNo[0]==1||YesNo[2]==1){
  Bugey_4_RatioN = x[0]*BugeyF[0]*Bugey_4_U235[0];
  Bugey_4_RatioN += x[1]*BugeyF[1]*Bugey_4_U238[0];
  Bugey_4_RatioN += x[2]*BugeyF[2]*Bugey_4_Pu239[0];
  Bugey_4_RatioN += x[3]*BugeyF[3]*Bugey_4_Pu241[0];

  Bugey_4_RatioD = BugeyF[0]*Bugey_4_U235_0[0];
  Bugey_4_RatioD += BugeyF[1]*Bugey_4_U238_0[0];
  Bugey_4_RatioD += BugeyF[2]*Bugey_4_Pu239_0[0];
  Bugey_4_RatioD += BugeyF[3]*Bugey_4_Pu241_0[0];
}
if(YesNo[0]==1){
  delta[0] = RexpSM[0] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/************
*  Rovno91  *
*************/

  double *Rovno91_U235 = glbGetSignalFitRatePtr(exp+1, 0);
  double *Rovno91_U238 = glbGetSignalFitRatePtr(exp+1, 1);
  double *Rovno91_Pu239 = glbGetSignalFitRatePtr(exp+1, 2);
  double *Rovno91_Pu241 = glbGetSignalFitRatePtr(exp+1, 3);

  double *Rovno91_U235_0 = glbGetRuleRatePtr(exp+1, 0);
  double *Rovno91_U238_0 = glbGetRuleRatePtr(exp+1, 1);
  double *Rovno91_Pu239_0 = glbGetRuleRatePtr(exp+1, 2);
  double *Rovno91_Pu241_0 = glbGetRuleRatePtr(exp+1, 3);

  double Rovno91_RatioN = 0.0;
  double Rovno91_RatioD = 0.0;

if(YesNo[1]==1){
  Rovno91_RatioN = x[0]*Rovno_91F[0]*Rovno91_U235[0];
  Rovno91_RatioN += x[1]*Rovno_91F[1]*Rovno91_U238[0];
  Rovno91_RatioN += x[2]*Rovno_91F[2]*Rovno91_Pu239[0];
  Rovno91_RatioN += x[3]*Rovno_91F[3]*Rovno91_Pu241[0];

  Rovno91_RatioD = Rovno_91F[0]*Rovno91_U235_0[0];
  Rovno91_RatioD += Rovno_91F[1]*Rovno91_U238_0[0];
  Rovno91_RatioD += Rovno_91F[2]*Rovno91_Pu239_0[0];
  Rovno91_RatioD += Rovno_91F[3]*Rovno91_Pu241_0[0];

  delta[1] = RexpSM[1] - Rovno91_RatioN/Rovno91_RatioD;
}

/*******************
*  Bugey-3 (15 m)  *
********************/

/*
  No need to recalculate rates -- reuse the output for Bugey-4!
*/

if(YesNo[2]==1){
  delta[2] = RexpSM[2] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/*******************
*  Bugey-3 (40 m)  *
********************/

  double *Bugey_3_40_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *Bugey_3_40_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);
   
  double *Bugey_3_40_U235_0 = glbGetRuleRatePtr(exp+2, 0);
  double *Bugey_3_40_U238_0 = glbGetRuleRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239_0 = glbGetRuleRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241_0 = glbGetRuleRatePtr(exp+2, 3);
   
  double Bugey_3_40_RatioN = 0.0;
  double Bugey_3_40_RatioD = 0.0;

if(YesNo[3]==1){   
  Bugey_3_40_RatioN = x[0]*BugeyF[0]*Bugey_3_40_U235[0];
  Bugey_3_40_RatioN += x[1]*BugeyF[1]*Bugey_3_40_U238[0];
  Bugey_3_40_RatioN += x[2]*BugeyF[2]*Bugey_3_40_Pu239[0];
  Bugey_3_40_RatioN += x[3]*BugeyF[3]*Bugey_3_40_Pu241[0];

  Bugey_3_40_RatioD = BugeyF[0]*Bugey_3_40_U235_0[0];
  Bugey_3_40_RatioD += BugeyF[1]*Bugey_3_40_U238_0[0];
  Bugey_3_40_RatioD += BugeyF[2]*Bugey_3_40_Pu239_0[0];
  Bugey_3_40_RatioD += BugeyF[3]*Bugey_3_40_Pu241_0[0];

  delta[3] = RexpSM[3] - Bugey_3_40_RatioN/Bugey_3_40_RatioD;
}

/*******************
*  Bugey-3 (95 m)  *
********************/

  double *Bugey_3_95_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *Bugey_3_95_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);
   
  double *Bugey_3_95_U235_0 = glbGetRuleRatePtr(exp+3, 0);
  double *Bugey_3_95_U238_0 = glbGetRuleRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239_0 = glbGetRuleRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241_0 = glbGetRuleRatePtr(exp+3, 3);
   
  double Bugey_3_95_RatioN = 0.0;
  double Bugey_3_95_RatioD = 0.0;

if(YesNo[4]==1){   
  Bugey_3_95_RatioN = x[0]*BugeyF[0]*Bugey_3_95_U235[0];
  Bugey_3_95_RatioN += x[1]*BugeyF[1]*Bugey_3_95_U238[0];
  Bugey_3_95_RatioN += x[2]*BugeyF[2]*Bugey_3_95_Pu239[0];
  Bugey_3_95_RatioN += x[3]*BugeyF[3]*Bugey_3_95_Pu241[0];

  Bugey_3_95_RatioD = BugeyF[0]*Bugey_3_95_U235_0[0];
  Bugey_3_95_RatioD += BugeyF[1]*Bugey_3_95_U238_0[0];
  Bugey_3_95_RatioD += BugeyF[2]*Bugey_3_95_Pu239_0[0];
  Bugey_3_95_RatioD += BugeyF[3]*Bugey_3_95_Pu241_0[0];

  delta[4] = RexpSM[4] - Bugey_3_95_RatioN/Bugey_3_95_RatioD;
}

/********************
*  Gosgen (37.9 m)  *
*********************/

  double *Gosgen_38_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *Gosgen_38_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
  double *Gosgen_38_U235_0 = glbGetRuleRatePtr(exp+4, 0);
  double *Gosgen_38_U238_0 = glbGetRuleRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239_0 = glbGetRuleRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241_0 = glbGetRuleRatePtr(exp+4, 3);
   
  double Gosgen_38_RatioN = 0.0;
  double Gosgen_38_RatioD = 0.0;

if(YesNo[5]==1){
   
  Gosgen_38_RatioN = x[0]*Gosgen38F[0]*Gosgen_38_U235[0];
  Gosgen_38_RatioN += x[1]*Gosgen38F[1]*Gosgen_38_U238[0];
  Gosgen_38_RatioN += x[2]*Gosgen38F[2]*Gosgen_38_Pu239[0];
  Gosgen_38_RatioN += x[3]*Gosgen38F[3]*Gosgen_38_Pu241[0];

  Gosgen_38_RatioD = Gosgen38F[0]*Gosgen_38_U235_0[0];
  Gosgen_38_RatioD += Gosgen38F[1]*Gosgen_38_U238_0[0];
  Gosgen_38_RatioD += Gosgen38F[2]*Gosgen_38_Pu239_0[0];
  Gosgen_38_RatioD += Gosgen38F[3]*Gosgen_38_Pu241_0[0];

  delta[5] = RexpSM[5] - Gosgen_38_RatioN/Gosgen_38_RatioD;
}

/********************
*  Gosgen (45.9 m)  *
*********************/

  double *Gosgen_46_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *Gosgen_46_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);
   
  double *Gosgen_46_U235_0 = glbGetRuleRatePtr(exp+5, 0);
  double *Gosgen_46_U238_0 = glbGetRuleRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239_0 = glbGetRuleRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241_0 = glbGetRuleRatePtr(exp+5, 3);
   
  double Gosgen_46_RatioN = 0.0;
  double Gosgen_46_RatioD = 0.0;

if(YesNo[6]==1){
   
  Gosgen_46_RatioN = x[0]*Gosgen46F[0]*Gosgen_46_U235[0];
  Gosgen_46_RatioN += x[1]*Gosgen46F[1]*Gosgen_46_U238[0];
  Gosgen_46_RatioN += x[2]*Gosgen46F[2]*Gosgen_46_Pu239[0];
  Gosgen_46_RatioN += x[3]*Gosgen46F[3]*Gosgen_46_Pu241[0];

  Gosgen_46_RatioD = Gosgen46F[0]*Gosgen_46_U235_0[0];
  Gosgen_46_RatioD += Gosgen46F[1]*Gosgen_46_U238_0[0];
  Gosgen_46_RatioD += Gosgen46F[2]*Gosgen_46_Pu239_0[0];
  Gosgen_46_RatioD += Gosgen46F[3]*Gosgen_46_Pu241_0[0];

  delta[6] = RexpSM[6] - Gosgen_46_RatioN/Gosgen_46_RatioD;
}

/********************
*  Gosgen (64.7 m)  *
*********************/

  double *Gosgen_65_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *Gosgen_65_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);

  double *Gosgen_65_U235_0 = glbGetRuleRatePtr(exp+6, 0);
  double *Gosgen_65_U238_0 = glbGetRuleRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239_0 = glbGetRuleRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241_0 = glbGetRuleRatePtr(exp+6, 3);
   
  double Gosgen_65_RatioN = 0.0;
  double Gosgen_65_RatioD = 0.0;

if(YesNo[7]==1){
   
  Gosgen_65_RatioN = x[0]*Gosgen65F[0]*Gosgen_65_U235[0];
  Gosgen_65_RatioN += x[1]*Gosgen65F[1]*Gosgen_65_U238[0];
  Gosgen_65_RatioN += x[2]*Gosgen65F[2]*Gosgen_65_Pu239[0];
  Gosgen_65_RatioN += x[3]*Gosgen65F[3]*Gosgen_65_Pu241[0];

  Gosgen_65_RatioD = Gosgen65F[0]*Gosgen_65_U235_0[0];
  Gosgen_65_RatioD += Gosgen65F[1]*Gosgen_65_U238_0[0];
  Gosgen_65_RatioD += Gosgen65F[2]*Gosgen_65_Pu239_0[0];
  Gosgen_65_RatioD += Gosgen65F[3]*Gosgen_65_Pu241_0[0];

  delta[7] = RexpSM[7] - Gosgen_65_RatioN/Gosgen_65_RatioD;
}

/*****************
*  ILL (8.76 m)  *
******************/

  double *ILL_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *ILL_U235_0 = glbGetRuleRatePtr(exp+7, 0);

if(YesNo[8]==1){
  delta[8] = RexpSM[8] - x[0]*ILL_U235[0]/ILL_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (32.8 m)  *
****************************/

  double *Krasnoyarsk_33_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *Krasnoyarsk_33_U235_0 = glbGetRuleRatePtr(exp+8, 0);

if(YesNo[9]==1){
  delta[9] = RexpSM[9] - x[0]* Krasnoyarsk_33_U235[0]/Krasnoyarsk_33_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (92.3 m)  *
****************************/

  double *Krasnoyarsk_92_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *Krasnoyarsk_92_U235_0 = glbGetRuleRatePtr(exp+9, 0);

if(YesNo[10]==1){
  delta[10] = RexpSM[10] - x[0]* Krasnoyarsk_92_U235[0]/Krasnoyarsk_92_U235_0[0];
}

/***************************
*  Krasnoyarsk94 (57.0 m)  *
****************************/

  double *Krasnoyarsk_57_U235 = glbGetSignalFitRatePtr(exp+10, 0);
  double *Krasnoyarsk_57_U235_0 = glbGetRuleRatePtr(exp+10, 0);

if(YesNo[11]==1){
  delta[11] = RexpSM[11] - x[0]* Krasnoyarsk_57_U235[0]/Krasnoyarsk_57_U235_0[0];
}

/***************************
*  Krasnoyarsk99 (34.0 m)  *
****************************/

  double *Krasnoyarsk_34_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *Krasnoyarsk_34_U235_0 = glbGetRuleRatePtr(exp+11, 0);

if(YesNo[12]==1){
  delta[12] = RexpSM[12] - x[0]* Krasnoyarsk_34_U235[0]/Krasnoyarsk_34_U235_0[0];
}

/********************
*  SRP-18 (18.2 m)  *
*********************/

  double *SRP_18_U235 = glbGetSignalFitRatePtr(exp+12, 0);
  double *SRP_18_U235_0 = glbGetRuleRatePtr(exp+12, 0);

if(YesNo[13]==1){
  delta[13] = RexpSM[13] - x[0]* SRP_18_U235[0]/SRP_18_U235_0[0];
}

/********************
*  SRP-24 (23.8 m)  *
*********************/

  double *SRP_24_U235 = glbGetSignalFitRatePtr(exp+13, 0);
  double *SRP_24_U235_0 = glbGetRuleRatePtr(exp+13, 0);

if(YesNo[14]==1){
  delta[14] = RexpSM[14] - x[0]*SRP_24_U235[0]/SRP_24_U235_0[0];
}

/***************
*  Rovno88-1I  *
****************/

  double *Rovno88_1I_U235 = glbGetSignalFitRatePtr(exp+14, 0);
  double *Rovno88_1I_U238 = glbGetSignalFitRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239 = glbGetSignalFitRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241 = glbGetSignalFitRatePtr(exp+14, 3);

  double *Rovno88_1I_U235_0 = glbGetRuleRatePtr(exp+14, 0);
  double *Rovno88_1I_U238_0 = glbGetRuleRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239_0 = glbGetRuleRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241_0 = glbGetRuleRatePtr(exp+14, 3);

  double Rovno88_1I_RatioN = 0.0;
  double Rovno88_1I_RatioD = 0.0;

if(YesNo[15]==1){
  Rovno88_1I_RatioN = x[0]* Rovno88_1IF[0]*Rovno88_1I_U235[0];
  Rovno88_1I_RatioN += x[1]* Rovno88_1IF[1]*Rovno88_1I_U238[0];
  Rovno88_1I_RatioN += x[2]* Rovno88_1IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1I_RatioN += x[3]* Rovno88_1IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1I_RatioD = Rovno88_1IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[3]*Rovno88_1I_Pu241_0[0];

  delta[15] = RexpSM[15] - Rovno88_1I_RatioN/Rovno88_1I_RatioD;
}

/***************
*  Rovno88-2I  *
****************/

  double Rovno88_2I_RatioN = 0.0;
  double Rovno88_2I_RatioD = 0.0;

if(YesNo[16]==1){   
  Rovno88_2I_RatioN = x[0]* Rovno88_2IF[0]*Rovno88_1I_U235[0];
  Rovno88_2I_RatioN += x[1]* Rovno88_2IF[1]*Rovno88_1I_U238[0];
  Rovno88_2I_RatioN += x[2]* Rovno88_2IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_2I_RatioN += x[3]* Rovno88_2IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_2I_RatioD = Rovno88_2IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[3]*Rovno88_1I_Pu241_0[0];

  delta[16] = RexpSM[16] - Rovno88_2I_RatioN/Rovno88_2I_RatioD;
}
/***************
*  Rovno88-1S  *
****************/

  double Rovno88_1S_RatioN = 0.0;
  double Rovno88_1S_RatioD = 0.0;

if(YesNo[17]==1){
  Rovno88_1S_RatioN = x[0]* Rovno88_1SF[0]*Rovno88_1I_U235[0];
  Rovno88_1S_RatioN += x[1]* Rovno88_1SF[1]*Rovno88_1I_U238[0];
  Rovno88_1S_RatioN += x[2]* Rovno88_1SF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1S_RatioN += x[3]* Rovno88_1SF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1S_RatioD = Rovno88_1SF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[3]*Rovno88_1I_Pu241_0[0];

  delta[17] = RexpSM[17] - Rovno88_1S_RatioN/Rovno88_1S_RatioD;
}

/**********************
*  Rovno88-2S (25 m)  *
***********************/
   
  double Rovno88_2S_25_RatioN = 0.0;
  double Rovno88_2S_25_RatioD = 0.0;

  double *Rovno88_2S_25_U235 = glbGetSignalFitRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238 = glbGetSignalFitRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239 = glbGetSignalFitRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241 = glbGetSignalFitRatePtr(exp+15, 3);

  double *Rovno88_2S_25_U235_0 = glbGetRuleRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238_0 = glbGetRuleRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239_0 = glbGetRuleRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241_0 = glbGetRuleRatePtr(exp+15, 3);

if(YesNo[18]==1){
  Rovno88_2S_25_RatioN = x[0]*Rovno88_2S_25F[0]*Rovno88_2S_25_U235[0];
  Rovno88_2S_25_RatioN += x[1]*Rovno88_2S_25F[1]*Rovno88_2S_25_U238[0];
  Rovno88_2S_25_RatioN += x[2]*Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239[0];
  Rovno88_2S_25_RatioN += x[3]*Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241[0];

  Rovno88_2S_25_RatioD = Rovno88_2S_25F[0]*Rovno88_2S_25_U235_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[1]*Rovno88_2S_25_U238_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241_0[0];

  delta[18] = RexpSM[18] - Rovno88_2S_25_RatioN/Rovno88_2S_25_RatioD;
}

/**********************
*  Rovno88-2S (18 m)  *
***********************/

  double Rovno88_2S_18_RatioN = 0.0;
  double Rovno88_2S_18_RatioD = 0.0;

if(YesNo[19]==1){
  Rovno88_2S_18_RatioN = x[0]*Rovno88_2S_18F[0]*Rovno88_1I_U235[0];
  Rovno88_2S_18_RatioN += x[1]*Rovno88_2S_18F[1]*Rovno88_1I_U238[0];
  Rovno88_2S_18_RatioN += x[2]*Rovno88_2S_18F[2]*Rovno88_1I_Pu239[0];
  Rovno88_2S_18_RatioN += x[3]*Rovno88_2S_18F[3]*Rovno88_1I_Pu241[0];

  Rovno88_2S_18_RatioD = Rovno88_2S_18F[0]*Rovno88_1I_U235_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[1]*Rovno88_1I_U238_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[3]*Rovno88_1I_Pu241_0[0];

  delta[19] = RexpSM[19] - Rovno88_2S_18_RatioN/Rovno88_2S_18_RatioD;
}

/************
*  Nucifer  *
*************/

  double *Nucifer_U235 = glbGetSignalFitRatePtr(exp+16, 0);
  double *Nucifer_U238 = glbGetSignalFitRatePtr(exp+16, 1);
  double *Nucifer_Pu239 = glbGetSignalFitRatePtr(exp+16, 2);
  double *Nucifer_Pu241 = glbGetSignalFitRatePtr(exp+16, 3);

  double *Nucifer_U235_0 = glbGetRuleRatePtr(exp+16, 0);
  double *Nucifer_U238_0 = glbGetRuleRatePtr(exp+16, 1);
  double *Nucifer_Pu239_0 = glbGetRuleRatePtr(exp+16, 2);
  double *Nucifer_Pu241_0 = glbGetRuleRatePtr(exp+16, 3);

  double Nucifer_RatioN=0.0;
  double Nucifer_RatioD=0.0;

if(YesNo[20]==1){
  Nucifer_RatioN = x[0]*NuciferF[0]*Nucifer_U235[0];
  Nucifer_RatioN += x[1]*NuciferF[1]*Nucifer_U238[0];
  Nucifer_RatioN += x[2]*NuciferF[2]*Nucifer_Pu239[0];
  Nucifer_RatioN += x[3]*NuciferF[3]*Nucifer_Pu241[0];

  Nucifer_RatioD = NuciferF[0]*Nucifer_U235_0[0];
  Nucifer_RatioD += NuciferF[1]*Nucifer_U238_0[0];
  Nucifer_RatioD += NuciferF[2]*Nucifer_Pu239_0[0];
  Nucifer_RatioD += NuciferF[3]*Nucifer_Pu241_0[0];

  delta[20] = Rexp[20] - Nucifer_RatioN/Nucifer_RatioD;
}

/**********************************
*  Palo Verde -- 750 m and 890 m  *
***********************************/

  double PV_RatioN=0.0;
  double PV_RatioD=0.0;

  double *PV_750_U235 = glbGetSignalFitRatePtr(exp+17, 0);
  double *PV_750_U238 = glbGetSignalFitRatePtr(exp+17, 1);
  double *PV_750_Pu239 = glbGetSignalFitRatePtr(exp+17, 2);
  double *PV_750_Pu241 = glbGetSignalFitRatePtr(exp+17, 3);

  double *PV_750_U235_0 = glbGetRuleRatePtr(exp+17, 0);
  double *PV_750_U238_0 = glbGetRuleRatePtr(exp+17, 1);
  double *PV_750_Pu239_0 = glbGetRuleRatePtr(exp+17, 2);
  double *PV_750_Pu241_0 = glbGetRuleRatePtr(exp+17, 3);

  double *PV_890_U235 = glbGetSignalFitRatePtr(exp+18, 0);
  double *PV_890_U238 = glbGetSignalFitRatePtr(exp+18, 1);
  double *PV_890_Pu239 = glbGetSignalFitRatePtr(exp+18, 2);
  double *PV_890_Pu241 = glbGetSignalFitRatePtr(exp+18, 3);

  double *PV_890_U235_0 = glbGetRuleRatePtr(exp+18, 0);
  double *PV_890_U238_0 = glbGetRuleRatePtr(exp+18, 1);
  double *PV_890_Pu239_0 = glbGetRuleRatePtr(exp+18, 2);
  double *PV_890_Pu241_0 = glbGetRuleRatePtr(exp+18, 3);

if(YesNo[21]==1){
  PV_RatioN = x[0]*PVF[0]*(PV_750_U235[0] + PV_890_U235[0]);
  PV_RatioN += x[1]*PVF[1]*(PV_750_U238[0] + PV_890_U238[0]);
  PV_RatioN += x[2]*PVF[2]*(PV_750_Pu239[0] + PV_890_Pu239[0]);
  PV_RatioN += x[3]*PVF[3]*(PV_750_Pu241[0] + PV_890_Pu241[0]);

  PV_RatioD = PVF[0]*(PV_750_U235_0[0]+PV_890_U235_0[0]);
  PV_RatioD += PVF[1]*(PV_750_U238_0[0]+PV_890_U238_0[0]);
  PV_RatioD += PVF[2]*(PV_750_Pu239_0[0]+PV_890_Pu239_0[0]);
  PV_RatioD += PVF[3]*(PV_750_Pu241_0[0]+PV_890_Pu241_0[0]);

  delta[21] = RexpSM[21] - PV_RatioN/PV_RatioD;
}

/************************************
*  Double Chooz -- 355 m and 469 m  *
*************************************/

  double DC_RatioN=0.0;
  double DC_RatioD=0.0;

  double *DC_355_U235 = glbGetSignalFitRatePtr(exp+19, 0);
  double *DC_355_U238 = glbGetSignalFitRatePtr(exp+19, 1);
  double *DC_355_Pu239 = glbGetSignalFitRatePtr(exp+19, 2);
  double *DC_355_Pu241 = glbGetSignalFitRatePtr(exp+19, 3);

  double *DC_355_U235_0 = glbGetRuleRatePtr(exp+19, 0);
  double *DC_355_U238_0 = glbGetRuleRatePtr(exp+19, 1);
  double *DC_355_Pu239_0 = glbGetRuleRatePtr(exp+19, 2);
  double *DC_355_Pu241_0 = glbGetRuleRatePtr(exp+19, 3);

  double *DC_469_U235 = glbGetSignalFitRatePtr(exp+20, 0);
  double *DC_469_U238 = glbGetSignalFitRatePtr(exp+20, 1);
  double *DC_469_Pu239 = glbGetSignalFitRatePtr(exp+20, 2);
  double *DC_469_Pu241 = glbGetSignalFitRatePtr(exp+20, 3);

  double *DC_469_U235_0 = glbGetRuleRatePtr(exp+20, 0);
  double *DC_469_U238_0 = glbGetRuleRatePtr(exp+20, 1);
  double *DC_469_Pu239_0 = glbGetRuleRatePtr(exp+20, 2);
  double *DC_469_Pu241_0 = glbGetRuleRatePtr(exp+20, 3);

if(YesNo[22]==1){
  DC_RatioN = x[0]*DCF[0]*(DC_355_U235[0] + DC_469_U235[0]);
  DC_RatioN += x[1]*DCF[1]*(DC_355_U238[0] + DC_469_U238[0]);
  DC_RatioN += x[2]*DCF[2]*(DC_355_Pu239[0] + DC_469_Pu239[0]);
  DC_RatioN += x[3]*DCF[3]*(DC_355_Pu241[0] + DC_469_Pu241[0]);

  DC_RatioD = DCF[0]*(DC_355_U235_0[0] + DC_469_U235_0[0]);
  DC_RatioD += DCF[1]*(DC_355_U238_0[0] + DC_469_U238_0[0]);
  DC_RatioD += DCF[2]*(DC_355_Pu239_0[0] + DC_469_Pu239_0[0]);
  DC_RatioD += DCF[3]*(DC_355_Pu241_0[0] + DC_469_Pu241_0[0]);

  delta[22] = RexpSM[22] - DC_RatioN/DC_RatioD;
}

/******************************
*  Chooz -- 998 m and 1115 m  *
*******************************/

  double Chooz_RatioN=0.0;
  double Chooz_RatioD=0.0;

  double *Chooz_998_U235 = glbGetSignalFitRatePtr(exp+21, 0);
  double *Chooz_998_U238 = glbGetSignalFitRatePtr(exp+21, 1);
  double *Chooz_998_Pu239 = glbGetSignalFitRatePtr(exp+21, 2);
  double *Chooz_998_Pu241 = glbGetSignalFitRatePtr(exp+21, 3);

  double *Chooz_998_U235_0 = glbGetRuleRatePtr(exp+21, 0);
  double *Chooz_998_U238_0 = glbGetRuleRatePtr(exp+21, 1);
  double *Chooz_998_Pu239_0 = glbGetRuleRatePtr(exp+21, 2);
  double *Chooz_998_Pu241_0 = glbGetRuleRatePtr(exp+21, 3);

  double *Chooz_1115_U235 = glbGetSignalFitRatePtr(exp+22, 0);
  double *Chooz_1115_U238 = glbGetSignalFitRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239 = glbGetSignalFitRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241 = glbGetSignalFitRatePtr(exp+22, 3);

  double *Chooz_1115_U235_0 = glbGetRuleRatePtr(exp+22, 0);
  double *Chooz_1115_U238_0 = glbGetRuleRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239_0 = glbGetRuleRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241_0 = glbGetRuleRatePtr(exp+22, 3);

if(YesNo[23]==1){
  Chooz_RatioN = x[0]*ChoozF[0]*(Chooz_998_U235[0] + Chooz_1115_U235[0]);
  Chooz_RatioN += x[1]*ChoozF[1]*(Chooz_998_U238[0] + Chooz_1115_U238[0]);
  Chooz_RatioN += x[2]*ChoozF[2]*(Chooz_998_Pu239[0] + Chooz_1115_Pu239[0]);
  Chooz_RatioN += x[3]*ChoozF[3]*(Chooz_998_Pu241[0] + Chooz_1115_Pu241[0]);

  Chooz_RatioD = ChoozF[0]*(Chooz_998_U235_0[0] + Chooz_1115_U235_0[0]);
  Chooz_RatioD += ChoozF[1]*(Chooz_998_U238_0[0] + Chooz_1115_U238_0[0]);
  Chooz_RatioD += ChoozF[2]*(Chooz_998_Pu239_0[0] + Chooz_1115_Pu239_0[0]);
  Chooz_RatioD += ChoozF[3]*(Chooz_998_Pu241_0[0] + Chooz_1115_Pu241_0[0]);

  delta[23] = RexpSM[23] - Chooz_RatioN/Chooz_RatioD;
}

/************
*  EH1 AD1  *
*************/

  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+23, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+23, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+23, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+23, 3);
   
  double EH1_AD1_Total = 0.0;
  double EH1_AD1_Total_0 = 0.0;

/************
*  EH1 AD2  *
*************/

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+24, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+24, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+24, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+24, 3);

  double EH1_AD2_Total = 0.0;
  double EH1_AD2_Total_0 = 0.0;

/************
*  EH2 AD3  *
*************/

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+25, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+25, 3);

  double *EH2_AD3_U235_0 = glbGetRuleRatePtr(exp+25, 0);
  double *EH2_AD3_U238_0 = glbGetRuleRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239_0 = glbGetRuleRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241_0 = glbGetRuleRatePtr(exp+25, 3);
   
  double EH2_AD3_Total = 0.0;
  double EH2_AD3_Total_0 = 0.0;

/************
*  EH2 AD8  *
*************/

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+26, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+26, 3);

  double *EH2_AD8_U235_0 = glbGetRuleRatePtr(exp+26, 0);
  double *EH2_AD8_U238_0 = glbGetRuleRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239_0 = glbGetRuleRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241_0 = glbGetRuleRatePtr(exp+26, 3);

  double EH2_AD8_Total = 0.0;
  double EH2_AD8_Total_0 = 0.0;

/*********
*  RENO  *
**********/

  double *RENO_U235 = glbGetSignalFitRatePtr(exp+27, 0);
  double *RENO_U238 = glbGetSignalFitRatePtr(exp+27, 1);
  double *RENO_Pu239 = glbGetSignalFitRatePtr(exp+27, 2);
  double *RENO_Pu241 = glbGetSignalFitRatePtr(exp+27, 3);

  double *RENO_U235_0 = glbGetRuleRatePtr(exp+27, 0);
  double *RENO_U238_0 = glbGetRuleRatePtr(exp+27, 1);
  double *RENO_Pu239_0 = glbGetRuleRatePtr(exp+27, 2);
  double *RENO_Pu241_0 = glbGetRuleRatePtr(exp+27, 3);

  double RENO_Total = 0.0;
  double RENO_Total_0 = 0.0;

/*********************************
*  Assembling Daya Bay and RENO  *
**********************************/

  double chi2 = 0.0;
  int i,j;

  double num, denom;

  for (i=0; i<8; i++){
  if (YesNo[24] == 1){
    EH1_AD1_Total = x[0]*DayaBayFs[i][0]*EH1_AD1_U235[0];
    EH1_AD1_Total += x[1]*DayaBayFs[i][1]*EH1_AD1_U238[0];
    EH1_AD1_Total += x[2]*DayaBayFs[i][2]*EH1_AD1_Pu239[0];
    EH1_AD1_Total += x[3]*DayaBayFs[i][3]*EH1_AD1_Pu241[0];

    EH1_AD2_Total = x[0]*DayaBayFs[i][0]*EH1_AD2_U235[0];
    EH1_AD2_Total += x[1]*DayaBayFs[i][1]*EH1_AD2_U238[0];
    EH1_AD2_Total += x[2]*DayaBayFs[i][2]*EH1_AD2_Pu239[0];
    EH1_AD2_Total += x[3]*DayaBayFs[i][3]*EH1_AD2_Pu241[0];

    EH2_AD3_Total = x[0]*DayaBayFs[i][0]*EH2_AD3_U235[0];
    EH2_AD3_Total += x[1]*DayaBayFs[i][1]*EH2_AD3_U238[0];
    EH2_AD3_Total += x[2]*DayaBayFs[i][2]*EH2_AD3_Pu239[0];
    EH2_AD3_Total += x[3]*DayaBayFs[i][3]*EH2_AD3_Pu241[0];

    EH2_AD8_Total = x[0]*DayaBayFs[i][0]*EH2_AD8_U235[0];
    EH2_AD8_Total += x[1]*DayaBayFs[i][1]*EH2_AD8_U238[0];
    EH2_AD8_Total += x[2]*DayaBayFs[i][2]*EH2_AD8_Pu239[0];
    EH2_AD8_Total += x[3]*DayaBayFs[i][3]*EH2_AD8_Pu241[0];

    EH1_AD1_Total_0 = DayaBayFs[i][0]*EH1_AD1_U235_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][1]*EH1_AD1_U238_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][2]*EH1_AD1_Pu239_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][3]*EH1_AD1_Pu241_0[0];

    EH1_AD2_Total_0 = DayaBayFs[i][0]*EH1_AD2_U235_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][1]*EH1_AD2_U238_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][2]*EH1_AD2_Pu239_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][3]*EH1_AD2_Pu241_0[0];

    EH2_AD3_Total_0 = DayaBayFs[i][0]*EH2_AD3_U235_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][1]*EH2_AD3_U238_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][2]*EH2_AD3_Pu239_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][3]*EH2_AD3_Pu241_0[0];

    EH2_AD8_Total_0 = DayaBayFs[i][0]*EH2_AD8_U235_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][1]*EH2_AD8_U238_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][2]*EH2_AD8_Pu239_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][3]*EH2_AD8_Pu241_0[0];

    num = EH1_AD1_Total+ EH1_AD2_Total+ EH2_AD3_Total+ EH2_AD8_Total;
    denom = EH1_AD1_Total_0+ EH1_AD2_Total_0+ EH2_AD3_Total_0+ EH2_AD8_Total_0;

    delta[i+24] = num/denom - DBdataSM[i];
  }

  if (YesNo[25] == 1){
    RENO_Total = x[0]*RENOFs[i][0]*RENO_U235[0];
    RENO_Total += x[1]*RENOFs[i][1]*RENO_U238[0];
    RENO_Total += x[2]*RENOFs[i][2]*RENO_Pu239[0];
    RENO_Total += x[3]*RENOFs[i][3]*RENO_Pu241[0];

    RENO_Total_0 = RENOFs[i][0]*RENO_U235_0[0];
    RENO_Total_0 += RENOFs[i][1]*RENO_U238_0[0];
    RENO_Total_0 += RENOFs[i][2]*RENO_Pu239_0[0];
    RENO_Total_0 += RENOFs[i][3]*RENO_Pu241_0[0];

    delta[i+32] = RENO_Total/RENO_Total_0 - RENOdataSM[i] ;
  }
  }

/********************************
*  Putting the pieces together  *
*********************************/

  /* Adding the contributions from each experiment */

  for (i=0; i<24; i++){
    for (j=0; j<24; j++){
      chi2 += delta[i]*cov_matrix[i][j]*delta[j];
    }
  }

  /* Now, adding Daya Bay... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+24]*cov_matrix_DB[i][j]*delta[j+24];
    }
  }

  /* Now, adding RENO... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+32]*cov_matrix_RENO[i][j]*delta[j+32];
    }
  }

  /* Adding in the systematic uncertainties on the ratios... */

  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      chi2 += (x[i]-1.0)*(x[j]-1.0)*V_inv_SM[i][j]/(errors[i]*errors[j]);
    }
  }

  for (i = 0; i<4; i++){ /* Export nuisance parameters to main file */
    systematic[i] = x[i];
  }

  return chi2;
}

/********************************************************************
*                                                                   *
*     THE CALCULATION OF A CHI-SQUARED WITH SYSTEMATICS USING       *
*             HKSS CORRECTION TO HM ANTINEUTRINO FLUXES             *
*                                                                   *
*********************************************************************/

double combo_rate_chi_HKSS(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following, for instance, 1703.00860

  "delta" is given by (Rexp - Rpred)
*/

  int k;

  double delta[40] = {0.0};

  int YesNo[26];
  for (k=0; k<26; k++){
    YesNo[k] = ((int *)user_data)[k];
  }

/************
*  Bugey-4  *
*************/

  double *Bugey_4_U235 = glbGetSignalFitRatePtr(exp, 0);
  double *Bugey_4_U238 = glbGetSignalFitRatePtr(exp, 1);
  double *Bugey_4_Pu239 = glbGetSignalFitRatePtr(exp, 2);
  double *Bugey_4_Pu241 = glbGetSignalFitRatePtr(exp, 3);

  double *Bugey_4_U235_0 = glbGetRuleRatePtr(exp, 0);
  double *Bugey_4_U238_0 = glbGetRuleRatePtr(exp, 1);
  double *Bugey_4_Pu239_0 = glbGetRuleRatePtr(exp, 2);
  double *Bugey_4_Pu241_0 = glbGetRuleRatePtr(exp, 3);
   
  double Bugey_4_RatioN = 0.0;
  double Bugey_4_RatioD = 0.0;

if(YesNo[0]==1||YesNo[2]==1){
  Bugey_4_RatioN = x[0]*BugeyF[0]*Bugey_4_U235[0];
  Bugey_4_RatioN += x[1]*BugeyF[1]*Bugey_4_U238[0];
  Bugey_4_RatioN += x[2]*BugeyF[2]*Bugey_4_Pu239[0];
  Bugey_4_RatioN += x[3]*BugeyF[3]*Bugey_4_Pu241[0];

  Bugey_4_RatioD = BugeyF[0]*Bugey_4_U235_0[0];
  Bugey_4_RatioD += BugeyF[1]*Bugey_4_U238_0[0];
  Bugey_4_RatioD += BugeyF[2]*Bugey_4_Pu239_0[0];
  Bugey_4_RatioD += BugeyF[3]*Bugey_4_Pu241_0[0];
}
if(YesNo[0]==1){
  delta[0] = RexpHKSS[0] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/************
*  Rovno91  *
*************/

  double *Rovno91_U235 = glbGetSignalFitRatePtr(exp+1, 0);
  double *Rovno91_U238 = glbGetSignalFitRatePtr(exp+1, 1);
  double *Rovno91_Pu239 = glbGetSignalFitRatePtr(exp+1, 2);
  double *Rovno91_Pu241 = glbGetSignalFitRatePtr(exp+1, 3);

  double *Rovno91_U235_0 = glbGetRuleRatePtr(exp+1, 0);
  double *Rovno91_U238_0 = glbGetRuleRatePtr(exp+1, 1);
  double *Rovno91_Pu239_0 = glbGetRuleRatePtr(exp+1, 2);
  double *Rovno91_Pu241_0 = glbGetRuleRatePtr(exp+1, 3);

  double Rovno91_RatioN = 0.0;
  double Rovno91_RatioD = 0.0;

if(YesNo[1]==1){
  Rovno91_RatioN = x[0]*Rovno_91F[0]*Rovno91_U235[0];
  Rovno91_RatioN += x[1]*Rovno_91F[1]*Rovno91_U238[0];
  Rovno91_RatioN += x[2]*Rovno_91F[2]*Rovno91_Pu239[0];
  Rovno91_RatioN += x[3]*Rovno_91F[3]*Rovno91_Pu241[0];

  Rovno91_RatioD = Rovno_91F[0]*Rovno91_U235_0[0];
  Rovno91_RatioD += Rovno_91F[1]*Rovno91_U238_0[0];
  Rovno91_RatioD += Rovno_91F[2]*Rovno91_Pu239_0[0];
  Rovno91_RatioD += Rovno_91F[3]*Rovno91_Pu241_0[0];

  delta[1] = RexpHKSS[1] - Rovno91_RatioN/Rovno91_RatioD;
}

/*******************
*  Bugey-3 (15 m)  *
********************/

/*
  No need to recalculate rates -- reuse the output for Bugey-4!
*/

if(YesNo[2]==1){
  delta[2] = RexpHKSS[2] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/*******************
*  Bugey-3 (40 m)  *
********************/

  double *Bugey_3_40_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *Bugey_3_40_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);
   
  double *Bugey_3_40_U235_0 = glbGetRuleRatePtr(exp+2, 0);
  double *Bugey_3_40_U238_0 = glbGetRuleRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239_0 = glbGetRuleRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241_0 = glbGetRuleRatePtr(exp+2, 3);
   
  double Bugey_3_40_RatioN = 0.0;
  double Bugey_3_40_RatioD = 0.0;

if(YesNo[3]==1){   
  Bugey_3_40_RatioN = x[0]*BugeyF[0]*Bugey_3_40_U235[0];
  Bugey_3_40_RatioN += x[1]*BugeyF[1]*Bugey_3_40_U238[0];
  Bugey_3_40_RatioN += x[2]*BugeyF[2]*Bugey_3_40_Pu239[0];
  Bugey_3_40_RatioN += x[3]*BugeyF[3]*Bugey_3_40_Pu241[0];

  Bugey_3_40_RatioD = BugeyF[0]*Bugey_3_40_U235_0[0];
  Bugey_3_40_RatioD += BugeyF[1]*Bugey_3_40_U238_0[0];
  Bugey_3_40_RatioD += BugeyF[2]*Bugey_3_40_Pu239_0[0];
  Bugey_3_40_RatioD += BugeyF[3]*Bugey_3_40_Pu241_0[0];

  delta[3] = RexpHKSS[3] - Bugey_3_40_RatioN/Bugey_3_40_RatioD;
}

/*******************
*  Bugey-3 (95 m)  *
********************/

  double *Bugey_3_95_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *Bugey_3_95_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);
   
  double *Bugey_3_95_U235_0 = glbGetRuleRatePtr(exp+3, 0);
  double *Bugey_3_95_U238_0 = glbGetRuleRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239_0 = glbGetRuleRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241_0 = glbGetRuleRatePtr(exp+3, 3);
   
  double Bugey_3_95_RatioN = 0.0;
  double Bugey_3_95_RatioD = 0.0;

if(YesNo[4]==1){   
  Bugey_3_95_RatioN = x[0]*BugeyF[0]*Bugey_3_95_U235[0];
  Bugey_3_95_RatioN += x[1]*BugeyF[1]*Bugey_3_95_U238[0];
  Bugey_3_95_RatioN += x[2]*BugeyF[2]*Bugey_3_95_Pu239[0];
  Bugey_3_95_RatioN += x[3]*BugeyF[3]*Bugey_3_95_Pu241[0];

  Bugey_3_95_RatioD = BugeyF[0]*Bugey_3_95_U235_0[0];
  Bugey_3_95_RatioD += BugeyF[1]*Bugey_3_95_U238_0[0];
  Bugey_3_95_RatioD += BugeyF[2]*Bugey_3_95_Pu239_0[0];
  Bugey_3_95_RatioD += BugeyF[3]*Bugey_3_95_Pu241_0[0];

  delta[4] = RexpHKSS[4] - Bugey_3_95_RatioN/Bugey_3_95_RatioD;
}

/********************
*  Gosgen (37.9 m)  *
*********************/

  double *Gosgen_38_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *Gosgen_38_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
  double *Gosgen_38_U235_0 = glbGetRuleRatePtr(exp+4, 0);
  double *Gosgen_38_U238_0 = glbGetRuleRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239_0 = glbGetRuleRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241_0 = glbGetRuleRatePtr(exp+4, 3);
   
  double Gosgen_38_RatioN = 0.0;
  double Gosgen_38_RatioD = 0.0;

if(YesNo[5]==1){
   
  Gosgen_38_RatioN = x[0]*Gosgen38F[0]*Gosgen_38_U235[0];
  Gosgen_38_RatioN += x[1]*Gosgen38F[1]*Gosgen_38_U238[0];
  Gosgen_38_RatioN += x[2]*Gosgen38F[2]*Gosgen_38_Pu239[0];
  Gosgen_38_RatioN += x[3]*Gosgen38F[3]*Gosgen_38_Pu241[0];

  Gosgen_38_RatioD = Gosgen38F[0]*Gosgen_38_U235_0[0];
  Gosgen_38_RatioD += Gosgen38F[1]*Gosgen_38_U238_0[0];
  Gosgen_38_RatioD += Gosgen38F[2]*Gosgen_38_Pu239_0[0];
  Gosgen_38_RatioD += Gosgen38F[3]*Gosgen_38_Pu241_0[0];

  delta[5] = RexpHKSS[5] - Gosgen_38_RatioN/Gosgen_38_RatioD;
}

/********************
*  Gosgen (45.9 m)  *
*********************/

  double *Gosgen_46_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *Gosgen_46_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);
   
  double *Gosgen_46_U235_0 = glbGetRuleRatePtr(exp+5, 0);
  double *Gosgen_46_U238_0 = glbGetRuleRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239_0 = glbGetRuleRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241_0 = glbGetRuleRatePtr(exp+5, 3);
   
  double Gosgen_46_RatioN = 0.0;
  double Gosgen_46_RatioD = 0.0;

if(YesNo[6]==1){
   
  Gosgen_46_RatioN = x[0]*Gosgen46F[0]*Gosgen_46_U235[0];
  Gosgen_46_RatioN += x[1]*Gosgen46F[1]*Gosgen_46_U238[0];
  Gosgen_46_RatioN += x[2]*Gosgen46F[2]*Gosgen_46_Pu239[0];
  Gosgen_46_RatioN += x[3]*Gosgen46F[3]*Gosgen_46_Pu241[0];

  Gosgen_46_RatioD = Gosgen46F[0]*Gosgen_46_U235_0[0];
  Gosgen_46_RatioD += Gosgen46F[1]*Gosgen_46_U238_0[0];
  Gosgen_46_RatioD += Gosgen46F[2]*Gosgen_46_Pu239_0[0];
  Gosgen_46_RatioD += Gosgen46F[3]*Gosgen_46_Pu241_0[0];

  delta[6] = RexpHKSS[6] - Gosgen_46_RatioN/Gosgen_46_RatioD;
}

/********************
*  Gosgen (64.7 m)  *
*********************/

  double *Gosgen_65_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *Gosgen_65_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);

  double *Gosgen_65_U235_0 = glbGetRuleRatePtr(exp+6, 0);
  double *Gosgen_65_U238_0 = glbGetRuleRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239_0 = glbGetRuleRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241_0 = glbGetRuleRatePtr(exp+6, 3);
   
  double Gosgen_65_RatioN = 0.0;
  double Gosgen_65_RatioD = 0.0;

if(YesNo[7]==1){
   
  Gosgen_65_RatioN = x[0]*Gosgen65F[0]*Gosgen_65_U235[0];
  Gosgen_65_RatioN += x[1]*Gosgen65F[1]*Gosgen_65_U238[0];
  Gosgen_65_RatioN += x[2]*Gosgen65F[2]*Gosgen_65_Pu239[0];
  Gosgen_65_RatioN += x[3]*Gosgen65F[3]*Gosgen_65_Pu241[0];

  Gosgen_65_RatioD = Gosgen65F[0]*Gosgen_65_U235_0[0];
  Gosgen_65_RatioD += Gosgen65F[1]*Gosgen_65_U238_0[0];
  Gosgen_65_RatioD += Gosgen65F[2]*Gosgen_65_Pu239_0[0];
  Gosgen_65_RatioD += Gosgen65F[3]*Gosgen_65_Pu241_0[0];

  delta[7] = RexpHKSS[7] - Gosgen_65_RatioN/Gosgen_65_RatioD;
}

/*****************
*  ILL (8.76 m)  *
******************/

  double *ILL_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *ILL_U235_0 = glbGetRuleRatePtr(exp+7, 0);

if(YesNo[8]==1){
  delta[8] = RexpHKSS[8] - x[0]*ILL_U235[0]/ILL_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (32.8 m)  *
****************************/

  double *Krasnoyarsk_33_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *Krasnoyarsk_33_U235_0 = glbGetRuleRatePtr(exp+8, 0);

if(YesNo[9]==1){
  delta[9] = RexpHKSS[9] - x[0]* Krasnoyarsk_33_U235[0]/Krasnoyarsk_33_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (92.3 m)  *
****************************/

  double *Krasnoyarsk_92_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *Krasnoyarsk_92_U235_0 = glbGetRuleRatePtr(exp+9, 0);

if(YesNo[10]==1){
  delta[10] = RexpHKSS[10] - x[0]* Krasnoyarsk_92_U235[0]/Krasnoyarsk_92_U235_0[0];
}

/***************************
*  Krasnoyarsk94 (57.0 m)  *
****************************/

  double *Krasnoyarsk_57_U235 = glbGetSignalFitRatePtr(exp+10, 0);
  double *Krasnoyarsk_57_U235_0 = glbGetRuleRatePtr(exp+10, 0);

if(YesNo[11]==1){
  delta[11] = RexpHKSS[11] - x[0]* Krasnoyarsk_57_U235[0]/Krasnoyarsk_57_U235_0[0];
}

/***************************
*  Krasnoyarsk99 (34.0 m)  *
****************************/

  double *Krasnoyarsk_34_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *Krasnoyarsk_34_U235_0 = glbGetRuleRatePtr(exp+11, 0);

if(YesNo[12]==1){
  delta[12] = RexpHKSS[12] - x[0]* Krasnoyarsk_34_U235[0]/Krasnoyarsk_34_U235_0[0];
}

/********************
*  SRP-18 (18.2 m)  *
*********************/

  double *SRP_18_U235 = glbGetSignalFitRatePtr(exp+12, 0);
  double *SRP_18_U235_0 = glbGetRuleRatePtr(exp+12, 0);

if(YesNo[13]==1){
  delta[13] = RexpHKSS[13] - x[0]* SRP_18_U235[0]/SRP_18_U235_0[0];
}

/********************
*  SRP-24 (23.8 m)  *
*********************/

  double *SRP_24_U235 = glbGetSignalFitRatePtr(exp+13, 0);
  double *SRP_24_U235_0 = glbGetRuleRatePtr(exp+13, 0);

if(YesNo[14]==1){
  delta[14] = RexpHKSS[14] - x[0]*SRP_24_U235[0]/SRP_24_U235_0[0];
}

/***************
*  Rovno88-1I  *
****************/

  double *Rovno88_1I_U235 = glbGetSignalFitRatePtr(exp+14, 0);
  double *Rovno88_1I_U238 = glbGetSignalFitRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239 = glbGetSignalFitRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241 = glbGetSignalFitRatePtr(exp+14, 3);

  double *Rovno88_1I_U235_0 = glbGetRuleRatePtr(exp+14, 0);
  double *Rovno88_1I_U238_0 = glbGetRuleRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239_0 = glbGetRuleRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241_0 = glbGetRuleRatePtr(exp+14, 3);

  double Rovno88_1I_RatioN = 0.0;
  double Rovno88_1I_RatioD = 0.0;

if(YesNo[15]==1){
  Rovno88_1I_RatioN = x[0]* Rovno88_1IF[0]*Rovno88_1I_U235[0];
  Rovno88_1I_RatioN += x[1]* Rovno88_1IF[1]*Rovno88_1I_U238[0];
  Rovno88_1I_RatioN += x[2]* Rovno88_1IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1I_RatioN += x[3]* Rovno88_1IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1I_RatioD = Rovno88_1IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[3]*Rovno88_1I_Pu241_0[0];

  delta[15] = RexpHKSS[15] - Rovno88_1I_RatioN/Rovno88_1I_RatioD;
}

/***************
*  Rovno88-2I  *
****************/

  double Rovno88_2I_RatioN = 0.0;
  double Rovno88_2I_RatioD = 0.0;

if(YesNo[16]==1){   
  Rovno88_2I_RatioN = x[0]* Rovno88_2IF[0]*Rovno88_1I_U235[0];
  Rovno88_2I_RatioN += x[1]* Rovno88_2IF[1]*Rovno88_1I_U238[0];
  Rovno88_2I_RatioN += x[2]* Rovno88_2IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_2I_RatioN += x[3]* Rovno88_2IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_2I_RatioD = Rovno88_2IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[3]*Rovno88_1I_Pu241_0[0];

  delta[16] = RexpHKSS[16] - Rovno88_2I_RatioN/Rovno88_2I_RatioD;
}
/***************
*  Rovno88-1S  *
****************/

  double Rovno88_1S_RatioN = 0.0;
  double Rovno88_1S_RatioD = 0.0;

if(YesNo[17]==1){
  Rovno88_1S_RatioN = x[0]* Rovno88_1SF[0]*Rovno88_1I_U235[0];
  Rovno88_1S_RatioN += x[1]* Rovno88_1SF[1]*Rovno88_1I_U238[0];
  Rovno88_1S_RatioN += x[2]* Rovno88_1SF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1S_RatioN += x[3]* Rovno88_1SF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1S_RatioD = Rovno88_1SF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[3]*Rovno88_1I_Pu241_0[0];

  delta[17] = RexpHKSS[17] - Rovno88_1S_RatioN/Rovno88_1S_RatioD;
}

/**********************
*  Rovno88-2S (25 m)  *
***********************/
   
  double Rovno88_2S_25_RatioN = 0.0;
  double Rovno88_2S_25_RatioD = 0.0;

  double *Rovno88_2S_25_U235 = glbGetSignalFitRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238 = glbGetSignalFitRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239 = glbGetSignalFitRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241 = glbGetSignalFitRatePtr(exp+15, 3);

  double *Rovno88_2S_25_U235_0 = glbGetRuleRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238_0 = glbGetRuleRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239_0 = glbGetRuleRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241_0 = glbGetRuleRatePtr(exp+15, 3);

if(YesNo[18]==1){
  Rovno88_2S_25_RatioN = x[0]*Rovno88_2S_25F[0]*Rovno88_2S_25_U235[0];
  Rovno88_2S_25_RatioN += x[1]*Rovno88_2S_25F[1]*Rovno88_2S_25_U238[0];
  Rovno88_2S_25_RatioN += x[2]*Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239[0];
  Rovno88_2S_25_RatioN += x[3]*Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241[0];

  Rovno88_2S_25_RatioD = Rovno88_2S_25F[0]*Rovno88_2S_25_U235_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[1]*Rovno88_2S_25_U238_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241_0[0];

  delta[18] = RexpHKSS[18] - Rovno88_2S_25_RatioN/Rovno88_2S_25_RatioD;
}

/**********************
*  Rovno88-2S (18 m)  *
***********************/

  double Rovno88_2S_18_RatioN = 0.0;
  double Rovno88_2S_18_RatioD = 0.0;

if(YesNo[19]==1){
  Rovno88_2S_18_RatioN = x[0]*Rovno88_2S_18F[0]*Rovno88_1I_U235[0];
  Rovno88_2S_18_RatioN += x[1]*Rovno88_2S_18F[1]*Rovno88_1I_U238[0];
  Rovno88_2S_18_RatioN += x[2]*Rovno88_2S_18F[2]*Rovno88_1I_Pu239[0];
  Rovno88_2S_18_RatioN += x[3]*Rovno88_2S_18F[3]*Rovno88_1I_Pu241[0];

  Rovno88_2S_18_RatioD = Rovno88_2S_18F[0]*Rovno88_1I_U235_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[1]*Rovno88_1I_U238_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[3]*Rovno88_1I_Pu241_0[0];

  delta[19] = RexpHKSS[19] - Rovno88_2S_18_RatioN/Rovno88_2S_18_RatioD;
}

/************
*  Nucifer  *
*************/

  double *Nucifer_U235 = glbGetSignalFitRatePtr(exp+16, 0);
  double *Nucifer_U238 = glbGetSignalFitRatePtr(exp+16, 1);
  double *Nucifer_Pu239 = glbGetSignalFitRatePtr(exp+16, 2);
  double *Nucifer_Pu241 = glbGetSignalFitRatePtr(exp+16, 3);

  double *Nucifer_U235_0 = glbGetRuleRatePtr(exp+16, 0);
  double *Nucifer_U238_0 = glbGetRuleRatePtr(exp+16, 1);
  double *Nucifer_Pu239_0 = glbGetRuleRatePtr(exp+16, 2);
  double *Nucifer_Pu241_0 = glbGetRuleRatePtr(exp+16, 3);

  double Nucifer_RatioN=0.0;
  double Nucifer_RatioD=0.0;

if(YesNo[20]==1){
  Nucifer_RatioN = x[0]*NuciferF[0]*Nucifer_U235[0];
  Nucifer_RatioN += x[1]*NuciferF[1]*Nucifer_U238[0];
  Nucifer_RatioN += x[2]*NuciferF[2]*Nucifer_Pu239[0];
  Nucifer_RatioN += x[3]*NuciferF[3]*Nucifer_Pu241[0];

  Nucifer_RatioD = NuciferF[0]*Nucifer_U235_0[0];
  Nucifer_RatioD += NuciferF[1]*Nucifer_U238_0[0];
  Nucifer_RatioD += NuciferF[2]*Nucifer_Pu239_0[0];
  Nucifer_RatioD += NuciferF[3]*Nucifer_Pu241_0[0];

  delta[20] = Rexp[20] - Nucifer_RatioN/Nucifer_RatioD;
}

/**********************************
*  Palo Verde -- 750 m and 890 m  *
***********************************/

  double PV_RatioN=0.0;
  double PV_RatioD=0.0;

  double *PV_750_U235 = glbGetSignalFitRatePtr(exp+17, 0);
  double *PV_750_U238 = glbGetSignalFitRatePtr(exp+17, 1);
  double *PV_750_Pu239 = glbGetSignalFitRatePtr(exp+17, 2);
  double *PV_750_Pu241 = glbGetSignalFitRatePtr(exp+17, 3);

  double *PV_750_U235_0 = glbGetRuleRatePtr(exp+17, 0);
  double *PV_750_U238_0 = glbGetRuleRatePtr(exp+17, 1);
  double *PV_750_Pu239_0 = glbGetRuleRatePtr(exp+17, 2);
  double *PV_750_Pu241_0 = glbGetRuleRatePtr(exp+17, 3);

  double *PV_890_U235 = glbGetSignalFitRatePtr(exp+18, 0);
  double *PV_890_U238 = glbGetSignalFitRatePtr(exp+18, 1);
  double *PV_890_Pu239 = glbGetSignalFitRatePtr(exp+18, 2);
  double *PV_890_Pu241 = glbGetSignalFitRatePtr(exp+18, 3);

  double *PV_890_U235_0 = glbGetRuleRatePtr(exp+18, 0);
  double *PV_890_U238_0 = glbGetRuleRatePtr(exp+18, 1);
  double *PV_890_Pu239_0 = glbGetRuleRatePtr(exp+18, 2);
  double *PV_890_Pu241_0 = glbGetRuleRatePtr(exp+18, 3);

if(YesNo[21]==1){
  PV_RatioN = x[0]*PVF[0]*(PV_750_U235[0] + PV_890_U235[0]);
  PV_RatioN += x[1]*PVF[1]*(PV_750_U238[0] + PV_890_U238[0]);
  PV_RatioN += x[2]*PVF[2]*(PV_750_Pu239[0] + PV_890_Pu239[0]);
  PV_RatioN += x[3]*PVF[3]*(PV_750_Pu241[0] + PV_890_Pu241[0]);

  PV_RatioD = PVF[0]*(PV_750_U235_0[0]+PV_890_U235_0[0]);
  PV_RatioD += PVF[1]*(PV_750_U238_0[0]+PV_890_U238_0[0]);
  PV_RatioD += PVF[2]*(PV_750_Pu239_0[0]+PV_890_Pu239_0[0]);
  PV_RatioD += PVF[3]*(PV_750_Pu241_0[0]+PV_890_Pu241_0[0]);

  delta[21] = RexpHKSS[21] - PV_RatioN/PV_RatioD;
}

/************************************
*  Double Chooz -- 355 m and 469 m  *
*************************************/

  double DC_RatioN=0.0;
  double DC_RatioD=0.0;

  double *DC_355_U235 = glbGetSignalFitRatePtr(exp+19, 0);
  double *DC_355_U238 = glbGetSignalFitRatePtr(exp+19, 1);
  double *DC_355_Pu239 = glbGetSignalFitRatePtr(exp+19, 2);
  double *DC_355_Pu241 = glbGetSignalFitRatePtr(exp+19, 3);

  double *DC_355_U235_0 = glbGetRuleRatePtr(exp+19, 0);
  double *DC_355_U238_0 = glbGetRuleRatePtr(exp+19, 1);
  double *DC_355_Pu239_0 = glbGetRuleRatePtr(exp+19, 2);
  double *DC_355_Pu241_0 = glbGetRuleRatePtr(exp+19, 3);

  double *DC_469_U235 = glbGetSignalFitRatePtr(exp+20, 0);
  double *DC_469_U238 = glbGetSignalFitRatePtr(exp+20, 1);
  double *DC_469_Pu239 = glbGetSignalFitRatePtr(exp+20, 2);
  double *DC_469_Pu241 = glbGetSignalFitRatePtr(exp+20, 3);

  double *DC_469_U235_0 = glbGetRuleRatePtr(exp+20, 0);
  double *DC_469_U238_0 = glbGetRuleRatePtr(exp+20, 1);
  double *DC_469_Pu239_0 = glbGetRuleRatePtr(exp+20, 2);
  double *DC_469_Pu241_0 = glbGetRuleRatePtr(exp+20, 3);

if(YesNo[22]==1){
  DC_RatioN = x[0]*DCF[0]*(DC_355_U235[0] + DC_469_U235[0]);
  DC_RatioN += x[1]*DCF[1]*(DC_355_U238[0] + DC_469_U238[0]);
  DC_RatioN += x[2]*DCF[2]*(DC_355_Pu239[0] + DC_469_Pu239[0]);
  DC_RatioN += x[3]*DCF[3]*(DC_355_Pu241[0] + DC_469_Pu241[0]);

  DC_RatioD = DCF[0]*(DC_355_U235_0[0] + DC_469_U235_0[0]);
  DC_RatioD += DCF[1]*(DC_355_U238_0[0] + DC_469_U238_0[0]);
  DC_RatioD += DCF[2]*(DC_355_Pu239_0[0] + DC_469_Pu239_0[0]);
  DC_RatioD += DCF[3]*(DC_355_Pu241_0[0] + DC_469_Pu241_0[0]);

  delta[22] = RexpHKSS[22] - DC_RatioN/DC_RatioD;
}

/******************************
*  Chooz -- 998 m and 1115 m  *
*******************************/

  double Chooz_RatioN=0.0;
  double Chooz_RatioD=0.0;

  double *Chooz_998_U235 = glbGetSignalFitRatePtr(exp+21, 0);
  double *Chooz_998_U238 = glbGetSignalFitRatePtr(exp+21, 1);
  double *Chooz_998_Pu239 = glbGetSignalFitRatePtr(exp+21, 2);
  double *Chooz_998_Pu241 = glbGetSignalFitRatePtr(exp+21, 3);

  double *Chooz_998_U235_0 = glbGetRuleRatePtr(exp+21, 0);
  double *Chooz_998_U238_0 = glbGetRuleRatePtr(exp+21, 1);
  double *Chooz_998_Pu239_0 = glbGetRuleRatePtr(exp+21, 2);
  double *Chooz_998_Pu241_0 = glbGetRuleRatePtr(exp+21, 3);

  double *Chooz_1115_U235 = glbGetSignalFitRatePtr(exp+22, 0);
  double *Chooz_1115_U238 = glbGetSignalFitRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239 = glbGetSignalFitRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241 = glbGetSignalFitRatePtr(exp+22, 3);

  double *Chooz_1115_U235_0 = glbGetRuleRatePtr(exp+22, 0);
  double *Chooz_1115_U238_0 = glbGetRuleRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239_0 = glbGetRuleRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241_0 = glbGetRuleRatePtr(exp+22, 3);

if(YesNo[23]==1){
  Chooz_RatioN = x[0]*ChoozF[0]*(Chooz_998_U235[0] + Chooz_1115_U235[0]);
  Chooz_RatioN += x[1]*ChoozF[1]*(Chooz_998_U238[0] + Chooz_1115_U238[0]);
  Chooz_RatioN += x[2]*ChoozF[2]*(Chooz_998_Pu239[0] + Chooz_1115_Pu239[0]);
  Chooz_RatioN += x[3]*ChoozF[3]*(Chooz_998_Pu241[0] + Chooz_1115_Pu241[0]);

  Chooz_RatioD = ChoozF[0]*(Chooz_998_U235_0[0] + Chooz_1115_U235_0[0]);
  Chooz_RatioD += ChoozF[1]*(Chooz_998_U238_0[0] + Chooz_1115_U238_0[0]);
  Chooz_RatioD += ChoozF[2]*(Chooz_998_Pu239_0[0] + Chooz_1115_Pu239_0[0]);
  Chooz_RatioD += ChoozF[3]*(Chooz_998_Pu241_0[0] + Chooz_1115_Pu241_0[0]);

  delta[23] = RexpHKSS[23] - Chooz_RatioN/Chooz_RatioD;
}

/************
*  EH1 AD1  *
*************/

  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+23, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+23, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+23, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+23, 3);
   
  double EH1_AD1_Total = 0.0;
  double EH1_AD1_Total_0 = 0.0;

/************
*  EH1 AD2  *
*************/

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+24, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+24, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+24, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+24, 3);

  double EH1_AD2_Total = 0.0;
  double EH1_AD2_Total_0 = 0.0;

/************
*  EH2 AD3  *
*************/

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+25, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+25, 3);

  double *EH2_AD3_U235_0 = glbGetRuleRatePtr(exp+25, 0);
  double *EH2_AD3_U238_0 = glbGetRuleRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239_0 = glbGetRuleRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241_0 = glbGetRuleRatePtr(exp+25, 3);
   
  double EH2_AD3_Total = 0.0;
  double EH2_AD3_Total_0 = 0.0;

/************
*  EH2 AD8  *
*************/

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+26, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+26, 3);

  double *EH2_AD8_U235_0 = glbGetRuleRatePtr(exp+26, 0);
  double *EH2_AD8_U238_0 = glbGetRuleRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239_0 = glbGetRuleRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241_0 = glbGetRuleRatePtr(exp+26, 3);

  double EH2_AD8_Total = 0.0;
  double EH2_AD8_Total_0 = 0.0;

/*********
*  RENO  *
**********/

  double *RENO_U235 = glbGetSignalFitRatePtr(exp+27, 0);
  double *RENO_U238 = glbGetSignalFitRatePtr(exp+27, 1);
  double *RENO_Pu239 = glbGetSignalFitRatePtr(exp+27, 2);
  double *RENO_Pu241 = glbGetSignalFitRatePtr(exp+27, 3);

  double *RENO_U235_0 = glbGetRuleRatePtr(exp+27, 0);
  double *RENO_U238_0 = glbGetRuleRatePtr(exp+27, 1);
  double *RENO_Pu239_0 = glbGetRuleRatePtr(exp+27, 2);
  double *RENO_Pu241_0 = glbGetRuleRatePtr(exp+27, 3);

  double RENO_Total = 0.0;
  double RENO_Total_0 = 0.0;

/*********************************
*  Assembling Daya Bay and RENO  *
**********************************/

  double chi2 = 0.0;
  int i,j;

  double num, denom;

  for (i=0; i<8; i++){
  if (YesNo[24] == 1){
    EH1_AD1_Total = x[0]*DayaBayFs[i][0]*EH1_AD1_U235[0];
    EH1_AD1_Total += x[1]*DayaBayFs[i][1]*EH1_AD1_U238[0];
    EH1_AD1_Total += x[2]*DayaBayFs[i][2]*EH1_AD1_Pu239[0];
    EH1_AD1_Total += x[3]*DayaBayFs[i][3]*EH1_AD1_Pu241[0];

    EH1_AD2_Total = x[0]*DayaBayFs[i][0]*EH1_AD2_U235[0];
    EH1_AD2_Total += x[1]*DayaBayFs[i][1]*EH1_AD2_U238[0];
    EH1_AD2_Total += x[2]*DayaBayFs[i][2]*EH1_AD2_Pu239[0];
    EH1_AD2_Total += x[3]*DayaBayFs[i][3]*EH1_AD2_Pu241[0];

    EH2_AD3_Total = x[0]*DayaBayFs[i][0]*EH2_AD3_U235[0];
    EH2_AD3_Total += x[1]*DayaBayFs[i][1]*EH2_AD3_U238[0];
    EH2_AD3_Total += x[2]*DayaBayFs[i][2]*EH2_AD3_Pu239[0];
    EH2_AD3_Total += x[3]*DayaBayFs[i][3]*EH2_AD3_Pu241[0];

    EH2_AD8_Total = x[0]*DayaBayFs[i][0]*EH2_AD8_U235[0];
    EH2_AD8_Total += x[1]*DayaBayFs[i][1]*EH2_AD8_U238[0];
    EH2_AD8_Total += x[2]*DayaBayFs[i][2]*EH2_AD8_Pu239[0];
    EH2_AD8_Total += x[3]*DayaBayFs[i][3]*EH2_AD8_Pu241[0];

    EH1_AD1_Total_0 = DayaBayFs[i][0]*EH1_AD1_U235_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][1]*EH1_AD1_U238_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][2]*EH1_AD1_Pu239_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][3]*EH1_AD1_Pu241_0[0];

    EH1_AD2_Total_0 = DayaBayFs[i][0]*EH1_AD2_U235_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][1]*EH1_AD2_U238_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][2]*EH1_AD2_Pu239_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][3]*EH1_AD2_Pu241_0[0];

    EH2_AD3_Total_0 = DayaBayFs[i][0]*EH2_AD3_U235_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][1]*EH2_AD3_U238_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][2]*EH2_AD3_Pu239_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][3]*EH2_AD3_Pu241_0[0];

    EH2_AD8_Total_0 = DayaBayFs[i][0]*EH2_AD8_U235_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][1]*EH2_AD8_U238_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][2]*EH2_AD8_Pu239_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][3]*EH2_AD8_Pu241_0[0];

    num = EH1_AD1_Total+ EH1_AD2_Total+ EH2_AD3_Total+ EH2_AD8_Total;
    denom = EH1_AD1_Total_0+ EH1_AD2_Total_0+ EH2_AD3_Total_0+ EH2_AD8_Total_0;

    delta[i+24] = num/denom - DBdataHKSS[i];
  }

  if (YesNo[25] == 1){
    RENO_Total = x[0]*RENOFs[i][0]*RENO_U235[0];
    RENO_Total += x[1]*RENOFs[i][1]*RENO_U238[0];
    RENO_Total += x[2]*RENOFs[i][2]*RENO_Pu239[0];
    RENO_Total += x[3]*RENOFs[i][3]*RENO_Pu241[0];

    RENO_Total_0 = RENOFs[i][0]*RENO_U235_0[0];
    RENO_Total_0 += RENOFs[i][1]*RENO_U238_0[0];
    RENO_Total_0 += RENOFs[i][2]*RENO_Pu239_0[0];
    RENO_Total_0 += RENOFs[i][3]*RENO_Pu241_0[0];

    delta[i+32] = RENO_Total/RENO_Total_0 - RENOdataHKSS[i] ;
  }
  }

/********************************
*  Putting the pieces together  *
*********************************/

  /* Adding the contributions from each experiment */

  for (i=0; i<24; i++){
    for (j=0; j<24; j++){
      chi2 += delta[i]*cov_matrix[i][j]*delta[j];
    }
  }

  /* Now, adding Daya Bay... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+24]*cov_matrix_DB[i][j]*delta[j+24];
    }
  }

  /* Now, adding RENO... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+32]*cov_matrix_RENO[i][j]*delta[j+32];
    }
  }

  /* Adding in the systematic uncertainties on the ratios... */

  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      chi2 += (x[i]-1.0)*(x[j]-1.0)*V_inv_HKSS[i][j]/(errors[i]*errors[j]);
    }
  }

  for (i = 0; i<4; i++){ /* Export nuisance parameters to main file */
    systematic[i] = x[i];
  }
  return chi2;
}

/*************************************************************************
*                                                                        *
*        THE CALCULATION OF A CHI-SQUARED WITH U235 AND Pu239 FREE       *
*                                                                        *
**************************************************************************/

double combo_rate_chi_unfix(int exp, int rule, int n_params, double *x, double *errors,
                 void *user_data)
{ 

/* 
  Here is where I calculate the chi-squared following, for instance, 1703.00860

  "delta" is given by (Rexp - Rpred)

  "user_data" contains the list of experiments being considered in a particular
  invocation of the SBL2 code. 

  NOTE: main.c (SBL v.1) is hard-coded to forcibly include every experiment all the time.
  main2.c (SBL v.2) is not so constrained; changes come from the command line

  Note that there are only two nuisance parameters now -- the rescaling factors for
  U238 and Pu241!

*/

  int k;

  double delta[40] = {0.0};

  double Ratios[2];
  for (k=0; k<2; k++){
    Ratios[k] = ((double *)user_data)[k];
  }

  int YesNo[26];
  for (k=0; k<26; k++){
    YesNo[k] = round( ((double *)user_data)[k+2] );
  }


/************
*  Bugey-4  *
*************/

  double *Bugey_4_U235 = glbGetSignalFitRatePtr(exp, 0);
  double *Bugey_4_U238 = glbGetSignalFitRatePtr(exp, 1);
  double *Bugey_4_Pu239 = glbGetSignalFitRatePtr(exp, 2);
  double *Bugey_4_Pu241 = glbGetSignalFitRatePtr(exp, 3);

  double *Bugey_4_U235_0 = glbGetRuleRatePtr(exp, 0);
  double *Bugey_4_U238_0 = glbGetRuleRatePtr(exp, 1);
  double *Bugey_4_Pu239_0 = glbGetRuleRatePtr(exp, 2);
  double *Bugey_4_Pu241_0 = glbGetRuleRatePtr(exp, 3);
   
  double Bugey_4_RatioN = 0.0;
  double Bugey_4_RatioD = 0.0;

if(YesNo[0]==1||YesNo[2]==1){
  Bugey_4_RatioN = Ratios[0]*BugeyF[0]*Bugey_4_U235[0];
  Bugey_4_RatioN += x[0]*BugeyF[1]*Bugey_4_U238[0];
  Bugey_4_RatioN += Ratios[1]*BugeyF[2]*Bugey_4_Pu239[0];
  Bugey_4_RatioN += x[1]*BugeyF[3]*Bugey_4_Pu241[0];

  Bugey_4_RatioD = BugeyF[0]*Bugey_4_U235_0[0];
  Bugey_4_RatioD += BugeyF[1]*Bugey_4_U238_0[0];
  Bugey_4_RatioD += BugeyF[2]*Bugey_4_Pu239_0[0];
  Bugey_4_RatioD += BugeyF[3]*Bugey_4_Pu241_0[0];
}
if(YesNo[0]==1){
  delta[0] = Rexp[0] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/************
*  Rovno91  *
*************/

  double *Rovno91_U235 = glbGetSignalFitRatePtr(exp+1, 0);
  double *Rovno91_U238 = glbGetSignalFitRatePtr(exp+1, 1);
  double *Rovno91_Pu239 = glbGetSignalFitRatePtr(exp+1, 2);
  double *Rovno91_Pu241 = glbGetSignalFitRatePtr(exp+1, 3);

  double *Rovno91_U235_0 = glbGetRuleRatePtr(exp+1, 0);
  double *Rovno91_U238_0 = glbGetRuleRatePtr(exp+1, 1);
  double *Rovno91_Pu239_0 = glbGetRuleRatePtr(exp+1, 2);
  double *Rovno91_Pu241_0 = glbGetRuleRatePtr(exp+1, 3);

  double Rovno91_RatioN = 0.0;
  double Rovno91_RatioD = 0.0;

if(YesNo[1]==1){
  Rovno91_RatioN = Ratios[0]*Rovno_91F[0]*Rovno91_U235[0];
  Rovno91_RatioN += x[0]*Rovno_91F[1]*Rovno91_U238[0];
  Rovno91_RatioN += Ratios[1]*Rovno_91F[2]*Rovno91_Pu239[0];
  Rovno91_RatioN += x[1]*Rovno_91F[3]*Rovno91_Pu241[0];

  Rovno91_RatioD = Rovno_91F[0]*Rovno91_U235_0[0];
  Rovno91_RatioD += Rovno_91F[1]*Rovno91_U238_0[0];
  Rovno91_RatioD += Rovno_91F[2]*Rovno91_Pu239_0[0];
  Rovno91_RatioD += Rovno_91F[3]*Rovno91_Pu241_0[0];

  delta[1] = Rexp[1] - Rovno91_RatioN/Rovno91_RatioD;
}

/*******************
*  Bugey-3 (15 m)  *
********************/

/*
  No need to recalculate rates -- reuse the output for Bugey-4!
*/

if(YesNo[2]==1){
  delta[2] = Rexp[2] - Bugey_4_RatioN/Bugey_4_RatioD;
}

/*******************
*  Bugey-3 (40 m)  *
********************/

  double *Bugey_3_40_U235 = glbGetSignalFitRatePtr(exp+2, 0);
  double *Bugey_3_40_U238 = glbGetSignalFitRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239 = glbGetSignalFitRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241 = glbGetSignalFitRatePtr(exp+2, 3);
   
  double *Bugey_3_40_U235_0 = glbGetRuleRatePtr(exp+2, 0);
  double *Bugey_3_40_U238_0 = glbGetRuleRatePtr(exp+2, 1);
  double *Bugey_3_40_Pu239_0 = glbGetRuleRatePtr(exp+2, 2);
  double *Bugey_3_40_Pu241_0 = glbGetRuleRatePtr(exp+2, 3);
   
  double Bugey_3_40_RatioN = 0.0;
  double Bugey_3_40_RatioD = 0.0;

if(YesNo[3]==1){   
  Bugey_3_40_RatioN = Ratios[0]*BugeyF[0]*Bugey_3_40_U235[0];
  Bugey_3_40_RatioN += x[0]*BugeyF[1]*Bugey_3_40_U238[0];
  Bugey_3_40_RatioN += Ratios[1]*BugeyF[2]*Bugey_3_40_Pu239[0];
  Bugey_3_40_RatioN += x[1]*BugeyF[3]*Bugey_3_40_Pu241[0];

  Bugey_3_40_RatioD = BugeyF[0]*Bugey_3_40_U235_0[0];
  Bugey_3_40_RatioD += BugeyF[1]*Bugey_3_40_U238_0[0];
  Bugey_3_40_RatioD += BugeyF[2]*Bugey_3_40_Pu239_0[0];
  Bugey_3_40_RatioD += BugeyF[3]*Bugey_3_40_Pu241_0[0];

  delta[3] = Rexp[3] - Bugey_3_40_RatioN/Bugey_3_40_RatioD;
}

/*******************
*  Bugey-3 (95 m)  *
********************/

  double *Bugey_3_95_U235 = glbGetSignalFitRatePtr(exp+3, 0);
  double *Bugey_3_95_U238 = glbGetSignalFitRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239 = glbGetSignalFitRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241 = glbGetSignalFitRatePtr(exp+3, 3);
   
  double *Bugey_3_95_U235_0 = glbGetRuleRatePtr(exp+3, 0);
  double *Bugey_3_95_U238_0 = glbGetRuleRatePtr(exp+3, 1);
  double *Bugey_3_95_Pu239_0 = glbGetRuleRatePtr(exp+3, 2);
  double *Bugey_3_95_Pu241_0 = glbGetRuleRatePtr(exp+3, 3);
   
  double Bugey_3_95_RatioN = 0.0;
  double Bugey_3_95_RatioD = 0.0;

if(YesNo[4]==1){   
  Bugey_3_95_RatioN = Ratios[0]*BugeyF[0]*Bugey_3_95_U235[0];
  Bugey_3_95_RatioN += x[0]*BugeyF[1]*Bugey_3_95_U238[0];
  Bugey_3_95_RatioN += Ratios[1]*BugeyF[2]*Bugey_3_95_Pu239[0];
  Bugey_3_95_RatioN += x[1]*BugeyF[3]*Bugey_3_95_Pu241[0];

  Bugey_3_95_RatioD = BugeyF[0]*Bugey_3_95_U235_0[0];
  Bugey_3_95_RatioD += BugeyF[1]*Bugey_3_95_U238_0[0];
  Bugey_3_95_RatioD += BugeyF[2]*Bugey_3_95_Pu239_0[0];
  Bugey_3_95_RatioD += BugeyF[3]*Bugey_3_95_Pu241_0[0];

  delta[4] = Rexp[4] - Bugey_3_95_RatioN/Bugey_3_95_RatioD;
}

/********************
*  Gosgen (37.9 m)  *
*********************/

  double *Gosgen_38_U235 = glbGetSignalFitRatePtr(exp+4, 0);
  double *Gosgen_38_U238 = glbGetSignalFitRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239 = glbGetSignalFitRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241 = glbGetSignalFitRatePtr(exp+4, 3);
   
  double *Gosgen_38_U235_0 = glbGetRuleRatePtr(exp+4, 0);
  double *Gosgen_38_U238_0 = glbGetRuleRatePtr(exp+4, 1);
  double *Gosgen_38_Pu239_0 = glbGetRuleRatePtr(exp+4, 2);
  double *Gosgen_38_Pu241_0 = glbGetRuleRatePtr(exp+4, 3);
   
  double Gosgen_38_RatioN = 0.0;
  double Gosgen_38_RatioD = 0.0;

if(YesNo[5]==1){
   
  Gosgen_38_RatioN = Ratios[0]*Gosgen38F[0]*Gosgen_38_U235[0];
  Gosgen_38_RatioN += x[0]*Gosgen38F[1]*Gosgen_38_U238[0];
  Gosgen_38_RatioN += Ratios[1]*Gosgen38F[2]*Gosgen_38_Pu239[0];
  Gosgen_38_RatioN += x[1]*Gosgen38F[3]*Gosgen_38_Pu241[0];

  Gosgen_38_RatioD = Gosgen38F[0]*Gosgen_38_U235_0[0];
  Gosgen_38_RatioD += Gosgen38F[1]*Gosgen_38_U238_0[0];
  Gosgen_38_RatioD += Gosgen38F[2]*Gosgen_38_Pu239_0[0];
  Gosgen_38_RatioD += Gosgen38F[3]*Gosgen_38_Pu241_0[0];

  delta[5] = Rexp[5] - Gosgen_38_RatioN/Gosgen_38_RatioD;
}

/********************
*  Gosgen (45.9 m)  *
*********************/

  double *Gosgen_46_U235 = glbGetSignalFitRatePtr(exp+5, 0);
  double *Gosgen_46_U238 = glbGetSignalFitRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239 = glbGetSignalFitRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241 = glbGetSignalFitRatePtr(exp+5, 3);
   
  double *Gosgen_46_U235_0 = glbGetRuleRatePtr(exp+5, 0);
  double *Gosgen_46_U238_0 = glbGetRuleRatePtr(exp+5, 1);
  double *Gosgen_46_Pu239_0 = glbGetRuleRatePtr(exp+5, 2);
  double *Gosgen_46_Pu241_0 = glbGetRuleRatePtr(exp+5, 3);
   
  double Gosgen_46_RatioN = 0.0;
  double Gosgen_46_RatioD = 0.0;

if(YesNo[6]==1){
   
  Gosgen_46_RatioN = Ratios[0]*Gosgen46F[0]*Gosgen_46_U235[0];
  Gosgen_46_RatioN += x[0]*Gosgen46F[1]*Gosgen_46_U238[0];
  Gosgen_46_RatioN += Ratios[1]*Gosgen46F[2]*Gosgen_46_Pu239[0];
  Gosgen_46_RatioN += x[1]*Gosgen46F[3]*Gosgen_46_Pu241[0];

  Gosgen_46_RatioD = Gosgen46F[0]*Gosgen_46_U235_0[0];
  Gosgen_46_RatioD += Gosgen46F[1]*Gosgen_46_U238_0[0];
  Gosgen_46_RatioD += Gosgen46F[2]*Gosgen_46_Pu239_0[0];
  Gosgen_46_RatioD += Gosgen46F[3]*Gosgen_46_Pu241_0[0];

  delta[6] = Rexp[6] - Gosgen_46_RatioN/Gosgen_46_RatioD;
}

/********************
*  Gosgen (64.7 m)  *
*********************/

  double *Gosgen_65_U235 = glbGetSignalFitRatePtr(exp+6, 0);
  double *Gosgen_65_U238 = glbGetSignalFitRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239 = glbGetSignalFitRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241 = glbGetSignalFitRatePtr(exp+6, 3);

  double *Gosgen_65_U235_0 = glbGetRuleRatePtr(exp+6, 0);
  double *Gosgen_65_U238_0 = glbGetRuleRatePtr(exp+6, 1);
  double *Gosgen_65_Pu239_0 = glbGetRuleRatePtr(exp+6, 2);
  double *Gosgen_65_Pu241_0 = glbGetRuleRatePtr(exp+6, 3);
   
  double Gosgen_65_RatioN = 0.0;
  double Gosgen_65_RatioD = 0.0;

if(YesNo[7]==1){
   
  Gosgen_65_RatioN = Ratios[0]*Gosgen65F[0]*Gosgen_65_U235[0];
  Gosgen_65_RatioN += x[0]*Gosgen65F[1]*Gosgen_65_U238[0];
  Gosgen_65_RatioN += Ratios[1]*Gosgen65F[2]*Gosgen_65_Pu239[0];
  Gosgen_65_RatioN += x[1]*Gosgen65F[3]*Gosgen_65_Pu241[0];

  Gosgen_65_RatioD = Gosgen65F[0]*Gosgen_65_U235_0[0];
  Gosgen_65_RatioD += Gosgen65F[1]*Gosgen_65_U238_0[0];
  Gosgen_65_RatioD += Gosgen65F[2]*Gosgen_65_Pu239_0[0];
  Gosgen_65_RatioD += Gosgen65F[3]*Gosgen_65_Pu241_0[0];

  delta[7] = Rexp[7] - Gosgen_65_RatioN/Gosgen_65_RatioD;
}

/*****************
*  ILL (8.76 m)  *
******************/

  double *ILL_U235 = glbGetSignalFitRatePtr(exp+7, 0);
  double *ILL_U235_0 = glbGetRuleRatePtr(exp+7, 0);

if(YesNo[8]==1){
  delta[8] = Rexp[8] - Ratios[0]*ILL_U235[0]/ILL_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (32.8 m)  *
****************************/

  double *Krasnoyarsk_33_U235 = glbGetSignalFitRatePtr(exp+8, 0);
  double *Krasnoyarsk_33_U235_0 = glbGetRuleRatePtr(exp+8, 0);

if(YesNo[9]==1){
  delta[9] = Rexp[9] - Ratios[0]* Krasnoyarsk_33_U235[0]/Krasnoyarsk_33_U235_0[0];
}

/***************************
*  Krasnoyarsk87 (92.3 m)  *
****************************/

  double *Krasnoyarsk_92_U235 = glbGetSignalFitRatePtr(exp+9, 0);
  double *Krasnoyarsk_92_U235_0 = glbGetRuleRatePtr(exp+9, 0);

if(YesNo[10]==1){
  delta[10] = Rexp[10] - Ratios[0]* Krasnoyarsk_92_U235[0]/Krasnoyarsk_92_U235_0[0];
}

/***************************
*  Krasnoyarsk94 (57.0 m)  *
****************************/

  double *Krasnoyarsk_57_U235 = glbGetSignalFitRatePtr(exp+10, 0);
  double *Krasnoyarsk_57_U235_0 = glbGetRuleRatePtr(exp+10, 0);

if(YesNo[11]==1){
  delta[11] = Rexp[11] - Ratios[0]* Krasnoyarsk_57_U235[0]/Krasnoyarsk_57_U235_0[0];
}

/***************************
*  Krasnoyarsk99 (34.0 m)  *
****************************/

  double *Krasnoyarsk_34_U235 = glbGetSignalFitRatePtr(exp+11, 0);
  double *Krasnoyarsk_34_U235_0 = glbGetRuleRatePtr(exp+11, 0);

if(YesNo[12]==1){
  delta[12] = Rexp[12] - Ratios[0]* Krasnoyarsk_34_U235[0]/Krasnoyarsk_34_U235_0[0];
}

/********************
*  SRP-18 (18.2 m)  *
*********************/

  double *SRP_18_U235 = glbGetSignalFitRatePtr(exp+12, 0);
  double *SRP_18_U235_0 = glbGetRuleRatePtr(exp+12, 0);

if(YesNo[13]==1){
  delta[13] = Rexp[13] - Ratios[0]* SRP_18_U235[0]/SRP_18_U235_0[0];
}

/********************
*  SRP-24 (23.8 m)  *
*********************/

  double *SRP_24_U235 = glbGetSignalFitRatePtr(exp+13, 0);
  double *SRP_24_U235_0 = glbGetRuleRatePtr(exp+13, 0);

if(YesNo[14]==1){
  delta[14] = Rexp[14] - Ratios[0]*SRP_24_U235[0]/SRP_24_U235_0[0];
}

/***************
*  Rovno88-1I  *
****************/

  double *Rovno88_1I_U235 = glbGetSignalFitRatePtr(exp+14, 0);
  double *Rovno88_1I_U238 = glbGetSignalFitRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239 = glbGetSignalFitRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241 = glbGetSignalFitRatePtr(exp+14, 3);

  double *Rovno88_1I_U235_0 = glbGetRuleRatePtr(exp+14, 0);
  double *Rovno88_1I_U238_0 = glbGetRuleRatePtr(exp+14, 1);
  double *Rovno88_1I_Pu239_0 = glbGetRuleRatePtr(exp+14, 2);
  double *Rovno88_1I_Pu241_0 = glbGetRuleRatePtr(exp+14, 3);

  double Rovno88_1I_RatioN = 0.0;
  double Rovno88_1I_RatioD = 0.0;

if(YesNo[15]==1){
  Rovno88_1I_RatioN = Ratios[0]* Rovno88_1IF[0]*Rovno88_1I_U235[0];
  Rovno88_1I_RatioN += x[0]* Rovno88_1IF[1]*Rovno88_1I_U238[0];
  Rovno88_1I_RatioN += Ratios[1]* Rovno88_1IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1I_RatioN += x[1]* Rovno88_1IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1I_RatioD = Rovno88_1IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1I_RatioD += Rovno88_1IF[3]*Rovno88_1I_Pu241_0[0];

  delta[15] = Rexp[15] - Rovno88_1I_RatioN/Rovno88_1I_RatioD;
}

/***************
*  Rovno88-2I  *
****************/

  double Rovno88_2I_RatioN = 0.0;
  double Rovno88_2I_RatioD = 0.0;

if(YesNo[16]==1){   
  Rovno88_2I_RatioN = Ratios[0]* Rovno88_2IF[0]*Rovno88_1I_U235[0];
  Rovno88_2I_RatioN += x[0]* Rovno88_2IF[1]*Rovno88_1I_U238[0];
  Rovno88_2I_RatioN += Ratios[1]* Rovno88_2IF[2]*Rovno88_1I_Pu239[0];
  Rovno88_2I_RatioN += x[1]* Rovno88_2IF[3]*Rovno88_1I_Pu241[0];

  Rovno88_2I_RatioD = Rovno88_2IF[0]*Rovno88_1I_U235_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[1]*Rovno88_1I_U238_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2I_RatioD += Rovno88_2IF[3]*Rovno88_1I_Pu241_0[0];

  delta[16] = Rexp[16] - Rovno88_2I_RatioN/Rovno88_2I_RatioD;
}
/***************
*  Rovno88-1S  *
****************/

  double Rovno88_1S_RatioN = 0.0;
  double Rovno88_1S_RatioD = 0.0;

if(YesNo[17]==1){
  Rovno88_1S_RatioN = Ratios[0]* Rovno88_1SF[0]*Rovno88_1I_U235[0];
  Rovno88_1S_RatioN += x[0]* Rovno88_1SF[1]*Rovno88_1I_U238[0];
  Rovno88_1S_RatioN += Ratios[1]* Rovno88_1SF[2]*Rovno88_1I_Pu239[0];
  Rovno88_1S_RatioN += x[1]* Rovno88_1SF[3]*Rovno88_1I_Pu241[0];

  Rovno88_1S_RatioD = Rovno88_1SF[0]*Rovno88_1I_U235_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[1]*Rovno88_1I_U238_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_1S_RatioD += Rovno88_1SF[3]*Rovno88_1I_Pu241_0[0];

  delta[17] = Rexp[17] - Rovno88_1S_RatioN/Rovno88_1S_RatioD;
}

/**********************
*  Rovno88-2S (25 m)  *
***********************/
   
  double Rovno88_2S_25_RatioN = 0.0;
  double Rovno88_2S_25_RatioD = 0.0;

  double *Rovno88_2S_25_U235 = glbGetSignalFitRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238 = glbGetSignalFitRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239 = glbGetSignalFitRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241 = glbGetSignalFitRatePtr(exp+15, 3);

  double *Rovno88_2S_25_U235_0 = glbGetRuleRatePtr(exp+15, 0);
  double *Rovno88_2S_25_U238_0 = glbGetRuleRatePtr(exp+15, 1);
  double *Rovno88_2S_25_Pu239_0 = glbGetRuleRatePtr(exp+15, 2);
  double *Rovno88_2S_25_Pu241_0 = glbGetRuleRatePtr(exp+15, 3);

if(YesNo[18]==1){
  Rovno88_2S_25_RatioN = Ratios[0]*Rovno88_2S_25F[0]*Rovno88_2S_25_U235[0];
  Rovno88_2S_25_RatioN += x[0]*Rovno88_2S_25F[1]*Rovno88_2S_25_U238[0];
  Rovno88_2S_25_RatioN += Ratios[1]*Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239[0];
  Rovno88_2S_25_RatioN += x[1]*Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241[0];

  Rovno88_2S_25_RatioD = Rovno88_2S_25F[0]*Rovno88_2S_25_U235_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[1]*Rovno88_2S_25_U238_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[2]*Rovno88_2S_25_Pu239_0[0];
  Rovno88_2S_25_RatioD += Rovno88_2S_25F[3]*Rovno88_2S_25_Pu241_0[0];

  delta[18] = Rexp[18] - Rovno88_2S_25_RatioN/Rovno88_2S_25_RatioD;
}

/**********************
*  Rovno88-2S (18 m)  *
***********************/

  double Rovno88_2S_18_RatioN = 0.0;
  double Rovno88_2S_18_RatioD = 0.0;

if(YesNo[19]==1){
  Rovno88_2S_18_RatioN = Ratios[0]*Rovno88_2S_18F[0]*Rovno88_1I_U235[0];
  Rovno88_2S_18_RatioN += x[0]*Rovno88_2S_18F[1]*Rovno88_1I_U238[0];
  Rovno88_2S_18_RatioN += Ratios[1]*Rovno88_2S_18F[2]*Rovno88_1I_Pu239[0];
  Rovno88_2S_18_RatioN += x[1]*Rovno88_2S_18F[3]*Rovno88_1I_Pu241[0];

  Rovno88_2S_18_RatioD = Rovno88_2S_18F[0]*Rovno88_1I_U235_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[1]*Rovno88_1I_U238_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[2]*Rovno88_1I_Pu239_0[0];
  Rovno88_2S_18_RatioD += Rovno88_2S_18F[3]*Rovno88_1I_Pu241_0[0];

  delta[19] = Rexp[19] - Rovno88_2S_18_RatioN/Rovno88_2S_18_RatioD;
}

/************
*  Nucifer  *
*************/

  double *Nucifer_U235 = glbGetSignalFitRatePtr(exp+16, 0);
  double *Nucifer_U238 = glbGetSignalFitRatePtr(exp+16, 1);
  double *Nucifer_Pu239 = glbGetSignalFitRatePtr(exp+16, 2);
  double *Nucifer_Pu241 = glbGetSignalFitRatePtr(exp+16, 3);

  double *Nucifer_U235_0 = glbGetRuleRatePtr(exp+16, 0);
  double *Nucifer_U238_0 = glbGetRuleRatePtr(exp+16, 1);
  double *Nucifer_Pu239_0 = glbGetRuleRatePtr(exp+16, 2);
  double *Nucifer_Pu241_0 = glbGetRuleRatePtr(exp+16, 3);

  double Nucifer_RatioN=0.0;
  double Nucifer_RatioD=0.0;

if(YesNo[20]==1){
  Nucifer_RatioN = Ratios[0]*NuciferF[0]*Nucifer_U235[0];
  Nucifer_RatioN += x[0]*NuciferF[1]*Nucifer_U238[0];
  Nucifer_RatioN += Ratios[1]*NuciferF[2]*Nucifer_Pu239[0];
  Nucifer_RatioN += x[1]*NuciferF[3]*Nucifer_Pu241[0];

  Nucifer_RatioD = NuciferF[0]*Nucifer_U235_0[0];
  Nucifer_RatioD += NuciferF[1]*Nucifer_U238_0[0];
  Nucifer_RatioD += NuciferF[2]*Nucifer_Pu239_0[0];
  Nucifer_RatioD += NuciferF[3]*Nucifer_Pu241_0[0];

  delta[20] = Rexp[20] - Nucifer_RatioN/Nucifer_RatioD;
}

/**********************************
*  Palo Verde -- 750 m and 890 m  *
***********************************/

  double PV_RatioN=0.0;
  double PV_RatioD=0.0;

  double *PV_750_U235 = glbGetSignalFitRatePtr(exp+17, 0);
  double *PV_750_U238 = glbGetSignalFitRatePtr(exp+17, 1);
  double *PV_750_Pu239 = glbGetSignalFitRatePtr(exp+17, 2);
  double *PV_750_Pu241 = glbGetSignalFitRatePtr(exp+17, 3);

  double *PV_750_U235_0 = glbGetRuleRatePtr(exp+17, 0);
  double *PV_750_U238_0 = glbGetRuleRatePtr(exp+17, 1);
  double *PV_750_Pu239_0 = glbGetRuleRatePtr(exp+17, 2);
  double *PV_750_Pu241_0 = glbGetRuleRatePtr(exp+17, 3);

  double *PV_890_U235 = glbGetSignalFitRatePtr(exp+18, 0);
  double *PV_890_U238 = glbGetSignalFitRatePtr(exp+18, 1);
  double *PV_890_Pu239 = glbGetSignalFitRatePtr(exp+18, 2);
  double *PV_890_Pu241 = glbGetSignalFitRatePtr(exp+18, 3);

  double *PV_890_U235_0 = glbGetRuleRatePtr(exp+18, 0);
  double *PV_890_U238_0 = glbGetRuleRatePtr(exp+18, 1);
  double *PV_890_Pu239_0 = glbGetRuleRatePtr(exp+18, 2);
  double *PV_890_Pu241_0 = glbGetRuleRatePtr(exp+18, 3);

if(YesNo[21]==1){
  PV_RatioN = Ratios[0]*PVF[0]*(PV_750_U235[0] + PV_890_U235[0]);
  PV_RatioN += x[0]*PVF[1]*(PV_750_U238[0] + PV_890_U238[0]);
  PV_RatioN += Ratios[1]*PVF[2]*(PV_750_Pu239[0] + PV_890_Pu239[0]);
  PV_RatioN += x[1]*PVF[3]*(PV_750_Pu241[0] + PV_890_Pu241[0]);

  PV_RatioD = PVF[0]*(PV_750_U235_0[0]+PV_890_U235_0[0]);
  PV_RatioD += PVF[1]*(PV_750_U238_0[0]+PV_890_U238_0[0]);
  PV_RatioD += PVF[2]*(PV_750_Pu239_0[0]+PV_890_Pu239_0[0]);
  PV_RatioD += PVF[3]*(PV_750_Pu241_0[0]+PV_890_Pu241_0[0]);

  delta[21] = Rexp[21] - PV_RatioN/PV_RatioD;
}

/************************************
*  Double Chooz -- 355 m and 469 m  *
*************************************/

  double DC_RatioN=0.0;
  double DC_RatioD=0.0;

  double *DC_355_U235 = glbGetSignalFitRatePtr(exp+19, 0);
  double *DC_355_U238 = glbGetSignalFitRatePtr(exp+19, 1);
  double *DC_355_Pu239 = glbGetSignalFitRatePtr(exp+19, 2);
  double *DC_355_Pu241 = glbGetSignalFitRatePtr(exp+19, 3);

  double *DC_355_U235_0 = glbGetRuleRatePtr(exp+19, 0);
  double *DC_355_U238_0 = glbGetRuleRatePtr(exp+19, 1);
  double *DC_355_Pu239_0 = glbGetRuleRatePtr(exp+19, 2);
  double *DC_355_Pu241_0 = glbGetRuleRatePtr(exp+19, 3);

  double *DC_469_U235 = glbGetSignalFitRatePtr(exp+20, 0);
  double *DC_469_U238 = glbGetSignalFitRatePtr(exp+20, 1);
  double *DC_469_Pu239 = glbGetSignalFitRatePtr(exp+20, 2);
  double *DC_469_Pu241 = glbGetSignalFitRatePtr(exp+20, 3);

  double *DC_469_U235_0 = glbGetRuleRatePtr(exp+20, 0);
  double *DC_469_U238_0 = glbGetRuleRatePtr(exp+20, 1);
  double *DC_469_Pu239_0 = glbGetRuleRatePtr(exp+20, 2);
  double *DC_469_Pu241_0 = glbGetRuleRatePtr(exp+20, 3);

if(YesNo[22]==1){
  DC_RatioN = Ratios[0]*DCF[0]*(DC_355_U235[0] + DC_469_U235[0]);
  DC_RatioN += x[0]*DCF[1]*(DC_355_U238[0] + DC_469_U238[0]);
  DC_RatioN += Ratios[1]*DCF[2]*(DC_355_Pu239[0] + DC_469_Pu239[0]);
  DC_RatioN += x[1]*DCF[3]*(DC_355_Pu241[0] + DC_469_Pu241[0]);

  DC_RatioD = DCF[0]*(DC_355_U235_0[0] + DC_469_U235_0[0]);
  DC_RatioD += DCF[1]*(DC_355_U238_0[0] + DC_469_U238_0[0]);
  DC_RatioD += DCF[2]*(DC_355_Pu239_0[0] + DC_469_Pu239_0[0]);
  DC_RatioD += DCF[3]*(DC_355_Pu241_0[0] + DC_469_Pu241_0[0]);

  delta[22] = Rexp[22] - DC_RatioN/DC_RatioD;
}

/******************************
*  Chooz -- 998 m and 1115 m  *
*******************************/

  double Chooz_RatioN=0.0;
  double Chooz_RatioD=0.0;

  double *Chooz_998_U235 = glbGetSignalFitRatePtr(exp+21, 0);
  double *Chooz_998_U238 = glbGetSignalFitRatePtr(exp+21, 1);
  double *Chooz_998_Pu239 = glbGetSignalFitRatePtr(exp+21, 2);
  double *Chooz_998_Pu241 = glbGetSignalFitRatePtr(exp+21, 3);

  double *Chooz_998_U235_0 = glbGetRuleRatePtr(exp+21, 0);
  double *Chooz_998_U238_0 = glbGetRuleRatePtr(exp+21, 1);
  double *Chooz_998_Pu239_0 = glbGetRuleRatePtr(exp+21, 2);
  double *Chooz_998_Pu241_0 = glbGetRuleRatePtr(exp+21, 3);

  double *Chooz_1115_U235 = glbGetSignalFitRatePtr(exp+22, 0);
  double *Chooz_1115_U238 = glbGetSignalFitRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239 = glbGetSignalFitRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241 = glbGetSignalFitRatePtr(exp+22, 3);

  double *Chooz_1115_U235_0 = glbGetRuleRatePtr(exp+22, 0);
  double *Chooz_1115_U238_0 = glbGetRuleRatePtr(exp+22, 1);
  double *Chooz_1115_Pu239_0 = glbGetRuleRatePtr(exp+22, 2);
  double *Chooz_1115_Pu241_0 = glbGetRuleRatePtr(exp+22, 3);

if(YesNo[23]==1){
  Chooz_RatioN = Ratios[0]*ChoozF[0]*(Chooz_998_U235[0] + Chooz_1115_U235[0]);
  Chooz_RatioN += x[0]*ChoozF[1]*(Chooz_998_U238[0] + Chooz_1115_U238[0]);
  Chooz_RatioN += Ratios[1]*ChoozF[2]*(Chooz_998_Pu239[0] + Chooz_1115_Pu239[0]);
  Chooz_RatioN += x[1]*ChoozF[3]*(Chooz_998_Pu241[0] + Chooz_1115_Pu241[0]);

  Chooz_RatioD = ChoozF[0]*(Chooz_998_U235_0[0] + Chooz_1115_U235_0[0]);
  Chooz_RatioD += ChoozF[1]*(Chooz_998_U238_0[0] + Chooz_1115_U238_0[0]);
  Chooz_RatioD += ChoozF[2]*(Chooz_998_Pu239_0[0] + Chooz_1115_Pu239_0[0]);
  Chooz_RatioD += ChoozF[3]*(Chooz_998_Pu241_0[0] + Chooz_1115_Pu241_0[0]);

  delta[23] = Rexp[23] - Chooz_RatioN/Chooz_RatioD;
}

/************
*  EH1 AD1  *
*************/

  double *EH1_AD1_U235 = glbGetSignalFitRatePtr(exp+23, 0);
  double *EH1_AD1_U238 = glbGetSignalFitRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239 = glbGetSignalFitRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241 = glbGetSignalFitRatePtr(exp+23, 3);

  double *EH1_AD1_U235_0 = glbGetRuleRatePtr(exp+23, 0);
  double *EH1_AD1_U238_0 = glbGetRuleRatePtr(exp+23, 1);
  double *EH1_AD1_Pu239_0 = glbGetRuleRatePtr(exp+23, 2);
  double *EH1_AD1_Pu241_0 = glbGetRuleRatePtr(exp+23, 3);
   
  double EH1_AD1_Total = 0.0;
  double EH1_AD1_Total_0 = 0.0;

/************
*  EH1 AD2  *
*************/

  double *EH1_AD2_U235 = glbGetSignalFitRatePtr(exp+24, 0);
  double *EH1_AD2_U238 = glbGetSignalFitRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239 = glbGetSignalFitRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241 = glbGetSignalFitRatePtr(exp+24, 3);

  double *EH1_AD2_U235_0 = glbGetRuleRatePtr(exp+24, 0);
  double *EH1_AD2_U238_0 = glbGetRuleRatePtr(exp+24, 1);
  double *EH1_AD2_Pu239_0 = glbGetRuleRatePtr(exp+24, 2);
  double *EH1_AD2_Pu241_0 = glbGetRuleRatePtr(exp+24, 3);

  double EH1_AD2_Total = 0.0;
  double EH1_AD2_Total_0 = 0.0;

/************
*  EH2 AD3  *
*************/

  double *EH2_AD3_U235 = glbGetSignalFitRatePtr(exp+25, 0);
  double *EH2_AD3_U238 = glbGetSignalFitRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239 = glbGetSignalFitRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241 = glbGetSignalFitRatePtr(exp+25, 3);

  double *EH2_AD3_U235_0 = glbGetRuleRatePtr(exp+25, 0);
  double *EH2_AD3_U238_0 = glbGetRuleRatePtr(exp+25, 1);
  double *EH2_AD3_Pu239_0 = glbGetRuleRatePtr(exp+25, 2);
  double *EH2_AD3_Pu241_0 = glbGetRuleRatePtr(exp+25, 3);
   
  double EH2_AD3_Total = 0.0;
  double EH2_AD3_Total_0 = 0.0;

/************
*  EH2 AD8  *
*************/

  double *EH2_AD8_U235 = glbGetSignalFitRatePtr(exp+26, 0);
  double *EH2_AD8_U238 = glbGetSignalFitRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239 = glbGetSignalFitRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241 = glbGetSignalFitRatePtr(exp+26, 3);

  double *EH2_AD8_U235_0 = glbGetRuleRatePtr(exp+26, 0);
  double *EH2_AD8_U238_0 = glbGetRuleRatePtr(exp+26, 1);
  double *EH2_AD8_Pu239_0 = glbGetRuleRatePtr(exp+26, 2);
  double *EH2_AD8_Pu241_0 = glbGetRuleRatePtr(exp+26, 3);

  double EH2_AD8_Total = 0.0;
  double EH2_AD8_Total_0 = 0.0;

/*********
*  RENO  *
**********/

  double *RENO_U235 = glbGetSignalFitRatePtr(exp+27, 0);
  double *RENO_U238 = glbGetSignalFitRatePtr(exp+27, 1);
  double *RENO_Pu239 = glbGetSignalFitRatePtr(exp+27, 2);
  double *RENO_Pu241 = glbGetSignalFitRatePtr(exp+27, 3);

  double *RENO_U235_0 = glbGetRuleRatePtr(exp+27, 0);
  double *RENO_U238_0 = glbGetRuleRatePtr(exp+27, 1);
  double *RENO_Pu239_0 = glbGetRuleRatePtr(exp+27, 2);
  double *RENO_Pu241_0 = glbGetRuleRatePtr(exp+27, 3);

  double RENO_Total = 0.0;
  double RENO_Total_0 = 0.0;

/*********************************
*  Assembling Daya Bay and RENO  *
**********************************/

  double chi2 = 0.0;
  int i,j;

  double num, denom;

  for (i=0; i<8; i++){
  if (YesNo[24] == 1){
    EH1_AD1_Total = Ratios[0]*DayaBayFs[i][0]*EH1_AD1_U235[0];
    EH1_AD1_Total += x[0]*DayaBayFs[i][1]*EH1_AD1_U238[0];
    EH1_AD1_Total += Ratios[1]*DayaBayFs[i][2]*EH1_AD1_Pu239[0];
    EH1_AD1_Total += x[1]*DayaBayFs[i][3]*EH1_AD1_Pu241[0];

    EH1_AD2_Total = Ratios[0]*DayaBayFs[i][0]*EH1_AD2_U235[0];
    EH1_AD2_Total += x[0]*DayaBayFs[i][1]*EH1_AD2_U238[0];
    EH1_AD2_Total += Ratios[1]*DayaBayFs[i][2]*EH1_AD2_Pu239[0];
    EH1_AD2_Total += x[1]*DayaBayFs[i][3]*EH1_AD2_Pu241[0];

    EH2_AD3_Total = Ratios[0]*DayaBayFs[i][0]*EH2_AD3_U235[0];
    EH2_AD3_Total += x[0]*DayaBayFs[i][1]*EH2_AD3_U238[0];
    EH2_AD3_Total += Ratios[1]*DayaBayFs[i][2]*EH2_AD3_Pu239[0];
    EH2_AD3_Total += x[1]*DayaBayFs[i][3]*EH2_AD3_Pu241[0];

    EH2_AD8_Total = Ratios[0]*DayaBayFs[i][0]*EH2_AD8_U235[0];
    EH2_AD8_Total += x[0]*DayaBayFs[i][1]*EH2_AD8_U238[0];
    EH2_AD8_Total += Ratios[1]*DayaBayFs[i][2]*EH2_AD8_Pu239[0];
    EH2_AD8_Total += x[1]*DayaBayFs[i][3]*EH2_AD8_Pu241[0];

    EH1_AD1_Total_0 = DayaBayFs[i][0]*EH1_AD1_U235_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][1]*EH1_AD1_U238_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][2]*EH1_AD1_Pu239_0[0];
    EH1_AD1_Total_0 += DayaBayFs[i][3]*EH1_AD1_Pu241_0[0];

    EH1_AD2_Total_0 = DayaBayFs[i][0]*EH1_AD2_U235_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][1]*EH1_AD2_U238_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][2]*EH1_AD2_Pu239_0[0];
    EH1_AD2_Total_0 += DayaBayFs[i][3]*EH1_AD2_Pu241_0[0];

    EH2_AD3_Total_0 = DayaBayFs[i][0]*EH2_AD3_U235_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][1]*EH2_AD3_U238_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][2]*EH2_AD3_Pu239_0[0];
    EH2_AD3_Total_0 += DayaBayFs[i][3]*EH2_AD3_Pu241_0[0];

    EH2_AD8_Total_0 = DayaBayFs[i][0]*EH2_AD8_U235_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][1]*EH2_AD8_U238_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][2]*EH2_AD8_Pu239_0[0];
    EH2_AD8_Total_0 += DayaBayFs[i][3]*EH2_AD8_Pu241_0[0];

    num = EH1_AD1_Total+ EH1_AD2_Total+ EH2_AD3_Total+ EH2_AD8_Total;
    denom = EH1_AD1_Total_0+ EH1_AD2_Total_0+ EH2_AD3_Total_0+ EH2_AD8_Total_0;

    delta[i+24] = num/denom - DBdata[i];
  }

  if (YesNo[25] == 1){
    RENO_Total = Ratios[0]*RENOFs[i][0]*RENO_U235[0];
    RENO_Total += x[0]*RENOFs[i][1]*RENO_U238[0];
    RENO_Total += Ratios[1]*RENOFs[i][2]*RENO_Pu239[0];
    RENO_Total += x[1]*RENOFs[i][3]*RENO_Pu241[0];

    RENO_Total_0 = RENOFs[i][0]*RENO_U235_0[0];
    RENO_Total_0 += RENOFs[i][1]*RENO_U238_0[0];
    RENO_Total_0 += RENOFs[i][2]*RENO_Pu239_0[0];
    RENO_Total_0 += RENOFs[i][3]*RENO_Pu241_0[0];

    delta[i+32] = RENO_Total/RENO_Total_0 - RENOdata[i] ;
  }
  }

/********************************
*  Putting the pieces together  *
*********************************/

  /* Adding the contributions from each experiment */

  for (i=0; i<24; i++){
    for (j=0; j<24; j++){
      chi2 += delta[i]*cov_matrix[i][j]*delta[j];
    }
  }

  /* Now, adding Daya Bay... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+24]*cov_matrix_DB[i][j]*delta[j+24];
    }
  }

  /* Now, adding RENO... */

  for (i=0; i<8; i++){
    for (j=0; j<8; j++){
      chi2 += delta[i+32]*cov_matrix_RENO[i][j]*delta[j+32];
    }
  }

/*
  Adding in the systematic uncertainties on the ratios...this is easy here, because
  the U238 and Pu241 fluxes are uncorrelated, per arXiv:1703.00860
*/

  for (i=0; i<2; i++){
      chi2 += (x[i]-1.0)*(x[i]-1.0) / (errors[i]*errors[i]);
  }

  return chi2;
}  


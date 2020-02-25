/* (C) 2019 Jeffrey M. Berryman */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_expint.h>
#include <globes/globes.h>   /* GLoBES library */

#include "rate_combo.h"
#include "rate_funcs.h"
#include "rate_data.h"

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


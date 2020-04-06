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

#ifndef GLF_SPECTRUM_AUX_2_H
#define GLF_SPECTRUM_AUX_2_H

/****************************************
*  INPUTS FOR THE SPECTRAL MEASUREMENTS *
*****************************************/

  /**********
  *  DANSS  *
  **********/

/*
  These are the inputs relevant for each experiment, including various ratios and
  uncertainties on these ratios; detail on each object are given below.
*/

/* Bottom-top ratio at DANSS */
static const double DANSSratio[24] = {0.7175, 0.7193, 0.7272, 0.7130, 0.7114, 0.7080, 0.7108, 0.7096, 0.6915, 0.7008, 0.7138, 0.6994, 0.7044, 0.6910, 0.6889, 0.7139, 0.7147, 0.7377, 0.7244, 0.7553, 0.7483, 0.6990, 0.6907, 0.6722};

/*
  These are arrays containing the fuel composition of each of the reactors, ordered:
  { U235, U238, Pu239, Pu241}

  DANSS: estimate from Kopeikin for VVER-1000; results shouldn't depend strongly on this
  assumption, given that we are taking a ratio.
*/

static const double DANSSF[4] = {0.56, 0.07, 0.31, 0.06};

/*
  To estimate the effect of updated DANSS data on these results, we include a pre-
  publication version of their results (see Danilov @ EPS-HEP 2019)
*/

static const double DANSSratio2[48] = {0.707, 0.722, 0.708, 0.705, 0.707, 0.718, 0.704, 
		0.712, 0.695, 0.706, 0.697, 0.689, 0.7, 0.695, 0.702, 0.708, 0.691, 0.701, 
		0.701, 0.711, 0.703, 0.699, 0.694, 0.701, 0.685, 0.678, 0.704, 0.685, 
		0.692, 0.686, 0.709, 0.718, 0.704, 0.697, 0.687, 0.717, 0.728, 0.726, 
		0.694, 0.669, 0.704, 0.712, 0.75, 0.712, 0.746, 0.681, 0.668, 0.615};

static const double DANSSF2[4] = {0.542, 0.072, 0.328, 0.058};

  /*************
  *  DAYA BAY  *
  **************/

/*
  The ratios of estimated number of events, N2/N1 and N3/N1. The binning is nonuniform:
	- The middle 24 bins are all 0.25 MeV from 1.3 to 7.3 MeV (prompt energy).
	- The first bin is 0.7-1.3 MeV.
	- The last bin is 7.3-12.0 MeV.
*/

static const double DayaBayRatio21[26] = {0.926932, 0.928954, 0.918739, 0.939571, 0.93593, 0.927909, 0.938499, 0.934976, 0.938284, 0.939889, 0.929927, 0.939565, 0.939061, 0.936539, 0.936932, 0.940709, 0.938614, 0.946836, 0.93802, 0.94908, 0.927143, 0.949198, 0.937487, 0.950238, 0.963763, 0.955083};

static const double DayaBayRatio31[26] = {0.288691, 0.287703, 0.276019, 0.277756, 0.275704, 0.273235, 0.271906, 0.271127, 0.274559, 0.272144, 0.271828, 0.275376, 0.27577, 0.276768, 0.275629, 0.280527, 0.279364, 0.281164, 0.283437, 0.281031, 0.268344, 0.27921, 0.276389, 0.293173, 0.271901, 0.273236};

/*
  The matrix "DayaBayCovariance" lives in "spectra_aux2.h"
*/


/*
  These are arrays containing the fuel fractions for each Daya Bay AD, ordered:
  { U235, U238, Pu239, Pu241}

  The number here is the AD number (8 is in EH2!)

  These are from Table 9 of arXiv:1607.05378
*/

static const double DayaBayF1[4] = {0.5678, 0.0761, 0.30075, 0.05545};
static const double DayaBayF2[4] = {0.56605, 0.07615, 0.3021, 0.0556};
static const double DayaBayF3[4] = {0.5618, 0.0761, 0.30665, 0.0553};
static const double DayaBayF8[4] = {0.56345, 0.076, 0.3052, 0.05555};

static const double DayaBayF4[4] = {0.559, 0.076, 0.310, 0.055};
static const double DayaBayF5[4] = {0.559, 0.076, 0.310, 0.055};
static const double DayaBayF6[4] = {0.559, 0.076, 0.310, 0.055};
static const double DayaBayF7[4] = {0.552, 0.076, 0.315, 0.057};

  /*********
  *  NEOS  *
  **********/

static const double NEOS_Ratio[60] = {0.978, 1.04, 1.031, 0.961, 0.958, 0.963, 0.963, 0.973, 0.989, 0.968, 0.986, 0.993, 1.003, 0.987, 1.015, 1.016, 0.983, 0.993, 0.993, 1.003, 1.015, 1.034, 1.023, 1.002, 1.002, 1., 0.995, 1.019, 1.013, 0.986, 0.994, 0.996, 1.004, 1.014, 1.015, 1.004, 1.031, 1.026, 1.011, 1.042, 1.013, 1.024, 1.009, 1.005, 0.975, 0.991, 0.978, 0.965, 1.039, 1.022, 0.979, 1.014, 1.002, 1.051, 1.04, 1.031, 1.017, 1.037, 0.997, 1.051};

/*
  The matrix "NEOS_Covariance" lives in "spectra_aux2.h"
*/

static const double NEOSF[4] = {0.655, 0.072, 0.235, 0.038};

  /*****************
  *  DOUBLE CHOOZ  *
  ******************/

static const double DC_Ratio[26] = {0.425358, 0.422108, 0.430931, 0.435575, 0.438825, 0.437896, 0.445791, 0.44254, 0.441147, 0.444397, 0.444862, 0.434181, 0.443933, 0.450899, 0.455542, 0.463436, 0.455078, 0.470402, 0.45322, 0.458793, 0.451363, 0.44997, 0.465294, 0.474581, 0.469938, 0.520089};

static const double DCF[4] = {0.511, 0.087, 0.340, 0.062};

  /**********
  *  BUGEY  *
  **********/

static const double Bugey_Ratio[25] = {0.137, 0.143, 0.139, 0.141, 0.141, 0.14, 0.138, 0.145, 0.14, 0.135, 0.138, 0.144, 0.144, 0.126, 0.139, 0.146, 0.133, 0.131, 0.144, 0.144, 0.149, 0.141, 0.141, 0.138, 0.131};

static const double BugeyF[4] = {0.538, 0.078, 0.328, 0.056};

  /*********
  *  RENO  *
  **********/

static const double RENO_Ratio[25] = {0.12592, 0.120914, 0.11629, 0.12038, 0.1154, 0.117704, 0.115853, 0.116933, 0.117296, 0.115695, 0.118126, 0.118856, 0.120047, 0.12019, 0.117643, 0.118033, 0.118237, 0.119066, 0.121666, 0.119507, 0.115003, 0.122545, 0.119143, 0.121214, 0.124873};

static const double RENO_ND_F[4] = {0.57264, 0.07309, 0.29911, 0.055161};
static const double RENO_FD_F[4] = {0.57447, 0.07340, 0.29742, 0.054752};

#endif /* !GLF_SPECTRUM_AUX_2_H */

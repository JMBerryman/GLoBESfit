#ifndef RATE_DATA_H
#define RATE_DATA_H

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

#endif /* !RATE_DATA_H */

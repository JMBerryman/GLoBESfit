#ifndef SBL_FUNCS_H
#define SBL_FUNCS_H

#define YES 1
#define NO -1

#define FAILURE -1
#define SUCCESS 0

#define MY_THETA_X 1
#define MY_DELTA_M 2

/**************************************************************************/

/* Functions for DANSS and NEOS */

int DANSS_up_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DANSS_down_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int NEOS_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

/**************************************************************************/

/* Functions for Bugey */

int bugey_15_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int bugey_40_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

/**************************************************************************/

#endif /* !SBL_FUNCS_H */

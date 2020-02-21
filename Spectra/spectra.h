#ifndef SPECTRA_H
#define SPECTRA_H

#define YES 1
#define NO -1

#define FAILURE -1
#define SUCCESS 0

#define MY_THETA_13 0
#define MY_THETA_X 1
#define MY_DELTA_M 2
#define MY_DELTA_ATM 3

/**************************************************************************/

/* Setting and getting parameters */

int spectra_set_oscillation_parameters(glb_params p, void *user_data);
int spectra_get_oscillation_parameters(glb_params p, void *user_data);

/**************************************************************************/

/* The overall chi-squared function */

double spectra_chi(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

#endif /* !SPECTRA_H */

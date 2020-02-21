#ifndef RATE_COMBO_H
#define RATE_COMBO_H

#define YES 1
#define NO -1

#define FAILURE -1
#define SUCCESS 0

#define MY_THETA_13 0
#define MY_THETA_X 1
#define MY_DELTA_M 2
#define MY_DELTA_ATM 3

/**************************************************************************/

/* Various chi-squared functions */

double combo_rate_chi(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

double combo_rate_chi_nosys(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

double combo_rate_chi_SM(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

double combo_rate_chi_HKSS(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

double combo_rate_chi_unfix(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

#endif /* !RATE_COMBO_H */

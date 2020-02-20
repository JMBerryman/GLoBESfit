#ifndef RATE_FUNCS_H
#define RATE_FUNCS_H

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

int combo_set_oscillation_parameters(glb_params p, void *user_data);
int combo_get_oscillation_parameters(glb_params p, void *user_data);

/**************************************************************************/

/* Generic probability matrices */

int combo_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int standard_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

/**************************************************************************/

/* Functions for Daya Bay near detectors */

int DB_EH1_AD1_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DB_EH1_AD2_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DB_EH2_AD3_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DB_EH2_AD8_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

/**************************************************************************/

/* Functions for RENO */

int RENO_ND_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int RENO_FD_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);


/**************************************************************************/

#endif /* !RATE_FUNCS_H */

#ifndef LBL_FUNCS
#define LBL_FUNCS

#define YES 1
#define NO -1

#define FAILURE -1
#define SUCCESS 0

#define MY_THETA_13 0
#define MY_THETA_X 1
#define MY_DELTA_M 2
#define MY_DELTA_ATM 3

/**************************************************************************/

/* Functions for Daya Bay */

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

int DB_EH3_AD4_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DB_EH3_AD5_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DB_EH3_AD6_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DB_EH3_AD7_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

/**************************************************************************/

/* Functions for Double Chooz */

int DC_ND_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int DC_FD_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
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

#endif /* !LBL_FUNCS */

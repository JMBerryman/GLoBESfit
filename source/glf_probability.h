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

#ifndef GLF_PROBABILITY_H
#define GLF_PROBABILITY_H

#define YES 1
#define NO -1

#define FAILURE -1
#define SUCCESS 0

#define MY_THETA_13 0
#define MY_THETA_14 1
#define MY_DELTA_M 2
#define MY_DELTA_ATM 3
#define MY_THETA_12 4
#define MY_DELTA_SOLAR 5

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

int glf_four_state_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);

int glf_two_state_probability_matrix(double P[3][3], int cp_sign, double E, int psteps,
                          const double *length, const double *density,
				 double filter_sigma, void *user_data);


/**************************************************************************/

#endif /* !GLF_PROBABILITY_H */

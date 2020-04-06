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

#ifndef GLF_RATE_CHI_H
#define GLF_RATE_CHI_H

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

#endif /* !GLF_RATE_CHI_H */

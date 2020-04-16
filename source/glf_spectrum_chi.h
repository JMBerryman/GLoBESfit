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

#ifndef GLF_SPECTRUM_H
#define GLF_SPECTRUM_H 1

/**************************************************************************/

/* The overall chi-squared function */

double glf_spectra_chi(int exp, int rule, int n_params, double *x, double *errors,
		void *user_data);

#endif /* !GLF_SPECTRUM_H */

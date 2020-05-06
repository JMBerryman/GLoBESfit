GLoBESfit - GLOBES fitting tools
===================================================

This is GLoBESfit, a collection of numerical routines based on GLoBES
to allow neutrino oscillation fits to global neutrino data.

GLoBESfit is free software, you can redistribute it and/or modify it
under the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Limitations
===========

GLoBESfit is *not* thread safe.
GLoBESfit is *not* designed to be used in *security critical* applications 
and environments, e.g. no care has been taken to avoid vulnerability
to heap corruption exploits etc.

Credit
======

GLoBESfit is mainly developed for academic use. Thus the authors would
appreciate being given academic credit for it. Whenever you use
GLoBESfit to produce a publication or a talk indicate that you have
used GLoBESfit and please cite the following references:

J.M. Berryman and P. Huber, Sterile Neutrinos and the Global Reactor
Antineutrino Dataset
arXiv: 2005.01756

Besides, many of the data which is used by GLoBESfit and distributed together
with it should be correctly referenced where possible. For details see the
documentation.

Installation
============

GLoBESfit follows the standard GLoBES installation procedure. That is,
it is assumed that GLoBES has been successfully installed on your
system and thus all dependencies of GLoBESfit should be satisfied. Go
into 'source/' and in 'Makefile' make sure that 'prefix' points to a
location where 'globes-config' can be found. Then use 'make'.

More information about GLoBES
=============================

The project homepage is http://www.globesfit.org

See the NEWS file for recent changes to the library.

Reporting Bugs
==============

A list of known bugs can be found in the BUGS file.  
Details of compilation problems can be found in the INSTALL file.

If you find a bug which is not listed in these files please report it
to <globesfit-g@vt.edu>.

All bug reports should include:

       The version number of GLoBESfit, and where you obtained it.
       The version number of GLoBES.	   
       The hardware and operating system
       The compiler used, including version number and compilation options
       A description of the bug behavior
       A short program which reproducibly exercises the bug

It is also useful if you can report whether the same problem occurs
when the library is compiled without optimization.

Thank you.

Contributing to GLoBESfit
=========================

If you are interested in participating in GLoBESfit development, please contact
us at <globesfit-g@vt.edu>
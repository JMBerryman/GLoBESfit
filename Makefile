# GLoBESfit -- GLoBES fitting tools
# (C) 2019-2020 The GLoBESfit Team
#
# GLoBESfit is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBESfit or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# This Makefile is adapted from the globes examples Makefile
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.

# Absolutely requires globes 3.2.17.0 or higher

prefix = /usr/local
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O3 -Wall -W -w
INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
lapack = -llapack -lblas -lf2c

BIN = glf_rate glf_spectrum
OBJ = source/glf_rate.o source/glf_spectrum.o source/glf_rate_chi.o source/glf_spectrum_chi.o source/glf_probability.o  

all: $(BIN)

glf_rate: source/glf_rate_chi.o source/glf_rate.o source/glf_probability.o
	g++ glf_rate_chi.o glf_rate.o glf_probability.o -o\
	 glf_rate $(LDFLAGS)  $(local_LDFLAGS)
glf_spectrum: source/glf_spectrum.o source/glf_spectrum_chi.o source/glf_probability.o
	g++ glf_spectrum.o glf_spectrum_chi.o  glf_probability.o -o\
	 glf_spectrum $(LDFLAGS)  $(local_LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS)  -c $< $(INCFLAGS) -DGITVERSION=\"$(gitversion)\"
.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)

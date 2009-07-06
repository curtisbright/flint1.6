#######################################################################
#   This file is part of FLINT.
#
#   FLINT is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   FLINT is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with FLINT; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

#******************************************************************************
#
#  default configuration file for compare.py
#
#  Please don't edit this file; make a copy and then tell compare.py about it
#  by running e.g.
#
#      python compare.py -c my_config.py
#
#******************************************************************************


#------------------------------------------------------------------------------
# dots per inch for output file
# larger values make a bigger graphics file

CONFIG_dpi = 96


#------------------------------------------------------------------------------
# radius of each dot (don't ask me what the units are, I have no idea)

CONFIG_dotsize = 48 


#------------------------------------------------------------------------------
# title at the top of the image

CONFIG_title = "FLINT vs Magma Z[x] GCD (Scale = 20)"


#------------------------------------------------------------------------------
# ignore all data points whose minimum and maximum times differ by a factor
# of more than CONFIG_tolerance

CONFIG_tolerance = 1000.0


#------------------------------------------------------------------------------
# xscale and yscale determine the scaling used on each axis.
# Possible values are None, "log10", "log2"

CONFIG_xscale = "log2" 
CONFIG_yscale = "log2"


#------------------------------------------------------------------------------
# text labels for each axis (None to leave blank)

CONFIG_xlabel = "Length" 
CONFIG_ylabel = "Bits" 


#------------------------------------------------------------------------------
# if the values from the two data sets differ by a factor of more than
# CONFIG_truncate, the maximum colour intensity is assigned

CONFIG_truncate = 20.0


#------------------------------------------------------------------------------
# The minimum colour intensity, in the range (0, 1).
# If this is zero, then negative through positive values are plotted on
# a continuous spectrum from blue through red. If this is 0.5, then a tiny
# positive value is already somewhat brown, and a tiny negative value is
# already somewhat blue. If this is 1.0, you get all flat blue and all flat
# red.

CONFIG_min_intensity = 0.2


################# end of config file

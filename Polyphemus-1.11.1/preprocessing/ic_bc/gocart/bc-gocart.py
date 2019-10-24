#!/usr/bin/env python

# Copyright (C) 2007, ENPC - INRIA - EDF R&D
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# Polyphemus is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Polyphemus is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the Polyphemus web site:
#      http://cerea.enpc.fr/polyphemus/

#
# This script processes Gocart raw data files into Polair3D input files.
# Gocart concentrations are expected to be daily averages.
#

import os
import sys
import datetime
from optparse import OptionParser
import ConfigParser

#################
# CONFIGURATION #

# List all species. Warning: the order is crucial; in case you need
# to change it, please contact polyphemus-help@lists.gforge.inria.fr.
species_list = ["CC.STD.tv15", "DU.STD.tv20", "SU.STD.tv15", "SS.STD.tv17"]

# Put your file format, where %Y is the year, %M is the month
# and %S is the species name.
file_name_format = "/u/cergrene/Y/fahey/CHIN-BCs/%Y%M.%S.g.day"

# Configuration files.
cfg_general_default = "../general.cfg"

# configuration files corresponding to species
# in "species_list" (same order !).
cfg_species_list = ["bc-gocart-CC.cfg", "bc-gocart-DU.cfg",
                    "bc-gocart-SU.cfg", "bc-gocart-SS.cfg"]

# configuration file for nh4.
cfg_nh4 = "bc-nh4.cfg"

# bc-gocart related program
bc_gocart_prog = "./bc-gocart"
bc_nh4_prog = "./bc-nh4"

# CONFIGURATION #
#################


##########
# SCRIPT #
usage = "%prog [options] general_config_file beg_date number_of_days \n\n" \
    + " Computes aerosol Polair3D boundary conditions" \
    + " from Gocart data.\n\n" \
    + " general_config_file : general configuration file, optional, default is " \
    + cfg_general_default + " .\n" \
    + " beg_date            : beginning date in the format YYYY-MM-DD.\n" \
    + " number_of_days      : number of days to perform."
version= "%prog 1.0"

parser = OptionParser(usage=usage,version=version)

parser.add_option("-n", "--dry-run", action="count", \
          dest="dry_run", default=0, help="says what it would do," \
          + " but actually does NOT do it.")
parser.add_option("--nh4", action="store_true", \
          dest="nh4", default=False, help="compute ammonia" \
          + " which is not directly provided by Gocart," \
          + " default is to not compute.")

(options, args) = parser.parse_args()

# Test number of options
if len(args) != 3 and len(args) != 2:
    os.system(sys.argv[0] + " -h")
    sys.exit(2)

# General config file
if len(args) == 3:
      cfg_general = args[0]
elif len(args) == 2:
    cfg_general = cfg_general_default

if os.path.exists(cfg_general):
    sys.stderr.write("using general config file " + cfg_general + "\n")
else:
    sys.stderr.write("general config file " + cfg_general + " NOT found !!\n")
    sys.exit(2)

# Dates.
if len(args) == 3:
    beg_date = [ int(x) for x in args[1].split("-") ]
    number_of_days = int(args[2])
elif len(args) == 2:
    beg_date = [ int(x) for x in args[0].split("-") ]
    number_of_days = int(args[1])

beg_date = datetime.date(beg_date[0], beg_date[1], beg_date[2])
number_of_days = datetime.timedelta(days = number_of_days)
end_date = beg_date + number_of_days

sys.stderr.write("beginning date=" + beg_date.strftime("%Y-%m-%d") + "\n")
sys.stderr.write("ending date=" + end_date.strftime("%Y-%m-%d") + "\n")

# Dry run option.
if options.dry_run > 0:
    sys.stderr.write("This is a DRY run, nothing actually performed.\n")

# Process month per month, because of Gocart files.
for year in range(beg_date.year, end_date.year + 1):
    str_year = str(year)
    sys.stderr.write("year " + str_year + "\n")

    month = [1, 12]
    if year == beg_date.year: month[0] = beg_date.month
    if year == end_date.year: month[1] = end_date.month

    for month in range(month[0], month[1] + 1):
        str_month = str(month)
        if month < 10: str_month = "0" + str_month

        sys.stderr.write("  month " + str_month + "\n")

        for species, cfg in zip(species_list, cfg_species_list):
            file_name = file_name_format.replace("%Y", str_year)
            file_name = file_name.replace("%M", str_month)
            file_name = file_name.replace("%S", species)

            command = bc_gocart_prog + " " + cfg_general + " " \
                  + cfg + " " + file_name + " " \
                  + str_year + str_month + " " + str(number_of_days.days)

            if options.dry_run > 0:
                sys.stderr.write("    " + command + "\n")
            else:
                sys.stderr.write("    compute species " + species + " ...\n")

                if os.system(command) != 0:
                    print " ERROR "
                    sys.exit(2)
                else:
                    sys.stderr.write("... done\n")

# Perform NH4 for all month range.
if options.nh4:
    command = bc_nh4_prog + " " + cfg_general \
          + " " + cfg_nh4

    if options.dry_run > 0:
        sys.stderr.write(command + "\n")
    else:
        sys.stderr.write("compute NH4 ...\n")
        if os.system(command) != 0:
            print " ERROR "
            sys.exit(2)
        else:
            sys.stderr.write("... done\n")
# SCRIPT #
##########

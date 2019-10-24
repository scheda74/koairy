#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2009, INERIS - INRIA
# Author(s): Édouard Debry, Vivien Mallet
#
# This file is part of the air quality modeling system Polyphemus.
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

# This script launched chimere.sh, with the necessary preparation steps. One
# is defining the command to be launched to run the interfaced Chimere.

import os, sys
import ops
from atmopy import talos
from datetime import datetime

if len(sys.argv) != 2:
    print "Please provide one configuration file to " \
        + sys.argv[0] + "."
    sys.exit(1)

configuration_file = os.path.abspath(sys.argv[1])
if not os.path.isfile(configuration_file):
    print "ERROR: unable to open the configuration file \"%s\"." \
        % configuration_file
    sys.exit(1)

# Reads the top configuration file.
configuration = ops.Ops(configuration_file)

run_list = configuration.GetEntryList("run")
configuration.SetPrefix("run.")
os.environ["polyphemus_path"] = configuration.GetString("polyphemus_path")
os.environ["verdandi_path"] = configuration.GetString("verdandi_path")
os.environ["main_program_path"] = configuration.GetString("main_program_path")
main_configuration_file = configuration.GetString("main_configuration_file")
os.environ["main_configuration_file"] = main_configuration_file
os.environ["compilation_option"] \
    = configuration.GetString("compilation_option")
chimere_script = configuration.GetString("chimere_script")
os.environ["LAMMPIF77"] = configuration.GetString("lammpi_fortran_compiler")

if not os.path.isfile(chimere_script):
    print "Error: the Chimere script \"" + chimere_script \
        + "\" does not exist."
    sys.exit(1)

if not os.path.isfile(main_configuration_file):
    print "Error: the main configuration file \"" + main_configuration_file \
        + "\" does not exist."
    sys.exit(1)

# Read beginning date from main configuration file.
def get_date_min(configuration_file):
    configuration = talos.Config(configuration_file,
                                 [("Date_min", "[domain]", "DateTime")])
    date_min = configuration.Date_min
    starting_day = datetime(date_min.year, date_min.month, date_min.day)
    if date_min != starting_day:
        print "Error: the starting date " \
            + date_min.strftime("%Y-%m-%d %h:%M") \
            + " is not supported: the starting date should be at midnight."
        sys.exit(1)
    return date_min

if main_configuration_file[-3:] == "cfg":
    date_min = get_date_min(main_configuration_file)
elif main_configuration_file[-3:] == "lua":
    main_configuration = ops.Ops(main_configuration_file)
    date_min = \
        get_date_min(main_configuration.GetString("model_configuration_file"))
else:
    print "ERROR: unknown suffix for configuration file,"
    + " must be either \"lua\" or \"cfg\"."
    sys.exit(1)


# Instead of launching the processing step of Chimere, the following command
# is launched. It compiles the C++ interface to Chimere and launches it.
compile_run_command = os.path.join(os.path.join(os.environ["polyphemus_path"],
                                   "include/models/chimere/compile_run.py"))
os.environ["chimere_interface_script"] = compile_run_command

os.execvp("bash", ("bash", chimere_script, date_min.strftime("%Y%m%d"),))

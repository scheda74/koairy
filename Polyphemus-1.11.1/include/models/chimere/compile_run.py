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

# This script compiles the Chimere C++ interface and runs it.

import os, sys

### Environment variables

main_program_path = os.environ["main_program_path"]
main_program_directory = os.path.dirname(main_program_path)
main_program_name = os.path.basename(main_program_path)
main_configuration_file = os.environ["main_configuration_file"]

compilation_option = os.environ["compilation_option"]

chimere_tmp_path = os.environ["chimere_tmp"]
golam = os.environ["golam"]
goopi = os.environ["goopi"]
nzdoms = int(os.environ["nzdoms"])
nmdoms = int(os.environ["nmdoms"])

### Compilation

os.chdir(main_program_directory)
status = os.system("scons %s %s" % (compilation_option, main_program_name))
if status != 0:
    sys.exit(status)

### Run

os.chdir(chimere_tmp_path)

print "*** Launching the main C++ program"

if golam == "1":
    mpi_arguments = "-ssi rpi tcp"
else:
    mpi_arguments = "-mca btl self,sm,tcp"

np = nzdoms * nmdoms + 1

os.execvp("mpirun", ("mpirun", "-np", "%d" % np,) +
          tuple(mpi_arguments.split()) + (main_program_path,
                                          main_configuration_file,))

# To debug.
#os.execvp("mpirun", ("mpirun", "-np", "%d" % np,) +
#          tuple(mpi_arguments.split())
#          + ("/usr/bin/urxvt -e /usr/bin/gdb " + main_program_path,))

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


import sys, os, optparse

from numpy import *
from atmopy import *

### Initializations.

parser = optparse.OptionParser(usage = "%prog configuration_file")
(options, args) = parser.parse_args()

if not args:
    parser.error("A configuration file is required.")

additional_content = [("Nc", "[input]", "Int"), \
                      ("File_landpar", "[input]", "String"), \
                      ("File_landuse", "[input]", "String"), \
                      ("Directory_out", "[output]", "String")]

config = talos.Config(sys.argv[1], additional_content = additional_content)

Nmonth = 12

file = open(config.File_landpar)
line = file.readline()

roughness = empty([Nmonth, config.Nc], dtype = 'f')
for m in range(Nmonth):
    line = file.readline()
    rough = line.split()[1:]
    for n in range(config.Nc):
        roughness[m, n] = float(rough[n])

file.close()

landuse = fromfile(config.File_landuse, dtype = 'f', sep = '  ')

landuse.shape = (config.Ny, config.Nx, config.Nc)
landuse_trans = landuse.transpose(2,0,1)
landuse_trans.tofile(os.path.join(config.Directory_out, 'LUC.bin'))

roughness_out = empty([Nmonth, config.Ny, config.Nx], dtype = 'f')
for m in range(Nmonth):
    for j in range(config.Ny):
        for i in range(config.Nx):
            roughness_out[m, j, i] = sum(multiply(landuse[j, i], roughness[m]))

roughness_out.tofile(os.path.join(config.Directory_out, 'Roughness.bin'))

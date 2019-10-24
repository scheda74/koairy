# Copyright (C) 2007-2008, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet
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

from numpy import *
from atmopy import *
from atmopy.ensemble import *

import sys

if len(sys.argv) == 1:
    configuration_file = "all.cfg"
elif len(sys.argv) == 2:
    configuration_file = sys.argv[1]
else:
    raise Exception, "No more than one argument!"
Nskip = 30

# Loads ensemble data.
ens = EnsembleData(configuration_file, verbose = True)

# Base combinations.
ebm = BestModel(ens, Nskip = Nskip)
em = EnsembleMean(ens)

# Least-square models.
elsd30 = ELSdN(ens, process = True, Nskip = 30, Nlearning = 30)
elsd30_5 = ELSdN(ens, process = True, Nskip = 5, Nlearning = 30,
                 Nlearning_min = 5)
elsd20_5 = ELSdN(ens, process = True,
                 Nskip = 5, Nlearning = 20, Nlearning_min = 5)
elsd10_5 = ELSdN(ens, process = True, Nskip = 5, Nlearning = 10,
                 Nlearning_min = 5)

# Statistics.
els = ELS(ens, configuration_file = configuration_file)
els.ComputeStatistics(els.all_dates[Nskip:])
print els.stat["rmse"]

# Adding new "models".
add_model(em, elsd30_5, ens)
add_model(em, elsd20_5, ens)
add_model(em, elsd10_5, ens)

# Learning algorithm.
for U in [.95, 1., 1.05]:
    eg = ExponentiatedGradient(ens,extended = True, U = U, Nskip = 1,
                               Nlearning = 1, learning_rate = 1.e-5 )
    eg.ComputeStatistics(eg.all_dates[(Nskip-1):])
    print U, eg.stat['rmse']

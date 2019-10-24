# Copyright (C) 2008, ENPC - INRIA - EDF R&D
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

# This program computes the performances of 'RidgeRegressionDiscounted'
# against the number of models in the ensemble.

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

# Inclusion order: default (no change), from_best (from best to worst model),
# from_worst (from worst to best model).
order = "default"

# Loads ensemble data, with all models.
ens = EnsembleData(configuration_file, verbose = True)
ens.ComputeStatistics(ens.all_dates[Nskip:])

Nmodel = ens.Nsim
# Ranks of the model, from the best to the worst model.
rmse_index = zip(*sorted(zip(ens.stat["rmse"], range(Nmodel))))[1]

rid = RidgeRegressionDiscounted(ens)
rid.ComputeStatistics(ens.all_dates[Nskip:])
print "Reference score (with all members):", rid.stat["rmse"]

restricted_ens = EnsembleData()
restricted_ens.DuplicateEnsemble(ens)

rmse = []
for i in range(Nmodel):
     if order is "default":
          ind = i
     elif order is "from_best":
          ind = rmse_index[i]
     elif order is "from_worst":
          ind = rmse_index[Nmodel - 1 - i]
     else:
          raise Exception, "Unrecognized option."
     restricted_ens.AddSimulation(ens.sim[ind])

     rid = RidgeRegressionDiscounted(restricted_ens)
     rid.ComputeStatistics(ens.all_dates[Nskip:])

     rmse.append(rid.stat["rmse"])
     print i, rmse[-1], ens.stat["rmse"][ind]

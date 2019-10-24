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

# This program computes the performances of the oracles. The ensemble 'ens'
# and the integer 'Nskip' must be available.

mode = "global"

# Best model.
ebm = BestModel(ens, Nskip = Nskip, option = mode)
print "Best model:", ebm.stat["rmse"]

# Best constant combination with weights in the simplex.
eccs = ELS(ens, process = True, Nskip = Nskip, constraint = "simplex",
           option = mode)
print "Best constant combination (simplex):", eccs.stat["rmse"]

# Best constant combination.
ecc = ELS(ens, process = True, Nskip = Nskip, option = mode)
print "Best constant combination:", ecc.stat["rmse"]

if mode == "global":
    # Best combination.
    eb = ELSd(ens, process = True, Nskip = Nskip)
    print "Best combination:", eb.stat["rmse"]

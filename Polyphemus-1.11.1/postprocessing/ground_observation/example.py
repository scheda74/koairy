# Copyright (C) 2011, INRIA
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


# Note that the SWIG module in 'atmopy.observation' must be compiled.
from atmopy import observation
# The Seldon SWIG interface is also needed.
from seldon import *

om = observation.GroundNetworkOM("observation.lua");

obs = VectorDouble()
om.SetDate("2001-06-15 09:00")
om.GetObservation(obs)
obs.Print()

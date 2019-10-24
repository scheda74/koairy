# Copyright (C) 2012, INRIA
#     Author(s): Vivien Mallet
#
# This file is part of AtmoPy library, a tool for data processing and
# visualization in atmospheric sciences.
#
# AtmoPy is developed in the INRIA - ENPC joint project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# AtmoPy is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# AtmoPy is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the AtmoPy home page:
#     http://cerea.enpc.fr/polyphemus/atmopy.html


"""
Observation is a module from AtmoPy package that provides an interface to
Polair3D.
"""


import os
current_directory = os.path.dirname(os.path.abspath(__file__))
if os.path.isfile(os.path.join(current_directory, "_polair3d.so")) \
        and os.path.isfile(os.path.join(current_directory, "polair3d.py")):
    from polair3d import *

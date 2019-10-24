# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vincent Picavet, Vivien Mallet
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
Copyright (C) 2005 - 2007, ENPC - INRIA - EDF R&D

This file is part of AtmoPy library, a tool for data processing and
visualization in atmospheric sciences.

AtmoPy is developed in the INRIA - ENPC joint project-team CLIME and in the
ENPC - EDF R&D joint laboratory CEREA.

AtmoPy is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

AtmoPy is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

For more information,visit the AtmoPy home page:
http://www.enpc.fr/cerea/polyphemus/atmopy.html
"""


import display
try:
    import ensemble
except:
    pass
import io
import observation
import stat
import talos


def atmopy_test():
    """
    Tests AtmoPy installation. If no error occurs while this function is
    executed, AtmoPy is properly installed.
    """
    import numpy
    import tempfile, os

    # Configuration file.
    o, name = tempfile.mkstemp()
    o = open(name, "w")
    o.write("x_min = -10.5 Delta_x = 0.5 Nx = 67\n"
            + "y_min = 35. Delta_y = 0.5 Ny = 48\n")
    o.close()

    # Generates data.
    def f(x_, y_):
        return (1. - x / 2. + x**5 + y**3) * numpy.exp(- x**2 - y**2)
    x = numpy.arange(67., dtype = 'f') * .1 - 3.
    y = numpy.arange(48., dtype = 'f') * .1 - 2.5
    x, y = numpy.meshgrid(x, y)
    d = f(x, y)
    d = numpy.array([d for i in range(3)], dtype = 'f') + 1.
    # Saves data in a binary file.
    o, name_data = tempfile.mkstemp()
    d.tofile(name_data)

    # Tests Talos and display capabilities.
    m = display.getm(name)
    d = display.getd(name, name_data, Nt = 0, Nz = 1)
    display.dispcf(m, stat.spatial_distribution(d, "mean")[0], V = 25)

    # Removes temporary files.
    os.remove(name)
    os.remove(name_data)

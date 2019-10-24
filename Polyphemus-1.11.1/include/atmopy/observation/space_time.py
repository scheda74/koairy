# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vincent Picavet
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


def remove_lowest(dates, data1, data2, mini):
    """
    Removes values of data1, data2 and dates where data1 values
    are lower than the given minimum.

    @type dates: sequence of datetime.datetime
    @param dates: Dates corresponding to the data.
    @type data1: 1D numpy.array
    @param data1: Data from which values lower than mini must be
    removed (observations most of the time).
    @type data2: Data where values corresponding to removed value in
    data1 are also removed (simulation data).
    @type mini: float
    @param mini: Minimum value. Dates, elements of data1 and data2
    corresponding to values lower than mini in data1 are removed.
    @rtype: sequence of datetime.datetime, 1D numpy.array,
    1D numpy.array
    @return: sequence of dates, data1 array and data2 array.
    Elements corresponding to values lower than mini in data1
    are not returned.
    """
    condition = (data1 > mini)
    for i in range(len(condition)-1, -1, -1):
        if condition[i] == 0:
            dates.pop(i)
    return dates, data1[numpy.where(condition)], \
           data2[numpy.where(condition)]


def remove_highest(dates, data1, data2, maxi):
    """
    Removes values of data1, data2 and dates where data1 values
    are higher than the given minimum.

    @type dates: sequence of datetime.datetime
    @param dates: Dates corresponding to the data.
    @type data1: 1D numpy.array
    @param data1: Data from which values higher than mini must be
    removed (observations most of the time).
    @type data2: Data where values corresponding to removed value in
    data1 are also removed (simulation data).
    @type maxi: float
    @param maxi: Maximum value. Dates, elements of data1 and data2
    corresponding to values higher than maxi in data1 are removed.
    @rtype: sequence of datetime.datetime, 1D numpy.array,
    1D numpy.array
    @return: sequence of dates, data1 array and data2 array.
    Elements corresponding to values higher than maxi in data1
    are not returned.
    """
    condition = (data1 < maxi)
    for i in range(len(condition)-1, -1, -1):
        if condition[i] == 0:
            dates.pop(i)
    return dates, data1[numpy.where(condition)], \
           data2[numpy.where(condition)]

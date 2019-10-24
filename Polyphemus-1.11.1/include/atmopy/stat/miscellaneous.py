# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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


from numpy import *


def spatial_distribution(data, function):
    """
    Applies a function to a time series in every cell. It therefore returns
    the spatial distribution of an indicator.

    @type data: numpy.array
    @param data: Data to be processed. Time is assumed to be the first
    dimension. There must be 1, 2 or 3 extra-dimensions.
    @type function: string or function
    @param function: The function to be applied to the time series. If
    'function' is a string, it is assumed to be a numpy.array method.
    """
    m = zeros(data.shape[1:], dtype = 'd')
    if data.ndim == 2:
        for i in range(data.shape[1]):
            if isinstance(function, str):
                m[i] = getattr(data[:, i], function)()
            else:
                m[i] = function(data[:, i])
    elif data.ndim == 3:
        for i in range(data.shape[1]):
            for j in range(data.shape[2]):
                if isinstance(function, str):
                    m[i, j] = getattr(data[:, i, j], function)()
                else:
                    m[i, j] = function(data[:, i, j])
    elif data.ndim == 4:
        for i in range(data.shape[1]):
            for j in range(data.shape[2]):
                for k in range(data.shape[3]):
                    if isinstance(function, str):
                        m[i, j, k] = getattr(data[:, i, j, k], function)()
                    else:
                        m[i, j, k] = function(data[:, i, j, k])
    elif data.ndim == 5:
        for i in range(data.shape[1]):
            for j in range(data.shape[2]):
                for k in range(data.shape[3]):
                    for l in range(data.shape[4]):
                        if isinstance(function, str):
                            m[i, j, k, l] = getattr(data[:, i, j, k, l],
                                                    function)()
                        else:
                            m[i, j, k, l] = function(data[:, i, j, k, l])
    elif data.ndim > 5:
        raise ValueError, "Too many dimensions (" + str(data.ndim) \
              + "). There should be 2, 3, 4 or 5 dimensions."
    else:
        raise ValueError, "Too few dimensions (" + str(data.ndim) \
              + "). There should be 2, 3, 4 or 5 dimensions."

    return m


def time_evolution(data, function):
    """
    Computes the time evolution of a given indicator on spatial fields.

    @type data: numpy.array
    @param data: Data to be processed. Time is assumed to be the first
    dimension.
    @type function: string or function
    @param function: The function to be applied to the fields. If 'function'
    is a string, it is assumed to be a numpy.array method.
    """
    if isinstance(function, str):
        return array([getattr(x, function)() for x in data])
    else:
        return array([function(x) for x in data])

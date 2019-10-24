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


import numpy
import datetime
import sys, os
sys.path.insert(0,
                os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
import observation


def get_filesize(filename):
    """
    Gets a file's size.

    @type filename: string
    @param filename: The name of the file.

    @rtype: int
    @return: File's size, 0 if no readable file found.
    """
    fileSize = 0
    try:
        f = open(filename, "r", 0)
        try:
            f.seek(0,2)
            fileSize = f.tell()
        finally:
            f.close()
    except IOError:
        pass
    return fileSize


def get_timesteps(filename, recordLength):
    """
    Gets the number of timesteps in a file given its
    record length (in Bytes).

    @type filename: string
    @param filename: The name of the file.

    @type recordLength: integer
    @param recordLength: Size of one record.

    @rtype: int
    @return: Number of timesteps.
    """
    ts = 0
    if (recordLength != 0):
        ts = get_filesize(filename) / recordLength
    return ts


def load_binary(filename, shape, type = 'f'):
    """
    Loads a binary file into an array using specified shape.
    Returns numpy.

    @type filename: string or Python file object.
    @param filename: The name of the file to load.

    @type shape: tuple
    @param shape: The shape of the array to load from file.

    @type type: string
    @param type: Type of data read. Default is 'f'

    @rtype: numpy.array
    @return: New array filled with binary data from specified file.
    """
    length = 1
    for l in shape:
        length *= l
    d = numpy.fromfile(filename, type, length)
    if d.shape[0] != length:
        raise Exception, "File \"" + filename \
              + "\" does not contain enough elements."
    d.shape = shape
    d = d.astype('d')
    return d


def load_binary_first_level(filename, shape, type = 'f'):
    """
    Loads a binary file into an array using specified 3D shape for
    X, Y and T dimensions (a time sequence of planes).
    If the given binary file is a 4D file (XYZT), the plane Z = 1 is
    extracted.

    @type filename: string or Python file object.
    @param filename: The name of the file to load.

    @type shape: tuple
    @param shape: The 3D shape of the array to load from file.

    @type type: string
    @param type: Type of data read. Default is 'f'

    @rtype: numpy.array
    @return: New 3D array os given shape filled with binary data
    from specified file.
    """
    res = []
    zsize = get_filesize(filename) \
            / (numerictypes.getType(type).bytes * shape[0] \
               * shape[1] * shape[2])
    if zsize != 1:
        res = load_binary(filename, shape, type)
    else:
        newshape = list(shape)
        newshape.insert(1, zsize)
        # Use temp array to be sure that memory is freed
        temp = load_binary(filename, newshape, type)
        res = temp[:,0,:,:]
        del temp
    return res


def save_binary(arrayToSave, filename, type = 'f'):
    """
    Saves a numpy in a binary file using specified type.

    @type arrayToSave: numpy.array
    @param arrayToSave: The array to save.

    @type filename: string or python file object
    @param filename: The name of the file to save the array into.

    @type type: string
    @param type: Format of data to save the array in file.
    """
    numpy.array(arrayToSave, dtype = type).tofile(filename)

def filter_config(config, data):
    """
    Filters data based on the cells and the days to be discarded according to
    a given configuration.

    @type config: Config
    @param config: The configuration associated with 'data'.
    @type data: numpy.array
    @param data: The data array to be filtered.

    @rtype: (list of datetime, numpy.array)
    @return: The dates and the output data. The output data array does not
    contain:
       0. the first and/or last days if they are missing data, except if
       'config.discarded_days' is negative;
       1. days to be removed at the beginning (if the first day is incomplete,
       it is removed, but not counted as a discarded day)
       'config.discarded_days'.
       2. the cells to be removed in the domain edges according to
       'config.discarded_cells'.
    """
    dates = observation.get_simulation_dates(config.t_min, config.Delta_t,
                                             config.Nt)
    if config.discarded_days >= 0:
        dates, data = observation.remove_incomplete_days(dates, data)
        dates, data \
               = observation.remove_days(dates, data, config.discarded_days)
    if config.discarded_cells <= 0:
        return dates, data
    if len(data.shape) == 4:   # Z included.
        return dates, \
               data[:, :, config.discarded_cells:-config.discarded_cells,
                    config.discarded_cells:-config.discarded_cells]
    else:   # Only X and Y.
        return dates, data[:, config.discarded_cells:-config.discarded_cells,
                           config.discarded_cells:-config.discarded_cells]

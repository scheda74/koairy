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
import os, sys
sys.path.insert(0,
                os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
import observation
sys.path.pop(0)


def load_stations(filename, type, origins = (0, 0), \
                  deltas = (0, 0), lengths = (0, 0)):
    """
    Loads stations description from text file. Removes stations
    outside the domain described by origins, deltas and
    lengths (if specified).

    @type filename: string
    @param filename: Name of the file that describes the stations.

    @type type: string
    @param type: The type of the initialization strings in the file. It could
    be Default, Pioneer, Emep or Bdqa.

    @type origins: (float, float)
    @param origins: (latitude, longitude) as coordinates of the first cell in
    the mesh.

    @type deltas: (float, float)
    @param deltas: mesh sizes on latitude and longitude.

    @type lengths: (int, int)
    @param lengths: lengths of the mesh on latitude and longitude.

    @rtype: sequence of Station
    @return: sequence of stations loaded from specified file.
    """
    stations = []
    try:
        f = open(filename)
        for i in f.readlines():
            station = observation.Station(i, type)
            if deltas == (0,0) or \
                   station.IsInsideGridBox(origins, deltas, lengths):
                stations.append(station)
        f.close()
    except IOError:
        pass
    return stations


def load_station(filename, type, station_name):
    """
    Loads a specific station description from text file containing
    stations descriptions.

    @type filename: string
    @param filename: Name of the file that describes the stations.
    @type type: string
    @param type: The type of the initialization strings in the file. It could
    be Default, Pioneer, Emep or Bdqa.
    @type station_name: string
    @param station_name: Name of the station to load from the file.

    @rtype: Station
    @return: Station loaded from specified file.
    """
    lines = open(filename).readlines()
    for line in lines:
        station = observation.Station(line, type)
        if station.name == station_name:
            return station

def load_file_observations(name, directory):
    """
    Loads observations data from a file and puts it in a sequence.

    @type name: string
    @param name: Name of the file to load data from (without directory).
    @type directory: string
    @param directory: location of the specified file.

    @rtype: datetime.datetime sequence, 1D numpy.array
    @return: Sequence of datetime and corresponding observation values in 1D
    array.
    """

    filename = os.path.normpath(directory) + '/' + name
    dates = []
    observations = []
    try:
        f = open(filename)
        for i in f.readlines():
            line = i.split()
            date = line[0]
            year = int(date[0:4])
            month = int(date[4:6])
            day = int(date[6:8])
            if len(date) == 10:
                hour = int(date[8:10])
                dates.append(datetime.datetime(year, month, day, hour))
            else:
                dates.append(datetime.datetime(year, month, day))
            observations.append(float(line[1]))
        f.close()
    except IOError:
        pass
    return dates, numpy.array(observations, 'Float32')


def load_observations(stations, directory):
    """
    Loads observations data from files for given stations
    File names are guessed from station.name.

    @type stations: list of Station
    @param stations: Stations for which observation data must be returned.
    @type directory: string
    @param directory: location of the observation files.

    @rtype: list of datetime.datetime lists, list of numpy.array.
    @return: A list of dates list (one list of dates per Station), and a
    list of observations arrays (one array per Station).
    """
    obs_dates_list = []
    obs_list = []
    for i in stations:
        dates, obs = load_file_observations(i.name, directory)
        obs_dates_list.append(dates)
        obs_list.append(obs)
    return obs_dates_list, obs_list

# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Irene Korsakissok
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


import os, sys
sys.path.insert(0,
                os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
import talos
sys.path.pop(0)


class Source:
    """
    Stores information about a point source
    """

    def __init__(self, content = []):
        """
        Initializes the instance in case 'content' is not empty.

        @type  content: list.
        @param content: list of values to initialize the source with.
        It contains the following fields:
           0. The source abscissa (float);
           1. The source ordinate (float);
           2. The source altitude (float);
           2. The emission beginning date (date);
           2. The emission ending date (date);
           3. The list of species emitted from the source;
           5. The list of species rates.
        """
        if len(content) == 7:
            self.abscissa = float(content[0])
            self.ordinate = float(content[1])
            self.altitude = float(content[2])
            self.date_beg = content[3]
            self.date_end = content[4]
            self.species = content[5]
            self.rate = content[6]
        else:
            print 'Error: list of attributes must be of length 7.'


    def SetSourceLocation(self, coordinates):
        """
        Sets the source coordinates.
        @type coordinates: list
        @param coordinates: The list of coordinates (abscissa, ordinate).
        """
        self.abscissa = float(coordinates[0])
        self.ordinate = float(coordinates[1])


def read_sources(config_file):
    """
    Loads sources attributes from a configuration file.

    @type config_file: string
    @param config_file: The configuration file in which source attributes are
    to be found.
    @return: The list of sources.
    """
    sources = []
    config_stream = talos.ConfigStream(config_file)
    attribute_list = ['Abscissa', 'Ordinate', 'Altitude', 'Date_beg', \
                      'Date_end', 'Species', 'Rate']
    section = config_stream.ListAll().split('[source]')
    for i in range(0,len(section)):
        content = []
        coord = []
        date = []
        lines = section[i].split('\n')
        if lines != ['']:
            for j in range(len(lines)):
                l = lines[j].split(':')
                for k in range(len(attribute_list)):
                    if l[0] == attribute_list[k]:
                        if k == 0 or k == 1 or k == 2:
                            coord.append(float(l[1]))
                        elif k == 3 or k == 4:
                            d = config_stream.StringToDateTime(l[1])
                            date.append(d)
                        elif k == 5:
                            species = l[1].split()
                        elif k == 6:
                            a = l[1].split()
                            rate = [float(a[s]) for s in range(len(a))]
            for k in range(3):
                content.append(coord[k])
            for k in range(2):
                content.append(date[k])
            content.append(species)
            content.append(rate)
            source = Source(content)
            sources.append(source)
    return sources


def get_sources_coordinates(config_file):
    """
    Loads sources coordinates from a configuration file.

    @type config_file: string
    @param config_file: The configuration file in which source attributes are
    to be found.
    @return: The lists of abscissa and ordinates of all sources.
    """
    sourcelist = read_sources(config_file)
    Nsource = len(sourcelist)
    xslist = [sourcelist[i].abscissa for i in range(Nsource)]
    yslist = [sourcelist[i].ordinate for i in range(Nsource)]
    return xslist, yslist


def get_species_sources_coordinates(config_file, species, min_rate = 0.):
    """
    Loads sources coordinates from a configuration file, for sources emitting
    a given species.

    @type config_file: string
    @type species: string
    @type min_rate: float
    @param config_file: The configuration file in which source attributes are
    to be found.
    @param species: The species name.
    @param min_rate: the minimal rate to be considered for the species.
    @return: The lists of abscissa and ordinates of all sources emitting the
    species 'species' with a higher rate than 'min_rate'.
    """
    sourcelist = read_sources(config_file)
    Nsource = len(sourcelist)
    xslist = []
    yslist = []
    for i in range(Nsource):
        s = sourcelist[i]
        species_index = -1
        for j in range(len(s.species)):
            if s.species[j] == species:
                species_index = j
        if species_index != -1:
            if s.rate[species_index] > min_rate:
                xslist.append(s.abscissa)
                yslist.append(s.ordinate)
    return xslist, yslist


def get_level_sources_coordinates(config_file, level_file):
    """
    Loads sources coordinates from a configuration file, sorted by vertical
    level.
    @type config_file: string
    @type level_file: string
    @param config_file: The configuration file in which source attributes are
    to be found.
    @param level_file: The file with the list of vertical level interfaces.
    @return: The lists of (list of) abscissa and ordinates of all sources,
    sorted by level.
    """
    f = open(level_file)
    levels = []
    for line in f.xreadlines():
        l = line.split()
        if l != []:
            levels.append(float(l[0]))
    f.close()
    xslist , yslist = get_level_sources_coordinates(config_file, levels)
    return xslist, yslist


def get_level_sources_coordinates(config_file, levels):
    """
    Loads sources coordinates from a configuration file, sorted by vertical
    level.
    @type config_file: string
    @type levels: list
    @param config_file: The configuration file in which source attributes are
    to be found.
    @param levels: The list of vertical level interfaces.
    @return: The lists of (list of) abscissa and ordinates of all sources,
    sorted by level.
    """
    sourcelist = read_sources(config_file)
    Nsource = len(sourcelist)
    xslist = []
    yslist = []
    zslist = []
    for k in range(len(levels)-1):
        x = []
        y = []
        z = []
        for i in range(Nsource):
            s = sourcelist[i]
            if s.altitude > levels[k] and s.altitude < levels[k+1]:
                x.append(s.abscissa)
                y.append(s.ordinate)
                z.append(s.altitude)
        xslist.append(x)
        yslist.append(y)
        zslist.append(z)
    return xslist, yslist, zslist

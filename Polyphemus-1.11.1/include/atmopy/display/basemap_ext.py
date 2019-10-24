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


from pylab import *
from numpy import *
try:
    from matplotlib.toolkits.basemap import Basemap
except:
    pass
# From version 0.98.0 of Matplotlib.
try:
    from mpl_toolkits.basemap import Basemap
except:
    pass
import sys, os
sys.path.insert(0,
                os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
import talos
sys.path.pop(0)


def getm(config = None, y_min = None, x_min = None,
         Delta_y = None, Delta_x = None, Ny = None, Nx = None, cbar = True,
         open_figure = True, resolution = 'l', area_thresh = 1000):
    """
    Generates a map with Basemap.

    @type config: Config or string
    @param config: The configuration or the configuration file.
    @type cbar: Boolean
    @param cbar: True if there is a colormap, false otherwise.
    @type open_figure: Boolean
    @param open_figure: Should a figure be opened?

    @rtype: Basemap
    @return: The map.
    """
    if isinstance(config, str):
        if os.path.exists(config):
            config = talos.Config(config)
        else:
            raise IOError, "Configuration file \"" + config + "\" not found."

    if x_min is None:
        x_min = config.x_min
    if y_min is None:
        y_min = config.y_min
    if Delta_x is None:
        Delta_x = config.Delta_x
    if Delta_y is None:
        Delta_y = config.Delta_y
    if Nx is None:
        Nx = config.Nx
    if Ny is None:
        Ny = config.Ny

    m = Basemap(projection = 'cyl',
                llcrnrlon = x_min - Delta_x / 2.,
                llcrnrlat = y_min - Delta_y / 2.,
                urcrnrlon = x_min + Delta_x / 2. + Delta_x * float(Nx - 1),
                urcrnrlat = y_min + Delta_y / 2. + Delta_y * float(Ny - 1),
                resolution = resolution, suppress_ticks = False,
                area_thresh = area_thresh)
    if open_figure:
        fig_num = get_current_fig_manager().num
        xsize = rcParams['figure.figsize'][0]
        fig = figure(num = fig_num)
        if cbar:
            ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
            axes(ax)
            axes([0.875, 0.1, 0.05, 0.75])
        else:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            axes(ax)
    return m


def getd(config = None, filename = "", Nt = None,
         Nz = None, Ny = None, Nx = None):
    """
    Reads data from a binary file.

    @type config: Config or string
    @param config: The configuration or the configuration file.
    @type filename: string
    @param filename: The file to be loaded. If filename is empty, then the
    file from 'config' is loaded.
    @type Nt: integer.
    @param Nt: The number of time steps in the file to be loaded. If it is not
    given, it is read in 'config'.
    @type Nz: integer.
    @param Nz: The number of levels in the file to be loaded. If it is not
    given, it is read in 'config'.
    @type Ny: integer.
    @param Ny: The number of space steps along y in the file to be loaded. If
    it is not given, it is read in 'config'.
    @type Nx: integer.
    @param Nx: The number of space steps along xin the file to be loaded. If
    it is not given, it is read in 'config'.

    @rtype: numpy.array
    @return: The data.
    """
    if isinstance(config, str):
        config = talos.Config(config)
    if filename == "":
        filename = config.input_file

    import os
    if Nx is None:
        Nx = config.Nx
    if Ny is None:
        Ny = config.Ny
    if Nz is None:
        Nz = config.Nz
    if Nt is None:
        Nt = config.Nt

    if Nx == 0:
        Nx = int(os.stat(filename)[6] / 4) / Ny / Nz / Nt
    if Ny == 0:
        Ny = int(os.stat(filename)[6] / 4) / Nx / Nz / Nt
    if Nz == 0:
        Nz = int(os.stat(filename)[6] / 4) / Nx / Ny / Nt
    if Nt == 0:
        Nt = int(os.stat(filename)[6] / 4) / Nx / Ny / Nz

    length = 1
    for l in [Nt, Nz, Ny, Nx]:
        length *= l
    d = fromfile(filename, 'f', length)
    d.shape = (Nt, Nz, Ny, Nx)
    return d.astype('f8')


def getdJ(config = None, filename = "", Ndays = None,
         Ntheta = None, Ny = None, Nz = None):
    """
    Reads data from a binary file.

    @type config: Config or string
    @param config: The configuration or the configuration file.
    @type filename: string
    @param filename: The file to be loaded. If filename is empty,
    then the file from 'config' is loaded.
    @type Nd: integer.
    @param Nd: The number of day steps in the file to be loaded.
    If it is not given, it is read in 'config'.
    @type Ntheta: integer.
    @param Ntheta: The number of time solar angle in the file
    to be loaded. If it is not given, it is read in 'config'.
    @type Ny: integer.
    @param Ny: The number of space steps along y in the file
    to be loaded. If it is not given, it is read in 'config'.
    @type Nz: integer.
    @param Nx: The number of levels in the file to be loaded. If
    it is not given, it is read in 'config'.

    @rtype: numpy.array
    @return: The data.
    """
    if isinstance(config, str):
        config = talos.Config(config)
    if filename == "":
        filename = config.input_file

    import os
    if Ndays is None:
        Ndays = config.Ndays
    if Ntheta is None:
        Ntheta = config.Ntheta
    if Ny is None:
        Ny = config.Ny
    if Nz is None:
        Nz = config.Nz

    if Ndays == 0:
        Ndays = int(os.stat(filename)[6] / 4) / Ntheta / Ny / Nz
    if Ntheta == 0:
        Ntheta = int(os.stat(filename)[6] / 4) / Ndays / Ny / Nz
    if Ny == 0:
        Ny = int(os.stat(filename)[6] / 4) / Ndays / Ntheta / Nz
    if Nz == 0:
        Nz = int(os.stat(filename)[6] / 4) / Ndays / Ntheta / Ny

    length = 1
    for l in [Ndays, Ntheta, Ny, Nz]:
        length *= l
    daux = fromfile(filename, 'f', length)
    daux.shape = (Ndays, Ntheta, Ny, Nz)

    d=zeros([Nz,Ny,Ntheta,Ndays])
    for i in range(Ndays):
        for j in range(Ntheta):
            for k in range(Ny):
                for l in range(Nz):
                    d[l,k,j,i] = daux[i,j,k,l]

    return d.astype('f8')


def getmd(config, cbar = True):
    """
    Reads data from a binary file and generates the corresponding map.

    @type config: Config or string
    @param config: The configuration or the configuration file.
    @type cbar: Boolean
    @param cbar: True is there is a colormap, false otherwise.

    @rtype: (Basemap, numpy.array)
    @return: The map and the data.
    """
    return (getm(config), getd(config))


def disp(map, data, **kwargs):
    """
    Displays a 2D array on a given map.

    @type map: Basemap
    @param map: The map on which data is displayed.
    @type data: 2D numpy.array
    @param data: Data (2D) to be displayed.
    """
    if data.ndim != 2:
        raise Exception, "Function \"disp\" proceeds with 2D data," \
              + " but input data has " + str(data.ndim) + " dimension(s)."

    # If the figure is empty, sets new axes.
    if len(gcf().axes) == 0:
        xsize = rcParams['figure.figsize'][0]
        fig_num = get_current_fig_manager().num
        fig = figure(num = fig_num)
        ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
        axes(ax)
        axes([0.875, 0.1, 0.05, 0.75])

    # Clears current image.
    with_states = kwargs.has_key("states") and kwargs["states"]
    try:
        kwargs.pop("states")
    except:
        pass
    gcf().axes[0].clear()
    axes(gcf().axes[0])
    map.imshow(data, **kwargs)
    map.drawcountries()
    map.drawcoastlines()
    if with_states:
        map.drawstates()

    # Colorbar.
    if len(gcf().axes) > 1:
        gcf().axes[1].clear()
        cax = gcf().axes[1]
        colorbar(cax = cax)


def dispcf(map, data, V = None, **kwargs):
    """
    Displays a 2D array on a given map with filled contours.

    @type map: Basemap
    @param map: The map on which data is displayed.
    @type data: 2D numpy.array
    @param data: Data (2D) to be displayed.
    @type V: integer, list or 1D numpy.array
    @param V: The number of levels or the list of thresholds for the contours.
    """
    if data.ndim != 2:
        raise Exception, "Function \"dispcf\" proceeds with 2D data," \
              + " but input data has " + str(data.ndim) + " dimension(s)."

    # If the figure is empty, sets new axes.
    if len(gcf().axes) == 0:
        xsize = rcParams['figure.figsize'][0]
        fig_num = get_current_fig_manager().num
        fig = figure(num = fig_num)
        ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
        axes(ax)
        axes([0.875, 0.1, 0.05, 0.75])

    # Clears current image.
    with_states = kwargs.has_key("states") and kwargs["states"]
    try:
        kwargs.pop("states")
    except:
        pass
    gcf().axes[0].clear()
    axes(gcf().axes[0])

    xrange, yrange = map.makegrid(data.shape[1], data.shape[0])
    if kwargs.has_key("colors"):
        V = len(kwargs["colors"]) - 1
    if V is None:
        map.contourf(xrange, yrange, data, **kwargs)
    else:
        map.contourf(xrange, yrange, data, V, **kwargs)
    map.drawcountries()
    map.drawcoastlines()
    if with_states:
        map.drawstates()

    # Colorbar.
    if len(gcf().axes) > 1:
        gcf().axes[1].clear()
        cax = gcf().axes[1]
        colorbar(cax = cax,format="%1.2f")


def gridprofile(map, nx = -1, ny = -1, zz = -1):
    """
    Returns arrays containing lon,lat of an equally spaced native projection
    horizontal grid (if zz = -1) or containing lon,height (if ny = -1)
    or lat,height (if nx = -1) of a non-equally spaced vertical grid.
    """
    if(nx != -1 and ny != -1 and zz != -1) \
              or (nx == -1 and ny == -1) or (nx == -1 and zz == -1) \
              or (ny == -1 and zz == -1):
        print "ERROR : cannot execute program \"gridprofile\"."

    if zz == -1:
        lons, lats = map.makegrid(nx, ny)
        return lons, lats

    if ny == -1:
        lons2, lats2 = map.makegrid(nx, 2)
        lons = empty([len(zz), nx], dtype = float32)
        heights = empty([len(zz), nx], dtype = float32)
        for i in range(len(zz)):
            lons[i] = lons2[0]
            for j in range(nx):
                heights[i, j] = zz[i]
        return lons, heights

    if nx == -1:
        lons2, lats2 = map.makegrid(2, ny)
        lats = empty([len(zz), ny], dtype = float32)
        heights = empty([len(zz), ny], dtype = float32)
        for i in range(len(zz)):
            for j in range(ny):
                lats[i, j] = lats2[j, 0]
                heights[i, j] = zz[i]
        return lats, heights


def profile2Dcf(config, data, xx = -1, yy = -1, zz = -1,
                lon = -1, lat = -1, alt = -1,
                NLevels = 10, V = None, **kwargs):
    """
    Displays a 2D section of a 3D array with filled contours.
    """
    # Some Errors' messages.
    if data.ndim != 3:
        raise Exception, "Function \"profile\" proceeds with 3D data," \
              + " but input data has " + str(data.ndim) + " dimension(s)."

    if(xx == -1 and yy == -1 and zz == -1 \
       and lon == -1 and lat == -1 and alt == -1):
        raise Exception, "Function \"profilecf\" cannot display 3D maps." \
              + " Missing one argument among xx, yy, zz, lon, lat, or alt."

    if xx != -1 and lon != -1:
        raise Exception, "You cannot specify a \"xx\" value and a \"lon\"" \
              + " value at the same time."

    if yy != -1 and lat != -1:
        raise Exception, "You cannot specify a \"yy\" value and a \"lat\"" \
              + " value at the same time."

    if zz != -1 and alt != -1:
        raise Exception, "You cannot specify a \"zz\" value and a \"alt\"" \
              + " value at the same time."

    arg_num = 0
    arg_list = [xx, yy, zz, lon, lat, alt]
    for item in arg_list:
        if item != -1:
            arg_num += 1

    if arg_num != 1:
        raise Exception, "Cannot execute \"profilecf\". Only one argument" \
              + " among \"xx\", \"yy\", \"zz\", \"lon\", \"lat\", \"alt\"" \
              + " should be specified."

    # Reads configuration file 'config'.
    if isinstance(config, str):
        if os.path.isfile(config):
            config = talos.Config(config)
        else:
            raise Exception, "Cannot read configuration file : \" " \
                  + str(config) + " \"."
    else:
        raise Exception, "First argument has to be a config file's name."

    x_min = config.x_min
    y_min = config.y_min
    Delta_x = config.Delta_x
    Delta_y = config.Delta_y
    Nx = config.Nx
    Ny = config.Ny
    Nz = config.Nz
    levels = config.levels

    x_max = x_min + Delta_x * Nx
    y_max = y_min + Delta_y * Ny
    z_max = levels[-1]

    # Some other Errors' messages.
    if(xx > data.shape[2] or (lon !=-1 and (lon < x_min or lon > x_max))):
        raise Exception, "Cannot display profile. data's maximum x_" \
              + "coordinate is " + str(data.shape[2]) + " and data's" \
              + " longitude is in [" + str(x_min) + "; " + str(x_max) + "]" \
              + " (degrees)."

    if(yy > data.shape[1] or (lat !=-1 and (lat < y_min or lat > y_max))):
        raise Exception, "Cannot display profile. data's maximum y_" \
              + "coordinate is " + str(data.shape[1]) + " and data's" \
              + " latitude is in [" + str(y_min) + "; " + str(y_max) + "]" \
              + " (degrees)."

    if(zz > data.shape[0] or alt > z_max):
        raise Exception, "Cannot display profile. data's maximum z_" \
              + "coordinate is " + str(data.shape[0]) + " and data's" \
              + " maximum altitude is " +  str(z_max) + " (m)."

    # Transformation of lon, lat, alt in xx, yy, zz if necessary.

    if lon != -1:
        xx = int((lon - x_min) / Delta_x)

    if lat != -1:
        yy = int((lat - y_min) / Delta_y)

    if alt != -1:
        zz = 0
        for z_index in range(len(levels)):
            if alt > levels[z_index]:
                zz = z_index + 1

    # Initialization of the map.

    mref = Basemap(projection = 'cyl',
                   llcrnrlon = x_min - Delta_x / 2.,
                   llcrnrlat = y_min - Delta_y / 2.,
                   urcrnrlon = x_min + Delta_x / 2. + Delta_x * float(Nx - 1),
                   urcrnrlat = y_min + Delta_y / 2. + Delta_y * float(Ny - 1),
                   resolution = 'l', suppress_ticks = False,
                   area_thresh = 1000)

    # If the figure is empty, sets new axes.
    if len(gcf().axes) == 0:
        xsize = rcParams['figure.figsize'][0]
        fig_num = get_current_fig_manager().num
        fig = figure(num = fig_num)
        ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
        axes(ax)
        axes([0.875, 0.1, 0.05, 0.75])

    # Clears current image.
    gcf().axes[0].clear()
    axes(gcf().axes[0])

    # Display of profiles.
    if(zz != -1):
        # Is equivalent to dispcf ( = horizontal section)
        data2 = data[zz]
        if V is None and NLevels != 1:
            V = arange(NLevels) / float(NLevels - 1) *\
                (data2.max() - data2.min()) + data2.min()
        xrange, yrange = mref.makegrid(data2.shape[1], data2.shape[0])
        if kwargs.has_key("colors"):
            V = len(kwargs["colors"]) - 1
        mref.contourf(xrange, yrange, data2, V, extend='both', **kwargs)
        mref.drawcountries()
        mref.drawcoastlines()

    if(yy != -1):
        data2 = data[:,yy,:]
        if data2.shape[0] == len(levels) + 1:
            # Skips the first level for Kz (Kz[0] = 0.2).
            data2 = data2[1:(len(levels) + 1),:]
        if V is None and NLevels != 1:
            V = arange(NLevels) / float(NLevels - 1) *\
                (data2.max() - data2.min()) + data2.min()
        xrange, zrange = gridprofile(mref, data2.shape[1], -1, levels)
        if kwargs.has_key("colors"):
            V = len(kwargs["colors"]) - 1
        contourf(xrange, zrange, data2, V, extend='both', **kwargs)
        gca().set_xlim([xrange[0][0], xrange[0][-1]])
        gca().set_ylim([zrange[0][0], zrange[-1][0]])

    if xx != -1:
        data2 = data[:,:,xx]
        if data2.shape[0] == len(levels) + 1:
            data2 = data2[1:(len(levels) + 1),:]
        if V is None and NLevels != 1:
            V = arange(NLevels) / float(NLevels - 1) *\
                (data2.max() - data2.min()) + data2.min()
        yrange, zrange = gridprofile(mref, -1, data2.shape[1], levels)
        print yrange[0][0], yrange[0][-1], zrange[0][0], zrange[-1][0]
        print zrange
        if kwargs.has_key("colors"):
            V = len(kwargs["colors"]) - 1
        contourf(yrange, zrange, data2, V, extend='both', **kwargs)
        gca().set_xlim([yrange[0][0], yrange[0][-1]])
        gca().set_ylim([zrange[0][0], zrange[-1][0]])

    # Colorbar.
    if len(gcf().axes) > 1:
        gcf().axes[1].clear()
        cax = gcf().axes[1]
        colorbar(cax = cax)



def profile2DJ(config, data, yy=-1, zz =-1, dd=-1, theta=-1,
               lat = -1, alt = -1, days = -1, time_angle=-1,
               NLevels = 10, V = None, **kwargs):
    """
    Displays an ALTITUDE-LATITUDE section of photolysis rate tabulation
    (4D array) with filled contours. dd=0 or days=0 :
    averaged values over the all year.
    """
    # Some Errors' messages.
    if data.ndim != 4:
        raise Exception, "Function \"profileJ\" proceeds with," \
            + " 4D data, but input data has " + \
            str(data.ndim) + " dimension(s)."

    if(yy == -1 and zz == -1 and dd == -1 and theta == -1 \
       and lat == -1 and alt == -1 and days == -1 and time_angle == -1):
        raise Exception, "Function \"profilecf\" cannot display" \
            + " 4D maps. Missing one argument among  yy, zz, " \
            + "dd, theta, lat, alt, days, time_angle."

    if yy != -1 and lat != -1:
        raise Exception, "You cannot specify a \"yy\" value" \
            + " and a \"lat\" value at the same time."

    if zz != -1 and alt != -1:
        raise Exception, "You cannot specify a \"zz\" value " \
            + "  and a \"alt\" value at the same time."

    if dd != -1 and days != -1:
        raise Exception, "You cannot specify a \"dd\" value " \
            + " and a \"days\" value at the same time."

    if theta != -1 and time_angle != -1:
        raise Exception, "You cannot specify a \"theta\" value " \
            + " and a \"time_angle\" value at the same time."

    arg_num = 0
    arg_list = [yy, zz, dd, theta, lat, alt, days, time_angle]
    for item in arg_list:
        if item != -1:
            arg_num += 1

    if arg_num != 2:
        raise Exception, "Cannot execute \"profilecf\"." \
            + " Two argument among \"yy\", \"zz\", \"dd\"," \
            + " \"theta\", \"lat\", \"alt\", \"days\", " \
            + " \"time_angle\" should be specified."

    # Reads configuration file 'config'.
    if isinstance(config, str):
        if os.path.isfile(config):
            config = talos.Config(config)
        else:
            raise Exception, "Cannot read configuration " \
                + "file : \" " \
                + str(config) + " \"."
    else:
        raise Exception, "First argument has to be a config file's name."

    y_min = config.y_min
    time_angle_min = config.time_angle_min
    day_min = config.day_min
    Delta_y = config.Delta_y
    Delta_time_angle = config.Delta_time_angle
    Delta_t = config.Delta_t
    Ny = config.Ny
    Nz = config.Nz
    Ntheta = config.Ntheta
    Delta_t = config.Delta_t
    Ndays = config.Ndays
    levels = config.Jlevels

    y_max = y_min + Delta_y * Ny
    time_angle_max = time_angle_min + Delta_time_angle * Ntheta
    z_max = levels[-1]
    day_max=Ndays

    # Errors messages

    if(yy > data.shape[1] or (lat !=-1 and (lat < y_min or lat > y_max))):
        raise Exception, "Cannot display profile. Data's maximum y_" \
            + "coordinate is " + str(data.shape[1]) + " and data's" \
            + " latitude is in [" + str(y_min) + "; " + str(y_max) + "]" \
            + " (degrees)."

    if(zz > data.shape[0] or alt > z_max):
        raise Exception, "Cannot display profile. Data's maximum z_" \
            + "coordinate is " + str(data.shape[0]) + " and data's" \
            + " maximum altitude is " +  str(z_max) + " (m)."

    if(theta > data.shape[2] or (time_angle !=-1 and
                                 (time_angle < time_angle_min or
                                  time_angle > time_angle_max))):
        raise Exception, "Cannot display profile. Data's maximum " \
            + "time_angle coordinate is " + str(data.shape[2]) \
            + " and data's maximum time angle is " + str(time_angle_max) \
            + " (degrees)."

    if(dd > data.shape[3] or days > day_max):
        raise Exception, "Cannot display profile. data's maximum dd_" \
            + "coordinate is " + str(data.shape[3]) + " and data's" \
            + " maximum days is " +  str(day_max) + "."

    # Transformation of lon, lat, alt in xx, yy, zz if necessary.

    if lat != -1:
        yy = int((lat - y_min) / Delta_y)

    if alt != -1:
        zz = 0
        for z_index in range(len(levels)):
            if alt > levels[z_index]:
                zz = z_index + 1

    if time_angle != -1:
        theta = int((time_angle - time_angle_min) / Delta_time_angle)

    if days != -1:
        dd = int((days - day_min) / Delta_t)

    if (days == 0 or dd ==0):
        dd_sum=0
        for i in range(data.shape[3]):
            dd_sum = dd_sum+data[:,:,:,i]
        data[:,:,:,0] = dd_sum / float(data.shape[3])

    # Initialization of the map.

    x_min = 0
    Nx = 1
    Delta_x = 1
    mref = Basemap(projection = 'cyl',
                   llcrnrlon = x_min - Delta_x / 2.,
                   llcrnrlat = y_min,
                   urcrnrlon = x_min + Delta_x / 2. + Delta_x * float(Nx - 1),
                   urcrnrlat = y_min + Delta_y * float(Ny - 1),
                   resolution = 'l', suppress_ticks = False,
                   area_thresh = 1000)

    # If the figure is empty, sets new axes.
    if len(gcf().axes) == 0:
        xsize = rcParams['figure.figsize'][0]
        fig_num = get_current_fig_manager().num
        fig = figure(num = fig_num)
        ax = fig.add_axes([0.1, 0.1, 0.75, 0.75])
        axes(ax)
        axes([0.875, 0.1, 0.05, 0.75])

    # Clears current image.
    gcf().axes[0].clear()
    axes(gcf().axes[0])

    # Display of profiles.

    if(dd != -1 and theta != -1):
        data2 = data[:,:,theta,dd]
        if V is None and NLevels != 1:
            V = arange(NLevels) / float(NLevels - 1) *\
                (data2.max() - data2.min()) + data2.min()
        yrange, zrange = gridprofile(mref, -1, data2.shape[1],  levels)

        if kwargs.has_key("colors"):
            V = len(kwargs["colors"]) - 1
        contourf(yrange, zrange, data2, V, extend='both', **kwargs)

    if(dd == 0 and theta != -1):
        data2 = data[:,:,theta,0]
        if V is None and NLevels != 1:
            V = arange(NLevels) / float(NLevels - 1) *\
                (data2.max() - data2.min()) + data2.min()
        yrange, zrange = gridprofile(mref, -1, data2.shape[1],  levels)

        if kwargs.has_key("colors"):
            V = len(kwargs["colors"]) - 1
    contourf(yrange, zrange, data2, V, extend='both', **kwargs)

    # Colorbar.
    if len(gcf().axes) > 1:
        gcf().axes[1].clear()
        cax = gcf().axes[1]
        colorbar(cax=cax)


def cbar():
    """
    Displays a colorbar.
    """
    if len(gcf().axes) > 1:
        gcf().axes[1].clear()
        cax = gcf().axes[1]
    else:
        cax = axes([0.875, 0.1, 0.05, 0.75])
    colorbar(cax = cax)

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


from pylab import *


def segplot(x, y, fmt, maxdelta, **kwargs):
    """
    Plots x versus y, breaking the plot at any point where x[i] -
    x[i-1] > maxdelta. kwargs are passed on to plot.

    @type x: sequence of float or int
    @param x: Data to plot in absciss.
    @type y: sequence of float or int
    @param y: Data to plot in ordinates
    @type fmt: string
    @param fmt: Line style and color, combined in a single format string, as
    in 'bo' for blue circles. See Matplotlib 'plot' command for more details.
    @type maxdelta: float or int
    @param maxdelta: Maximum delta between two consecutive x values for
    which a line should be drawn.
    @type kwargs: keyword argument list
    @param kwargs: The **kwargs can be used to set line properties
    (any property that has a set_* method).  You can use this to set
    a line label (for auto legends), linewidth, anitialising, marker
    face color, etc. See Matplotlib 'plot' command documentation for
    more details.

    @rtype: matplotlib.lines
    @return: The list of plotted lines.
    """
    x = asarray(x)
    y = asarray(y)
    d = diff(x)
    lines = []
    ind = nonzero(greater(d, maxdelta))[0]
    ind = ind+1
    if not len(ind):
        lines.extend( plot(x,y,fmt,**kwargs) )
    else:
        allind = [0]
        allind.extend(ind)
        allind.append(len(x))
        for i1,i2 in zip(allind[:-1], allind[1:]):
            lines.extend( plot(x[i1:i2], y[i1:i2], fmt, **kwargs) )
    return lines


def segplot_logx(x, y, fmt, maxdelta, **kwargs):
    """
    Plots log(x) versus y, breaking the plot at any point where x[i] -
    x[i-1] > maxdelta. kwargs are passed on to plot.

    @type x: sequence of float or int
    @param x: Data to plot in abscissa.
    @type y: sequence of float or int
    @param y: Data to plot in ordinates
    @type fmt: string
    @param fmt: Line style and color, combined in a single format string, as
    in 'bo' for blue circles. See Matplotlib 'plot' command for more details.
    @type maxdelta: float or int
    @param maxdelta: Maximum delta between two consecutive x values for
    which a line should be drawn.
    @type kwargs: keyword argument list
    @param kwargs: The **kwargs can be used to set line properties
    (any property that has a set_* method).  You can use this to set
    a line label (for auto legends), linewidth, anitialising, marker
    face color, etc. See Matplotlib 'plot' command documentation for
    more details.

    @rtype: matplotlib.lines
    @return: The list of plotted lines.
    """
    x = asarray(x)
    y = asarray(y)
    d = diff(x)
    lines = []
    ind = nonzero(greater(d, maxdelta))[0]
    ind = ind+1
    if not len(ind):
        lines.extend( semilogx(x,y,fmt,**kwargs) )
    else:
        allind = [0]
        allind.extend(ind)
        allind.append(len(x))
        for i1,i2 in zip(allind[:-1], allind[1:]):
            lines.extend( semilogx(x[i1:i2], y[i1:i2], fmt, **kwargs) )
    return lines

def segplot_date(x, y, fmt, maxdelta, **kwargs):
    """
    Plots x versus y with dates, breaking the plot at any point where
    x[i] - x[i-1] > maxdelta. kwargs are passed on to plot

    @type x: sequence of dates represented as float days
    @param x: Dates to plot in absciss.
    @type y: sequence of float or int
    @param y: y values at those dates.
    @type fmt: string
    @param fmt: Line style and color, combined in a single format string, as
    in 'bo' for blue circles. See Matplotlib 'plot' command for more details.
    @type maxdelta: float or int
    @param maxdelta: Maximum delta between two consecutive x values for
    which a line should be drawn.
    @type kwargs: keyword argument list
    @param kwargs: The **kwargs can be used to set line properties
    (any property that has a set_* method).  You can use this to set
    a line label (for auto legends), linewidth, anitialising, marker
    face color, etc. See Matplotlib 'plot' command documentation for
    more details.

    @rtype: matplotlib.lines
    @return: The list of plotted lines.
    """
    x = asarray(x)
    y = asarray(y)
    d = diff(x)
    lines = []
    ind = nonzero(greater(d, maxdelta))[0]
    ind = ind+1
    if not len(ind):
        lines.extend(plot_date(x,y,fmt,**kwargs) )
    else:
        allind = [0]
        allind.extend(ind)
        allind.append(len(x))
        for i1,i2 in zip(allind[:-1], allind[1:]):
            lines.extend(plot_date(x[i1:i2], y[i1:i2], fmt, **kwargs) )
    return lines


def set_style1(lines):
    """
    Sets parameters for specified lines, style #1.
    Style is red, continuous line, 1.5 width, antialiased.

    @type lines: matplotlib.lines
    @param lines: Lines to set this predefined style to.
    """
    set(lines, \
        antialiased = True, \
        color = 'r', \
        linestyle = '-', \
        linewidth = 1.5)
    grid(True)
    yearsFmt = DateFormatter('%d/%m/%y')
    axes().xaxis.set_major_locator(WeekdayLocator())
    axes().xaxis.set_minor_locator(DayLocator())
    axes().xaxis.set_major_formatter(yearsFmt)
    labels = axes().get_xticklabels()
    set(labels, rotation=60)
    draw()


def set_style2(lines):
    """
    Sets parameters for specified lines, style #2.
    Style is blue, discontinuous line, 1.0 width, antialiased.

    @type lines: matplotlib.lines
    @param lines: Lines to set this predefined style to.
    """
    set(lines, \
        antialiased = True, \
        color = 'b', \
        linestyle = '--', \
        linewidth = 1.0)
    grid(True)
    yearsFmt = DateFormatter('%d/%m/%y')
    axes().xaxis.set_major_locator(WeekdayLocator())
    axes().xaxis.set_minor_locator(DayLocator())
    axes().xaxis.set_major_formatter(yearsFmt)
    labels = axes().get_xticklabels()
    set(labels, rotation=60)
    draw()


def set_style_fromconfig(config, section, lines):
    """
    Searches for style variables in specified section of configstream,
    and applies the style to the given lines.
    Style must have following entries :
    antialiased (True/False)
    color (r,y,...)
    linestyle (-, --, ...)
    linewidth (float)
    grid (True/False)
    date_format (date format, ie %d/%m/%y)
    labels_rotation (integer)
    For more details about style options, see matplotlib reference page.

    @type config: ConfigStream
    @param config: The configstream containing the style definition.
    @type section: string
    @param section: Name of the section in ConfigStream in which the style is
    defined.
    @type lines: matplotlib.lines
    @param lines: Lines to apply the loaded style to.
    """
    antialiased_arg = config.GetElement("antialiased", section)
    color_arg = config.GetElement("color", section)
    linestyle_arg = config.GetElement("linestyle", section)
    linewidth_arg = config.GetElement("linewidth", section)
    grid_arg = config.GetElement("grid", section)
    date_format = config.GetElement("date_format", section)
    labels_rotation = config.GetElement("labels_rotation", section)
    set(lines, \
        antialiased = antialiased_arg, \
        color = color_arg, \
        linestyle = linestyle_arg, \
        linewidth = linewidth_arg)
    grid(grid_arg)
    yearsFmt = DateFormatter(date_format)
    axes().xaxis.set_major_locator(WeekdayLocator())
    axes().xaxis.set_minor_locator(DayLocator())
    axes().xaxis.set_major_formatter(yearsFmt)
    labels = axes().get_xticklabels()
    set(labels, rotation=labels_rotation)
    draw()


def dstep(config, date):
    """
    Returns the time index of a given date in a data array defined by a
    configuration file.

    @type config: Config or string
    @param config: The configuration or the configuration file.
    @type date: string or datetime.datetime

    @rtype: int
    @return: The time index corresponding to 'date' in the data array defined
    in 'config'.
    """
    if isinstance(config, str):
        config = talos.Config(config)
    if isinstance(date, (float, int)):
        date = str(int(date))
    if isinstance(date, str):
        date = config.stream.StringToDateTime(date)

    delta = date - config.t_min

    if date >= config.t_min:
        return int((delta.days * 24 * 3600 + delta.seconds
                    + config.Delta_t * 1800) / (config.Delta_t * 3600))
    else:
        return int((delta.days * 24 * 3600 + delta.seconds
                    - config.Delta_t * 1800) / (config.Delta_t * 3600))

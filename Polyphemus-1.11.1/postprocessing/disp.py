#!/usr/bin/env python

# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet, Vincent Picavet
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

"""
Displays Polair3D output concentrations and observations at a given station.
"""


import sys, os, datetime, optparse
from numpy import *
from atmopy import *


### Initializations.

parser = optparse.OptionParser(usage = "%prog configuration_file")
(options, args) = parser.parse_args()

if not args:
    parser.error("A configuration file is required.")

additional_content = [("concentrations", "[output]", "String"), \
                      ("paired", "[output]", "Bool"), \
                      ("daily_basis", "[output]", "Bool"), \
                      ("scatter", "[output]", "Int"), \
                      ("multiple_files", "[input]", "Bool"), \
                      ("meas_style", "[output]", "StringList"), \
                      ("sim_style", "[output]", "StringList"), \
                      ("[legend]", "", "legend", "StringSection")]
config = talos.Config(sys.argv[1], additional_content = additional_content)

origin = (config.t_min, config.y_min, config.x_min)
delta = (config.Delta_t, config.Delta_y, config.Delta_x)
shape = (config.Nt, config.Ny, config.Nx)


### Reads input data.

# Observations.
station = io.load_station(config.station_file, \
                          config.station_file_type, config.station)
obs_dates, obs = io.load_file_observations(station.name, config.obs_dir)
obs_dates, obs = observation.remove_missing(obs_dates, obs, (-999, 0))

if config.concentrations == "daily" \
       and not config.daily_basis \
       and len(obs) != 0:
    dates_daily, obs_daily  = observation.split_into_days(obs_dates, obs)
    obs_dates = []
    obs_d = []
    for i in range(len(dates_daily)):
        # Complete day.
        if len(dates_daily[i]) == int(24 / delta[0]):
            obs_dates.append(dates_daily[i][0])
            obs_d.append(obs_daily[i].mean())
    obs = array(obs_d, 'd')

# Computed.

if config.multiple_files:
    Nsim = len(config.file_list)
    files = config.file_list
else:
    Nsim = 1
    files = [config.input_file]
sim = []
for i in range(Nsim):
    # Reads the concentrations and extracts the first level.
    tmp = io.load_binary(files[i], \
                         [config.Nt, config.Nz, config.Ny, config.Nx])[:,0,:]
    # Extracts simulated data.
    sim_tmp = observation.get_simulated_at_station(origin, \
                                                   delta, tmp, station)
    if config.concentrations == "daily":
        sim_dates_ref = \
                      observation.get_simulation_dates(origin[0], delta[0], shape[0])
        dates_daily, sim_daily  = \
                     observation.split_into_days(sim_dates_ref, sim_tmp)
        sim_d = []
        sim_dates_ref = []
        for i in range(len(dates_daily)):
            # Complete day.
            if len(dates_daily[i]) == int(24 / delta[0]):
                sim_dates_ref.append(dates_daily[i][0])
                sim_d.append(sim_daily[i].mean())
        sim_tmp = array(sim_d, 'd')
        sim_dates = sim_dates_ref
    else:
        sim_dates = observation.get_simulation_dates(origin[0], delta[0], shape[0])
    sim.append(sim_tmp)

# Removes concentrations outside the considered period.
obs_dates, obs = \
           observation.restrict_to_period(obs_dates, obs, config.t_range)
for i in range(Nsim):
    tmp, sim[i] = \
         observation.restrict_to_period(sim_dates, sim[i], config.t_range)
sim_dates = tmp

if config.concentrations == "peak" and not config.paired:
    obs_dates, obs = observation.get_daily_peaks(obs_dates, obs, [11, 17])
    for i in range(Nsim):
        tmp, sim[i] = observation.get_daily_peaks(sim_dates, sim[i], [11, 17])
    sim_dates = tmp
if config.concentrations == "peak" and config.paired:
    sim_mask, obs_mask \
              = observation.masks_to_common_dates(sim_dates, obs_dates)
    obs = obs[obs_mask]
    obs_dates = [d for d, b in zip(obs_dates, obs_mask) if b]
    sim_dates = [d for d, b in zip(sim_dates, sim_mask) if b]
    for i in range(Nsim): sim[i] = sim[i][obs_mask]
    # Here is the restriction that prevents the use of multiple files with
    # paired peaks.
    dates, sim[0], obs \
           = observation.get_daily_obs_peaks(dates, sim[0], obs, \
                                             [11, 17], paired = True)
    sim_dates = dates[:]
    obs_dates = dates[:]
if config.concentrations == "peak":
    sim_dates = [observation.midnight(x) for x in sim_dates]
    obs_dates = [observation.midnight(x) for x in obs_dates]


### Visualization.

from pylab import *
from matplotlib.dates import WeekdayLocator, DayLocator, DateFormatter

clf()

# Sets styles.
if hasattr(config, "meas_style"):
    obs_fmt = config.meas_style[0]
    if len(config.meas_style) == 2:
        obs_width = float(config.meas_style[1])
    else:
        obs_width = 0.5
else:
    obs_fmt = 'r-'
    obs_width = 0.5
if hasattr(config, "sim_style"):
    line = " ".join(config.sim_style).split('&')
    options = [x.split() for x in line]
    # If there is not enough options.
    if len(line) < Nsim or min([len(x) for x in options[:Nsim]]) == 0:
        sim_fmt = [options[0][0] for i in range(Nsim)]
        if len(options[0]) == 1:
            sim_width = [obs_width for i in range(Nsim)]
        else:
            sim_width = [float(options[0][1]) for i in range(Nsim)]
        # A single style leads to a single entry in the legend.
        if hasattr(config, "legend"):
            config.legend = []
    else: # Enough options.
        sim_fmt = [x[0] for x in options[:Nsim]]
        sim_width = []
        for x in options[:Nsim]:
            if len(x) > 1:
                sim_width.append(float(x[1]))
            else:
                sim_width.append(obs_width)
else:
    sim_fmt = ['b-' for i in range(Nsim)]
    sim_width = [obs_width for i in range(Nsim)]
    # A single style leads to a single entry in the legend.
    if hasattr(config, "legend"):
        config.legend = []

# Multiple files are not supported in scatter plots.
if config.scatter != 0:
    dates, sim[0], obs = \
           observation.restrict_to_common_dates(sim_dates, sim[0], \
                                                obs_dates, obs)
    plot(obs, sim[0], obs_fmt, linewidth = obs_width)

    if len(config.y_range) != 1:
        lim_min = config.y_range[0]
        lim_max = config.y_range[1]
    else:
        limits = axis()
        lim_min = min(limits[0], limits[2])
        lim_max = min(limits[1], limits[3])

    if config.scatter == 3:
        plot([lim_min, lim_max], [lim_min, lim_max], 'k--')
        plot([lim_min, lim_max], \
             [lim_min, lim_min + (lim_max - lim_min) / 2.], 'k:')
        plot([lim_min, lim_min + (lim_max - lim_min) / 2.], \
             [lim_min, lim_max], 'k:')

    if len(config.y_range) != 1:
        axis([config.y_range[0], config.y_range[1], \
              config.y_range[0], config.y_range[1]])
    elif config.scatter == 2 or config.scatter == 3:
        axis([lim_min, lim_max, lim_min, lim_max])

    xlabel(r'$\rm{Measurements,}\ \ \mu g\ \cdot\ m^{-3}$')
    ylabel(r'$\rm{Polyphemus,}\ \ \mu g\ \cdot\ m^{-3}$')

# Not a scatter plot.
else:
    ax = subplot(1, 1, 1)
    sim_lines = []
    for i in range(Nsim):
        sim_lines.append(display.segplot_date(date2num(sim_dates), sim[i], \
                                              sim_fmt[i], config.Delta_t, \
                                              linewidth = sim_width[i]))
    obs_lines = display.segplot_date(date2num(obs_dates), obs, \
                                     obs_fmt, config.Delta_t, \
                                     linewidth = obs_width)
    if not config.multiple_files:
        legend((sim_lines[0][0], obs_lines[0]), \
               ("Polyphemus", "Measurements"))
    else:
        if not hasattr(config, "legend") or len(config.legend) != Nsim:
            legend([sim_lines[0], obs_lines[0]], \
                   ["Polyphemus", "Measurements"])
        else:
            legend([x[0] for x in sim_lines] + [obs_lines[0]], \
                   config.legend + ["Measurements"])

    # Formats the date axis.
    xlim(date2num(config.t_range[0]), date2num(config.t_range[1]))
    if len(config.y_range) != 1:
        ylim(config.y_range[0], config.y_range[1])
    Ndays = (config.t_range[1] - config.t_range[0]).days
    if Ndays < 62:
        locator = DayLocator(interval = int(ceil(float(Ndays) / 7.)))
        format = DateFormatter("%d %b")
    else:
        locator = MonthLocator()
        format = DateFormatter("%b")
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(format)

    ylabel(r'$\mu g\ \cdot\ m^{-3}$')

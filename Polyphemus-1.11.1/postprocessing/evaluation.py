#!/usr/bin/env python

# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet
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
Computes performances of a chemistry-transport model.
"""


import sys, os, datetime, optparse, inspect
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
                      ("measure", "[output]", "StringList"), \
                      ("cutoff", "[output]", "Float"), \
                      ("select_station", "[output]", "StringList"), \
                      ("output", "[output]", "StringList"), \
                      ("ratio", "[output]", "Float")]
config = talos.Config(sys.argv[1], additional_content)

origin = (config.t_min, config.y_min, config.x_min)
delta = (config.Delta_t, config.Delta_y, config.Delta_x)
shape = (config.Nt, config.Ny, config.Nx)

# Airbase files are preprocessed in the same way as EMEP files.
if config.station_file_type == "airbase":
    config.station_file_type = "emep"

# Statistical measures.
measures = []
config_measure = config.measure[:]
epa_ratio = None
if "epa" in config_measure:
    i = config_measure.index("epa")
    # Checks that "epa" is followed by three numbers in [0, 1].
    try:
        epa_ratio = {"mnbe": float(config_measure[i + 1]),
                     "mnge": float(config_measure[i + 2]),
                     "upa": float(config_measure[i + 3])}
        for j in range(4):
            config_measure.pop(i)
    except:
        raise Exception, \
              "In measures, \"epa\" must be followed by three numbers."
    if min(epa_ratio.values()) < 0. or max(epa_ratio.values()) > 1.:
        raise Exception, \
              "In measures, \"epa\" must be followed by numbers in [0, 1]."
if "all" in config_measure:
    measures = talos.get_module_functions(stat.measure, (1, 2, 3))
    measures = ["meas_mean", "sim_mean"] + measures
else:
    for meas in config_measure:
        if meas in ["meas_mean", "sim_mean"]:
            measures.append(meas)
            continue
        if meas.lower() not in dir(stat) \
               or not inspect.isfunction(getattr(stat, meas.lower())):
            raise Exception, meas + " is not a valid statistical measure."
        Nargs = len(inspect.getargspec(getattr(stat, meas.lower()))[0])
        if Nargs in (1, 2, 3):
            measures.append(meas)
        else:
            raise Exception, meas + " cannot be handled: it takes " \
                  + str(Nargs) + " arguments instead of 1, 2 or 3."

# Formats.
formats_nb = {"mnbe": "%.1f", "nmb": "%.1f", "nme": "%.1f", "mnge": "%.1f",
              "upa": "%.1f", "correlation": "%.1f", "determination": "%.1f",
              "epa": "%.1f", "rmse": "%.1f", "meas_mean": "%.1f",
              "sim_mean": "%.1f"}
formats_ratio = {"mnbe": 100., "nmb": 100., "nme": 100., "mnge": 100.,
                 "upa": 100., "correlation": 100., "determination": 100.,
                 "epa": 100.}
formats_end = {"mnbe": '%', "nmb": '%', "nme": '%', "mnge": '%', "upa": '%',
               "correlation": '%', "determination": '%', "epa": '%'}

def format(meas, nb):
    format_nb = "%.2f"
    format_ratio = 1.
    format_end = ''
    meas = meas.lower()
    if meas in globals()["formats_nb"].keys():
        format_nb = globals()["formats_nb"][meas]
    if meas in globals()["formats_ratio"].keys():
        format_ratio = globals()["formats_ratio"][meas]
    if meas in globals()["formats_end"].keys():
        format_end = globals()["formats_end"][meas]
    return format_nb % (format_ratio * nb) + format_end

# Reads computed concentrations.
sim_dates_ref = \
              observation.get_simulation_dates(origin[0], delta[0], shape[0])
# Reads the computed concentrations and extracts the first level.
sim_ref = io.load_binary(config.input_file, \
                         [config.Nt, config.Nz, config.Ny, config.Nx])[:,0,:]

if len(config.select_station) == 1 and config.select_station[0] == "single":
    stations = [io.load_station(config.station_file, \
                                config.station_file_type, config.station)]
else:
    stations = \
             io.load_stations(config.station_file, config.station_file_type, \
                              origin[1:], delta[1:], shape[1:])
if len(config.select_station) == 2:
    stations = [x for x in stations if getattr(x, config.select_station[0]) \
                == config.select_station[1]]


### Main loop.

# Results.
output = []
functions = []
# Output file.
if "file" in config.output:
    i = config.output.index("file")
    output_file = open(config.output[i+1], 'w')
else:
    output_file = None

for station in stations:

    ### Formats observations and computed concentrations.

    # Observations.
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

    # Extracts simulated data.
    sim = \
        observation.get_simulated_at_station(origin, delta, sim_ref, station)

    if config.concentrations == "daily":
        sim_dates_ref = \
                      observation.get_simulation_dates(origin[0], delta[0], shape[0])
        dates_daily, sim_daily  = \
                     observation.split_into_days(sim_dates_ref, sim)
        sim_d = []
        sim_dates_ref = []
        for i in range(len(dates_daily)):
            # Complete day.
            if len(dates_daily[i]) == int(24 / delta[0]):
                sim_dates_ref.append(dates_daily[i][0])
                sim_d.append(sim_daily[i].mean())
        sim = array(sim_d, 'd')

    # Removes concentrations outside the considered period.
    obs_dates, obs = \
               observation.restrict_to_period(obs_dates, obs, config.t_range)
    sim_dates, sim = observation.restrict_to_period(sim_dates_ref, \
                                                    sim, config.t_range)

    # For peaks.
    if config.concentrations == "peak" and not config.paired:
        sim_dates, sim = observation.get_daily_peaks(sim_dates, sim, [11, 17])
        obs_dates, obs = observation.get_daily_peaks(obs_dates, obs, [11, 17])
        sim_dates = [observation.midnight(x) for x in sim_dates]
        obs_dates = [observation.midnight(x) for x in obs_dates]

    dates, sim, obs = observation.restrict_to_common_dates(sim_dates, sim, \
                                                           obs_dates, obs)

    if config.concentrations == "peak" and config.paired:
        dates, sim, obs = \
               observation.get_daily_obs_peaks(dates, sim, obs, \
                                               [11, 17], paired = True)


    ### Statistical measures.

    # Are there enough observations?
    # Number of hours or days:
    diff = config.t_range[1] - config.t_range[0]
    if config.concentrations == "peak"  or config.concentrations == "daily":
        length = diff.days
    else:
        length = float(diff.days) * 24. + float(diff.seconds) / 3600.
    if float(len(obs)) / float(length) < config.ratio:
        continue

    results = []
    for meas in measures:
        if meas == "upa":   # Special case: performed per day.
            if config.concentrations == "peak":
                dsim, dobs = sim, obs
            else:   # Hourly concentrations.
                dsim = observation.get_daily_peaks(dates, sim, [11, 17])[1]
                dobs = observation.get_daily_peaks(dates, obs, [11, 17])[1]
            tmp = [(x - y) / y for x, y \
                   in zip(dsim, dobs) if y > config.cutoff]
            if len(tmp) == 0:
                raise Exception, "There are not enough observations at" \
                      + " station " + station.name + " to compute UPA."
            else:
                results.append(array(tmp).mean())
        elif meas == "meas_mean":   # Not in module 'stat'.
            results.append(obs.mean())
        elif meas == "sim_mean":   # Not in module 'stat'.
            results.append(sim.mean())
        else:
            Nargs = len(inspect.getargspec(getattr(stat, meas.lower()))[0])
            if Nargs == 1:
                results.append(getattr(stat, meas.lower())(obs))
            elif Nargs == 2:
                results.append(getattr(stat, meas.lower())(sim, obs))
            else:   # 3 arguments.
                results.append(getattr(stat, meas.lower())(sim, obs,
                                                           config.cutoff))

    output.append(results)

    if "station_names" in config.output or "all_stats" in config.output:
        if station.GetRealName() != station.name:
            name = "\n[" + station.name + "]   " + station.GetRealName()
        else:
            name = "\n" + station.GetRealName()
        talos.print_stdout_file(name, output_file)
    if "all_stats" in config.output:
        for f, r in zip(measures, results):
            talos.print_stdout_file("   " + f + ": " + format(f, r),
                                    output_file)

results = array(output)

if "summary" in config.output:
    talos.print_stdout_file("\n* SUMMARY *", output_file)
    talos.print_stdout_file("\n \ number of stations: %s" % len(results),
                            output_file)
    for i in range(len(measures)):
        talos.print_stdout_file(" \ " + measures[i] + ": " \
                                + format(measures[i], results[:, i].mean()),
                                output_file)
        if (measures[i] == "mnbe" or measures[i] == "mnge" \
            or measures[i] == "upa") and epa_ratio != None:
            res = abs(results[:, i])
            ratio = float(len(res[res <= epa_ratio[measures[i]]])) \
                    / float(len(res))
            talos.print_stdout_file("   EPA agreement ratio: " \
                                    + format("epa", ratio), output_file)

if "file" in config.output:
    output_file.write('\n')
    output_file.close()

print

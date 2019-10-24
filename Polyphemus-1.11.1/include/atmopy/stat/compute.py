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
import datetime
import os, sys
sys.path.insert(0,
                os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
import talos, observation, measure
sys.path.pop(0)


def collect(sim, obs, dates = None, stations = None, period = None,
            stations_out = None):
    """
    Collects data (observations and simulated concentrations) over a given
    period and at a given set of stations.

    @type sim: list of list of 1D-array, or list of 1D-array.
    @param sim: The list (indexed by simulations) of lists (indexed by
    stations) of simulated concentrations, or the list (indexed by stations)
    of simulated concentrations.
    @type obs: list of 1D-array
    @param obs: The list (indexed by stations) of observed concentrations.
    @type dates: list of list of datetime
    @param dates: The list (indexed by stations) of list of dates at which the
    data is defined. Both observations and simulated data (of every ensemble)
    are assumed to be designed at the same dates.
    @type stations: list of Station
    @param stations: The stations at which the concentrations are given.
    @type period: 2-tuple of datetime, or datetime, or list of datetime
    @param period: The period where to select the concentrations (bounds
    included). A single date may be provided. If a list is provided, the
    period is defined by the first and the last dates in the list.
    @type stations_out: list of Station, or Station
    @param stations_out: The station(s) at which the concentrations are
    selected.

    @rtype: (2D-array, 1D-array)
    @return: The simulated concentrations in a 2D-array (simulations x
    concentrations) and the corresponding observed concentrations in a
    1D-array.
    """
    # Initializations.
    if isinstance(sim[0], ndarray):
        sim = (sim, )

    if dates == None:
        dates = [range(len(x)) for x in obs]
        period = None
    elif isinstance(dates[0], datetime.datetime) \
             or isinstance(dates[0], datetime.date):
        dates = (dates, )
    if period == None:
        period = (min([x[0] for x in dates]), max([x[-1] for x in dates]))
    elif isinstance(period, datetime.datetime) \
             or isinstance(period, datetime.date):
        period = (period, period)
    elif len(period) == 1:
        period = (period[0], period[0])
    elif len(period) > 2:
        period = (period[0], period[-1])

    if stations == None:
        stations = range(len(obs))
        stations_out = None
    elif isinstance(stations, observation.Station) \
             or isinstance(stations, str) \
             or isinstance(stations, int):
        stations = (stations, )
    if stations_out == None:
        stations_out = stations
    elif isinstance(stations_out, observation.Station) \
           or isinstance(stations_out, str) \
           or isinstance(stations_out, int):
        stations_out = (stations_out, )

    # Output arrays.
    out_obs = []
    out_sim = [[] for i in range(len(sim))]

    for istation in range(len(stations)):
        if stations[istation] in stations_out:
            # Searches for the first date in the considered period.
            i = 0
            while i < len(dates[istation]) and dates[istation][i] < period[0]:
                i += 1
            while i < len(dates[istation]) \
                      and dates[istation][i] <= period[1]:
                # Observations.
                out_obs.append(obs[istation][i])
                # Simulations.
                for isim in range(len(sim)):
                    out_sim[isim].append(sim[isim][istation][i])
                i += 1

    return array(out_sim), array(out_obs)


def compute_stat(sim, obs, measures, dates = None, stations = None, period =
                 None, stations_out = None, cutoff = None):
    """
    Computes a set of statistical measures for one simulation or for a set of
    simulations.

    @type sim: list of list of 1D-array, or list of 1D-array.
    @param sim: The list (indexed by simulations) of lists (indexed by
    stations) of simulated concentrations, or the list (indexed by stations)
    of simulated concentrations.
    @type obs: list of 1D-array
    @param obs: The list (indexed by stations) of observed concentrations.
    @type dates: list of list of datetime, or list of datetime
    @param dates: The list (indexed by stations) of list of dates at which the
    data is defined. Both observations and simulated data (of every ensemble)
    are assumed to be designed at the same dates. 'dates' may also be a list
    of datetime in case there is a single station. If 'dates' is set to None,
    all dates are included.
    @type stations: list of Station, or Station
    @param stations: The station(s) at which the concentrations are given. If
    it is set to None, all stations will be included.
    @type period: 2-tuple of datetime, or datetime, or list of datetime
    @param period: The period where to select the concentrations (bounds
    included). A single date may be provided. If 'period' is set to None, all
    input dates are included. If a list is provided, the period is defined by
    the first and the last dates in the list.
    @type stations_out: list of Station, or Station, or string
    @param stations_out: The station(s) at which the concentrations are
    selected. If 'stations_out' is set to None, all input dates are included.
    @type cutoff: float, or None
    @param cutoff: The value below (or equal) which data is discarded. This
    filters 'obs' and corresponding 'sim' values. Nothing is filtered if
    'cutoff' is set to None.

    @rtype: dict of array
    @return: The statistical measures are a key of the output dictionary. Each
    value is a 1D-array (indexed by simulations).
    """

    ### Initializations.

    import inspect

    if isinstance(sim[0], ndarray):
        sim = (sim, )

    if dates == None:
        dates = [range(len(x)) for x in obs]
        period = None
    elif isinstance(dates[0], datetime.datetime) \
             or isinstance(dates[0], datetime.date):
        dates = (dates, )
    if isinstance(period, datetime.datetime) \
           or isinstance(period, datetime.date):
        period = (period, period)
    elif period == None:
        period = (min([x[0] for x in dates]), max([x[-1] for x in dates]))

    if isinstance(stations, observation.Station):
        stations = (stations, )
    elif stations == None:
        stations = range(len(obs))
        stations_out = None
    if stations_out == None:
        stations_out = stations

    Nsim = len(sim)

    # Functions to be applied.
    if cutoff == None:
        functions = talos.get_module_functions(measure, (1, 2), measures)
    else:
        functions = talos.get_module_functions(measure, (1, 2, 3), measures)

    ### Statistics.

    stat_all = dict(zip(functions, [[] for f in functions]))

    s, o = collect(sim, obs, dates, stations, period, stations_out)
    for i in range(Nsim):
        for f in functions:
            Nargs = len(inspect.getargspec(getattr(measure, f))[0])
            if Nargs == 1:
                stat_all[f].append(getattr(measure, f)(o))
            elif Nargs == 2:
                stat_all[f].append(getattr(measure, f)(s[i], o))
            else:   # Nargs == 3
                stat_all[f].append(getattr(measure, f)(s[i], o, cutoff))

    # To arrays.
    for k in stat_all.keys():
        if Nsim == 1:
            stat_all[k] = stat_all[k][0]
        else:
            stat_all[k] = array(stat_all[k])

    return stat_all


def compute_stat_step(dates, sim, obs, obs_type, measures, stations = None,
                      period = None, stations_out = None, ratio = 0.,
                      cutoff = None):
    """
    Computes a set of statistical measures for one simulation or for a set of
    simulations, and for all time step.

    @type dates: list of list of datetime, or list of datetime
    @param dates: The list (indexed by stations) of list of dates at which the
    data is defined. Both observations and simulated data (of every ensemble)
    are assumed to be designed at the same dates. 'dates' may also be a list
    of datetime in case there is a single station.
    @type sim: list of list of 1D-array, or list of 1D-array.
    @param sim: The list (indexed by simulations) of lists (indexed by
    stations) of simulated concentrations, or the list (indexed by stations)
    of simulated concentrations.
    @type obs: list of 1D-array
    @param obs: The list (indexed by stations) of observed concentrations.
    @type obs_type: string
    @param obs_type: The type of the concentrations: "hourly" or "peak".
    @type stations: list of Station, or Station
    @param stations: The station(s) at which the concentrations are given. If
    it is set to None, all stations will be included.
    @type period: 2-tuple of datetime, or datetime, or list of datetime
    @param period: The period where to select the concentrations (bounds
    included). A single date may be provided. If 'period' is set to None, all
    input dates are included. If a list is provided, the period is defined by
    the first and the last dates in the list.
    @type stations_out: list of Station, or Station, or string
    @param stations_out: The station(s) at which the concentrations are
    selected. If 'stations_out' is set to None, all input dates are included.
    @type ratio: float
    @param ratio: Minimum ratio of the number of available observations (per
    step) and the number of stations. A step at which the actual ratio is
    below this minimum is discarded.
    @type cutoff: float, or None
    @param cutoff: The value below (or equal) which data is discarded. This
    filters 'obs' and corresponding 'sim' values. Nothing is filtered if
    'cutoff' is set to None.

    @rtype: (list of datetime, dict of array)
    @return: The statistical measures are a key of the output dictionary. Each
    value is a (simulation x step)-array. The dates associated with the steps
    with enough measurements are returned in a list.
    """

    ### Initializations.

    import inspect

    if isinstance(sim[0], ndarray):
        sim = (sim, )

    if obs_type != "hourly" and obs_type != "peak":
        raise Exception, "Concentrations must be hourly concentrations" \
              + " or peaks."

    if isinstance(period, datetime.datetime) \
           or isinstance(period, datetime.date):
        period = (period, period)
    elif period == None:
        period = (min([x[0] for x in dates]), max([x[-1] for x in dates]))

    Nsim = len(sim)
    Nstations = len(sim[0])

    # Functions to be applied.
    if cutoff == None:
        functions = talos.get_module_functions(measure, (1, 2), measures)
    else:
        functions = talos.get_module_functions(measure, (1, 2, 3), measures)

    ### Statistics.

    # List of all steps.
    start_date = period[0]
    end_date = period[-1]
    if obs_type == "hourly":
        range_delta = datetime.timedelta(0, 3600)
        Nsteps = (end_date - start_date).days * 24 \
                 + (end_date - start_date).seconds / 3600 + 1
    else:
        start_date = observation.midnight(start_date)
        end_date = observation.midnight(end_date)
        range_delta = datetime.timedelta(1)
        Nsteps = (end_date - start_date).days + 1
    range_dates = [start_date + x * range_delta for x in range(Nsteps)]

    stat_step = dict(zip(functions,
                         [[[] for i in range(Nsim)] for f in functions]))

    output_dates = []
    for date in range_dates:
        s, o = collect(sim, obs, dates, stations, date, stations_out)
        # Enough observations?
        if float(len(o)) / float(Nstations) < ratio:
            continue
        output_dates.append(date)
        for i in range(Nsim):
            for f in functions:
                Nargs = len(inspect.getargspec(getattr(measure, f))[0])
                if Nargs == 1:
                    stat_step[f][i].append(getattr(measure, f)(o))
                elif Nargs == 2:
                    stat_step[f][i].append(getattr(measure, f)(s[i], o))
                else:   # Nargs == 3
                    stat_step[f][i].append(getattr(measure, f)(s[i], o,
                                                               cutoff))

    # To arrays.
    for k in stat_step.keys():
        if Nsim == 1:
            stat_step[k] = array(stat_step[k][0])
        else:
            stat_step[k] = array(stat_step[k])

    return output_dates, stat_step


def compute_stat_station(sim, obs, measures, dates = None, stations = None,
                         period = None, stations_out = None, cutoff = None):
    """
    Computes a set of statistical measures for one simulation or for a set of
    simulations, at given stations.

    @type sim: list of list of 1D-array, or list of 1D-array.
    @param sim: The list (indexed by simulations) of lists (indexed by
    stations) of simulated concentrations, or the list (indexed by stations)
    of simulated concentrations.
    @type obs: list of 1D-array
    @param obs: The list (indexed by stations) of observed concentrations.
    @type dates: list of list of datetime, or list of datetime
    @param dates: The list (indexed by stations) of list of dates at which the
    data is defined. Both observations and simulated data (of every ensemble)
    are assumed to be designed at the same dates. 'dates' may also be a list
    of datetime in case there is a single station. If 'dates' is set to None,
    all dates are included.
    @type stations: list of Station, or Station
    @param stations: The station(s) at which the concentrations are given. If
    it is set to None, all stations will be included.
    @type period: 2-tuple of datetime, or datetime, or list of datetime
    @param period: The period where to select the concentrations (bounds
    included). A single date may be provided. If 'period' is set to None, all
    input dates are included. If a list is provided, the period is defined by
    the first and the last dates in the list.
    @type stations_out: list of Station, or Station, or string
    @param stations_out: The station(s) at which the concentrations are
    selected. If 'stations_out' is set to None, all input dates are included.
    @type cutoff: float, or None
    @param cutoff: The value below (or equal) which data is discarded. This
    filters 'obs' and corresponding 'sim' values. Nothing is filtered if
    'cutoff' is set to None.

    @rtype: dict of array
    @return: The statistical measures are a key of the output dictionary. Each
    value is a (simulation x station)-array.
    """

    ### Initializations.

    import inspect

    if isinstance(sim[0], ndarray):
        sim = (sim, )

    if isinstance(stations, observation.Station):
        stations = (stations, )
    elif stations == None:
        stations = range(len(obs))
        stations_out = None
    if stations_out == None:
        stations_out = stations

    Nsim = len(sim)

    # Functions to be applied.
    if cutoff == None:
        functions = talos.get_module_functions(measure, (1, 2), measures)
    else:
        functions = talos.get_module_functions(measure, (1, 2, 3), measures)

    ### Statistics.

    stat_station = dict(zip(functions,
                            [[[] for i in range(Nsim)] for f in functions]))

    for station in stations_out:
        s, o = \
           collect(sim, obs, dates, stations, period, station)
        for i in range(Nsim):
            for f in functions:
                Nargs = len(inspect.getargspec(getattr(measure, f))[0])
                if Nargs == 1:
                    stat_station[f][i].append(getattr(measure, f)(o))
                elif Nargs == 2:
                    stat_station[f][i].append(getattr(measure, f)(s[i], o))
                else:   # Nargs == 3
                    stat_station[f][i].append(getattr(measure, f)(s[i], o,
                                                                  cutoff))

    # To arrays.
    for k in stat_station.keys():
        if Nsim == 1:
            stat_station[k] = array(stat_station[k][0])
        else:
            stat_station[k] = array(stat_station[k])

    return stat_station

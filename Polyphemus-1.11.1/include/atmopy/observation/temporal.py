# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vincent Picavet, Vivien Mallet
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


class Period:
    """Stores a time period."""
    start = datetime.datetime(1,1,1)
    end = datetime.datetime(1,1,1,1)

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        """ Returns Period's presentation."""
        return str(self.start) + " -> " + str(self.end)


def get_period(dates):
    """
    Returns the period englobing all dates in dates. dates
    must be sorted.

    @type dates: sequence of datetime.datetime
    @param dates: Sequence of dates for which to find the englobing period.
    @rtype: Period
    @return: Period englobing all dates in dates.
    """
    return Period(dates[0], dates[-1])


def get_periods(start, end, length=datetime.timedelta(1), \
                interperiod = datetime.timedelta(0), \
                fit_last = False):
    """
    Returns a list of periods of given length (except last period
    if fit_last is True), beginning at given start date, ending
    at given end date, and with an interperiod time between them.
    By default, interperiod is null and length of period is one day.

    @type start: datetime.datetime
    @param start: Start of desired periods sequence.
    @type end: datetime.datetime
    @param end: End of desired periods sequence.
    @type length: datetime.timedelta
    @param length: Length of one period (1 day by default).
    @type interperiod: datetime.timedelta
    @param interperiod: Interval of time between two consecutive
    periods (0 by default).
    @type fit_last: Boolean
    @param fit_last: If last period does not match the end date and
    fit_last is True, then the last period is shorter than the others.
    Otherwise the last period is of given length and is inside the
    start-end time interval (a time interval not covered with any period
    of the result sequence is therefore present)
    @rtype: sequence of Period
    @return: A sequence of periods, covering the given start-end interval.
    """
    periods = []
    if length != datetime.timedelta(0):
        in_period = False
        last_current = start
        current = start + length
        while current <= end:
            if in_period == True:

                last_current = current
                current = current + length
                in_period = not in_period
            else:
                periods.append(Period(last_current, current))
                last_current = current
                current = current + interperiod
                in_period = not in_period
        if not in_period and fit_last and last_current != end:
            periods.append(Period(last_current,end))
    return periods


def split_into_days(dates, data):
    """
    Gets a sequence of arrays which store the values for each day.

    @type dates: sequence of datetime.datetime
    @param dates: Sequence of dates corresponding to the given data.
    @type data: 1D numpy.array
    @param data: Array of data to split into arrays for days.
    @rtype: sequence of dates sequences, sequence of numpy.array
    @return: A list of dates sequences, splitted by days, and the
    corresponding sequence of numpy.array.
    """
    if len(data) != len(dates):
        raise ValueError, "Data and dates are not of the same length."

    if len(data) == 0:
        return [numpy.array([])], [[]]

    output_data = [[data[0]]]
    output_dates = [[dates[0]]]

    day = dates[0].date()
    for i in range(1, len(dates)):
        if dates[i].date() == day:
            output_data[-1].append(data[i])
            output_dates[-1].append(dates[i])
        else:
            output_data.append([data[i]])
            output_dates.append([dates[i]])
            day = dates[i].date()

    return output_dates, map(lambda x: numpy.array(x), output_data)


def get_daily_obs_peaks(dates, sim, obs, hour_range = [0, 23], \
                        nb_range_min = 24, nb_min = 0, paired = False):
    """
    Returns the daily peaks for both observations and computed
    concentrations.

    @type dates: list of datetime
    @param dates: The dates at which the concentrations are provided.
    @type sim: numpy.array
    @param sim: The simulated concentrations.
    @type obs: numpy.array
    @param obs: Observations.
    @type hour_range: list or tuple with two elements
    @param hour_range: Range of hours over which the peak is to be sought. If
    not all hours of 'hour_range' are available in the observations, the day
    is discarded. The peak in the simulated concentrations is also taken in
    'hour_range'.
    @type nb_range_min: int
    @param nb_range_min: The minimum number of available observations in the
    'hour_range' so that the daily peak should be included.
    @type nb_min: int
    @param nb_min: The minimum number of available observations in the day so
    that the daily peak should be included.
    @type paired: Boolean
    @param paired: True if observations and simulated peaks are assumed to be
    paired (i.e., occuring at the same hour), False otherwise.

    @rtype: (list of datetime, numpy.array, list of datetime,
    numpy.array) or (list of datetime, numpy.array, numpy.array)
    @return: The new dates for computed concentrations, the computed peaks,
    the corresponding measured peaks preceded by their own dates.
    """
    nb_range_min = min(nb_range_min, hour_range[1] - hour_range[0] + 1)
    tmp_dates = dates[:]
    dates, sim = split_into_days(tmp_dates, sim)
    dates, obs = split_into_days(tmp_dates, obs)
    output_dates = []
    if not paired:
        output_dates_sim = []
    output_sim = []
    output_obs = []
    for i in range(len(dates)):
        if len(dates[i]) < nb_min:
            continue
        j = 0
        while j < len(dates[i]) and dates[i][j].hour < hour_range[0]:
            j += 1
        if j == len(dates[i]):
            continue
        tmp_dates = []
        tmp_sim = []
        tmp_obs = []
        while j < len(dates[i]) and dates[i][j].hour <= hour_range[1]:
            tmp_dates.append(dates[i][j])
            tmp_sim.append(sim[i][j])
            tmp_obs.append(obs[i][j])
            j += 1
        if len(tmp_dates) < nb_range_min:
            continue
        j = numpy.array(tmp_obs).argmax()
        output_dates.append(tmp_dates[j])
        output_obs.append(tmp_obs[j])
        if paired:
            output_sim.append(tmp_sim[j])
        else:
            j = numpy.array(tmp_sim).argmax()
            output_dates_sim.append(tmp_dates[j])
            output_sim.append(tmp_sim[j])
    if paired:
        return output_dates, numpy.array(output_sim), \
               numpy.array(output_obs)
    else:
        return output_dates_sim, numpy.array(output_sim), \
               output_dates, numpy.array(output_obs)


def get_daily_peaks(dates, conc, hour_range = [0, 23], \
                    nb_range_min = 24, nb_min = 0):
    """
    Returns the daily peaks.

    @type dates: list of datetime
    @param dates: The dates at which the concentrations are provided.
    @type conc: numpy.array
    @param conc: The simulated concentrations.
    @type hour_range: list or tuple with two elements
    @param hour_range: Range of hours over which the peak is to be sought. If
    not all hours of 'hour_range' are available in the observations, the day
    is discarded. The peak in the simulated concentrations is also taken in
    'hour_range'.
    @type nb_range_min: int
    @param nb_range_min: The minimum number of available observations in the
    'hour_range' so that the daily peak should be included.
    @type nb_min: int
    @param nb_min: The minimum number of available observations in the day so
    that the daily peak should be included.

    @rtype: (list of datetime, numpy.array)
    @return: The concentration peaks preceded by their dates.
    """
    nb_range_min = min(nb_range_min, hour_range[1] - hour_range[0] + 1)
    dates, conc = split_into_days(dates, conc)
    output_dates = []
    output_conc = []
    for i in range(len(dates)):
        if len(dates[i]) < nb_min:
            continue
        j = 0
        while j < len(dates[i]) and dates[i][j].hour < hour_range[0]:
            j += 1
        if j == len(dates[i]):
            continue
        tmp_dates = []
        tmp_conc = []
        while j < len(dates[i]) and dates[i][j].hour <= hour_range[1]:
            tmp_dates.append(dates[i][j])
            tmp_conc.append(conc[i][j])
            j += 1
        if len(tmp_dates) < nb_range_min:
            continue
        j = numpy.array(tmp_conc).argmax()
        output_dates.append(tmp_dates[j])
        output_conc.append(tmp_conc[j])
    return output_dates, numpy.array(output_conc)


def mask_for_common_days(sim_dates, simulated, obs_dates, obs):
    """
    Creates masks for simulated data and observation data
    corresponding to data of common dates.

    @type sim_dates: sequence of datetime.datetime
    @param sim_dates: Dates for simulation data.
    @type simulated: 1D numpy.array
    @param simulated: Simulation data, 1D array
    @type obs_dates: sequence of datetime.datetime
    @param obs_dates: Dates for observation data.
    @type obs: 1D numpy.array
    @param obs: Observation data, 1D array.
    @rtype: numpy.array, numpy.array
    @return: Array masks for simulated and observation data. In each array,
    null values indicates that the corresponding data value has no common date
    in the other array, whereas a 1 value indicates the opposite.
    """
    if len(simulated) != len(sim_dates) or len(obs) != len(obs_dates):
        print len(simulated), len(sim_dates), len(obs), len(obs_dates)
        raise ValueError, "Incompatible dimensions!"

    # Dates.
    sim_dates_float = list(sim_dates)
    map(lambda x: float(x.toordinal()) + float(x.hour) / 24., sim_dates_float)
    obs_dates_float = list(obs_dates)
    map(lambda x: float(x.toordinal()) + float(x.hour) / 24., obs_dates_float)

    # Selection.
    sim_condition = numpy.zeros(len(simulated))
    obs_condition = numpy.zeros(len(obs))

    # Selects the common days.
    j = 0
    for i in range(len(sim_dates_float)):
        while j < len(obs_dates_float) \
                  and obs_dates_float[j] < sim_dates_float[i]:
            j += 1
        if j == len(obs_dates_float):
            break
        if obs_dates_float[j] == sim_dates_float[i]:
           sim_condition[i] = 1
           obs_condition[j] = 1

    return sim_condition, obs_condition


def apply_mask_for_common_days(sim_dates, simulated, obs_dates, obs, \
                               mask_sim, mask_obs):
    """ Applies a mask returned by mask_for_common_days on
    simulated and observation data, and gets corresponding dates.

    @type sim_dates: sequence of datetime.datetime
    @param sim_dates: Dates for simulation data.
    @type simulated: 1D numpy.array
    @param simulated: Simulation data, 1D array
    @type obs_dates: sequence of datetime.datetime
    @param obs_dates: Dates for observation data.
    @type obs: 1D numpy.array
    @param obs: Observation data, 1D array.
    @type mask_sim: 1D numpy.array
    @param mask_sim: Mask array for simulation data corresponding to dates
    common to simulation and observation. This array can be obtain thanks to
    mask_for_common_days.
    @type mask_obs: 1D numpy.array
    @param mask_obs: Mask array for observation data corresponding to dates
    common to simulation and observation. This array can be obtain thanks to
    mask_for_common_days.
    @rtype: sequence of datetime.datetime, numpy.array, numpy.array
    @return: List of dates common to observation and simulation data, and
    corresponding data arrays for simulation and observation.
    """
    if len(simulated) != len(sim_dates) or len(obs) != len(obs_dates):
        print len(simulated), len(sim_dates), len(obs), len(obs_dates)
        raise ValueError, "Incompatible dimensions!"

    common_dates = []
    # Selects the common days.
    for i in mask_sim:
        common_dates.append(sim_dates[i])

    return common_dates, simulated[numpy.where(mask_sim)], \
           obs[numpy.where(mask_obs)]


def restrict_to_common_dates(sim_dates, simulated, obs_dates, obs):
    """
    Gets items from data and dates so as to keep only
    dates and corresponding data which are both in observations and
    simulated data.

    @type sim_dates: sequence of datetime.datetime
    @param sim_dates: Dates for simulation data.
    @type simulated: 1D numpy.array
    @param simulated: Simulation data, 1D array.
    @type obs_dates: sequence of datetime.datetime
    @param obs_dates: Dates for observation data.
    @type obs: 1D numpy.array
    @param obs: Observation data, 1D array.
    @rtype: sequence of datetime.datetime, numpy.array, numpy.array
    @return: List of dates common to observation and simulation data, and
    corresponding data arrays for simulation and observation.
    """

    if len(simulated) != len(sim_dates) or len(obs) != len(obs_dates):
        print len(simulated), len(sim_dates), len(obs), len(obs_dates)
        raise ValueError, "Incompatible dimensions!"

    dates = list(sim_dates)
    sim_condition = numpy.zeros(len(simulated))
    obs_condition = numpy.zeros(len(obs))

    for i in range(len(dates)):
       try:
           ind = obs_dates.index(dates[i])
           sim_condition[i] = 1
           obs_condition[ind] = 1
       except ValueError:
           pass

    for i in range(len(sim_condition)-1, -1, -1):
        if sim_condition[i] == 0:
            dates.pop(i)

    return dates, simulated[numpy.where(sim_condition)], \
           obs[numpy.where(obs_condition)]


def masks_for_common_dates(dates0, dates1):
    """
    Computes masks to be applied so that data sets available at 'dates0' and
    'dates1' may be defined at the same dates.

    @type dates0: list of datetime
    @param dates0: The first list of dates.
    @type dates1: list of datetime
    @param dates1: The second list of dates.

    @rtype: (numpy.array(type=Bool), numpy.array(type=Bool))
    @return: The masks are returned for both lists in Boolean arrays. There
    are common dates wherever a Boolean is True.
    """
    N0 = len(dates0)
    N1 = len(dates1)

    # Masks: lists of Booleans.
    mask0 = numpy.zeros(N0, "Bool")
    mask1 = numpy.zeros(N1, "Bool")

    i0 = 0
    i1 = 0

    while i0 < N0:
        while i1 < N1 and dates1[i1] < dates0[i0]:
            i1 += 1
        if i1 == N1:
            break
        if dates1[i1] == dates0[i0]:
            mask0[i0] = True
            mask1[i1] = True
        i0 += 1

    return mask0, mask1


def restrict_to_common_days(sim_dates, simulated, obs_dates, obs):
    """
    Gets items from data and dates so as to keep only
    daily data which are both in observations and
    simulated data. Dates lists and data have to be sorted (by increasing
    time).

    @type sim_dates: sequence of datetime.datetime
    @param sim_dates: Dates for simulation data.
    @type simulated: 1D numpy.array
    @param simulated: Simulation data, 1D array.
    @type obs_dates: sequence of datetime.datetime
    @param obs_dates: Dates for observation data.
    @type obs: 1D numpy.array
    @param obs: Observation data, 1D array.
    @rtype: sequence of datetime.datetime, numpy.array, numpy.array
    @return: List of dates common to observation and simulation data, and
    corresponding data arrays for simulation and observation.
    """

    if len(simulated) != len(sim_dates) or len(obs) != len(obs_dates):
        print len(simulated), len(sim_dates), len(obs), len(obs_dates)
        raise ValueError, "Incompatible dimensions!"

    # Dates.
    sim_dates_iso = list(sim_dates)
    map(lambda x: x.date().isoformat(), sim_dates_iso)
    obs_dates_iso = list(obs_dates)
    map(lambda x: x.date().isoformat(), obs_dates_iso)

    # Selection.
    sim_condition = numpy.zeros(len(simulated))
    obs_condition = numpy.zeros(len(obs))

   # Gets the start date index for observations.
    obs_start = 0
    while(obs_dates_iso[obs_start] < sim_dates_iso[0]):
        obs_start = obs_start + 1

    common_dates = []
    old_ind = 0
    # Loops over observation dates and mark indices of matching dates.
    for i in range(obs_start, len(obs_dates_iso)):
        ind = old_ind
        while(sim_dates_iso[ind] < obs_dates_iso[i] \
              and ind < len(sim_dates_iso) - 1 ):
            ind = ind + 1
        if sim_dates_iso[ind] == obs_dates_iso[i]:
            sim_condition[ind] = 1
            obs_condition[i] = 1
            old_ind = ind
            common_dates.append(sim_dates[ind])
        else:
            old_ind = ind - 1

    return common_dates, simulated[numpy.where(sim_condition)], \
           obs[numpy.where(obs_condition)]


def restrict_to_common_days2(sim_dates, simulated, obs_dates, obs):
    """
    Gets items from data and dates so as to keep only
    daily data which are both in observations and
    simulated.

    @type sim_dates: sequence of datetime.datetime
    @param sim_dates: Dates for simulation data.
    @type simulated: 1D numpy.array
    @param simulated: Simulation data, 1D array.
    @type obs_dates: sequence of datetime.datetime
    @param obs_dates: Dates for observation data.
    @type obs: 1D numpy.array
    @param obs: Observation data, 1D array.
    @rtype: sequence of datetime.datetime, numpy.array, numpy.array
    @return: List of dates common to observation and simulation data, and
    corresponding data arrays for simulation and observation.
    """

    if len(simulated) != len(sim_dates) or len(obs) != len(obs_dates):
        print len(simulated), len(sim_dates), len(obs), len(obs_dates)
        raise ValueError, "Incompatible dimensions!"

    # Dates.
    sim_dates_iso = list(sim_dates)
    map(lambda x: x.date().isoformat(), sim_dates_iso)
    obs_dates_iso = list(obs_dates)
    map(lambda x: x.date().isoformat(), obs_dates_iso)

    # Selection.
    sim_condition = numpy.zeros(len(simulated))
    obs_condition = numpy.zeros(len(obs))

    # Selects the common days.
    for i in range(len(sim_dates_iso)):
       try:
           ind = obs_dates_iso.index(sim_dates_iso[i])
           sim_condition[i] = 1
           obs_condition[ind] = 1
       except ValueError:
           pass

    # Extracts the common dates.
    common_dates = sim_dates
    for i in range(len(sim_condition)-1, -1, -1):
        if sim_condition[i] == 0:
            common_dates.pop(i)

    return common_dates, simulated[numpy.where(sim_condition)], \
           obs[numpy.where(obs_condition)]


def restrict_to_period(dates, data, period_date, end_date = None):
    """
    Returns data and associated dates within a given period.

    @type dates: list of datetime
    @param dates: The dates associated with data.
    @type data: array
    @param data: Array of data.
    @type period_date: Period, list of datetime, or datetime
    @param period_date: Defines:
       0. a period (Period object);
       1. a period through its bounds (list of datetime);
       2. the first date of the selected period.
    @type end_date: datetime
    @param end_date: the last date of the selected period (if not provided by
    'period_date').

    @rtype: (list of datetime, array)
    @return: The dates and data over the selected period.
    """
    condition = numpy.zeros(len(dates))
    if isinstance(period_date, Period):
        start_date = period_date.start
        end_date = period_date.end
    elif isinstance(period_date, (list, tuple)):
        start_date = period_date[0]
        end_date = period_date[-1]
    else:
        start_date = period_date
    istart = 0
    while istart < len(dates) and dates[istart] < start_date: istart += 1
    if istart == len(dates) or dates[istart] > end_date :
        return [], numpy.array([])
    iend = istart + 1
    while iend < len(dates) and dates[iend] <= end_date: iend += 1
    return dates[istart:iend], data[istart:iend]


def mask_for_series(dates, delta, Ndates):
    """
    Removes the dates if there are not enough preceding contiguous dates.

    @type dates: list of datetime
    @param dates: Input dates.
    @type delta: timedelta
    @param delta: The minimum time-delta between two dates.
    @type Ndates: integer
    @param Ndates: The number of contiguous dates that must be available
    before a given date so that the latter should be kept. Two dates are
    contiguous if the time delta between them is less than (or equal to)
    'delta'.

    @rtype: array
    @return: An array of Boolean filled with True for all dates to be kept.
    """

    mask = numpy.ones(len(dates), "Bool")

    if len(dates) == 0:
        return mask

    # Number of preceding contiguous dates.
    count = 0
    # Previous date.
    prev_date = dates[0]
    if Ndates > 0:   # The first date is then removed.
        mask[0] = False
    for idate in range(1, len(dates)):
        if dates[idate] - prev_date <= delta:
            count += 1
        else:
            count = 0
        if count < Ndates:
            mask[idate] = False
        prev_date = dates[idate]

    return mask


def remove_incomplete_days(dates, data):
    """
    Removes dates and data from the first and/or last days of the period if
    there is missing data in these days. There is missing data in a day if
    there is not as many timesteps as possible. The timestep is assumed to be
    constant and is computed with the number of hours between the first two
    timesteps.

    @type dates: list of datetime
    @param dates: The dates associated with data.
    @type data: array
    @param data: Array of data.

    @rtype: (list of datetime, array)
    @return: The dates and data with the first and/or last days removed.
    """
    if len(dates) < 2:
        return dates, data

    # Time step.
    delta = dates[1] - dates[0]

    # If the time step is greater than one day, there cannot be missing data.
    if delta.days != 0:
        return dates, data

    # Timestep in hours.
    delta = delta.seconds / 3600
    # Number of steps per day.
    steps = int(24. / float(delta))

    # Number of steps in the first day.
    i = 1
    while i < len(dates) and dates[i].date() == dates[0].date():
        i += 1
    steps_first = i
    # In case there is only one day...
    if steps_first == len(dates):
        if steps_first != steps:
            return [], numpy.array([])
        return dates, data
    # First valid output-step.
    ind_first = 0
    if steps_first != steps:
        ind_first = steps_first

    # Number of steps in the last day.
    i = len(dates) - 2
    while i >= 0 and dates[i].date() == dates[-1].date():
        i -= 1
    steps_last = len(dates) - i - 1
    # Last valid output-step.
    ind_last = len(dates)
    if steps_last != steps:
        ind_last = -steps_last

    return dates[ind_first:ind_last], data[ind_first:ind_last]


def remove_days(dates, data, days):
    """
    Removes data in the first days.

    @type dates: list of datetime
    @param dates: The dates associated with 'data'.
    @type data: array
    @param data: The data to be filtered.
    @type days: integer
    @param days: The number of days to be removed at the beginning.

    @rtype: (list of datetime, array)
    @return: The dates and associated data without the first days.
    """
    i = 0
    while i < len(dates) and (dates[i] - dates[0]).days < days:
        i += 1
    return dates[i:], data[i:]


def midnight(date):
    """
    Move to midnight in the current day. Midnight is assumed to be
    the beginning of the day.

    @type date: datetime.datetime
    @param date: The date in the current day.
    @rtype: datetime.datetime
    @return: Midnight in the current day specified by date.
    """
    return date - datetime.timedelta(0, 3600 * date.hour \
                                     + 60 * date.minute + date.second, \
                                     date.microsecond)


def timedelta2num(delta):
    """
    Converts datetime.timedelta to float day number as used in Matplotlib.

    @type delta: datetime.timedelta
    @param delta: The time-delta to be converted in days.

    @rtype: float
    @return: The number of days in 'delta'.
    """
    if delta < datetime.timedelta(0):
        num = -(date2num(datetime.datetime(1,1,1) - delta) \
              - date2num(datetime.datetime(1,1,1)))
    else:
        num = date2num(datetime.datetime(1,1,1) + delta) \
              - date2num(datetime.datetime(1,1,1))
    return num


def get_simulation_dates(t_min, delta_t, Nt):
    """
    Gets a list of dates corresponding to the simulation data.

    @type t_min: datetime.datetime
    @param t_min: t_min is the time reference for simulation data.
    @type delta_t: int
    @param delta_t: The time delta between two consecutive time values,
    in hours.
    @type Nt: int
    @param Nt: Number of time values.
    @rtype: sequence of datetime.datetime
    @return: A sequence of time values corresponding to the simulation data.
    """
    sim_dates = []
    for i in range(Nt):
        sim_dates.append(t_min + datetime.timedelta(hours = i * delta_t))
    return sim_dates


def remove_missing(dates, data, rm_value = -999):
    """
    Removes given values from a data array and removes the corresponding
    dates.

    @type dates: list of datetime
    @param dates: The dates at which data is provided.
    @type data: 1D numpy.array
    @param data: The data array to be filtered.
    @type rm_value: float or list of floats
    @param rm_value: The value(s) to be removed from 'data'.

    @rtype: (list of datetime, 2D numpy.array)
    @return: The data array and its dates without the specified values.
    """
    if isinstance(rm_value, (list, tuple)):
        for x in rm_value:
            dates, data = remove_missing(dates, data, x)
        return dates, data
    condition = data != rm_value
    data = data[condition]
    dates = [d for d, c in zip(dates, condition) if c]
    return dates, data

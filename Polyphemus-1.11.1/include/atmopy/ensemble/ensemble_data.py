# Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
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


import sys, os
current_file = os.path.abspath(__file__)
atmopy_path = os.path.split(os.path.dirname(current_file))[0]
atmopy_path = os.path.split(os.path.dirname(atmopy_path))[0]
sys.path.insert(0, atmopy_path)
from atmopy import talos, io, observation, stat
sys.path.pop(0)

from numpy import *
import datetime
import scipy


################
# ENSEMBLEDATA #
################


class EnsembleData:
    """
    This class stores ensemble simulations and observations, at given
    stations. When all is loaded, it mainly contains:
       0. Nsim: number of simulations in the ensemble;
       1. sim: list (for simulations) of list (for stations) of arrays
       (concentrations at a given station);
       2. Nstation: number of stations;
       3. station: list of Station instances;
       4. obs: observations at stations;
       5. date: dates of observations;
       6. all_dates: list of all dates in the time period starting with the
       first observation and ending with the last observation;
       7. stat (possibly): global statistics;
       8. stat_step (possibly): statistics per time step.
       9. stat_station (possibly): statistics per station.
    """


    def __init__(self, configuration_file = None, verbose = False):
        """
        If a configuration file is provided, ensemble simulations and
        observations are read.

        @type configuration_file: string
        @param configuration_file: The path to the configuration file that
        describes ensemble simulations and observations.
        @type verbose: Boolean
        @param verbose: Should information be displayed on screen?
        """
        self.verbose = verbose
        self.prt = talos.PrintInPlace()
        self.stat = {}
        self.stat_step = {}
        self.stat_station = {}
        self.Nsim = 0
        self.sim = []

        if configuration_file != None:
            self.LoadConfiguration(configuration_file)
            self.LoadStation()
            self.LoadObservation()
            self.LoadSimulation()


    def LoadConfiguration(self, configuration_file):
        """
        Loads the configuration and checks it. No simulation data or
        observations are read.

        @type configuration_file: string
        @param configuration_file: The path to the configuration file that
        describes ensemble simulations and observations.
        """
        import os
        if not os.path.isfile(configuration_file):
            raise Exception, "Unable to open configuration file " \
                  + "\"" + configuration_file + "\"."

        self.configuration_file = configuration_file

        add_content = [("discarded_cells", "[input]", "Int"),
                       ("discarded_days", "[input]", "Int"),
                       ("concentrations", "[output]", "String"),
                       ("measure", "[output]", "StringList"),
                       ("cutoff", "[output]", "Float"),
                       ("select_station", "[output]", "StringList"),
                       ("paired", "[output]", "Bool"),
                       ("ratio", "[output]", "Float")]

        self.config = talos.Config(configuration_file,
                                   additional_content = add_content)

        # Filters config.measure if it exists.
        try:
            if "all" in self.config.measure:
                self.config.measure = "all"
        except:
            pass

        self.Nsim = len(self.config.file_list)

        self.CheckConfiguration()


    def CheckConfiguration(self):
        """
        Checks that the configuration is valid.
        """
        # Checks that the considered period is included in the simulated
        # period. This is a restrictive constraint to be sure that the user
        # knows what he is doing.
        if self.config.t_range[0] < self.config.origin[0] \
               or self.config.t_range[1] > self.config.origin[0] \
               + datetime.timedelta(0, 3600 * self.config.Delta_t
                                    * self.config.Nt):
            raise Exception, "The period considered for computations must " \
                  + "be included in the simulated period."

        # Checks that the required concentrations are supported.
        if self.config.concentrations == "peak":
            # Checks that peaks are not paired.
            if self.config.paired:
                raise Exception, "Unable to deal with paired peaks."
        elif self.config.concentrations != "hourly":
            raise Exception, "Field \"concentrations\" is set to \"" \
                  + self.config.concentrations \
                  + "\" but should be \"hourly\" or \"peak\"."


    def LoadStation(self):
        """
        Loads information about possibly involved stations.
        """
        if len(self.config.select_station) == 1 \
               and self.config.select_station[0] == "single":
            self.station = [io.load_station(self.config.station_file,
                                            self.config.station_file_type,
                                            self.config.station)]
        else:
            self.station = io.load_stations(self.config.station_file,
                                            self.config.station_file_type,
                                            self.config.origin[1:],
                                            self.config.Delta[1:],
                                            self.config.shape[1:])
        # Filters stations.
        if len(self.config.select_station) == 2:
            self.station = [x for x in self.station \
                            if getattr(x, self.config.select_station[0]) \
                            == self.config.select_station[1]]

        self.Nstation = len(self.station)


    def LoadObservation(self):
        """
        Loads observations.
        """
        # Output attributes.
        self.obs = []
        self.date = []

        station_restricted = []

        for istation in range(self.Nstation):
            station = self.station[istation]

            if self.verbose:
                self.prt(str(istation) + "/" + str(len(self.station))
                         + " " + station.real_name)

            ### Observations (temporary vectors).

            # Initial import.
            obs_date, obs = io.load_file_observations(station.name,
                                                      self.config.obs_dir)
            obs_date, obs = observation.remove_missing(obs_date,
                                                       obs, (-999, 0))
            # Removes concentrations outside the considered period.
            obs_date, obs = \
                      observation.restrict_to_period(obs_date, obs,
                                                     self.config.t_range)
            # Peaks.
            if self.config.concentrations == "peak":
                obs_date, obs = observation.get_daily_peaks(obs_date,
                                                            obs, [11, 17], 2)
                obs_date = [observation.midnight(x) for x in obs_date]

            ### Enough measurements?

            diff = self.config.t_range[1] - self.config.t_range[0]
            if self.config.concentrations == "peak":
                # Maximum number of peaks.
                length = float(diff.days)
            else:
                # Maximum number of hourly observations.
                length = float(diff.days) * 24. + float(diff.seconds) / 3600.
            if float(len(obs)) / length < self.config.ratio :
                continue   # Not enough observations...
            # Adds the stations in the output-station list.
            station_restricted.append(station)

            ### Stores the new concentrations.
            self.obs.append(array(obs))
            self.date.append(obs_date)

        if self.verbose:
            # Displays information.
            self.prt.Clear()
            print "Number of stations:", \
                  str(len(station_restricted)) + "/" + str(len(self.station))
            print "Number of observations:", \
                  array([len(x) for x in self.obs]).sum()

        # Update the station list.
        self.station = station_restricted
        self.Nstation = len(self.station)

        self.GetAllDates()


    def LoadSimulation(self):
        """
        Loads ensemble results.
        """
        # Output attribute.
        self.sim = [[] for i in range(self.Nsim)]

        mask = [None for i in range(self.Nstation)]

        # Simulation dates.
        date = observation.get_simulation_dates(self.config.t_min,
                                                self.config.Delta_t,
                                                self.config.Nt)

        for ifile in range(self.Nsim):

            if self.verbose:
                name = self.config.file_list[ifile].split("/")
                if len(name) < 2:
                    name = self.config.file_list[ifile].upper()
                else:
                    name = name[-2].upper()
                self.prt(str(ifile) + "/" + str(self.Nsim) + " " + name)

            # Reads computed concentrations and extracts the first level.
            sim_ref = io.load_binary(self.config.file_list[ifile],
                                     [self.config.Nt, self.config.Nz,
                                      self.config.Ny, self.config.Nx])[:, 0]

            for istation in range(self.Nstation):
                station = self.station[istation]
                # Extracts simulated data.
                sim = observation.get_simulated_at_station(self.config.origin,
                                                           self.config.Delta,
                                                           sim_ref, station)

                # Removes concentrations outside the considered period.
                sim_date, sim = \
                          observation.restrict_to_period(date, sim,
                                                         self.config.t_range)
                if self.config.concentrations == "peak":
                    # Peaks.
                    sim_date, sim = observation.get_daily_peaks(sim_date,
                                                                sim, [11, 17],
                                                                2)
                    sim_date = [observation.midnight(x) for x in sim_date]

                # Time selections.
                if mask[istation] == None:   # Masks are not computed.
                    # Hourly concentrations.
                    function = observation.masks_for_common_dates
                    mask[istation], tmp = function(sim_date,
                                                   self.date[istation])

                # Adds the new concentrations in the dedicated arrays.
                self.sim[ifile].append(array(sim[mask[istation]]))

        self.prt.Clear()


    def AddSimulation(self, sim, date = None, position = None,
                      duplicate = False):
        """
        Adds a new simulation.

        @type sim: list of 1D-array
        @param sim: The list (indexed by stations) of simulation data to be
        added.
        @type date: None or list of list of datetime
        @param date: If not None, the list (indexed by stations) of the list
        of dates at which the values in 'data' are defined. If 'date' is not
        None, consistency with dates in the ensemble is checked.
        @type position: None or integer
        @param position: index of the simulation in the ensemble. If
        'position' is None, the simulation is appended.
        @type duplicate: Boolean
        @param duplicate: True if the concentrations should be duplicated,
        False otherwise (copy by reference).
        """
        if date != None:
            if len(date) != len(self.date):
                raise Exception, "Inconsistent dates."
            for istation in range(len(date)):
                if len(date[istation]) != len(self.date[istation]):
                    raise Exception, "Inconsistent dates."
                for idate in range(len(date[istation])):
                    if date[istation][idate] != self.date[istation][idate]:
                        raise Exception, "Inconsistent dates."

        if duplicate:
            add_sim = [x.copy() for x in sim]
        else:
            add_sim = sim

        self.Nsim += 1
        if position is None:
            self.sim.append(add_sim)
        else:
            self.sim.insert(position, add_sim)


    def RemoveSimulation(self, index = -1):
        """
        Removes a simulation.

        @type index: integer
        @param index: The index of the model to be removed. The default index
        is -1 and corresponds to the last simulation.

        @rtype: list of 1D-array
        @return: The list (indexed by stations) of simulated data removed from
        the ensemble.
        """
        self.Nsim -= 1
        return self.sim.pop(index)


    def ComputeStatistics(self, period = None):
        """
        Computes statistics (simulations against measurements) over all
        stations and over the whole time period. It updates attribute "stat".
        """
        self.stat = stat.compute_stat(self.sim, self.obs,
                                      self.config.measure,
                                      cutoff = self.config.cutoff,
                                      period = period, dates = self.date)


    def ComputeStepStatistics(self, period = None):
        """
        Computes statistics (simulations against measurements) over all
        stations and for every time step. It updates attribute "stat_step".
        """
        self.stat_step = \
                       stat.compute_stat_step(self.date, self.sim, self.obs,
                                              self.config.concentrations,
                                              self.config.measure,
                                              cutoff = self.config.cutoff,
                                              period = period)[1]


    def ComputeStationStatistics(self, period = None):
        """
        Computes statistics (simulations against measurements) over all time
        steps and for every station. It updates attribute "stat_station".
        """
        self.stat_station = \
                          stat.compute_stat_station(self.sim, self.obs,
                                                    self.config.measure,
                                                    dates = self.date,
                                                    cutoff
                                                    = self.config.cutoff,
                                                    period = period)


    def GetAllDates(self):
        """
        Finds out all dates within the considered period. These dates are put
        in attribute "all_dates".
        """
        # Period covered by observations.
        self.period = [min([x[0] for x in self.date]),
                       max([x[-1] for x in self.date])]
        if self.config.concentrations == "hourly":
            delta = datetime.timedelta(0, 3600)
        else:   # Peaks.
            delta = datetime.timedelta(1)
        self.all_dates = [self.period[0]]
        while self.all_dates[-1] < self.period[-1]:
            self.all_dates.append(self.all_dates[-1] + delta)


    def RestrictToPeriod(self, period):
        """
        Removes simulated concentrations and observations that are outside a
        given period.

        @type period: list of datetime
        @param period: It defines the selected period through its bounds.
        """
        obs_out = []
        sim_out = [[] for i in range(self.Nsim)]
        date_out = []
        restrict = observation.restrict_to_period
        for i in range(self.Nstation):
            tmp = restrict(self.date[i], self.obs[i], period)
            date_out.append(tmp[0])
            obs_out.append(tmp[1])
        for i in range(self.Nsim):
            for j in range(self.Nstation):
                sim_out[i].append(restrict(self.date[j], self.sim[i][j],
                                           period)[1])
        self.sim = sim_out
        self.obs = obs_out
        self.date = date_out
        self.GetAllDates()


    def DuplicateEnsemble(self, ensemble):
        """
        Duplicates an EnsembleData instance, except its simulation data. The
        ensembles will share the same observations, stations, dates and
        configuration file. They will have their own independent simulation
        data and statistics.

        @type ensemble: EnsembleData
        @param ensemble: The ensemble to be duplicated.
        """
        self.Nstation = ensemble.Nstation
        self.station = ensemble.station
        self.obs = ensemble.obs
        self.date = ensemble.date
        self.all_dates = ensemble.all_dates
        self.configuration_file = ensemble.configuration_file
        self.config = ensemble.config


    def RankArray(self):
        """
        Computes the rank array to be used in a rank diagram.

        @rtype: 1D-array
        @return: The rank array of length Nsim+1.
        """
        result = zeros(self.Nsim + 1, 'd')
        for s in range(self.Nstation):
            for t in range(len(self.obs[s])):
                obs = self.obs[s][t]
                sim = array([x[s][t] for x in self.sim], 'd')
                sim.sort()

                rank = 0
                while rank != self.Nsim and obs > sim[rank]:
                    rank += 1
                result[rank] += 1.

        return result


####################
# USEFUL FUNCTIONS #
####################


def merge(date, data, date_update, data_update):
    """
    Merges two data sets assumed to be defined for the same stations, but not
    for the same dates.

    @type date: list of list of datetime
    @param date: The list (indexed by stations) of the list of dates at which
    the values in 'data' are defined.
    @type data: list of 1D-array.
    @param data: The list (indexed by stations) of data assumed to be the base
    of the output data.
    @type date_update: list of list of datetime
    @param date_update: The list (indexed by stations) of the list of dates
    at which the values in 'data_update' are defined.
    @type data_update: list of 1D-array.
    @param data_update: The list (indexed by stations) of data that should be
    added to 'data'. If data is found in both 'data' and 'data_update' for a
    given date, the latter is used.

    @rtype: list of list of datetime, list of 1D-array
    @return: Updated dates and associated data.
    """
    Nstation = len(data)

    output_date = [[] for i in range(Nstation)]
    output_data = [[] for i in range(Nstation)]

    for istation in range(Nstation):
        i = 0
        i_update = 0
        while i < len(date[istation]) \
                  and i_update < len(date_update[istation]):
            if date[istation][i] < date_update[istation][i_update]:
                output_date[istation].append(date[istation][i])
                output_data[istation].append(data[istation][i])
                i += 1
            else:
                output_date[istation].append(date_update[istation][i_update])
                output_data[istation].append(data_update[istation][i_update])
                if date[istation][i] == date_update[istation][i_update]:
                    i += 1
                i_update += 1
        if i == len(date[istation]):
            if i_update != len(date_update[istation]):
                output_date[istation] += date_update[istation][i_update:]
                output_data[istation] += data_update[istation][i_update:]
        else:
            output_date[istation] += date[istation][i:]
            output_data[istation] += data[istation][i:]

    return output_date, [array(x) for x in output_data]


def add_model(ref, model, ens):
    """
    Adds a model to an ensemble.

    @type ref: EnsembleMethod
    @param ref: The reference model that is used when 'model' has no data.
    @type model: EnsembleMethod
    @param model: The model to be added to the ensemble.
    @type ens: EnsembleData
    @param ens: The ensemble to which 'model' should be added. The model is
    appended in the model list; it therefore becomes the last model.
    """
    date, sim = merge(ref.date, ref.sim, model.date, model.sim)
    ens.AddSimulation(sim, date = date)

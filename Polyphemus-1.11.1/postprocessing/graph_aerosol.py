#!/usr/bin/env python

# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Marilyne Tombette
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
Displays Polair3D output concentrations for aerosols.
Should be launched after init_aerosol.py.
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
                      ("Nbins", "[input]", "Int"), \
                      ("computed", "[input]", "Bool"), \
                      ("Dmin", "[input]", "Float"), \
                      ("Dmax", "[input]", "Float"), \
                      ("with_organics", "[output]", "Bool"), \
                      ("primary", "[output]", "StringList"), \
                      ("organics", "[output]", "StringList"), \
                      ("inorganics", "[output]", "StringList"), \
                      ("primary_names", "[output]", "StringList"), \
                      ("inorganics_names", "[output]", "StringList"), \
                      ("file_bounds", "[input]", "String"), \
                      ("log_plot", "[output]", "Bool"), \
                      ("graph_type", "[output]", "StringList"), \
                      ("sim_style", "[output]", "StringList"), \
                      ("graphs_at_station", "[output]", "Bool"), \
                      ("i_range", "[output]", "IntList"), \
                      ("j_range", "[output]", "IntList"), \
                      ("[directory_list]", "", "directory_list",
                       "StringSection"), \
                      ("bin_index_shift", "[input]", "Int"), \
                      ("[legend]", "", "legend", "StringSection")]
config = talos.Config(sys.argv[1], additional_content = additional_content)

origin = (config.t_min, config.y_min, config.x_min)
delta = (config.Delta_t, config.Delta_y, config.Delta_x)
shape_sim = (config.Nt, config.Ny, config.Nx)

if config.graphs_at_station:
    station = io.load_station(config.station_file, \
                              config.station_file_type, config.station)
else:
    i_range = config.i_range
    j_range = config.j_range

# Aerosols.

if config.with_organics:
    aerosol_species = config.primary + config.inorganics + ["organics"]
    species_names = config.primary_names + config.inorganics_names + ["organics"]
else:
    aerosol_species = config.primary + config.inorganics
    species_names = config.primary_names + config.inorganics_names

Nspecies = len(aerosol_species)

Nsim = len(config.directory_list)
bin_index_shift = config.bin_index_shift

### Visualization.

from pylab import *
from matplotlib.dates import WeekdayLocator, DayLocator, DateFormatter

ccolors = ( "b", "g", "c", "r", "m", "y",
           "#3D0099","#D9BF01","#B280FF", "#E6BFFF", "#CE80FF",
           "#FFBFFF", "#FF80FF", "#FFBFEF", "#FFFF80")


if hasattr(config, "sim_style"):
    line = " ".join(config.sim_style).split('&')
    options = [x.split() for x in line]
    # If there is not enough options.
    if len(line) < Nsim or min([len(x) for x in options[:Nsim]]) == 0:
        sim_fmt = [options[0][0] for i in range(Nsim)]
        if len(options[0]) == 1:
            sim_width = [1 for i in range(Nsim)]
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
                sim_width.append(1)
else:
    sim_fmt = ['b-' for i in range(Nsim)]
    sim_width = [1 for i in range(Nsim)]
    # A single style leads to a single entry in the legend.
    if hasattr(config, "legend"):
        config.legend = []

sim_dates = observation.get_simulation_dates(origin[0], delta[0], shape_sim[0])
Ndates = len(sim_dates)
legend_sim = []
if len(config.legend) == Nsim:
    legend_sim = config.legend
else:
    for i in range(Nsim):
        legend_sim.append("Simulation " + str(i+1))


### Computes data to be displayed.

# Chemical species.
if "pie_chart" in config.graph_type or \
   "bar_chart" in config.graph_type or \
   "temporal" in config.graph_type:
    PM_species = []
    for i in range(Nsim):
        try:
            for ispecies in range(len(aerosol_species)):
                if aerosol_species[ispecies] == "organics":
                    if config.graphs_at_station:
                        total = zeros([config.Nt], dtype = 'd')
                    else:
                        total = zeros([config.Nt,j_range[1]-j_range[0],
                                       i_range[1]-i_range[0]],
                                      dtype = 'd')

                    for species in config.organics:
                        file_in = os.path.join(config.directory_list[i],
                                               species + ".bin")
                        tmp = io.load_binary(file_in, [config.Nt,
                                                       config.Nz,
                                                       config.Ny,
                                                       config.Nx])[:,0,:]
                        if config.graphs_at_station:
                            total += \
                                  observation.get_simulated_at_station(origin,
                                                                       delta,
                                                                       tmp,
                                                                       station)

                        else:
                            total += tmp[:,j_range[0]:j_range[1],
                                         i_range[0]:i_range[1]]
                    PM_species.append(total)
                else:
                    file_in = os.path.join(config.directory_list[i],
                                           aerosol_species[ispecies] + ".bin")
                    tmp = io.load_binary(file_in,
                                         [config.Nt, config.Nz,
                                          config.Ny,config.Nx])[:,0,:]
                    if config.graphs_at_station:
                        PM_species.append(
                            observation.get_simulated_at_station(origin,
                                                                 delta,
                                                                 tmp,
                                                                 station))
                    else:
                        PM_species.append(tmp[:,j_range[0]:j_range[1],
                                              i_range[0]:i_range[1]])
        except:
            raise Exception, \
                  "You have to launch init_aerosol.py with computation "+ \
                  "of total species activated (option total_species)."

# Total mass per bin.
if "mass_distribution" in config.graph_type:
    PM_bins = []
    try:
        for j in range(Nsim):
            for i in range(config.Nbins):
                file_in = os.path.join(config.directory_list[j],
                                       "total_mass_" + str(i+bin_index_shift) + ".bin")
                tmp = io.load_binary(file_in,
                                     [config.Nt, config.Nz,
                                      config.Ny, config.Nx])[:,0,:]
                if config.graphs_at_station:
                    PM_bins.append(
                        observation.get_simulated_at_station(origin,
                                                             delta,
                                                             tmp,
                                                             station))
                else:
                    PM_bins.append(tmp[:,j_range[0]:j_range[1],
                                       i_range[0]:i_range[1]])
    except:
        raise Exception, \
              "You have to launch init_aerosol.py with computation "+ \
              "of total bins activated (option total_bins)."

# Granulometry.
if "number_distribution" in config.graph_type:
    granulo = []
    try:
        for j in range(Nsim):
            for i in range(config.Nbins):
                file_in = os.path.join(config.directory_list[j],
                                       "number_" + str(i+bin_index_shift) + ".bin")
                tmp = io.load_binary(file_in,
                                     [config.Nt, config.Nz,
                                      config.Ny, config.Nx])[:,0,:]
                if config.graphs_at_station:
                    granulo.append(
                        observation.get_simulated_at_station(origin,
                                                             delta,
                                                             tmp,
                                                             station))
                else:
                    granulo.append(tmp[:,j_range[0]:j_range[1],
                                       i_range[0]:i_range[1]])
    except:
        raise Exception, \
              "You have to launch init_aerosol.py with computation "+ \
              "of granulometry activated (option granulometry)."

# Diameters.
if "mass_distribution" in config.graph_type or \
       "number_distribution" in config.graph_type:
    dry_bound_diameter = [] # Dry bound diameters
    dry_section_diameter = [] # Dry section diameters
    if config.computed:
        hx = log(config.Dmax / config.Dmin) / float(config.Nbins)
        # Computes bounds diameters.
        for i in range(config.Nbins + 1):
            dry_bound_diameter.append(config.Dmin * exp(float(i) * hx))
    else:
        f = open(config.file_bounds,'r')
        for i in range(config.Nbins + 1):
            dry_bound_diameter.append(float(f.readline()))

    for i in range(config.Nbins):
        dry_section_diameter.append(sqrt(dry_bound_diameter[i]
                                         * dry_bound_diameter[i+1]))


### Graphs.

if "pie_chart" in config.graph_type:
    for j in range(Nsim):
        figure(figsize = (8,8))
        ax = axes([0.1, 0.1, 0.8, 0.8])
        pm_composition = []
        for i in range(Nspecies):
            sim_species = PM_species[j * Nspecies + i]
            pm_composition.append(sim_species.mean())
        pie(pm_composition, labels = species_names, colors = ccolors,\
            autopct='%0.1f')
        title(legend_sim[j])

if "bar_chart" in config.graph_type:
    figure()
    indbar = arange(Nspecies)
    widthbar = 0.1
    pbar = []
    color_legend = []
    for j in range(Nsim):
        pm_composition = []
        for i in range(Nspecies):
            sim_species = PM_species[j * Nspecies + i]
            pm_composition.append(sim_species.mean())
        pbar = bar(indbar + j * widthbar, pm_composition,
                   widthbar, color = ccolors[j])
        color_legend.append(pbar[0])
    title('Concentrations for each species')
    xlabel('Species')
    xticks(indbar + widthbar, species_names)
    xlim(-widthbar, len(indbar))
    ylabel('Concentrations')
    legend(color_legend, legend_sim)

if "temporal" in config.graph_type:
    for j in range(Nsim):
        color_legend = []
        sim_lines = []
        figure()
        ax = gca()
#        sim_to_plot = zeros(shape = (Ndates), dtype = 'd')
        sim_to_plot = []
        sim_to_plot.append(zeros(shape = (Ndates), dtype = 'd'))
        for i in range(Nspecies):
            sim_species = PM_species[j * Nspecies + i]
            if config.graphs_at_station:
                sim_to_plot.append(sim_to_plot[len(sim_to_plot)-1] \
                                   + sim_species)
            else:
                sim_to_plot.append(sim_to_plot[len(sim_to_plot)-1] \
                                   + mean(mean(sim_species, 2), 1))

            sim_lines.append(display.segplot_date(date2num(sim_dates), \
                                                  sim_to_plot[len(sim_to_plot)-1], \
                                                  '-', 0.5))
            setp(sim_lines[len(sim_lines)-1], color = ccolors[i])
            color_legend.append(sim_lines[len(sim_lines)-1][0])
        legend(color_legend, aerosol_species)
        title(legend_sim[j])

        # Formats the date axis.
        xlim(date2num(config.t_range[0]), date2num(config.t_range[1]))
        if len(config.y_range) != 1:
            ylim(config.y_range[0], config.y_range[1])
        Ndays = (config.t_range[1] - config.t_range[0]).days
        if Ndays < 62:
            locator = DayLocator(interval = ceil(float(Ndays) / 7.))
            format = DateFormatter("%d %b")
        else:
            locator = MonthLocator()
            format = DateFormatter("%b")
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(format)

        ylabel(r'$\mu g\ \cdot\ m^{-3}$')

if "mass_distribution" in config.graph_type:
    figure()
    mass_distr_lines = []
    for j in range(Nsim):
        mass_distribution = []
        for i in range(config.Nbins):
            sim_bin = PM_bins[j * config.Nbins + i]
            mass_distribution.append(sim_bin.mean())
        if config.log_plot:
            mass_distr_line = display.segplot_logx(dry_section_diameter,
                                                   mass_distribution,
                                                   sim_fmt[j], 10)
            mass_distr_lines.append(mass_distr_line)
        else:
            mass_distr_line = display.segplot(dry_section_diameter,
                                              mass_distribution,
                                              sim_fmt[j], 10)
            mass_distr_lines.append(mass_distr_line)
    legend([x[0] for x in mass_distr_lines], legend_sim)
    xlabel('Diameter')
    ylabel('Concentration')

if "number_distribution" in config.graph_type:
    figure()
    num_distr_lines = []
    for j in range(Nsim):
        num_distribution = []
        for i in range(config.Nbins):
            sim_bin = granulo[j * config.Nbins + i]
            num_distribution.append(sim_bin.mean())
        if config.log_plot:
            num_distr_line = display.segplot_logx(dry_section_diameter,
                                                  num_distribution,
                                                  sim_fmt[j], 10)
            num_distr_lines.append(num_distr_line)
        else:
            num_distr_line = display.segplot(dry_section_diameter,
                                             num_distribution, sim_fmt[j], 10)
            num_distr_lines.append(num_distr_line)
    legend([x[0] for x in num_distr_lines], legend_sim)
    xlabel('Diameter')
    ylabel('Number')

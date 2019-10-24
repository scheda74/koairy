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
Initializes fields for aerosols.
"""


import sys, os, datetime, optparse, inspect
from math import *
from numpy import *
from atmopy import *


### Initializations.

parser = optparse.OptionParser(usage = "%prog [configuration file]")
(options, args) = parser.parse_args()

if not args:
    parser.error("A configuration file is required.")

additional_content = [("Nbins", "[input]", "Int"), \
                      ("computed", "[input]", "Bool"), \
                      ("Dmin", "[input]", "Float"), \
                      ("Dmax", "[input]", "Float"), \
                      ("file_bounds", "[input]", "String"), \
                      ("with_organics", "[output]", "Bool"), \
                      ("primary", "[output]", "StringList"), \
                      ("organics", "[output]", "StringList"), \
                      ("inorganics", "[output]", "StringList"), \
                      ("output_species", "[output]", "StringList"), \
                      ("[directory_list]", "", "directory_list",
                       "StringSection"), \
                      ("bin_index_shift", "[input]", "Int")]

config = talos.Config(sys.argv[1], additional_content)

# Dry bound diameters.
dry_bound_diameter = empty(shape = [config.Nbins + 1], dtype = 'd')
# Dry bound masses.
dry_bound_masses = empty(shape = [config.Nbins + 1], dtype = 'd')
# Dry section diameters.
dry_section_diameter = empty(shape = [config.Nbins], dtype = 'd')
# Dry section masses.
dry_section_masses = empty(shape = [config.Nbins], dtype = 'd')

# Aerosol density.
rho_aerosol = 1.4e-06
# Other constants.
pi = 3.14159
cst = pi / 6.

if config.with_organics:
    aerosol_species = config.primary + config.inorganics + config.organics
else:
    aerosol_species = config.primary + config.inorganics

if config.computed:
    delta_diameter = (log(config.Dmax / config.Dmin)) / float(config.Nbins)
    # Computes bounds diameters.
    for i in range(config.Nbins + 1):
        dry_bound_diameter[i] = config.Dmin * exp(float(i) * delta_diameter)
        dry_bound_masses[i] = rho_aerosol * cst * dry_bound_diameter[i] ** 3
else:
    f = open(config.file_bounds, 'r')
    for i in range(config.Nbins + 1):
       dry_bound_diameter[i] = float(f.readline())
       dry_bound_masses[i] = rho_aerosol * cst * dry_bound_diameter[i] ** 3

for i in range(config.Nbins):
    dry_section_diameter[i] = sqrt(dry_bound_diameter[i]
                                   * dry_bound_diameter[i+1])
    dry_section_masses[i] = sqrt(dry_bound_masses[i] * dry_bound_masses[i+1])


# PM limits.
Npm25 = argmin(abs(dry_bound_diameter - 2.5))
Npm10 = argmin(abs(dry_bound_diameter - 10.))

# Options.
compute_pm10 = compute_pm25 = compute_tot_species \
               = compute_bins = compute_granulo = False

compute_pm10 = "PM10" in config.output_species
compute_pm25 = "PM2.5" in config.output_species
compute_tot_species = "total_species" in config.output_species
compute_bins = "total_bins" in config.output_species
compute_granulo = "granulometry" in config.output_species


### Main computations.

directory_list = config.directory_list
Nsim = len(directory_list)
bin_index_shift = config.bin_index_shift

sum_species = zeros([config.Nt, config.Nz, config.Ny, config.Nx], dtype = 'd')

if compute_pm10 or compute_pm25 or compute_tot_species \
       or compute_bins or compute_granulo:
    for j in range(Nsim):
        if compute_pm10 or compute_pm25 or compute_bins \
               or compute_granulo:
            sum_bin = zeros([config.Nbins, config.Nt, config.Nz,
                             config.Ny, config.Nx], dtype = 'd')
        for ispecies in aerosol_species:
            sum_species = zeros([config.Nt, config.Nz, config.Ny, config.Nx],
                                dtype = 'd')
            for i in range(config.Nbins):
                file_s = os.path.join(directory_list[j],
                                      ispecies + "_"+ str(i+bin_index_shift) + ".bin")
                # Reads the concentrations.
                sim = io.load_binary(file_s,
                                     [config.Nt, config.Nz,
                                      config.Ny, config.Nx])
                sum_species += sim
                if compute_pm10 or compute_pm25 or compute_bins \
                       or compute_granulo:
                    sum_bin[i] += sim

            io.save_binary(sum_species,
                           os.path.join(directory_list[j], ispecies + ".bin"))

        # Saves aggregated data.

        if compute_pm10:
            sum_pm10 = sum(sum_bin[0:Npm10], 0)
            io.save_binary(sum_pm10,
                           os.path.join(directory_list[j], "PM10.bin"))

        if compute_pm25:
            sum_pm25 = sum(sum_bin[0:Npm25], 0)
            io.save_binary(sum_pm25,
                           os.path.join(directory_list[j], "PM2.5.bin"))

        if compute_bins:
            for i  in range(config.Nbins):
                io.save_binary(sum_bin[i],
                               os.path.join(directory_list[j], "total_mass_"
                                            + str(i+bin_index_shift) + ".bin"))
        if compute_granulo:
            for i in range(config.Nbins):
                io.save_binary(sum_bin[i] / dry_section_masses[i],
                               os.path.join(directory_list[j],
                                            "number_" + str(i+bin_index_shift) + ".bin"))

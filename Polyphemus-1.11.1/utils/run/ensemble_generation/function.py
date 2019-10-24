# Copyright (C) 2008-2009 EDF R&D
# Author: Damien Garaud
#
# This file is part of the air quality modeling system Polyphemus. It is used
# to generate ensembles.
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

"""\package ensemble_generation.function

A few useful functions.

The functions list:
 - create_output_dir
 - get_bc_init_date
 - check_program
 - get_perturbed_field
"""


def create_output_dir(program, dependency_list, parameter_dict):
    """Creates the path where the binary files will be saved for a given
    program.

    \param program the name of a progam.
    \param dependency_list a list of names of parameter.
    \param parameter_dict the dictionary of the description of a model.
    @return A string.
    """
    result = program + "/"
    general_list = ["luc", "vertical_resolution", "first_layer_height", \
                    "cloud_attenuation", "min_height_cloud", \
                    "critical_relative_humidity", "min_kz", "min_urban_kz", \
                    "apply_vert", "tm_stable", "exponent_p_tm", "sbl_ratio", \
                    "boundary_layer", "deposition_velocity", "Ra", "Rb", \
                    "ground_emission", "mozart_bc"]
    if "None" in dependency_list:
        return result
    # Checks if each dependency is in general_list.
    for d in dependency_list:
        if d not in general_list:
            raise Exception, "The dependency: \"" + d + "\" does not appear" \
                  + " in this function. You must change its name " \
                  + " or add it in this function."
    # Luc
    if "luc" in dependency_list:
        result += parameter_dict['luc'] + "/"
    # Levels
    if "vertical_resolution" in dependency_list:
        result += "Nz-" + parameter_dict["vertical_resolution"] + "/"
    if "first_layer_height" in dependency_list:
        result += "flh-" + parameter_dict["first_layer_height"] + "/"
    # Meteo
    if "cloud_attenuation" in dependency_list:
        result += parameter_dict["cloud_attenuation"] + "/"
    if "min_height_cloud" in dependency_list:
        result += "cloud-" + parameter_dict["min_height_cloud"] + "/"
    if "critical_relative_humidity" in dependency_list:
        result += parameter_dict["critical_relative_humidity"] + "/"
    if "min_kz" in dependency_list:
        result += "min_kz-" + parameter_dict["min_kz"] + "/"
    if "min_urban_kz" in dependency_list:
        result += "min_urban_kz-" + parameter_dict["min_urban_kz"] + "/"
    if "apply_vert" in dependency_list:
        result += "apply_vert-" + parameter_dict["apply_vert"] + "/"
    if "tm_stable" in dependency_list:
        result += "TM_stable-" + parameter_dict["tm_stable"] + "/"
    if "exponent_p_tm" in dependency_list:
        result += "pTM-" + parameter_dict["exponent_p_tm"] + "/"
    if "sbl_ratio" in dependency_list:
        result += "sbl-" + parameter_dict["sbl_ratio"] + "/"
    if "boundary_layer" in dependency_list:
        result += "BL-" + parameter_dict["boundary_layer"] + "/"
    # Deposition
    if "deposition_velocity" in dependency_list:
        result += parameter_dict["deposition_velocity"] + "/"
    if "Ra" in dependency_list:
        result += "Ra-" + parameter_dict["Ra"] + "/"
    if "Rb" in dependency_list:
        result += "Rb-" + parameter_dict["Rb"] + "/"
    # Emissions
    if "ground_emission" in dependency_list:
        result += parameter_dict["ground_emission"] + "/"
    # Mozart bc
    if "mozart_bc" in dependency_list:
        if parameter_dict['mozart_bc'] == "perturb":
            result += parameter_dict['mozart_bc'] + "/"

    return result[0:-1]


def get_bc_init_date(date_min):
    """Returns the minimal date for the program 'bc'.
    @return An instance 'datetime.datetime'.
    """

    import datetime
    init_date = datetime.datetime(date_min.year, date_min.month, date_min.day)
    if init_date.day == 1 and init_date.month == 1:
        init_date = datetime.datetime(init_date.year - 1, 12, 31)
    else:
        year = init_date.year
        Nday = (init_date - datetime.datetime(year, 01, 01)).days
        while ((Nday + 6) % 10 != 0):
            init_date = init_date - datetime.timedelta(1)
            Nday = (init_date - datetime.datetime(year, 01, 01)).days

    return init_date


def check_program(polyphemus_dir, program):
    """Checks if a program exists.

    Throws an Exception if the program is not in the subdirectories
    'processing' and 'preprocessing' in the Polyphemus directory.
    \param polyphemus_dir the Polyphemus directory.
    \param program the name of the program.
    """

    import os
    list_file = os.listdir(polyphemus_dir + "/processing/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/bc/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/bio/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/dep/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/emissions/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/ground/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/ic/")
    list_file += os.listdir(polyphemus_dir + "/preprocessing/meteo/")
    if program not in list_file:
        raise Exception, "The program \"" + program \
              + "\" does not exist in " + polyphemus_dir \
              + "/driver and /preprocessing."


def get_perturbed_field(field, parameter_dict, name):
    """Returns strings lines.

    It is used to replace several variables from the generic perturbation file
    by perturbed fields and lists of species.
    \param field the name of a field which will be perturbed.
    \param parameter_dict the dictionary of the description of a model.
    \param name strings or a list of strings which corresponds to the name of
    the field or the species.
    @return Strings lines.
    """
    if not isinstance(name, list) and not isinstance(name, str):
        raise Exception, "In function::get_perturbed_field, the name: \"" \
              + str(name) + "\" must be a list or strings."
    if not isinstance(parameter_dict, dict):
        raise Exception, "In function::get_perturbed_field, the" \
              + " paramter_dict: \"" + str(parameter_dict) \
              + "\" must be a dictionary."
    if field not in parameter_dict.keys():
        return ""
    else:
        result = ""
        if not isinstance(parameter_dict[field], tuple):
            raise Exception, "The value of the key \"" + field + "\" must " \
                  + "be a tuple."
        value = parameter_dict[field][0]
        type_perturbation = parameter_dict[field][-1]
        if isinstance(name, str):
            name = [name]
        if len(name) == 1:
            return name[0] + " \t " + type_perturbation + " \t " + str(value)
        for s in name:
            result += s + " \t " + type_perturbation + " \t " + str(value) \
                      + "\n"
        return result


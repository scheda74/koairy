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

"""\package ensemble_generation.config.ConfigPolair3D

This module contains the class 'ConfigPolair3D'.

It is used to generate an ensemble for Polair3D.
"""

from run.ensemble_generation import function


class ConfigPolair3D:
    """This class contains a few methods which return dictionaries.

    These dictionaries contain variables which appear in the generic
    configuration files and will be replaced by the values of the parameters.
    """


    def __init__(self):
        """An empty constructor.
        """
        pass


    def GetConfigVariable(self, model_index, ensemble_program):
        """Returns the dictionary with as keys the variables which will be
        replaced in the generic configuration files and as values the values
        of the parameters.

        \param model_index the index of a model.
        \param ensemble_program an instance 'EnsembleProgram'.
        @return A dictionary.
        """

        param_dict = ensemble_program.GetParameterDict(model_index)

        # Vertical distribution and first layer height.
        Nz = param_dict["vertical_resolution"]
        flh = param_dict["first_layer_height"]
        levels = "levels-" + Nz  + "-" + flh  + ".dat"
        # Ground
        luc = param_dict["luc"]
        if luc == "usgs":
            sea_index = 15
            urban_index = 0
        elif luc == "glcf":
            sea_index = 0
            urban_index = 13
        # Meteo
        cloud_attenuation = param_dict["cloud_attenuation"]
        if cloud_attenuation == "radm":
            attenuation_type = 1
        elif cloud_attenuation == "esquif":
            attenuation_type = 2
        min_height_cloud = param_dict["min_height_cloud"]
        critical_relative_humidity = param_dict["critical_relative_humidity"]
        if critical_relative_humidity == "sigma":
            critical_relative_humidity = "1"
        elif critical_relative_humidity == "two_layers":
            critical_relative_humidity = "2"
        # Meteo - Kz
        min_kz = param_dict["min_kz"]
        min_urban_kz = param_dict["min_urban_kz"]
        apply_vert = param_dict["apply_vert"]
        sbl_ratio = param_dict["sbl_ratio"]
        if param_dict["kz"] == "tm":
            kz_file = \
                function.create_output_dir('Kz_TM',
                                           ensemble_program.dependency['Kz_TM'],
                                           param_dict) \
                                           + "/Kz_TM.bin"
        elif param_dict["kz"] == "louis":
            kz_file = \
                function.create_output_dir('Kz_TM',
                                           ensemble_program.dependency['Kz_TM'],
                                           param_dict) \
                                           + "/Kz_Louis.bin"
        bool_tm_stable = param_dict["tm_stable"]
        exponent_p_tm = param_dict["exponent_p_tm"]
        boundary_layer = param_dict["boundary_layer"]
        # Dep
        dep = param_dict["deposition_velocity"]
        Ra = param_dict["Ra"]
        Rb = param_dict["Rb"]
        if Ra == "heat":
            Ra = "fh"
        elif Ra == "moment":
            Ra = "fm"
        elif Ra == "diag":
            Ra = "diag"
        # Emissions
        ground_emission = param_dict["ground_emission"]
        vertical_distribution = "vertical_distribution-" \
                                + ground_emission + "-" + Nz + ".dat"
        if ground_emission == "ground":
            bool_volume_emission = "no"
            Ncells = 0
        elif ground_emission == "2_cells":
            bool_volume_emission = "yes"
            Ncells = 2
        elif ground_emission in ["low", "medium", "high"]:
            bool_volume_emission = "yes"
            vertical_distribution = "vertical_distribution-" \
                                    + ground_emission + "-" + Nz \
                                    + "-" + flh + ".dat"
            if Nz == "5":
                Ncells = 2
            elif Nz == "9":
                if ground_emission == "low":
                    Ncells = 2
                elif ground_emission == "medium":
                    Ncells = 3
                else:
                    Ncells = 4
            else:
                raise Exception, "The level number: \"" \
                      + Nz + "\" is unknown."
        else:
            raise Exception, "The ground emission parameter: \"" \
                  + ground_emission + "\" is unknown."
        Nz_in = str(int(Nz) + 1)
        # Photolysis
        photolytic_constant = param_dict["photolytic_constant"]

        # General dictionary.
        result = {"%levels%": levels,
                  "%vertical_resolution%": Nz,
                  "%vertical_distribution%": vertical_distribution,
                  "%luc%": luc,
                  "%sea_index%": sea_index,
                  "%urban_index%": urban_index,
                  "%attenuation_type%": attenuation_type,
                  "%critical_relative_humidity%": \
                  critical_relative_humidity,
                  "%min_height_cloud%": min_height_cloud,
                  "%min_kz%": min_kz,
                  "%min_urban_kz%": min_urban_kz,
                  "%apply_vert%": apply_vert,
                  "%sbl_ratio%": sbl_ratio,
                  "%boundary_layer%": boundary_layer,
                  "%kz_file%": kz_file,
                  "%bool_tm_stable%": bool_tm_stable,
                  "%exponent_p_tm%": exponent_p_tm,
                  "%dep%": dep,
                  "%Ra%": Ra,
                  "%Rb%": Rb,
                  "%photolytic_constant%": photolytic_constant,
                  "%bool_volume_emission%": bool_volume_emission,
                  "%Ncells%": str(Ncells),
                  "%Nz_in%": Nz_in}

        return result


    def GetDefaultDict(self):
        """Returns the default dictionary.

        Returns the values of each parameter by default in the case where the
        parameters do not appear in the parameter configuration file.
        @return A dictionary.
        """
        result = {"luc": "usgs",
                  "chemistry": "racm",
                  "cloud_attenuation": "radm",
                  "min_height_cloud": "500.",
                  "critical_relative_humidity": "sigma",
                  "kz": "tm",
                  "min_kz": "0.2",
                  "min_urban_kz": "0.2",
                  "apply_vert": "yes",
                  "sbl_ratio": "0.1",
                  "tm_stable": "no",
                  "exponent_p_tm": "2",
                  "boundary_layer": "1.",
                  "deposition_velocity": "zhang",
                  "Ra": "heat",
                  "Rb": "friction",
                  "ground_emission": "ground",
                  "photolytic_constant": "jproc",
                  "vertical_resolution": "5",
                  "time_step": "600.",
                  "first_layer_height": "50",
                  "with_air_density": "yes",
                  "with_source_splitting": "yes",
                  "splitting_method": "first_order",
                  "with_forced_concentration": "no"}
        return result


    def GetBinaryFile(self, model_index, ensemble_program):
        """Returns the dictionary with the path of the preprocessing binary files.

        \param model_index the index of a model.
        \param ensemble_program an instance 'EnsembleProgram'.
        @return A dictionary.
        """
        param_dict = ensemble_program.GetParameterDict(model_index)
        luc = param_dict["luc"]
        result = {"luc": "LUC-" + luc + ".bin",
                  "roughness": "Roughness-" + luc + ".bin",
                  "luc-convert": "LUC-" + luc + "-zhang.bin",
                  "extract-glcf": "LUC.bin",
                  "meteo": "ZonalWind.bin",
                  "attenuation": "Rain.bin",
                  "Kz": "Kz_Louis.bin",
                  "Kz_TM": "Kz_TM.bin",
                  "emissions": "NO.bin",
                  "bio": "NO.bin",
                  "dep": "NO.bin",
                  "ic": "NO.bin",
                  "bc-dates": "NO_z.bin"}
        return result


    def GetPerturbedFieldDict(self, model_index, ensemble_program):
        """Returns the name of fields and the associated perturbation values.

        \param model_index the index of a model.
        \param ensemble_program an instance 'EnsembleProgram'.
        @return A dictionary.
        """

        field_list = ["DepositionVelocity",  "BoundaryCondition_x",
                      "BoundaryCondition_y",  "BoundaryCondition_z",
                      "SurfaceEmission", "VolumeEmission", "PhotolysisRate"]
        param_dict = ensemble_program.GetParameterDict(model_index)
        # List of species.
        dep_species = ["O3", "NO", "NO2", "H2O2", "HCHO", "ALD", "PAN",
                       "HONO", "SO2", "HNO3", "OP1", "PAA", "ORA1", "CO"]
        NOx = ["NO", "NO2"]
        VOC_species_bc = ["API", "ETE", "ETH", "HC3", "HC8", "HCHO",
                          "ISO", "KET", "MO2", "OLT", "ONIT", "OP1", "PAN"]
        VOC_species_emission = ["ALD", "ETE", "ETH", "HC3", "HC5", "HC8",
                                "HCHO", "KET", "OLI", "OLT", "ORA2", "TOL",
                                "XYL"]
        bio_species_emission = ["API", "ISO", "LIM"]
        photolysis_species = ["NO2", "O3O1D", "O3O3P", "HONO", "HNO3",
                              "HNO4", "NO3NO", "NO3NO2", "H2O2",
                              "HCHOmol", "HCHOrad", "ALD", "MHP",
                              "HOP", "PAA", "GLYform", "GLYmol",
                              "MGLY", "UDC", "ORGNIT", "MACR", "HKET"]
        # If chemical mechanism is RADM.
        if param_dict["chemistry"] == "radm":
            VOC_species_bc.remove("ETE")
            VOC_species_bc.remove("API")
            VOC_species_bc.append("ALD")
            VOC_species_bc.append("OL2")
            VOC_species_emission.remove("ETE")
            VOC_species_emission.append("OL2")
            bio_species_emission.remove("API")
            bio_species_emission.remove("LIM")
            photolysis_species.remove("MACR")
            photolysis_species.remove("HKET")
        # Checks if perturbed field exists.
        if "bc_O3" not in param_dict.keys() and \
           "bc_NOx" not in param_dict.keys() and \
           "bc_VOC" not in param_dict.keys():
            field_list.remove("BoundaryCondition_x")
            field_list.remove("BoundaryCondition_y")
            field_list.remove("BoundaryCondition_z")
        if "dep" not in param_dict.keys():
            field_list.remove("DepositionVelocity")
        if "photolysis" not in param_dict.keys():
            field_list.remove("PhotolysisRate")
        if ("emission_NOx" not in param_dict.keys() and \
            "emission_VOC" not in param_dict.keys() and \
            "bio" in param_dict.keys()) or \
            param_dict["ground_emission"] == "ground":
            field_list.remove("VolumeEmission")
        if "bio" not in param_dict.keys() and \
           "emission_NOx" not in param_dict.keys() and \
           "emission_VOC" not in param_dict.keys():
            field_list.remove("SurfaceEmission")

        # Values for the dictionary.
        if len(field_list) == 0:
            string_field_list = "---"
        else:
            string_field_list = ""
            for field in field_list:
                string_field_list += field + " "
        temperature = function.get_perturbed_field("temperature", param_dict,
                                                   "Temperature")
        wind_angle = function.get_perturbed_field("wind_angle", param_dict,
                                                  "WindAngle")
        wind_velocity = function.get_perturbed_field("wind_velocity",
                                                     param_dict,
                                                     "WindModule")
        vertical_diff = \
            function.get_perturbed_field("vertical_diff",
                                         param_dict,
                                         "VerticalDiffusionCoefficient")
        bc_O3 = function.get_perturbed_field("bc_O3", param_dict, "O3")
        bc_NOx = function.get_perturbed_field("bc_NOx", param_dict, NOx)
        bc_VOC = function.get_perturbed_field("bc_VOC", param_dict,
                                              VOC_species_bc)
        bio = function.get_perturbed_field("bio", param_dict,
                                           bio_species_emission)
        emission_NOx = function.get_perturbed_field("emission_NOx",
                                                    param_dict, NOx)
        emission_VOC = function.get_perturbed_field("emission_NOx",
                                                    param_dict,
                                                    VOC_species_emission)
        dep = function.get_perturbed_field("dep", param_dict, dep_species)
        photolysis = function.get_perturbed_field("photolysis", param_dict,
                                                  photolysis_species)
        emission = emission_NOx + emission_VOC

        result = {"%field_list%": string_field_list,
                  "%temperature_perturbation%": temperature,
                  "%wind_velocity_perturbation%": wind_velocity,
                  "%wind_angle_perturbation%": wind_angle,
                  "%vertical_diffusion_coefficient_perturbation%": \
                  vertical_diff,
                  "%bc_NOx_perturbation%": bc_NOx,
                  "%bc_VOC_perturbation%": bc_VOC,
                  "%bc_O3_perturbation%": bc_O3,
                  "%surface_emission_perturbation%": emission,
                  "%volume_emission_perturbation%": emission,
                  "%bio_emission_perturbation%": bio,
                  "%dep_perturbation%": dep,
                  "%photolysis_rate_perturbation%": photolysis}

        return result

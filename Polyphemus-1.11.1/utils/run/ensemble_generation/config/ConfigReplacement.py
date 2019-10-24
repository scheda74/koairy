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

"""\package ensemble_generation.config.ConfigReplacement

This module contains the class 'ConfigReplacement'.

It is just the prototype.
"""


class ConfigReplacement:
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
        return {}


    def GetDefaultDict(self):
        """Returns the default dictionary.

        Returns the values of each parameter by default in the case where the
        parameters do not appear in the parameter configuration file.
        @return A dictionary.
        """
        return {}


    def GetBinaryFile(self, model_index, ensemble_program):
        """Returns the dictionary with the path of the preprocessing binary files.

        \param model_index the index of a model.
        \param ensemble_program an instance 'EnsembleProgram'.
        @return A dictionary.
        """
        return {}


    def GetPerturbedFieldDict(self, model_index, ensemble_program):
        """Gets the name of fields and the associated perturbation values.

        \param model_index the index of a model.
        \param ensemble_program an instance 'EnsembleProgram'.
        @return A dictionary.
        """
        return {}

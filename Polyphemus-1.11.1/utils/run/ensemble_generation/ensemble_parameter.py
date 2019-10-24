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

"""\package ensemble_generation.ensemble_parameter

This module contains the class 'EnsembleParameter'.

It is used to generate an ensemble of models.
"""

from os import path
from atmopy import talos


class EnsembleParameter:
    """This class builds an ensemble of models.
    Each model depends on a several parameters described in a configuration
    file. Each choice of parameter is described by an integer. A model has its
    own identity which is a list of integers. An automatic random generation
    of an ensemble can be done. Moreover, a weight can be given to each
    parameter which corresponds to the occurrence frequency at the time of the
    random generation.
    """


    ###############
    # CONSTRUCTOR #
    ###############


    def __init__(self, configuration_file = None, id_ensemble_file = None):
        """The constructor.
        Initializes the attributes and reads the configuration file if it is
        provided.
        \param configuration_file the path to the configuration file which
        describes the parameters.
        \param id_ensemble_file the path for the ensemble identity file.
        """

        ## Number of models.
        self.Nmodel = 0

        ## \brief A list of dictionaries.
        ## \details Each dictionary describes a model with the names of
        ## parameter and the associated values.
        self.dict = []

        ## \brief A list of a list of integers.
        ## \details Each list of integers describes the identity of one
        ## model.
        self.id_ensemble = []

        ## Configuration file name.
        self.config_filename = ""

        ## \brief Identity ensemble file.
        ## \details File where is stored the identity of each model. Each line
        ## describes one model.
        self.id_ens_file = ""

        if configuration_file != None:
            self.LoadConfiguration(configuration_file)
        if id_ensemble_file != None:
            self.ReadIdEnsembleFile(id_ensemble_file, copy = True)


    def LoadConfiguration(self, configuration_file):
        """Loads the configuration file.
        \param configuration_file the path to the configuration
        file which describes the parameters.
        """

        ## Number of parameters.
        self.Nparameter = 0

        ## List of parameter names.
        self.name = []

        ## List of the choices of parameters.
        self.value = []

        ## Occurrence frequency for each value of a parameter.
        self.occurrence_frequency = []

        ## A list of names of parameters which appear in preprocessing.
        self.preprocessing = []

        ## The number of values for each parameter.
        self.parameter_size = []

        # Does the configuration file exist?
        if not path.isfile(configuration_file):
            raise Exception, "Unable to open configuration file " \
                  + "\"" + configuration_file + "\"."

        self.config_filename = path.abspath(configuration_file)
        content = [("Nmodel", "[general]", "Int"), \
                   ("date_min", "[general]", "DateTime"), \
                   ("date_max", "[general]", "DateTime"), \
                   ("x_min", "[general]", "Float"), \
                   ("x_max", "[general]", "Float"), \
                   ("y_min", "[general]", "Float"), \
                   ("y_max", "[general]", "Float"), \
                   ("Delta_x", "[general]", "Float"), \
                   ("Delta_y", "[general]", "Float"), \
                   ("Nx", "[general]", "Int"), \
                   ("Ny", "[general]", "Int"), \
                   ("data_dir", "[general]", "String"), \
                   ("model_dir", "[general]", "String"), \
                   ("[physics]", "", "physics", "StringSection"), \
                   ("[numerics]", "", "numerics", "StringSection"), \
                   ("[input_data]", "", "input_data", "StringSection")]
        config = talos.Config(configuration_file, new_content = content)

        ## A atmopy.talos.Config instance.
        self.config = config

        # Do the directories exist?
        if not path.isdir(self.config.data_dir):
            raise Exception, "Unable to find the directory \"" \
                  + self.config.data_dir + "\"."
        if not path.isdir(self.config.model_dir):
            raise Exception, "Unable to find the directory \"" \
                  + self.config.model_dir + "\"."

        # Extracts the section 'physics'.
        self.ExtractSection('physics')

        # Extracts the section 'numerics'.
        self.ExtractSection('numerics')

        # Extracts the section 'input_data'.
        self.ExtractSection('input_data')

        self.Nmodel = config.Nmodel
        self.Nparameter = len(self.value)

        # Initializes the identity of the ensemble.
        self.InitIdEnsemble()


    def InitIdEnsemble(self):
        """Initializes the identity of the ensemble.
        Creates an empty list.
        """

        self.id_ensemble = []


    ##################
    # ACCESS METHODS #
    ##################


    def GetNparameter(self):
        """Returns the number of parameters.

        @return The number of parameters.
        """

        return self.Nparameter


    def GetNmodel(self):
        """Returns the number of models.

        @return The number of models.
        """

        return self.Nmodel


    def GetDict(self, model_index = None):
        """Returns the dictionary of a model.
        \param model_index the index of the model.
        @return A dictionary or a list of dictionaries.
        """

        if len(self.dict) == 0:
            self.SetDict()
        if len(self.id_ensemble) != self.Nmodel:
            raise Exception, "The length of your models identity list (" \
                + str(len(self.id_ensemble)) + ") must be equal to " \
                + str(self.Nmodel) + ". Please load or generate an ensemble."
        if model_index != None:
            if model_index >= self.Nmodel:
                raise Exception, "The model index \"" + str(model_index) + \
                      "\" is too high. The model number must be less than " \
                      + str(self.Nmodel) + "."
            # Returns the dictionary of the model.
            return self.dict[model_index]
        else:
            # Returns the list of dictionaries.
            return self.dict


    def GetStringModel(self, model_index):
        """Returns the string identity of a model.
        \param model_index index of a model.
        @return A string.
        """

        if model_index >= self.Nmodel:
            raise Exception, "The model index \"" + str(model_index) + \
                  "\" is too high. The model number must be less than " \
                  + str(self.Nmodel) + "."
        return self.IdModelToString(self.id_ensemble[model_index])


    def GetParameterFrequency(self):
        """Returns the occurrence frequency of all parameters.
        @return A dictionary.
        """

        frequency = []
        result = {}
        for i in range(self.Nparameter):
            frequency.append([])
            tmp = [x[i] for x in self.id_ensemble]
            for j in range(self.parameter_size[i]):
                count = 0
                for n in tmp:
                    if j == n:
                        count += 1
                frequency[i].append(float(count) / float(self.Nmodel))
            result[self.name[i]] = frequency[i]

        return result


    def GetIdProgram(self, id_model, name):
        """Returns the suite of integers which represent the values
        of a list of parameters.

        \param id_model a string of integers which describes the model
        identity.
        \param name a list of parameter names.
        @return A string.
        """

        result = ""
        if name == ["None"]:
            return "0"
        for n in name:
            if n not in self.name:
                raise Exception, "The parameter name \"" + n + "\" is not" \
                      +" in the parameters list of the file \"" \
                      + self.config_filename + "\"."
            result += str(id_model[self.name.index(n)])
        return result


    def GetNormalUncertainty(self, value):
        """Returns the uncertainty for a normal law.
        \param value value.
        @return A list of 3 floats.
        """

        return [0., -value / 2., +value / 2.]


    def GetLogNormalUncertainty(self, value):
        """Returns the uncertainty for a log-normal law.
        \param value value.
        @return A list of 3 floats.
        """

        from math import sqrt, exp, log
        if value <= 0:
            raise Exception, "The value: \"" + str(value) + "\" must be" \
                  + " higher than 0 for a log-normal law."
        beta = exp(0.5 * (0.5 * log(value))**2)
        delta = 4. * beta**4  - 7. * beta**2 + 6. * beta - 3.
        if 3. * beta - 1. - sqrt(delta) <= 0:
            raise Exception, "The value :\"" + str(value) \
                  + "\" is not permitted."
        a = 2. * log((3. * beta - 1 - sqrt(delta)) / 2.) / log(value)
        b = 2. * log((3. * beta - 1 + sqrt(delta)) / 2.) / log(value)

#         mean = 1./ 3. * (value**(a / 2.) + value**(b / 2.) + 1.)
#         sigma = 0.5 * ((value**(a / 2.) - mean)**2 \
#                        + (value**(b / 2.) - mean)**2 \
#                        + (1. - mean)**2)
#         Var_x = beta**2 * (beta**2 - 1.)

        return [1., value**(a / 2.), value**(b / 2.)]


    ##############
    # PROCESSING #
    ##############


    def ExtractSection(self, section):
        """Extracts the parameters names and their values from a section of
        the configuration file.

        \param section the name of the section which must be 'physics',
        'numerics' or 'input_data'.
        """

        all_section = ['physics', 'numerics', 'input_data']
        if section not in all_section:
            raise Exception, "The section \"" + section \
                + "\" is not supported."
        # A list of strings.
        parameter_list = ['']
        if hasattr(self.config, section):
            parameter_list = getattr(self.config, section)
        else:
            raise Exception, "Unable to find the section \"" \
                + section + "\" in the configuration file \"" \
                + self.config_filename + "\"."
        N = len(self.name)
        if len(parameter_list[0]) != 0:
            # Loop on the parameters.
            for i in range(len(parameter_list)):
                tmp = parameter_list[i].split()
                # Removes sepration characters.
                for s in ',:;-=':
                    while s in tmp:
                        tmp.remove(s)
                # Name of the parameter.
                self.name.append(tmp[0].strip(":=,;"))
                if section in ['physics', 'numerics']:
                    self.value.append([])
                    self.occurrence_frequency.append([])
                    self.parameter_size.append(0)
                    # Loop on the choices and the occurrence frequencies.
                    for s in tmp[1:]:
                        if s == "preprocessing":
                            self.preprocessing.append\
                                (self.name[i + N])
                            continue
                        # If there is an occurrence frequency.
                        if '(' in s:
                            number = talos.to_num(s.strip(",()-;"))
                            self.occurrence_frequency[i + N]\
                                .append(number)
                        # The name of a choice of parameter.
                        else:
                            self.value[i + N].append(s.strip(":=,;-"))
                            self.parameter_size[i + N] += 1
                    # If there is any occurrence frequency.
                    if len(self.occurrence_frequency[i + N]) == 0:
                        for _ in range(self.parameter_size[i + N]):
                            self.occurrence_frequency[i + N].append\
                                (100. / float(self.parameter_size[i + N]))
                    # If the number of occurrence frequencies is wrong.
                    elif len(self.occurrence_frequency[i + N]) \
                            != self.parameter_size[i + N]:
                        raise Exception, "The number of occurrence " \
                            + "frequencies (" \
                            + str(len(self.occurrence_frequency[i + N])) \
                            + ") for the parameter \"" + self.name[i + N] \
                            + "\" must be equal to "\
                            + str(self.parameter_size[i + N]) + "."
                # If the section is 'input_data'.
                else:
                    value = float(tmp[1].strip(":=,;"))
                    type_perturbation = tmp[-1].strip(":=,;")
                    if type_perturbation == "normal":
                        uncertainty = self.GetNormalUncertainty(value)
                    elif type_perturbation == "log-normal":
                        uncertainty = self.GetLogNormalUncertainty(value)
                    else:
                        raise Exception, "The type: \"" + type_perturbation \
                            + "\" in \"" + self.config_filename + "\"" \
                            + " is not supported."
                    self.value.append((uncertainty, type_perturbation))
                    self.parameter_size.append(len(uncertainty))
        # If there are any parameters.
        else:
            pass


    def SetDict(self):
        """Sets the dictionary of each model.

        From the attributes 'name' and 'value', it sets the dictionaries of
        each model.
        """

        self.dict = []
        for id_model in self.id_ensemble:
            tmp = {}
            for i in range(self.Nparameter):
                # For the tuple ([uncertainty], type) for input data.
                if isinstance(self.value[i], tuple):
                    if self.value[i][-1] == "log-normal":
                        type_perturbation = "multiply"
                    elif self.value[i][-1] == "normal":
                        type_perturbation = "add"
                    else:
                        raise Exception, "The type: \"" + type_perturbation \
                              + "\" in \"" + self.config_filename + "\"" \
                              + " is not supported."
                    tmp[self.name[i]] = (self.value[i][0][id_model[i]], \
                                         type_perturbation)
                # For the list of choices
                # (for the sections 'physics' and 'numerics').
                else:
                    tmp[self.name[i]] = self.value[i][id_model[i]]
            self.dict.append(tmp)


    def ParameterChoice(self, index_parameter):
        """Chooses randomly the parameter value.
        Depending on the occurrence frequency, it returns an integer that
        represents a value of a parameter.

        \param index_parameter index of the parameter.
        @return An integer.
        """

        import random
        # Number of input data.
        if len(self.config.input_data[0]) == 0:
            input_data_length = 0
        else:
            input_data_length = len(self.config.input_data)
        # If 'index_parameter' is wrong.
        # The input data do not have occurrence frequencies.
        if index_parameter >= self.Nparameter - input_data_length:
            raise Exception, "The index of the parameter must " \
                + "equal or less than " \
                + str(self.Nparameter - input_data_length) + "."
        # The sum of occurrence frequencies must be equal to 100.
        if sum(self.occurrence_frequency[index_parameter]) != 100:
            raise Exception, "The sum of the occurrence frequencies \"" \
                  + str(self.occurrence_frequency[index_parameter]) \
                  + "\" must be equal to 100 for the parameter \"" \
                  + self.name[index_parameter] + "\"."
        length = self.parameter_size[index_parameter]
        random_number = random.random() * 100.
        # The cumulative sum of occurrence frequencies.
        occurrence_cumsum = \
            [sum(self.occurrence_frequency[index_parameter][0 : x + 1]) \
                 for x in range(length)]
        occurrence_cumsum.append(0.)
        occurrence_cumsum.sort()
        # Selection of the parameter choice in function to the random number.
        for i in range(length - 1):
            if random_number >= occurrence_cumsum[i] \
                    and random_number < occurrence_cumsum[i + 1]:
                return i

        return length - 1


    def GenerateIdModel(self):
        """Generates randomly the identity of a model.

        @return A list of integers.
        """

        import random
        id_model = []
        count = 0
        if len(self.config.input_data[0]) == 0:
            input_data_length = 0
        else:
            input_data_length = len(self.config.input_data)
        # For the parameters of the sections 'physics', 'numerics'.
        for i in range(self.Nparameter - input_data_length):
            parameter = self.ParameterChoice(i)
            id_model.append(parameter)
            count += 1
        # For the pertubed input data.
        for i in range(count, self.Nparameter, 1):
            random_number = random.randint(0, self.parameter_size[i] - 1)
            id_model.append(random_number)

        # If there is Kz Louis parameterization, these parameters (for Kz_TM):
        # sbl_ratio, tm_stable and exponent_tm must stay by default.
        if "kz" in self.name:
            index_kz = self.name.index("kz")
            if self.value[index_kz][id_model[index_kz]] == "louis":
                if "sbl_ratio" in self.name:
                    index = self.name.index("sbl_ratio")
                    id_model[index] = 0
                if "tm_stable" in self.name:
                    index = self.name.index("tm_stable")
                    id_model[index] = 0
                if "exponent_p_tm" in self.name:
                    index = self.name.index("exponent_p_tm")
                    id_model[index] = 0
        # If there is RADM parameterization for the chemical mechanism,
        # there will not be zenith_angle photolytic rates.
        if "chemistry" in self.name:
            index_chemistry = self.name.index("chemistry")
            if self.value[index_chemistry][id_model[index_chemistry]] \
                    == "radm":
                if "photolytic_constant" in self.name:
                    index = self.name.index("photolytic_constant")
                    id_model[index] = 0;

        return id_model


    def GenerateIdEnsemble(self):
        """Generates randomly an ensemble.
        """
        self.InitIdEnsemble()
        Nmodel_possible = 1
        # Checks an approximation of number of possible models.
        for i in self.parameter_size:
            Nmodel_possible *= i
        if self.Nmodel > Nmodel_possible:
            raise Exception, "Unable to generate an ensemble of " \
                  + str(self.Nmodel) + " members. There are only " \
                  + str(Nmodel_possible) + " possibilities."

        # Generates an ensemble of models identities.
        for i in range(self.Nmodel):
            self.id_ensemble.append([])
            self.id_ensemble[i] = self.GenerateIdModel()
            while not self.CheckIdEnsemble():
                self.id_ensemble[i] = self.GenerateIdModel()
        self.CheckParameterFrequency()


    def CheckIdEnsemble(self):
        """Checks the identity of the latest added model to the ensemble.

        Checks if the latest model added to the ensemble has not been in the
        ensemble yet.

        @return 0 or 1.
        """

        if len(self.id_ensemble) == 0 or len(self.id_ensemble) == 1:
            return 1
        for i in range(len(self.id_ensemble) - 1):
            if self.id_ensemble[i] == self.id_ensemble[-1]:
                return 0
        return 1


    def CheckParameterFrequency(self):
        """Checks the occurrence frequency of each choice of parameters.

        Prints a message about the parameters values which do not appear in
        the ensemble
        """

        for i in range(self.Nparameter):
            count = 0
            index = []
            tmp = [x[i] for x in self.id_ensemble]
            for j in range(self.parameter_size[i]):
                if j in tmp:
                    count += 1
                else:
                    index.append(j)
            if count != self.parameter_size[i]:
                print ("Checks parameter frequency : the index values \"" \
                       + str(index) + "\" is not present in the id_ensemble" \
                       + " for the parameter \"" + self.name[i] + "\".")


    def AddModel(self, id_model):
        """Adds a model.
        """
        id_model = self.IdModelToList(id_model)
        if id_model in self.id_ensemble:
            raise Exception, "The model what you want to add already exists."
        if (len(self.id_ensemble)) >= self.Nmodel:
            raise Exception, "Unable to add a model, the Nmodel is " \
                  + "already equal to " + str(self.Nmodel) \
                  + ". Just modify the Nmodel."
        self.id_ensemble.append(id_model)
        self.SetDict()


    def IdModelToString(self, id_model):
        """Converts a list of integers to a string.

        \param id_model a list of integers which describes the model identity.
        @return A string.
        """

        if len(id_model) != self.Nparameter:
            raise Exception, "The number of parameters is wrong (" \
                      + str(len(id_model)) + "). It must be equal to: " \
                      + str(self.Nparameter) + "."
        result = ""
        if isinstance(id_model, str):
            return id_model
        for i in id_model:
            result += str(i)
        return result


    def IdModelToList(self, id_model):
        """Converts a string of integers to a list of integers.

        \param id_model a string of integers which describes the model
        identity.
        @return A list of integers.
        """

        if len(id_model) != self.Nparameter:
            raise Exception, "The number of parameters is wrong (" \
                  + str(len(id_model)) + "). It must be equal to: " \
                  + str(self.Nparameter) + "."
        if isinstance(id_model, list):
            return id_model
        result = [int(x) for x in id_model]
        for i in range(self.Nparameter):
            if result[i] >= self.parameter_size[i]:
                raise Exception, "The number of values in the parameter \"" \
                      + self.name[i] + "\" must be less than " \
                      + str(self.parameter_size[i]) + " (instead of " \
                      + str(result[i]) + " )."
        return result


    def WriteIdEnsembleInFile(self, filename):
        """Writes the identity of each model in a file.

        \param filename the name of the file.
        """

        f = open(filename, 'w')
        for m in range(self.Nmodel):
            id_model = self.IdModelToString(self.id_ensemble[m])
            f.write(id_model + "\n")
        f.close()


    def ReadIdEnsembleFile(self, filename, copy = False):
        """Reads the identity of each model from a file.

        \param filename the name of the file.
        \param copy will the identity of the ensemble be copied?
        @return A list of model identities.
        """

        result = []
        if not path.isfile(filename):
            raise Exception, "Unable to open file \"" + filename \
                + "\" and read id_ensemble."
        self.id_ens_file = filename
        f = open(filename, 'r')
        tmp = f.readline()
        result.append(self.IdModelToList(tmp.strip()))
        count = 1
        while len(tmp) != 0:
            tmp = f.readline()
            if len(tmp) == 0:
                continue
            result.append(self.IdModelToList(tmp.strip()))
            count += 1
        f.close()
        print "Reads \"" + filename + "\"..."
        if not copy:
            return result
        else:
            self.id_ensemble = result
            self.Nmodel = count


    def ShowModelParameter(self, id_model):
        """Prints the choices of each parameter for a given model.

        \param id_model a string of integers.
        """
        model = self.IdModelToList(id_model)
        for i in range(self.Nparameter):
            if isinstance(self.value[i], tuple):
                print(self.name[i] + " : " \
                      + str(self.value[i][0][model[i]])) \
                      + " - " + self.value[i][-1]
            else:
                print(self.name[i] + " : " + str(self.value[i][model[i]]))


# -*- coding: utf-8 -*-
# Copyright (C) 2008-2014 EDF R&D
# Author: Damien Garaud, Sylvain Doré
#
# This file is part of the air quality modeling system Polyphemus. It is used
# to check ensemble generation.
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
from run.ensemble_generation import EnsembleParameter

import os, sys, unittest


class EnsembleParameterTest(unittest.TestCase):

    def setUp(self):
        # The loop number for the random number generation.
        self.Nloop = 1000
        self.cwd = os.getcwd()
        os.chdir(os.path.join(os.path.dirname(__file__), "../../example"))
        self.config_parameter = "parameter.cfg"

    def tearDown(self):
        os.chdir(self.cwd)

    def testLoadEnsembleParameter(self):
        # Declarations.
        ep = EnsembleParameter()

        # Raises an Exception for a wrong file.
        self.assertRaises(Exception, ep.LoadConfiguration, 'nofile_sorry')
        self.assertRaises(Exception, ep.ReadIdEnsembleFile, 'nofile_sorry',
                          True)
        # An empty ensemble.
        self.assert_(len(ep.id_ensemble) == 0)

        # Equal lengths.
        self.assert_(ep.Nmodel == 0)
        self.assert_(len(ep.value) == len(ep.name))
        self.assert_(len(ep.value) == ep.Nparameter)
        self.assert_(len(ep.name) == ep.Nparameter)
        self.assert_(len(ep.parameter_size) == ep.Nparameter)
        for i in range(ep.Nparameter):
            self.assert_(len(ep.value[i]) == ep.parameter_size[i])
            self.assert_(len(ep.occurrence_frequency[i]) == ep.parameter_size[i])

        # Loading.
        ep.LoadConfiguration(self.config_parameter)
        if len(ep.config.input_data[0]) != 0:
            Ntmp = ep.Nparameter - len(ep.config.input_data)
        else:
            Ntmp = ep.Nparameter

        # Equal lengths.
        self.assert_(len(ep.value) == len(ep.name))
        self.assert_(len(ep.value) == ep.Nparameter)
        self.assert_(len(ep.name) == ep.Nparameter)
        self.assert_(len(ep.parameter_size) == ep.Nparameter)
        self.assert_(len(ep.occurrence_frequency) == Ntmp)
        for i in range(Ntmp):
            self.assert_(len(ep.value[i]) == ep.parameter_size[i])
            self.assert_(len(ep.occurrence_frequency[i]) \
                             == ep.parameter_size[i])

    def testGenerateEnsemble(self):
        import random

        # Declarations.
        ep = EnsembleParameter(self.config_parameter)
        if len(ep.config.input_data[0]) != 0:
            Ntmp = ep.Nparameter - len(ep.config.input_data)
        else:
            Ntmp = ep.Nparameter

        # Equal lengths.
        self.assert_(len(ep.value) == len(ep.name))
        self.assert_(len(ep.value) == ep.Nparameter)
        self.assert_(len(ep.name) == ep.Nparameter)
        self.assert_(len(ep.parameter_size) == ep.Nparameter)
        self.assert_(len(ep.occurrence_frequency) == Ntmp)
        for i in range(Ntmp):
            self.assert_(len(ep.value[i]) == ep.parameter_size[i])
            self.assert_(len(ep.occurrence_frequency[i]) \
                             == ep.parameter_size[i])

        # Generates an ensemble.
        ep.GenerateIdEnsemble()
        self.assert_(len(ep.id_ensemble) == ep.Nmodel)
        for i in range(ep.Nmodel):
            self.assert_(len(ep.id_ensemble[i]) == ep.Nparameter)

        # Checks the ensemble.
        self.assert_(ep.CheckIdEnsemble())

        # For the random number generation.
        for i in range(self.Nloop):
            index_parameter = random.randint(0, Ntmp - 1)
            choice_list = range(ep.parameter_size[index_parameter])
            index_choice = ep.ParameterChoice(index_parameter)
            self.assert_(index_choice in choice_list)


if __name__ == '__main__':
    unittest.main(argv=sys.argv)

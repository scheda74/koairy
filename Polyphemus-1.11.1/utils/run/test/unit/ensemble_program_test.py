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
from run.ensemble_generation import config, EnsembleProgram

import os, sys, unittest


class EnsembleProgramTest(unittest.TestCase):

    def setUp(self):
        test_dir = os.path.dirname(__file__)
        self.cwd = os.getcwd()
        os.chdir(os.path.join(test_dir, "../../example"))
        self.config_parameter = "parameter.cfg"
        self.config_program = "program.cfg"
        self.polyphemus_dir = \
          os.path.abspath(os.path.join(test_dir, "../../../.."))

    def tearDown(self):
        os.chdir(self.cwd)

    def testLoadEnsembleProgram(self):
        # Declaration.
        p = EnsembleProgram()

        # Raises an Exception for a wrong file.
        self.assertRaises(Exception, p.Init, 'nofile_sorry', 'nofile_sorry')
        self.assertRaises(Exception, p.ReadIdEnsembleFile, 'nofile_sorry')

        # Loading.
        p.Init(self.config_parameter, self.config_program)

        # Raises an Exception because the ensemble is empty.
        self.assertRaises(Exception, p.GetGeneralDict, 0)
        p.parameter.GenerateIdEnsemble()

        # Raises an Exception because the 'config_replacement' is not set.
        self.assertRaises(Exception, p.GetGeneralDict, 0)

        # Sets the Polyphemus directory.
        p.SetPolyphemusDirectory(self.polyphemus_dir)

        # Declaration of the ConfigReplacement object.
        config_replacement = config.ConfigPolair3D()
        p.SetConfigReplacement(config_replacement)

        # Gets the general dictionary for the first model.
        _ = p.GetGeneralDict(0)


if __name__ == '__main__':
    unittest.main(argv=sys.argv)

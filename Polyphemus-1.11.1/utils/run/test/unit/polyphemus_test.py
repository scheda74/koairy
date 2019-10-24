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
from run.polyphemus import Polyphemus, Program

import io, unittest, sys


class PolyphemusTest(unittest.TestCase):

    def testProgram(self):
        msg = "Hello World!"
        program = Program("/bin/echo", arguments_format = msg)

        log = io.BytesIO()
        program.Run(log)
        self.assertTrue(log.getvalue().startswith(msg))

    def testPolyphemus(self):
        simulation = Polyphemus()

        msg_list = [ "first message", "second message" ]

        for msg in msg_list:
            simulation.AddProgram(Program("/bin/echo",
                                          arguments_format = msg))
        log = io.BytesIO()
        simulation.Run(log)
        log_str = log.getvalue()
        for msg in msg_list:
            self.assertTrue(msg in log_str)


if __name__ == '__main__':
    unittest.main(argv=sys.argv)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2016, ENPC
#     Author(s): Sylvain Doré
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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
# Polyphemus is developed in the CEREA, a joint laboratory between ENPC and
# EDF R&D.

import multiprocessing as mp
import sys
from functools import partial


class ParallelTaskRunner:
    """Runs in parallel a list of function calls.

    Function calls are added with 'append' and ran as soon as possible in a
    subprocess. The results are returned with 'retrieve'.

    Rational: For running tests, parallel execution is needed, and, generally,
    results should be kept in the order the tasks were submitted (so that test
    results can be compared). The reason of this class is the lack of such a
    tool in the standard Python libraries.
    """
    def __init__(self, process_count=mp.cpu_count()):
        self.process_count = process_count
        self.pool = None
        self.task_list = []
        self.callback_list = []
        self.result_list = []

    def _create_pool(self):
        self._destruct_pool()
        self.pool = mp.Pool(self.process_count)

    def _destruct_pool(self):
        if self.pool:
            self.pool.terminate()
            self.pool = None

    def append(self, func, args=(), callback=lambda x : None):
        """Runs 'func' with 'args' in a subprocess.

        A call to 'append' is non-blocking (it returns immediately).

        @type func: function
        @param func: The function to run in parallel.

        @type args: tupple
        @param func: The arguments to be given to 'func'.

        @type callback: function
        @param callback: An optional function that will be ran on the result
        when 'retrieve' is called.
        """
        if self.process_count > 1:
            if not self.pool:
                self._create_pool()
            task = self.pool.apply_async(partial(func, *args))
            self.task_list.append(task)
            self.callback_list.append(callback)
        else:
            self.result_list.append(callback(func(*args)))

    def retrieve(self):
        """Runs the optional callbacks and returns the list of result.

        Each callback is ran on its result, in the order the tasks were added,
        waiting for the corresponding process to end as necessary.

        A call to 'retrieve' is blocking, the list of tasks and their
        results is cleared.
        """
        try:
            for task, callback in zip(self.task_list, self.callback_list):
                result = task.get()
                callback(result)
                self.result_list.append(result)
            return self.result_list
        finally:
            self._destruct_pool()
            self.task_list = []
            self.callback_list = []
            self.result_list = []

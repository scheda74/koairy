#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2015, ENPC
#     Author(s): Sylvain Doré
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

import numpy as np
import numpy.ma as ma

from pprint import pprint

from test_utils.pytest_ext import recursive_equal
from test_utils.compare_array import compare_array


def test_compare_array():
    def test(lhs, rhs, expected):
        result = compare_array(lhs, rhs)
        result
        if not recursive_equal(result, expected, comp=np.equal):
            print "Here is the actual comparison result:"
            pprint(result)
            return False
        return True

    # Tests 'compare_array'.
    a = [[1., 2.], [3., 4.]]
    m=ma.masked_array
    F = False
    T = True

    assert test(a, a, {})

    assert test(a, [[1., 2.], [3., np.NaN]],
                {'NaN': {'cur': 'no NaN', 'ref': 'has 1 NaN'},
                 'diff_array' : m([[0., 0.], [0., -1]], [[F, F], [F, T]])})

    assert test(a, [[np.PINF, 2.], [np.NINF, 4.]],
                {'negative infinity': {'cur': 'no negative infinity',
                                       'ref': 'has 1 negative infinity'},
                 'positive infinity': {'cur': 'no positive infinity',
                                       'ref': 'has 1 positive infinity'},
                 'diff_array' : m([[-1, 0.], [-1, 0.]], [[T, F], [T, F]])})

    assert test([[1., 2.], [3., -8.]], [[np.NaN, 2*2.], [3., 5.]],
                {'NaN': {'cur': 'no NaN', 'ref': 'has 1 NaN'},
                 'different sign': {'cur': 'contains some strictly negative values',
                                    'ref': 'all values are positive or null'},
                 'abs max (rtol=5.00e-02)': {'cur': 8.0, 'diff': 3.,'ref': 5.0},
                 'mean (rtol=1.00e-02)': {'cur': -1.0, 'diff': -5., 'ref': 4.0},
                 'diff_array' : m([[-1, -2.], [0., -13.]], [[T, F], [F, F]])})

    assert test([[np.PINF, 2.], [3., 4.]], [[1., 2.], [3., 4.]],
                {'positive infinity': {'cur': 'has 1 positive infinity',
                                       'ref': 'no positive infinity'},
                 'diff_array' : m([[-1, 0.], [0., 0.]], [[T, F], [F, F]])})

    assert test([[1., np.NaN], [3., 4.]], [[np.NaN, 2.], [3., 4.]],
                {'NaN': {'cur': "1 NaN not in 'ref' (and 0 also in 'ref')",
                         'diff': '2 NaN differences, representing a rate of 50.00% on the total number of values (tolerance=1.00%)',
                         'ref': "1 NaN not in 'cur' (and 0 also in 'cur')"},
                 'diff_array' : m([[-1, -1], [0., 0.]], [[T, T], [F, F]])})

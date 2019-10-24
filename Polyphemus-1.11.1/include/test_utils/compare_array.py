#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2015, ENPC
#     Author(s): Sylvain Doré
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the INRIA project-team CLIME and in
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


def compare_array(array, ref_array,
                  mean_rtol=1e-02, max_rtol=5e-02,
                  nan_rtol=1e-02, inf_rtol=1e-02):
    """Returns true if 'array' is an approximation of 'ref_array'.

    A difference is detected in the following cases:
    - Array shape (dimensions and size) are different
    - The quotient of the means exceeds the relative tolerance 'mean_rtol'.
    - The quotient of the absolute maximums exceeds the relative tolerance
        'max_rtol'.
    - Negative value in 'array' when 'ref_array' is all positive or null. Idem
        for positive value.
    - NaN in 'array' when 'ref_array' has none. Idem for positive infinity and
        negative infinity.
    - NaN in 'ref_array' when 'array' has none. Idem for positive infinity
        and negative infinity.
    - When both 'array' and 'ref_array' have NaN: The number of NaN only in
        'array' or 'ref_array' divided by the total number of 'ref_array' values
        exceeds the relative tolerance 'nan_rtol'. Idem for positive and
        negative infinity. The number of differences in positive and negative
        infinity are summed together before comparing with 'inf_rtol'.

    @type array: array like
    @param array: The array that is to be compared.

    @type ref_array: array like
    @param ref_array: The reference of the comparison.

    @type mean_rtol: float
    @param mean_rtol: The relative tolerance for the mean of the array.

    @type max_rtol: float
    @param max_rtol: The relative tolerance for the absolute maximum of the
        array.

    @type nan_rtol: float
    @param nan_rtol: The tolerance for NaN proportion change, relative to the
        total number of value.

    @type inf_rtol: float
    @param inf_rtol: The tolerance for infinity proportion change, relative to
        the total number of value.

    @rtype: boolean
    @return: A self-describing dictionary reporting the comparison results. If
    a field is called 'diff_array', it is used to plot a map of difference.
    """
    if np.shape(array) != np.shape(ref_array):
        return{'different dimensions' : { 'cur' : str(np.shape(array)),
                                          'ref' : str(np.shape(ref_array)) } }
    if np.size(ref_array) == 0:
        return {}

    difference = {}

    def hide_special_values(raw_array):
        """Hides NaN and infinity, and returns a report on those values."""
        fixed_array = ma.masked_invalid(raw_array, copy=False)
        hidding = ma.is_masked(fixed_array)
        return (fixed_array if hidding else raw_array, {
        'positive infinity': np.equal(raw_array, np.PINF) if hidding else [False],
        'negative infinity': np.equal(raw_array, np.NINF) if hidding else [False],
        'NaN': np.isnan(raw_array) if hidding else [False]
        } )

    array, special_values = hide_special_values(array)
    ref_array, ref_special_values = hide_special_values(ref_array)

    for value_type, value_mask in special_values.iteritems():
        ref_value_mask = ref_special_values[value_type]
        value_count = np.sum(value_mask)
        ref_value_count = np.sum(ref_value_mask)

        if ref_value_count == 0 and value_count == 0:
            continue

        diff = np.logical_xor(value_mask, ref_value_mask)
        only_in_cur = np.sum(np.logical_and(value_mask, diff))
        only_in_ref = np.sum(np.logical_and(ref_value_mask, diff))

        if ref_value_count == 0:
            difference[value_type] = {
                'cur': "has {} {}".format(value_count, value_type),
                'ref': "no {}".format(value_type) }
        elif value_count == 0:
            difference[value_type] = {
                'cur': "no {}".format(value_type),
                'ref': "has {} {}".format(ref_value_count, value_type) }
        else:
            rtol = inf_rtol
            if value_type == "Nan":
                rtol = nan_rtol
            ref_size = np.size(ref_array)
            tol = rtol * ref_size
            common_count = value_count - only_in_cur
            diff_count = only_in_cur + only_in_ref
            if diff_count > tol:
                difference[value_type] = {
                    'cur': "{} {} not in 'ref' (and {} also in 'ref')"
                            .format(only_in_cur, value_type, common_count),
                    'diff': "{} {} differences, representing a \
rate of {:.2f}% on the total number of values (tolerance={:.2f}%)"
                            .format(diff_count, value_type,
                                    (100. * diff_count / ref_size), 100. * rtol),
                    'ref': "{} {} not in 'cur' (and {} also in 'cur')"
                            .format(only_in_ref, value_type, common_count) }

    # Puts a mask so that the numerical comparison does not take into account\
    # invalid values:
    if special_values or ref_special_values:
        mask_all = np.logical_or(ma.getmaskarray(array),
                                 ma.getmaskarray(ref_array))
        if special_values:
            ref_array = ma.masked_array(ref_array, mask_all)
        if ref_special_values:
            array = ma.masked_array(array, mask_all)

    array_min = np.min(array)
    ref_array_min = np.min(ref_array)
    if np.sign(ref_array_min) >= 0 and np.sign(array_min) == -1:
        difference['different sign'] = {
                            'cur' : "contains some strictly negative values",
                            'ref' : "all values are positive or null"
                            }

    array_max = np.max(array)
    ref_array_max = np.max(ref_array)
    if np.sign(ref_array_max) <= 0 and np.sign(array_max) == 1:
        difference['different sign'] = {
                            'cur' : "contains some strictly positive values",
                            'ref' : "all values are negative or null"
                            }

    def is_error(value, ref_value, rtol):
        return np.abs(value - ref_value) > rtol * np.abs(ref_value)

    array_max_abs = max(abs(array_max), abs(array_min))
    ref_array_max_abs = max(abs(ref_array_max), abs(ref_array_min))

    if is_error(array_max_abs, ref_array_max_abs, max_rtol):
        difference['abs max (rtol={:.2e})'.format(max_rtol)] = {
                                        'cur' : array_max_abs,
                                        'ref' : ref_array_max_abs,
                                        'diff' : array_max_abs - ref_array_max_abs,
                                        }

    array_mean = np.mean(array)
    ref_array_mean = np.mean(ref_array)
    if is_error(np.mean(array), np.mean(ref_array), mean_rtol):
        difference['mean (rtol={:.2e})'.format(mean_rtol)] = {
                                        'cur' : array_mean,
                                        'ref' : ref_array_mean,
                                        'diff' : array_mean - ref_array_mean,
                                        }

    # Puts an array of the differences, not accounting for special values.
    if difference:
        difference['diff_array'] = array - ref_array
    return difference

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

import filecmp
import fnmatch
import os

from array import array
from functools import partial
from numpy import fromfile
from pprint import pformat

from atmopy import talos
from atmopy.display import getd, getm

from config_helpers import as_config, get_result_dir_list
from file_helpers import (check_path,
                          normalize_path,
                          remove_duplicate_path,
                          file_list)
from test_utils.compare_array import compare_array
from test_utils.display.matplotlib_ext import *
from test_utils.ParallelTaskRunner import ParallelTaskRunner


def compare_result(config,
                   base_dir="$POLYPHEMUS",
                   ref_base_dir="$REF_POLYPHEMUS",
                   task_runner=ParallelTaskRunner(),
                   **kwargs):
    """Compares the data files between two result directories.

    The comparison itself is done by the function 'compare_path' called
    with the result directory read from 'config'.

    The full result directory must begins with 'base_dir', the reference
    result directory is then deduced by replacing 'base_dir' with
    'ref_base_dir'.

    @type config: Talos.Config or string or a list of them
    @param config: The configuration object or the configuration filename that
    gives the simulation result directory. It can be relative to the current
    directory. A list of config can also be given.

    @type base_dir: string
    @param base_dir: A path that:
    - is the beginning of the absolute result directory.
    - is replaced with 'ref_base_dir' to get the reference result directory.

    @type ref_base_dir: string
    @param ref_base_dir: The base path of the reference data.

    @type kwargs: named arguments
    @param kwargs: Arguments to be transmitted to the function
    'compare_bin_file'.

    @rtype: dictionary
    @return: A self-describing dictionary reporting the comparison results.
    It is actually the one returned by 'compare_path', but with 'data_dir'
    as top-level key. If there is no difference, an empty dictionary is
    returned.
    """
    if type(config) not in (str, unicode, talos.Config):
        difference = {}
        for c in config:
            difference.update(compare_result(c, base_dir, ref_base_dir,
                                             task_runner, **kwargs))
        return difference

    config = as_config(config)
    result_dir_list = get_result_dir_list(config)
    return compare_dir(result_dir_list, config, base_dir, ref_base_dir,
                       task_runner, **kwargs)


def compare_dir(data_dir=None,
                config=None,
                base_dir="$POLYPHEMUS",
                ref_base_dir="$REF_POLYPHEMUS",
                task_runner=ParallelTaskRunner(),
                **kwargs):
    """Compares the data files between two directories.

    The full 'data_dir' must begins with 'base_dir', the reference data
    directory is then deduced by replacing 'base_dir' with 'ref_base_dir'.

    The comparison itself is done by  the function 'compare_path'.

    @type data_dir: string or a list of strings
    @param data_dir: The data directory to be compared, that can be relative
    to the current directory. A list of data directory can also be given.

    @type base_dir: string
    @param base_dir: A path that:
    - is the beginning of the full 'data_dir'.
    - is replaced with 'ref_base_dir' to get the reference data directory.

    @type ref_base_dir: string
    @param ref_base_dir: The base path of the reference data.

    @type task_runner: ParallelTaskRunner
    @param task_runner: An instance of ParallelTaskRunner to run the
    file comparisons in parallel.

    @type kwargs: named arguments
    @param kwargs: Arguments to be transmitted to the function
    'compare_bin_file'.

    @rtype: dictionary
    @return: A self-describing dictionary reporting the comparison results.
    It is actually the one returned by 'compare_path', but with 'data_dir'
    as top-level key. If there is no difference, an empty dictionary is
    returned.
    """
    if type(data_dir) not in (str, unicode):
        difference = {}
        for d in remove_duplicate_path(data_dir):
            difference.update(compare_dir(d, config, base_dir, ref_base_dir,
                                          task_runner, **kwargs))
        return difference

    base_dir = normalize_path(base_dir)
    ref_base_dir = normalize_path(ref_base_dir)
    check_path(base_dir)
    check_path(ref_base_dir)

    abs_data_dir = normalize_path(data_dir)

    if not abs_data_dir.startswith(base_dir):
        raise ValueError, "The compared data directory 'data_dir' (\"%s\") "\
                    "must be under the base directory 'base_dir' (\"%s\")." %\
                    (abs_data_dir, base_dir)

    ref_data_dir = abs_data_dir.replace(base_dir, ref_base_dir)

    data_key = data_dir
    if not data_key:
        data_key = abs_data_dir

    if not os.path.exists(data_dir):
        if not os.path.exists(ref_data_dir):
            return {}
        else:
            return {data_key: 'lacking directory'}
    elif not os.path.exists(ref_data_dir):
        return {data_key: 'unexpected directory'}

    difference = compare_path(data_dir, ref_data_dir, task_runner,
                              config=config, **kwargs)
    task_runner.retrieve()

    if difference:
        return {data_key: difference}

    # No difference.
    return {}


def compare_path(data_dir, ref_data_dir, task_runner, show_figure=True,
                 **kwargs):
    """Compares the data files between the two input directories.

    A self-describing dictionary reporting the differences is generated.
    When a difference is detected between two ".bin" files, a plot
    representing this difference is generated along with some statistics.

    @type data_dir
    @param data_dir: The path whose content is to be compared.

    @type ref_data_dir: string
    @param ref_data_dir: The path whose content is the comparison reference.

    @type show_figure: boolean
    @param show_figure: True if the generated image files are also displayed
    to the user.

    @type task_runner: ParallelTaskRunner
    @param task_runner: An instance of ParallelTaskRunner to run the
    file comparisons in parallel.

    @type kwargs: named arguments
    @param kwargs: Arguments to be transmitted to the function
    'compare_bin_file'.

    @rtype: dictionary
    @return: A self-describing dictionary reporting the comparison results.
    It can contains path to plot image files generated by the comparison,
    along with some statistics.
    If there is no difference, an empty dictionary is returned.
    """
    full_data_dir = os.path.realpath(data_dir)
    full_ref_data_dir = os.path.realpath(ref_data_dir)
    print "Comparing {}\nwith {}\n".format(full_data_dir, full_ref_data_dir)

    data_file_list = set(file_list(data_dir, ))
    ref_data_file_list = set(file_list(ref_data_dir))
    difference = {}

    different_content = {}
    for data_file in data_file_list & ref_data_file_list:
        data_path = os.path.join(data_dir, data_file)
        ref_data_path = os.path.join(ref_data_dir, data_file)

        def callback(name, diff):
            if diff:
                different_content[name] = diff
                if 'different' not in difference:
                    difference['different'] = different_content
                if show_figure and 'plot' in diff:
                    show_img(diff['plot'])

        task_runner.append(compare_file,
                           (data_file, data_path, ref_data_path, kwargs),
                           partial(callback, data_file))

    unexpected_list = list(data_file_list - ref_data_file_list)
    lacking_list = list(ref_data_file_list - data_file_list)

    if unexpected_list:
        unexpected_dict = {}
        for unexpected_file in unexpected_list:
            unexpected_size = os.path.getsize(os.path.join(data_dir,
                                                           unexpected_file))
            unexpected_dict[unexpected_file] = 'Unexpected "{}" ({:,} bytes).'\
                                               .format(unexpected_file,
                                                       unexpected_size)
        difference['unexpected'] = unexpected_dict

    if lacking_list:
        lacking_dict = {}
        for lacking_file in lacking_list:
            lacking_size = os.path.getsize(os.path.join(ref_data_dir,
                                                        lacking_file))
            lacking_dict[lacking_file] = 'Cannot find "{}" ({:,} bytes).'\
                                    .format(lacking_file,
                                            lacking_size)
        difference['lacking'] = lacking_dict

    return difference


def compare_file(name, data_path, ref_data_path, kwargs):
    if fnmatch.fnmatch(name, "*.cfg"):
        # "*.cfg" files are not compared, they should be able to change
        # without breaking the tests, provided that the final results
        # are similar.
        return None
    elif fnmatch.fnmatch(name, "*.bin"):
        return compare_bin_file(name, data_path, ref_data_path, **kwargs)
    else:
        return compare_raw_file(data_path, ref_data_path)


def compare_bin_file(name, bin_path, ref_bin_path,
                     config=None,
                     show_figure=False,
                     with_map=True,
                     nb_image_t=5,
                     nb_image_z=5,
                     comparison_func=compare_array):
    """Compares two ".bin" files.

    A self-describing dictionary reporting the differences is generated.
    When a difference is detected, a plot representing
    this difference is generated along with some statistics.

    @type name
    @param name: The name of the data that is to be compared.

    @type bin_path
    @param bin_path: The '.bin' file path whose content is to be compared.

    @type ref_bin_path: string
    @param ref_bin_path: The '.bin' file path whose content is the comparison
    reference.

    @type config: Talos.Config or string
    @param config: The configuration that describes the domain of the data.
    If no configuration is given, ".bin" files are just compared for size and
    content perfect equality, no plot nor statistics are generated.

    @type show_figure: boolean
    @param show_figure: True if 'show' is called on the Matplotlib figures,
    besides the generation of image files.

    @type with_map: boolean
    @param with_map: True if a map is drawn on the Matplotlib figures.
    The map will be drawn only if there is also a 'config' in argument.

    @type nb_image_t: positive integer
    @param nb_image_t: The number of instants displayed in the plot grid.

    @type nb_image_z: positive integer
    @param nb_image_z: The number of levels displayed in the plot grid.

    @rtype: dictionary
    @return: A self-describing dictionary reporting the comparison results.
    It can contains path to plot image files generated by the comparison, along
    with some statistics.
    If there is no difference, an empty dictionary is returned.
"""
    ## Compares the size.
    size_difference = compare_file_size(bin_path, ref_bin_path)
    if size_difference:
        return size_difference

    ## Compares the values.
    # The data is assumed to be (x,y,z) or (x,y,z,t)
    # Other array shapes are not correctly interpreted (concretely, the plot is
    # distorted) but it is a limitation of the current data encoding,
    # When HDF5 metadata will be used, this will be solved.
    if config:
        config = as_config(config)
        array = getd(filename=bin_path, config=config, Nt=0)
        ref_array = getd(filename=ref_bin_path, config=config, Nt=0)
        Nt = ref_array.shape[0]
        Nz = ref_array.shape[1]
    else:
        # Without config file, binary files are interpreted as mono dimensional.
        array = fromfile(bin_path, 'f')
        ref_array = fromfile(ref_bin_path, 'f')
        Nz = len(array)
        Nt = 1

    difference = comparison_func(array, ref_array)
    if not difference:
        # No difference found.
        return {}

    diff_array = None
    if 'diff_array' in difference:
        diff_array = difference['diff_array']
        del difference['diff_array']

    if diff_array is not None and config:
        value_map = normalize_colormap([array, ref_array], cmap='PiYG',
                                       log_normed=True)
        diff_map = normalize_colormap(diff_array, cmap='seismic',
                                      log_normed=True)
        cmap_kwargs = [value_map, value_map, diff_map]
        map_error = [None]

        def plot_colormap(fig, ax, img, n, t, z):
            cs = None
            if with_map and config:
                try:
                    m = getm(config, open_figure=False)
                    cs = contourf_map(ax, m, img, **cmap_kwargs[n])
                except Exception, e:
                    map_error[0] = e
            if cs is None:
                cs = ax.contourf(img, **cmap_kwargs[n])
            fig.colorbar(cs, ax=ax)

        png_filename = os.path.splitext(bin_path)[0] + ".png"
        try:
            plot_array_list(array_list=[array, ref_array, diff_array],
                            title=os.path.splitext(name)[0],
                            array_name_list=["current", "reference", "difference"],
                            nb_image_t=min(nb_image_t, Nt),
                            nb_image_z=min(nb_image_z, Nz),
                            plot_func=plot_colormap,
                            png_filename=png_filename,
                            show_figure=show_figure)
            difference['plot'] = png_filename
        except Exception, e:
            print "[ERROR] Cannot plot a comparison for \"{}\":\n{}\n\n"\
                    .format(name, e)

        if map_error[0]:
            print "[WARNING] Cannot display a map layer for \"{}\":\n{}\n\n"\
                    .format(name, map_error[0])

    return difference


def compare_raw_file(filepath, ref_filepath):
    ## Compares the size.
    size_difference = compare_file_size(filepath, ref_filepath)
    if size_difference:
        return size_difference

    ## Compares the values.
    if not filecmp.cmp(filepath, ref_filepath, False):
        return "Same size but different value.".format(filepath, ref_filepath)

    ## No difference.
    return None


def compare_file_size(filepath, ref_filepath):
    byte_size = os.path.getsize(filepath)
    ref_byte_size = os.path.getsize(ref_filepath)
    if byte_size != ref_byte_size:
        return "Size of {:,} bytes whereas reference file has {:,} bytes."\
                                             .format(byte_size, ref_byte_size)
    ## No difference.
    return None


def assert_no_difference(*difference_list):
    difference = {}
    for d in difference_list:
        difference.update(d)
    assert not difference, "The following differences were detected:\n" + \
                           pformat(difference)

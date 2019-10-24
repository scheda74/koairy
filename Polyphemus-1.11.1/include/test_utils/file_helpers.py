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

import fnmatch
import os


def file_list(root_dir='.',
              file_pattern=None,
              file_exclude_pattern="*.png",
              dir_exclude_pattern=".*"):
    """Returns the list of the filename that matches the arguments."""
    root_dir = os.path.normpath(root_dir)
    path_begin = len(root_dir) + 1
    for path, dirname_list, filename_list in os.walk(root_dir):
        path = path[path_begin:]
        if file_exclude_pattern:
            for filename in fnmatch.filter(filename_list,
                                           file_exclude_pattern):
                filename_list.remove(filename)
        if file_pattern:
            filename_list = fnmatch.filter(filename_list, file_pattern)
        for filename in filename_list:
            yield os.path.normpath(os.path.join(path, filename))
        if dir_exclude_pattern:
            for dirname in fnmatch.filter(dirname_list, dir_exclude_pattern):
                dirname_list.remove(dirname)


def normalize_path(path):
    path = os.path.expanduser(os.path.expandvars(path))
    path = os.path.realpath(path)
    return path


def check_path(path):
    if not os.path.exists(path):
        raise IOError, "'%s' not found." % path


def remove_duplicate_path(path_list):
    normalized_path_list = [ normalize_path(d) + "/" for d in path_list]
    normalized_path_list.sort(key=len, reverse=True)
    root_path_list = []
    for i, p in enumerate(normalized_path_list):
        for p_shorter in normalized_path_list[i+1:]:
            if p.startswith(p_shorter):
                p = None
            break
        if p is not None:
            root_path_list.append(p[:-1])
    return root_path_list

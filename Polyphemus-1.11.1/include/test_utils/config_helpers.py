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

from glob import glob
import os
import re

from atmopy import talos

# Patterns used with the multi-line mode (re.M) of the regexp module:
patterns = dict(sp=r"[ \t\r\f\v]*",  # \s would include line returns.
                eq=r"[:=]",
                val=r"[^#$]*",  # values include spaces
                rest_of_line=r"[^$]*")


def as_config(config):
    """Returns a 'talos.Config' instance from 'config'.

    - If 'config' is a string, it is considered as a config filename,
    - If 'config' is already a 'talos.Config' object, it is returned as is.
    """
    if isinstance(config, talos.Config):
        return config

    if not os.path.isfile(config):
        raise ValueError, "Cannot open config file \"{}\".".format(config)

    additional_content = [
        ("Configuration_file", "[output]", "saver_config", "String"),
        ("Output_file", "[save]", "String")
    ]
    return talos.Config(config, additional_content)


def get_result_dir_list(config):
    """Returns the simulation output directories from 'config'.

    Because of a Talos Python wrapper limitation with multiple section of the
    same name, a regex is used."""
    config = as_config(config)

    with open(config.saver_config, "r") as f:
        config = f.read()

    match_list = re.findall(
        r"{sp}Output_file{sp}{eq}{sp}({val})"\
        .format(**patterns), config, re.M)

    return set([os.path.normpath(os.path.dirname(v.strip())) for v in match_list])


def replace_value(config_filename, variable_name, new_value):
    """Replaces a variable value, putting the former value in comment.

    The current value is put in a comment of the form '# was: ...' unless there
    was already a former value commented this way. If there was already a
    commented value, the current value is simply overwritten and forgotten.

    The argument 'config_filename' can be a normal filename, or a globing
    expression using wild cards to match one or several files.

    The argument 'new_value' does not need to be a string, since what is
    written is the result of 'str(new_value)'.
    """

    new_value = str(new_value)
    for filename in glob(config_filename):
        with open(filename, "r+") as f:
            config = f.read()

            # The case where there is already a commented value:
            config, n = re.subn(
                r"({sp}{variable_name}{sp}{eq})({val})(# was: {rest_of_line})$"\
            .format(variable_name=variable_name, **patterns),
                r"\1 %s \3" % new_value, config, re.M)
            if n == 0:
                # The case there is no commented value yet:
                config, n = re.subn(
                    r"({sp}{variable_name}{sp}{eq})({rest_of_line})$"\
                        .format(variable_name=variable_name, **patterns),
                    r"\1 %s # was:\2" % new_value, config, re.M)

            if n == 0:
                print "No '%s' variable in \"%s\"." % (variable_name, filename)
            else:
                f.seek(0)
                f.write(config)
                f.truncate()
                print "'%s' variable has been set to '%s' in \"%s\"." \
                        % (variable_name, new_value, filename)


def put_back_value(config_filename, variable_name):
    """Replaces a variable value by the value commented alongside.

    If there is a value in a comment of the form '# was: ...', it is
    uncommented and the current value is removed.

    The argument 'config_filename' can be a normal filename, or a globing
    expression using wild cards to match one or several files.
    The values commented by 'replace_value' are well recognized and put back.
    """
    for filename in glob(config_filename):
        with open(filename, "r+") as f:
            config = f.read()
            config, n = re.subn(
                r"({sp}{variable_name}{sp}{eq}){val}# was:({rest_of_line})$"\
                .format(variable_name=variable_name, **patterns),
                r"\1\2", config, re.M)
            if n == 0:
                print "No '%s' variable in \"%s\"." % (variable_name, filename)
            else:
                f.seek(0)
                f.write(config)
                f.truncate()
                print "'%s' variable put back in \"%s\"." \
                        % (variable_name, filename)

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

import py

from pprint import pprint

from test_utils.pytest_ext import tmprundir, recursive_equal
from test_utils.file_helpers import file_list
from test_utils.display.matplotlib_ext import dummy_data

from test_utils.compare_data import compare_result
from test_utils.ParallelTaskRunner import ParallelTaskRunner


def create_file(path, content="dummy file"):
    path = py.path.local(path)
    if isinstance(content, str):
        path.write(content, ensure=True)
    else:
        path.write_binary(content.tobytes(), ensure=True)


def test_file_list(tmprundir):
    # Some valid files:
    valid_file_list = ['a.bin', 'b.bin', 'A/c.bin', 'A/B/d.bin']
    valid_file_list.sort()

    # Some garbage files:
    garbage_file_list = ['A/B/bad.txt', 'A/.git/bad.bin', 'A/B/bad.bin2',\
                         'A/B/badbin']

    for f in valid_file_list + garbage_file_list:
        create_file(f)

    retrieved_file_list = list(file_list(file_pattern="*.bin"))
    retrieved_file_list.sort()

    assert retrieved_file_list == valid_file_list


def create_config(config_dir, result_dir, **domain):
    saver_config = config_dir + "/saver_config.cfg"
    saver_config_content = """
[save]

Output_file: {result_dir}/&f.bin
""".format(result_dir=result_dir)
    create_file(saver_config, saver_config_content)

    config = config_dir + "/config.cfg"
    config_content = """
[domain]
Nx = {Nx}
Ny = {Ny}
Nz = {Nz}
Nt = {Nt}
[output]
Configuration_file: {saver_config}
""".format(saver_config=saver_config, **domain)

    create_file(config, config_content)
    return config


def test_compare_result(tmprundir):

    ref_base_path = "ref"
    base_path = "."
    result_subdir = "results"

    ref_result_path = ref_base_path + "/" + result_subdir + "/"
    result_path = base_path + "/" + result_subdir + "/"

    ## Creates a dummy config.
    domain = dict(Nx=64, Ny=128, Nz=1, Nt=32)
    config = create_config(base_path + "/config/", result_subdir, **domain)

    ## Creates dummy results.
    default_data = dummy_data(**domain)

    # Same data:
    dataA = dummy_data(speed=0.2, **domain)
    create_file(ref_result_path + "sameA.bin", dataA)
    create_file(result_path + "sameA.bin", dataA)

    dataB = dummy_data(amplitude=2, **domain)
    create_file(ref_result_path + "sameB.bin", dataB)
    create_file(result_path + "sameB.bin", dataB)

    # Should not detect any difference at this stage:
    difference = compare_result(config, base_path, ref_base_path,
                                show_figure=False, task_runner=ParallelTaskRunner(1))

    if difference:
        print "Difference found:"
        pprint(difference)
    assert not difference

    # Lacking data:
    create_file(ref_result_path + "lacking.bin", default_data)

    # Excess data:
    create_file(result_path + "unexpected.bin", default_data)

    # Different domain:
    bad_domain = domain.copy()
    bad_domain['Nt'] = domain['Nt'] * 2
    create_file(ref_result_path + "bad_domain.bin", default_data)
    create_file(result_path + "bad_domain.bin", dummy_data(**bad_domain))

    # Different values:
    create_file(ref_result_path + "diffA.bin", dummy_data(speed=0.3, **domain))
    create_file(result_path + "diffA.bin", dummy_data(speed=0.2, **domain))

    create_file(ref_result_path + "diffB.bin", dummy_data(amplitude=-1,
                                                          **domain))
    create_file(result_path + "diffB.bin", dummy_data(amplitude=1, **domain))

    run_dir = str(tmprundir) + '/'
    expected_report = { run_dir + 'results': {
        'different': {'bad_domain.bin': 'Size of 2,097,152 bytes whereas reference file has 1,048,576 bytes.',
                      'diffB.bin': {'different sign': {'cur': 'contains some strictly positive values',
                                                       'ref': 'all values are negative or null'},
                                    'mean (rtol=1.00e-02)': {'cur': 0.12500029871238896,
                                                             'quotient': -1.0,
                                                             'ref': -0.12500029871238896},
                                                             'plot': run_dir + 'results/diffB.png'}},
        'lacking': {'lacking.bin': 'Cannot find "lacking.bin" (1,048,576 bytes).'},
        'unexpected': {'unexpected.bin': 'Unexpected "unexpected.bin" (1,048,576 bytes).'}}}

    report = compare_result(config, base_path, ref_base_path,
                            show_figure=False, task_runner=ParallelTaskRunner(1))
    pprint(report)
    pprint(expected_report)
    assert recursive_equal(report, expected_report)

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
import os, sys
import pytest

from test_utils.notebook.NotebookKernel import (NotebookKernel,
                                                NotebookException)


def pytest_addoption(parser):
    group = parser.getgroup("general")
    group._addoption(
        '--ipynb-incremental', action="store_true",
        dest="ipynb_incremental", default=False,
        help="Only runs notebooks that have some cells without output or "\
             "having an error."
    )
    group._addoption(
        '--ipynb-verbose', action="store_true",
        dest="ipynb_verbose", default=False,
        help="Displays information on success and failure as soon as a test "\
             "terminates."
    )


def pytest_collect_file(path, parent):
    """Collection hook for Pytest.

    - It is assumed that any notebook in the crawled directories is a test.
    - When there are several notebooks in the same directory, they are run
      serially, in lexicographic order, as a single big test.
    """
    if path.ext != ".ipynb":
        return None

    notebooks = glob(path.dirname + "/*.ipynb")
    notebooks.sort()
    if not path.samefile(notebooks[0]):
        return None

    if parent.config.option.ipynb_incremental:
        already_successful = True
        for nb_path in notebooks:
            with NotebookKernel(nb_path) as nb:
                if not nb.has_successfully_run():
                    already_successful = False
                    break
        if already_successful:
            return None
    return NotebookCollection(path, parent, notebooks)


def pytest_ignore_collect(path):
    """Collection ignore hook for Pytest.

        - It is assumed that a ".py" script with the same base name as a
          ".ipynb" notebook cannot be a test (it is likely the notebook
          converted to a script).
    """
    # Checks path is not the translation of a IPython Notebook.
    if path.ext == ".py":
        if os.path.exists(os.path.splitext(str(path))[0] + ".ipynb"):
            return True
    return path.check(link=1)


class NotebookCollection(pytest.File):
    """A collection of IPython notebooks.

    The notebooks are run serially in the order given at construction.
    The execution stops at the first error: remaining cells and notebooks are
    not executed.
    """

    def __init__(self, path, parent, notebooks):
        super(NotebookCollection, self).__init__(path, parent)

        self.notebooks = notebooks

        if len(self.notebooks) > 1:
            basenames = [os.path.splitext(os.path.basename(p))[0]
                         for p in self.notebooks]
            self.name = path.dirname + "/{" + ','.join(basenames) + "}.ipynb"

    def collect(self):
        yield NotebookTestScenario(self.name, self)


class NotebookTestScenario(pytest.Item):
    """A test scenario composed of one or multiple IPython notebooks."""

    def __init__(self, name, parent):
        super(NotebookTestScenario, self).__init__(name, parent)
        self.notebooks = parent.notebooks
        self.verbose = parent.config.option.ipynb_verbose
        self.current_nb_path = self.fspath
        self.current_nb_name = name

    def runtest(self):
        """Reads and runs notebooks' cells one by one."""
        # It is assumed that more than a week of run duration is overkilled.
        timeout_in_sec = 7 * 24 * 60 * 60
        for nb_path in self.notebooks:
            self.current_nb_path = nb_path
            self.current_nb_name = os.path.basename(nb_path)
            with NotebookKernel(nb_path) as nb:
                # Blanks the notebook output to avoid confusion on what
                # has effectively been run.
                nb.strip_output()

                for cell_index, cell in enumerate(nb.cells()):
                    if cell.cell_type != "code":
                        continue
                    try:
                        nb.run(cell_index, timeout_in_sec)
                    except Exception, e:
                        self.verbose_print("[FAILURE]", str(e), cell_index)
                        raise
                self.verbose_print("[SUCCESS]")

    def repr_failure(self, excinfo):
        """Called when self.runtest() raises an exception."""
        if isinstance(excinfo.value, NotebookException):
            return excinfo.value.pretty_str()
        else:
            tb = excinfo.getrepr(showlocals=True, style='long', abspath=False)
            return "*** Internal unexpected exception:\n%s\n"\
                   "Exception message: %s" % (tb, str(excinfo.value))

    def reportinfo(self):
        """Returns the test location, which is displayed in verbose mode."""
        return self.current_nb_path, 0, \
            "test scenario: %s" % self.current_nb_name

    def verbose_print(self, status="[INFO]", message=None, cell_index=None):
        if not self.verbose:
            return
        # By default, sys.stdout is wrapped by py.test,
        # the trick is to use the hidden attribute as
        # py.test does.
        write = sys.__stdout__.write

        write("\n" + status + " in " + self.current_nb_name)
        if cell_index:
            write(", cell index " + str(cell_index))
        write(" (at \""+ self.current_nb_path + ")")
        if message:
            write(":\n" + message + "\n")

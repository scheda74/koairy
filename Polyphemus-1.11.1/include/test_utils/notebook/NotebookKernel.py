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

# Some resources:
#
# runipy implementation for running a notebook:
# https://github.com/paulgb/runipy/blob/master/runipy/notebook_runner.py
#
# running a notebook through pytest
# https://gist.github.com/timo/2621679
#
# running a notebook and compare its result with a reference
# (whose data is the current notebook data):
# https://gist.github.com/minrk/2620735
# https://gist.github.com/akimd/f60b2f26ea70c72db175
#
# running notebook where the notebook is in another process:
# http://pydoc.net/Python/ipython/0.13/IPython.zmq.tests.test_embed_kernel/
#
# if "Too many open files" error comes up,
# see https://gist.github.com/minrk/2620735#comment-918564
# in substance, one should use only one kernel, and
# reset it before each use with km.reset_kernel(now=True)

import os
import py
import sys

try:
    # Python 3
    from queue import Empty
except ImportError:
    # Python 2
    from Queue import Empty

# An ad hoc compatibility layer for the different ipython versions.
# Should be removed when ready to switch to ipython 3.0 with the notebook
# format v4.
# See https://github.com/ipython/ipython/wiki/IPEP-17:-Notebook-Format-4
#     http://ipython.org/ipython-doc/dev/notebook/nbformat.html
try:
    try:
        # jupyter (the successor and generalization of ipython to multiple
        # languages)
        from nbformat import reads, writes, NO_CONVERT
        from nbformat.v3 import NotebookNode
    except ImportError:
        # ipython 3+ (don't confuse with python 3, it is unrelated)
        from IPython.nbformat import reads, writes, NO_CONVERT
        from IPython.nbformat.v3 import NotebookNode
    nb_format_version = 3  # Will be 4 when ready to drop ipython 2.
    read_notebook = lambda s: reads(s, nb_format_version)
    write_notebook = lambda s: writes(s, NO_CONVERT)
except:
    # ipython 2
    from IPython.nbformat.current import reads, writes, NotebookNode
    read_notebook = lambda s: reads(s, 'json')
    write_notebook = lambda s: writes(s, 'json')

# Some libraries, eg. Pytest, set a special syd.stdin that would interfere
# with the Notebook manager. The real stdin must be put back during the
# IPython kernel import.
wrapped_stdin = sys.stdin
sys.stdin = sys.__stdin__
try:
    try:
        # jupyter (the successor and generalization of ipython to multiple
        # languages)
        from jupyter_client.manager import KernelManager
        from jupyter_client.manager import start_new_kernel
    except ImportError:
        # ipython 3+ (don't confuse with python 3, it is unrelated)
        from IPython.kernel.manager import KernelManager
        from IPython.kernel.manager import start_new_kernel
except:
    # ipython 2 does not have 'start_new_kernel', here is a backport.
    def start_new_kernel(startup_timeout=60, kernel_name='python', **kwargs):
        """Start a new kernel, and return its Manager and Client"""
        km = KernelManager(kernel_name=kernel_name)
        km.start_kernel(**kwargs)
        kc = km.client()
        kc.start_channels()
        kc.kernel_info()

        # Wait for kernel info reply on shell channel
        while True:
            msg = kc.get_shell_msg(block=True)
            if msg['msg_type'] == 'kernel_info_reply':
                break

        # Flush IOPub channel
        while True:
            try:
                msg = kc.get_iopub_msg(block=True, timeout=0.2)
                print(msg['msg_type'])
            except Empty:
                break

        return km, kc
finally:
    sys.stdin = wrapped_stdin


def highlighted_source(ipython_source, style='monokai'):
    """Returns an ANSI colored version of the IPython source code in argument.

    If coloring service is not available, the source code is unchanged.
    """
    try:
        try:
            from IPython.lib.lexers import IPythonLexer
        except:
            # Is deprecated with IPython 4+.
            from IPython.nbconvert.utils.lexers import IPythonLexer
        from pygments import highlight
        from pygments.formatters.terminal256 import Terminal256Formatter
    except:
        return ipython_source

    return highlight(ipython_source, IPythonLexer(),
                     Terminal256Formatter(style=style))


class NotebookException(Exception):
    """Custom exception for an IPython notebook kernel."""

    def __init__(self, message,
                 notebook_path="?",
                 cell_index=-1,
                 code=None,
                 text_output=None,
                 trace_back=None):
        super(Exception, self).__init__(message)
        self.cell_index = cell_index
        self.notebook_path = notebook_path
        self.code = code
        self.text_output = text_output
        self.trace_back = trace_back

    def pretty_str(self):
        msg = self.args[0]
        if self.code is not None:
            code = highlighted_source(self.code)
        else:
            code = "No code"
        if self.trace_back:
            trace_back = self.trace_back
        else:
            trace_back = "No trace back"
        if self.text_output:
            text_output = self.text_output
        else:
            text_output = "No text output"
        delim = '-' * 64

        return "\n".join(
            ["==== Notebook execution failed:", msg, "==== Path: %s" %
             self.notebook_path, "==== Cell number: %d" % self.cell_index,
             "==== Cell code:", delim, code, delim, "",
             "==== Cell text output:", delim, text_output, delim, "",
             "==== Traceback:", str(trace_back)])

    pass


class NotebookKernel(object):
    """Context manager for an IPython notebook kernel and its notebook.

    In addition of the kernel resource management, some helpers are provided
    notably for running code.
    """

    MIME_MAP = {
        'text/plain': 'text',
        'text/html': 'html',
        'image/svg+xml': 'svg',
        'image/png': 'png',
        'image/jpeg': 'jpeg',
        'text/latex': 'latex',
        'application/json': 'json',
        'application/javascript': 'javascript',
    }

    CONV_TABLE = {
        'prompt_number': 'execution_count',
        'pyin': 'execute_input',
        'pyout': 'execute_result',
        'pyerr': 'error',
        'text': 'data'
    }

    FALLBACK_TABLE = dict((v, k) for k, v in CONV_TABLE.iteritems())

    def __init__(self, notebook_path, **kernel_kwargs):
        self.notebook_path = notebook_path
        self.kernel_kwargs = kernel_kwargs

    def __enter__(self):
        """Starts a kernel and loads the associated notebook."""
        with open(self.notebook_path, 'r') as f:
            self.nb = read_notebook(f.read())

        notebook_dirname = os.path.dirname(self.notebook_path)
        self.km, self.kc = start_new_kernel(stderr=open(os.devnull, 'w'),
                                            cwd=notebook_dirname,
                                            **self.kernel_kwargs)

        return self

    def cells(self):
        """Returns the notebook cell list."""
        # self.nb.cells in v4+
        return self.nb.worksheets[0].cells

    def has_successfully_run(self):
        """True if the notebook was entirely and successfully run."""
        for cell in self.cells():
            if cell.cell_type != "code":
                continue
            # 'execution_count' in v4
            if 'prompt_number' not in cell or cell['prompt_number'] == None:
                return False
            if 'outputs' in cell:
                for out in cell['outputs']:
                    # 'error' in v4+
                    if out.output_type == 'pyerr':
                        return False
        return True

    def strip_output(self, start_index=0):
        """Removes any output from the notebook."""
        # nb.cells in v4
        for cell in self.cells()[start_index:]:
            # 'execution_count' in v4
            if 'prompt_number' in cell:
                cell['prompt_number'] = None
            if 'outputs' in cell:
                cell['outputs'] = []
            # 'autoscroll' in v4
            cell.metadata.pop('scrolled', 0)
            # cell.metadata.pop('collapsed', 0) in v4
            if cell.cell_type == "code":
                cell['collapsed'] = False

    def run(self, cell_index, timeout_in_sec, update_notebook=True):
        """Runs the cell given in argument, raising an exception on error."""
        try:
            cell = self.cells()[cell_index]
            # cell['source'] in v4+
            code = cell['input']
            self.kc.execute(code, allow_stdin=False)
            msg = self.kc.get_shell_msg(timeout=timeout_in_sec)
        except Exception as e:
            excinfo = py.code.ExceptionInfo()
            trace_back = excinfo.getrepr(showlocals=True,
                                         style='short',
                                         abspath=False)
            msg = '%s' % excinfo.exconly()
            if type(e) == Empty:
                # Gives an hint when the mysterious 'Empty' exception is raised.
                msg += "\nEmpty message from the ipython kernel after running "\
                "a cell, likely to be a timeout, did your cell run duration "\
                "exceeded %d minutes ?" % (timeout_in_sec / 60)
            raise NotebookException("Internal IPython notebook error: %s" % msg,
                                    self.notebook_path, cell_index, code,
                                    trace_back=trace_back)

        reply = msg['content']

        outputs = []
        text_node = None
        while True:
            msg = self.kc.get_iopub_msg(timeout=5)
            msg_type = msg['msg_type']
            content = msg['content']

            # Some helpers to deal with the different protocol versions.
            # Will be removed in v4
            conv = lambda k: self.CONV_TABLE.get(k, k)
            fallback = lambda k: self.FALLBACK_TABLE.get(k, k)
            get = lambda k: content.get(k, content.get(fallback(k)))

            msg_type = conv(msg_type)
            if msg_type == 'status':
                if content['execution_state'] == 'idle':
                    break
                continue
            if msg_type == 'clear_output':
                outputs = []
                continue
            if msg_type == 'execute_input':
                execution_count = get('execution_count')
                if execution_count:
                    # cell['execution_count'] in v4
                    cell['prompt_number'] = int(execution_count)
                continue

            out = NotebookNode(output_type=fallback(msg_type))

            if msg_type == 'stream':
                # out.name in v4
                out.stream = content['name']
                text = get('text')
                if text_node is None:
                    text_node = out
                else:
                    # Already have a text node:
                    out = None
                if text is not None:
                    # out.data in v4
                    if 'text' in text_node:
                        # Already have text in the text node:
                        text_node.text += text
                    else:
                        text_node.text = text
            elif msg_type in ('display_data', 'execute_result'):

                # Seems mandatory, not sure why since it is redundant with
                # the execution count of the 'execute_input' message.
                if msg_type == 'execute_result':
                    out.prompt_number = get('execution_count')

                for mime, encoded_data in content['data'].items():
                    # Short MIME name will be removed in v4
                    try:
                        mime_short_name = self.MIME_MAP[mime]
                    except KeyError as e:
                        raise NotebookException("Unhandled mime type: %s" % mime,
                                                self.notebook_path, cell_index,
                                                code)
                    setattr(out, mime_short_name, encoded_data)
            elif msg_type == 'error':
                out.ename = content['ename']
                out.evalue = content['evalue']
                out.traceback = content['traceback']
            else:
                raise NotebookException("Unhandled iopub message: %s" % msg_type,
                                        self.notebook_path, cell_index, code)

            if out:
                outputs.append(out)

        cell['outputs'] = outputs

        text_output = None
        if text_node and 'text' in text_node:
            text_output = text_node.text

        if reply['status'] == 'error':
            raise NotebookException("Error while running a cell",
                                    self.notebook_path, cell_index, code,
                                    text_output, '\n'.join(reply['traceback']))

    def save(self):
        with open(self.notebook_path, 'w') as f:
            f.write(write_notebook(self.nb))

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """Shuts down the kernel and frees associated resources."""
        self.save()
        if hasattr(self, 'kc'):
            self.kc.stop_channels()
            del self.kc
        if hasattr(self, 'km'):
            self.km.shutdown_kernel(now=True)
            del self.km

# -*- coding: utf-8 -*-
r"""
============================
Parser for Jupyter notebooks
============================

Class that holds the Jupyter notebook information

"""
# Author: Óscar Nájera
# License: 3-clause BSD

from __future__ import division, absolute_import, print_function
from functools import partial
import json
import os
import re
import sys


def ipy_notebook_skeleton():
    """Returns a dictionary with the elements of a Jupyter notebook"""
    py_version = sys.version_info
    notebook_skeleton = {
        "cells": [],
        "metadata": {
            "kernelspec": {
                "display_name": "Python " + str(py_version[0]),
                "language": "python",
                "name": "python" + str(py_version[0])
            },
            "language_info": {
                "codemirror_mode": {
                    "name": "ipython",
                    "version": py_version[0]
                },
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython" + str(py_version[0]),
                "version": '{0}.{1}.{2}'.format(*sys.version_info[:3])
            }
        },
        "nbformat": 4,
        "nbformat_minor": 0
    }
    return notebook_skeleton


def directive_fun(match, directive):
    """Helper to fill in directives"""
    directive_to_alert = dict(note="info", warning="danger")
    return ('<div class="alert alert-{0}"><h4>{1}</h4><p>{2}</p></div>'
            .format(directive_to_alert[directive], directive.capitalize(),
                    match.group(1).strip()))


def rst2md(text):
    """Converts the RST text from the examples docstrigs and comments
    into markdown text for the Jupyter notebooks"""

    top_heading = re.compile(r'^=+$\s^([\w\s-]+)^=+$', flags=re.M)
    text = re.sub(top_heading, r'# \1', text)

    math_eq = re.compile(r'^\.\. math::((?:.+)?(?:\n+^  .+)*)', flags=re.M)
    text = re.sub(math_eq,
                  lambda match: r'\begin{{align}}{0}\end{{align}}'.format(
                      match.group(1).strip()),
                  text)
    inline_math = re.compile(r':math:`(.+?)`', re.DOTALL)
    text = re.sub(inline_math, r'$\1$', text)

    directives = ('warning', 'note')
    for directive in directives:
        directive_re = re.compile(r'^\.\. %s::((?:.+)?(?:\n+^  .+)*)'
                                  % directive, flags=re.M)
        text = re.sub(directive_re,
                      partial(directive_fun, directive=directive), text)

    links = re.compile(r'^ *\.\. _.*:.*$\n', flags=re.M)
    text = re.sub(links, '', text)

    refs = re.compile(r':ref:`')
    text = re.sub(refs, '`', text)

    contents = re.compile(r'^\s*\.\. contents::.*$(\n +:\S+: *$)*\n',
                          flags=re.M)
    text = re.sub(contents, '', text)

    images = re.compile(
        r'^\.\. image::(.*$)(?:\n *:alt:(.*$)\n)?(?: +:\S+:.*$\n)*',
        flags=re.M)
    text = re.sub(
        images, lambda match: '![{1}]({0})\n'.format(
            match.group(1).strip(), (match.group(2) or '').strip()), text)

    return text


class Notebook(object):
    """Jupyter notebook object

    Constructs the file cell-by-cell and writes it at the end"""

    def __init__(self, file_name, target_dir):
        """Declare the skeleton of the notebook

        Parameters
        ----------
        file_name : str
            original script file name, .py extension will be renamed
        target_dir: str
            directory where notebook file is to be saved
        """

        self.file_name = file_name.replace('.py', '.ipynb')
        self.write_file = os.path.join(target_dir, self.file_name)
        self.work_notebook = ipy_notebook_skeleton()
        self.add_code_cell("%matplotlib inline")

    def add_code_cell(self, code):
        """Add a code cell to the notebook

        Parameters
        ----------
        code : str
            Cell content
        """

        code_cell = {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {"collapsed": False},
            "outputs": [],
            "source": [code.strip()]
            }
        self.work_notebook["cells"].append(code_cell)

    def add_markdown_cell(self, text):
        """Add a markdown cell to the notebook

        Parameters
        ----------
        code : str
            Cell content
        """
        markdown_cell = {
            "cell_type": "markdown",
            "metadata": {},
            "source": [rst2md(text)]
        }
        self.work_notebook["cells"].append(markdown_cell)

    def save_file(self):
        """Saves the notebook to a file"""
        with open(self.write_file, 'w') as out_nb:
            json.dump(self.work_notebook, out_nb, indent=2)

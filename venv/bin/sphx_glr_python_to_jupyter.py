#!/home/gordon/development/scikit-learn/venv/bin/python3
# -*- coding: utf-8 -*-
r"""
Sphinx Gallery Notebook converter
=================================

Exposes the Sphinx-Gallery Notebook renderer to directly convert Python
scripts into Jupyter Notebooks.

"""
# Author: Óscar Nájera
# License: 3-clause BSD

from __future__ import division, absolute_import, print_function

from sphinx_gallery.notebook import python_to_jupyter_cli


if __name__ == '__main__':
    python_to_jupyter_cli()

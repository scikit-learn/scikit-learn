"""
    This module contains pythran frontend
"""

from pythran.openmp import GatherOMPData
from pythran.syntax import check_syntax
from pythran.transformations import ExtractDocStrings, HandleImport

import gast as ast
import re


def raw_parse(code):
    # hacky way to turn OpenMP comments into strings
    code = re.sub(r'(\s*)#\s*(omp\s[^\n]+)', r'\1"\2"', code)

    return ast.parse(code)


def parse(pm, code):

    # front end
    ir = raw_parse(code)

    # Handle user-defined import
    pm.apply(HandleImport, ir)

    # parse openmp directive
    pm.apply(GatherOMPData, ir)

    # extract docstrings
    _, docstrings = pm.apply(ExtractDocStrings, ir)

    # avoid conflicts with cxx keywords
    check_syntax(ir)
    return ir, docstrings

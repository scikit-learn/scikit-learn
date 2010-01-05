#! /usr/bin/env python

"""Some basic tools for common operations in datasets."""

import textwrap

def dumpvar(var, varname, wrap = 79):
    """return a list of string representing the variable var with name varname.

    For example, if var is [1, 2, 3, 4] and varname is l, this will return the
    list ["l = [1, 2, 3, 4]"]. Each item in the list is a wrapped line, using
    the value in wrap."""
    strvar = varname + " = " + str(var)
    l = [i + '\n' for i in textwrap.wrap(strvar, wrap)]
    l.append('\n')
    return l

# encoding: utf-8
"""Tests for IPython.utils.text"""
from __future__ import print_function

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import os
import math
import random
import sys

import nose.tools as nt

from .. import text


def test_columnize():
    """Basic columnize tests."""
    size = 5
    items = [l*size for l in 'abc']
    out = text.columnize(items, displaywidth=80)
    nt.assert_equal(out, 'aaaaa  bbbbb  ccccc\n')
    out = text.columnize(items, displaywidth=12)
    nt.assert_equal(out, 'aaaaa  ccccc\nbbbbb\n')
    out = text.columnize(items, displaywidth=10)
    nt.assert_equal(out, 'aaaaa\nbbbbb\nccccc\n')

def test_columnize_random():
    """Test with random input to hopfully catch edge case """
    for nitems in [random.randint(2,70) for i in range(2,20)]:
        displaywidth = random.randint(20,200)
        rand_len = [random.randint(2,displaywidth) for i in range(nitems)]
        items = ['x'*l for l in rand_len]
        out = text.columnize(items, displaywidth=displaywidth)
        longer_line = max([len(x) for x in out.split('\n')])
        longer_element = max(rand_len)
        if longer_line > displaywidth:
            print("Columnize displayed something lager than displaywidth : %s " % longer_line)
            print("longer element : %s " % longer_element)
            print("displaywidth : %s " % displaywidth)
            print("number of element : %s " % nitems)
            print("size of each element :\n %s" % rand_len)
            assert False

def test_columnize_medium():
    """Test with inputs than shouldn't be wider tahn 80 """
    size = 40
    items = [l*size for l in 'abc']
    out = text.columnize(items, displaywidth=80)
    nt.assert_equal(out, '\n'.join(items+['']))

def test_columnize_long():
    """Test columnize with inputs longer than the display window"""
    size = 11
    items = [l*size for l in 'abc']
    out = text.columnize(items, displaywidth=size-1)
    nt.assert_equal(out, '\n'.join(items+['']))


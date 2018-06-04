# encoding: utf-8
"""Tests for IPython.utils.text"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os
import math
import random
import sys

import nose.tools as nt
try:
    from pathlib import Path
except ImportError:
    # for Python 3.3
    from pathlib2 import Path

from IPython.utils import text

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

def test_columnize():
    """Basic columnize tests."""
    size = 5
    items = [l*size for l in 'abcd']

    out = text.columnize(items, displaywidth=80)
    nt.assert_equal(out, 'aaaaa  bbbbb  ccccc  ddddd\n')
    out = text.columnize(items, displaywidth=25)
    nt.assert_equal(out, 'aaaaa  ccccc\nbbbbb  ddddd\n')
    out = text.columnize(items, displaywidth=12)
    nt.assert_equal(out, 'aaaaa  ccccc\nbbbbb  ddddd\n')
    out = text.columnize(items, displaywidth=10)
    nt.assert_equal(out, 'aaaaa\nbbbbb\nccccc\nddddd\n')

    out = text.columnize(items, row_first=True, displaywidth=80)
    nt.assert_equal(out, 'aaaaa  bbbbb  ccccc  ddddd\n')
    out = text.columnize(items, row_first=True, displaywidth=25)
    nt.assert_equal(out, 'aaaaa  bbbbb\nccccc  ddddd\n')
    out = text.columnize(items, row_first=True, displaywidth=12)
    nt.assert_equal(out, 'aaaaa  bbbbb\nccccc  ddddd\n')
    out = text.columnize(items, row_first=True, displaywidth=10)
    nt.assert_equal(out, 'aaaaa\nbbbbb\nccccc\nddddd\n')

    out = text.columnize(items, displaywidth=40, spread=True)
    nt.assert_equal(out, 'aaaaa      bbbbb      ccccc      ddddd\n')
    out = text.columnize(items, displaywidth=20, spread=True)
    nt.assert_equal(out, 'aaaaa          ccccc\nbbbbb          ddddd\n')
    out = text.columnize(items, displaywidth=12, spread=True)
    nt.assert_equal(out, 'aaaaa  ccccc\nbbbbb  ddddd\n')
    out = text.columnize(items, displaywidth=10, spread=True)
    nt.assert_equal(out, 'aaaaa\nbbbbb\nccccc\nddddd\n')


def test_columnize_random():
    """Test with random input to hopefully catch edge case """
    for row_first in [True, False]:
        for nitems in [random.randint(2,70) for i in range(2,20)]:
            displaywidth = random.randint(20,200)
            rand_len = [random.randint(2,displaywidth) for i in range(nitems)]
            items = ['x'*l for l in rand_len]
            out = text.columnize(items, row_first=row_first, displaywidth=displaywidth)
            longer_line = max([len(x) for x in out.split('\n')])
            longer_element = max(rand_len)
            if longer_line > displaywidth:
                print("Columnize displayed something lager than displaywidth : %s " % longer_line)
                print("longer element : %s " % longer_element)
                print("displaywidth : %s " % displaywidth)
                print("number of element : %s " % nitems)
                print("size of each element :\n %s" % rand_len)
                assert False, "row_first={0}".format(row_first)

def test_columnize_medium():
    """Test with inputs than shouldn't be wider than 80"""
    size = 40
    items = [l*size for l in 'abc']
    for row_first in [True, False]:
        out = text.columnize(items, row_first=row_first, displaywidth=80)
        nt.assert_equal(out, '\n'.join(items+['']), "row_first={0}".format(row_first))

def test_columnize_long():
    """Test columnize with inputs longer than the display window"""
    size = 11
    items = [l*size for l in 'abc']
    for row_first in [True, False]:
        out = text.columnize(items, row_first=row_first, displaywidth=size-1)
        nt.assert_equal(out, '\n'.join(items+['']), "row_first={0}".format(row_first))

def eval_formatter_check(f):
    ns = dict(n=12, pi=math.pi, stuff='hello there', os=os, u=u"café", b="café")
    s = f.format("{n} {n//4} {stuff.split()[0]}", **ns)
    nt.assert_equal(s, "12 3 hello")
    s = f.format(' '.join(['{n//%i}'%i for i in range(1,8)]), **ns)
    nt.assert_equal(s, "12 6 4 3 2 2 1")
    s = f.format('{[n//i for i in range(1,8)]}', **ns)
    nt.assert_equal(s, "[12, 6, 4, 3, 2, 2, 1]")
    s = f.format("{stuff!s}", **ns)
    nt.assert_equal(s, ns['stuff'])
    s = f.format("{stuff!r}", **ns)
    nt.assert_equal(s, repr(ns['stuff']))
    
    # Check with unicode:
    s = f.format("{u}", **ns)
    nt.assert_equal(s, ns['u'])
    # This decodes in a platform dependent manner, but it shouldn't error out
    s = f.format("{b}", **ns)
        
    nt.assert_raises(NameError, f.format, '{dne}', **ns)

def eval_formatter_slicing_check(f):
    ns = dict(n=12, pi=math.pi, stuff='hello there', os=os)
    s = f.format(" {stuff.split()[:]} ", **ns)
    nt.assert_equal(s, " ['hello', 'there'] ")
    s = f.format(" {stuff.split()[::-1]} ", **ns)
    nt.assert_equal(s, " ['there', 'hello'] ")
    s = f.format("{stuff[::2]}", **ns)
    nt.assert_equal(s, ns['stuff'][::2])
    
    nt.assert_raises(SyntaxError, f.format, "{n:x}", **ns)

def eval_formatter_no_slicing_check(f):
    ns = dict(n=12, pi=math.pi, stuff='hello there', os=os)
    
    s = f.format('{n:x} {pi**2:+f}', **ns)
    nt.assert_equal(s, "c +9.869604")
    
    s = f.format('{stuff[slice(1,4)]}', **ns)
    nt.assert_equal(s, 'ell')
    
    if sys.version_info >= (3, 4):
        # String formatting has changed in Python 3.4, so this now works.
        s = f.format("{a[:]}", a=[1, 2])
        nt.assert_equal(s, "[1, 2]")
    else:
        nt.assert_raises(SyntaxError, f.format, "{a[:]}")

def test_eval_formatter():
    f = text.EvalFormatter()
    eval_formatter_check(f)
    eval_formatter_no_slicing_check(f)

def test_full_eval_formatter():
    f = text.FullEvalFormatter()
    eval_formatter_check(f)
    eval_formatter_slicing_check(f)

def test_dollar_formatter():
    f = text.DollarFormatter()
    eval_formatter_check(f)
    eval_formatter_slicing_check(f)
    
    ns = dict(n=12, pi=math.pi, stuff='hello there', os=os)
    s = f.format("$n", **ns)
    nt.assert_equal(s, "12")
    s = f.format("$n.real", **ns)
    nt.assert_equal(s, "12")
    s = f.format("$n/{stuff[:5]}", **ns)
    nt.assert_equal(s, "12/hello")
    s = f.format("$n $$HOME", **ns)
    nt.assert_equal(s, "12 $HOME")
    s = f.format("${foo}", foo="HOME")
    nt.assert_equal(s, "$HOME")


def test_long_substr():
    data = ['hi']
    nt.assert_equal(text.long_substr(data), 'hi')


def test_long_substr2():
    data = ['abc', 'abd', 'abf', 'ab']
    nt.assert_equal(text.long_substr(data), 'ab')

def test_long_substr_empty():
    data = []
    nt.assert_equal(text.long_substr(data), '')

def test_strip_email():
    src = """\
        >> >>> def f(x):
        >> ...   return x+1
        >> ... 
        >> >>> zz = f(2.5)"""
    cln = """\
>>> def f(x):
...   return x+1
... 
>>> zz = f(2.5)"""
    nt.assert_equal(text.strip_email_quotes(src), cln)


def test_strip_email2():
    src = '> > > list()'
    cln = 'list()'
    nt.assert_equal(text.strip_email_quotes(src), cln)

def test_LSString():
    lss = text.LSString("abc\ndef")
    nt.assert_equal(lss.l, ['abc', 'def'])
    nt.assert_equal(lss.s, 'abc def')
    lss = text.LSString(os.getcwd())
    nt.assert_is_instance(lss.p[0], Path)

def test_SList():
    sl = text.SList(['a 11', 'b 1', 'a 2'])
    nt.assert_equal(sl.n, 'a 11\nb 1\na 2')
    nt.assert_equal(sl.s, 'a 11 b 1 a 2')
    nt.assert_equal(sl.grep(lambda x: x.startswith('a')), text.SList(['a 11', 'a 2']))
    nt.assert_equal(sl.fields(0), text.SList(['a', 'b', 'a']))
    nt.assert_equal(sl.sort(field=1, nums=True), text.SList(['b 1', 'a 2', 'a 11']))

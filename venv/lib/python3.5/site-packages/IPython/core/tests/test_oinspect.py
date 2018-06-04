"""Tests for the object inspection functionality.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.


from inspect import Signature, Parameter
import os
import re
import sys

import nose.tools as nt

from .. import oinspect
from IPython.core.magic import (Magics, magics_class, line_magic,
                                cell_magic, line_cell_magic,
                                register_line_magic, register_cell_magic,
                                register_line_cell_magic)
from decorator import decorator
from IPython import get_ipython
from IPython.testing.tools import AssertPrints, AssertNotPrints
from IPython.utils.path import compress_user


#-----------------------------------------------------------------------------
# Globals and constants
#-----------------------------------------------------------------------------

inspector = oinspect.Inspector()
ip = get_ipython()

#-----------------------------------------------------------------------------
# Local utilities
#-----------------------------------------------------------------------------

# WARNING: since this test checks the line number where a function is
# defined, if any code is inserted above, the following line will need to be
# updated.  Do NOT insert any whitespace between the next line and the function
# definition below.
THIS_LINE_NUMBER = 41  # Put here the actual number of this line

from unittest import TestCase

class Test(TestCase):

    def test_find_source_lines(self):
        self.assertEqual(oinspect.find_source_lines(Test.test_find_source_lines),
                    THIS_LINE_NUMBER+6)


# A couple of utilities to ensure these tests work the same from a source or a
# binary install
def pyfile(fname):
    return os.path.normcase(re.sub('.py[co]$', '.py', fname))


def match_pyfiles(f1, f2):
    nt.assert_equal(pyfile(f1), pyfile(f2))


def test_find_file():
    match_pyfiles(oinspect.find_file(test_find_file), os.path.abspath(__file__))


def test_find_file_decorated1():

    @decorator
    def noop1(f):
        def wrapper(*a, **kw):
            return f(*a, **kw)
        return wrapper

    @noop1
    def f(x):
        "My docstring"
    
    match_pyfiles(oinspect.find_file(f), os.path.abspath(__file__))
    nt.assert_equal(f.__doc__, "My docstring")


def test_find_file_decorated2():

    @decorator
    def noop2(f, *a, **kw):
        return f(*a, **kw)

    @noop2
    @noop2
    @noop2
    def f(x):
        "My docstring 2"
    
    match_pyfiles(oinspect.find_file(f), os.path.abspath(__file__))
    nt.assert_equal(f.__doc__, "My docstring 2")
    

def test_find_file_magic():
    run = ip.find_line_magic('run')
    nt.assert_not_equal(oinspect.find_file(run), None)


# A few generic objects we can then inspect in the tests below

class Call(object):
    """This is the class docstring."""

    def __init__(self, x, y=1):
        """This is the constructor docstring."""

    def __call__(self, *a, **kw):
        """This is the call docstring."""

    def method(self, x, z=2):
        """Some method's docstring"""

class HasSignature(object):
    """This is the class docstring."""
    __signature__ = Signature([Parameter('test', Parameter.POSITIONAL_OR_KEYWORD)])

    def __init__(self, *args):
        """This is the init docstring"""


class SimpleClass(object):
    def method(self, x, z=2):
        """Some method's docstring"""


class OldStyle:
    """An old-style class for testing."""
    pass


def f(x, y=2, *a, **kw):
    """A simple function."""


def g(y, z=3, *a, **kw):
    pass  # no docstring


@register_line_magic
def lmagic(line):
    "A line magic"


@register_cell_magic
def cmagic(line, cell):
    "A cell magic"


@register_line_cell_magic
def lcmagic(line, cell=None):
    "A line/cell magic"


@magics_class
class SimpleMagics(Magics):
    @line_magic
    def Clmagic(self, cline):
        "A class-based line magic"
        
    @cell_magic
    def Ccmagic(self, cline, ccell):
        "A class-based cell magic"
        
    @line_cell_magic
    def Clcmagic(self, cline, ccell=None):
        "A class-based line/cell magic"


class Awkward(object):
    def __getattr__(self, name):
        raise Exception(name)

class NoBoolCall:
    """
    callable with `__bool__` raising should still be inspect-able.
    """

    def __call__(self):
        """does nothing"""
        pass

    def __bool__(self):
        """just raise NotImplemented"""
        raise NotImplementedError('Must be implemented')


class SerialLiar(object):
    """Attribute accesses always get another copy of the same class.

    unittest.mock.call does something similar, but it's not ideal for testing
    as the failure mode is to eat all your RAM. This gives up after 10k levels.
    """
    def __init__(self, max_fibbing_twig, lies_told=0):
        if lies_told > 10000:
            raise RuntimeError('Nose too long, honesty is the best policy')
        self.max_fibbing_twig = max_fibbing_twig
        self.lies_told = lies_told
        max_fibbing_twig[0] = max(max_fibbing_twig[0], lies_told)

    def __getattr__(self, item):
        return SerialLiar(self.max_fibbing_twig, self.lies_told + 1)

#-----------------------------------------------------------------------------
# Tests
#-----------------------------------------------------------------------------

def test_info():
    "Check that Inspector.info fills out various fields as expected."
    i = inspector.info(Call, oname='Call')
    nt.assert_equal(i['type_name'], 'type')
    expted_class = str(type(type))  # <class 'type'> (Python 3) or <type 'type'>
    nt.assert_equal(i['base_class'], expted_class)
    nt.assert_regex(i['string_form'], "<class 'IPython.core.tests.test_oinspect.Call'( at 0x[0-9a-f]{1,9})?>")
    fname = __file__
    if fname.endswith(".pyc"):
        fname = fname[:-1]
    # case-insensitive comparison needed on some filesystems
    # e.g. Windows:
    nt.assert_equal(i['file'].lower(), compress_user(fname).lower())
    nt.assert_equal(i['definition'], None)
    nt.assert_equal(i['docstring'], Call.__doc__)
    nt.assert_equal(i['source'], None)
    nt.assert_true(i['isclass'])
    nt.assert_equal(i['init_definition'], "Call(x, y=1)")
    nt.assert_equal(i['init_docstring'], Call.__init__.__doc__)

    i = inspector.info(Call, detail_level=1)
    nt.assert_not_equal(i['source'], None)
    nt.assert_equal(i['docstring'], None)

    c = Call(1)
    c.__doc__ = "Modified instance docstring"
    i = inspector.info(c)
    nt.assert_equal(i['type_name'], 'Call')
    nt.assert_equal(i['docstring'], "Modified instance docstring")
    nt.assert_equal(i['class_docstring'], Call.__doc__)
    nt.assert_equal(i['init_docstring'], Call.__init__.__doc__)
    nt.assert_equal(i['call_docstring'], Call.__call__.__doc__)

def test_class_signature():
    info = inspector.info(HasSignature, 'HasSignature')
    nt.assert_equal(info['init_definition'], "HasSignature(test)")
    nt.assert_equal(info['init_docstring'], HasSignature.__init__.__doc__)

def test_info_awkward():
    # Just test that this doesn't throw an error.
    inspector.info(Awkward())

def test_bool_raise():
    inspector.info(NoBoolCall())

def test_info_serialliar():
    fib_tracker = [0]
    inspector.info(SerialLiar(fib_tracker))

    # Nested attribute access should be cut off at 100 levels deep to avoid
    # infinite loops: https://github.com/ipython/ipython/issues/9122
    nt.assert_less(fib_tracker[0], 9000)

def test_calldef_none():
    # We should ignore __call__ for all of these.
    for obj in [f, SimpleClass().method, any, str.upper]:
        print(obj)
        i = inspector.info(obj)
        nt.assert_is(i['call_def'], None)

def f_kwarg(pos, *, kwonly):
    pass

def test_definition_kwonlyargs():
    i = inspector.info(f_kwarg, oname='f_kwarg')  # analysis:ignore
    nt.assert_equal(i['definition'], "f_kwarg(pos, *, kwonly)")

def test_getdoc():
    class A(object):
        """standard docstring"""
        pass

    class B(object):
        """standard docstring"""
        def getdoc(self):
            return "custom docstring"
    
    class C(object):
        """standard docstring"""
        def getdoc(self):
            return None
    
    a = A()
    b = B()
    c = C()
    
    nt.assert_equal(oinspect.getdoc(a), "standard docstring")
    nt.assert_equal(oinspect.getdoc(b), "custom docstring")
    nt.assert_equal(oinspect.getdoc(c), "standard docstring")


def test_empty_property_has_no_source():
    i = inspector.info(property(), detail_level=1)
    nt.assert_is(i['source'], None)


def test_property_sources():
    import posixpath 
    # A simple adder whose source and signature stays
    # the same across Python distributions
    def simple_add(a, b):
        "Adds two numbers"
        return a + b
    
    class A(object):
        @property
        def foo(self):
            return 'bar'

        foo = foo.setter(lambda self, v: setattr(self, 'bar', v))

        dname = property(posixpath.dirname)
        adder = property(simple_add) 

    i = inspector.info(A.foo, detail_level=1)
    nt.assert_in('def foo(self):', i['source'])
    nt.assert_in('lambda self, v:', i['source'])

    i = inspector.info(A.dname, detail_level=1)
    nt.assert_in('def dirname(p)', i['source'])
    
    i = inspector.info(A.adder, detail_level=1)
    nt.assert_in('def simple_add(a, b)', i['source'])


def test_property_docstring_is_in_info_for_detail_level_0():
    class A(object):
        @property
        def foobar(self):
            """This is `foobar` property."""
            pass

    ip.user_ns['a_obj'] = A()
    nt.assert_equal(
        'This is `foobar` property.',
        ip.object_inspect('a_obj.foobar', detail_level=0)['docstring'])

    ip.user_ns['a_cls'] = A
    nt.assert_equal(
        'This is `foobar` property.',
        ip.object_inspect('a_cls.foobar', detail_level=0)['docstring'])


def test_pdef():
    # See gh-1914
    def foo(): pass
    inspector.pdef(foo, 'foo')


def test_pinfo_nonascii():
    # See gh-1177
    from . import nonascii2
    ip.user_ns['nonascii2'] = nonascii2
    ip._inspect('pinfo', 'nonascii2', detail_level=1)


def test_pinfo_docstring_no_source():
    """Docstring should be included with detail_level=1 if there is no source"""
    with AssertPrints('Docstring:'):
        ip._inspect('pinfo', 'str.format', detail_level=0)
    with AssertPrints('Docstring:'):
        ip._inspect('pinfo', 'str.format', detail_level=1)


def test_pinfo_no_docstring_if_source():
    """Docstring should not be included with detail_level=1 if source is found"""
    def foo():
        """foo has a docstring"""

    ip.user_ns['foo'] = foo

    with AssertPrints('Docstring:'):
        ip._inspect('pinfo', 'foo', detail_level=0)
    with AssertPrints('Source:'):
        ip._inspect('pinfo', 'foo', detail_level=1)
    with AssertNotPrints('Docstring:'):
        ip._inspect('pinfo', 'foo', detail_level=1)


def test_pinfo_docstring_if_detail_and_no_source():
    """ Docstring should be displayed if source info not available """
    obj_def = '''class Foo(object):
                  """ This is a docstring for Foo """
                  def bar(self):
                      """ This is a docstring for Foo.bar """
                      pass
              ''' 
    
    ip.run_cell(obj_def)
    ip.run_cell('foo = Foo()')
    
    with AssertNotPrints("Source:"):
        with AssertPrints('Docstring:'):
            ip._inspect('pinfo', 'foo', detail_level=0)
        with AssertPrints('Docstring:'):
            ip._inspect('pinfo', 'foo', detail_level=1)
        with AssertPrints('Docstring:'):
            ip._inspect('pinfo', 'foo.bar', detail_level=0)

    with AssertNotPrints('Docstring:'):
        with AssertPrints('Source:'):
            ip._inspect('pinfo', 'foo.bar', detail_level=1)


def test_pinfo_magic():
    with AssertPrints('Docstring:'):
        ip._inspect('pinfo', 'lsmagic', detail_level=0)

    with AssertPrints('Source:'):
        ip._inspect('pinfo', 'lsmagic', detail_level=1)


def test_init_colors():
    # ensure colors are not present in signature info
    info = inspector.info(HasSignature)
    init_def = info['init_definition']
    nt.assert_not_in('[0m', init_def)


def test_builtin_init():
    info = inspector.info(list)
    init_def = info['init_definition']
    # Python < 3.4 can't get init definition from builtins,
    # but still exercise the inspection in case of error-raising bugs.
    if sys.version_info >= (3,4):
        nt.assert_is_not_none(init_def)


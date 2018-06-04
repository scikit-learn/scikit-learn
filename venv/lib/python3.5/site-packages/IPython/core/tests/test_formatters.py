"""Tests for the Formatters."""

import warnings
from math import pi

try:
    import numpy
except:
    numpy = None
import nose.tools as nt

from IPython import get_ipython
from traitlets.config import Config
from IPython.core.formatters import (
    PlainTextFormatter, HTMLFormatter, PDFFormatter, _mod_name_key,
    DisplayFormatter, JSONFormatter,
)
from IPython.utils.io import capture_output

class A(object):
    def __repr__(self):
        return 'A()'

class B(A):
    def __repr__(self):
        return 'B()'

class C:
    pass

class BadRepr(object):
    def __repr__(self):
        raise ValueError("bad repr")

class BadPretty(object):
    _repr_pretty_ = None

class GoodPretty(object):
    def _repr_pretty_(self, pp, cycle):
        pp.text('foo')

    def __repr__(self):
        return 'GoodPretty()'

def foo_printer(obj, pp, cycle):
    pp.text('foo')

def test_pretty():
    f = PlainTextFormatter()
    f.for_type(A, foo_printer)
    nt.assert_equal(f(A()), 'foo')
    nt.assert_equal(f(B()), 'B()')
    nt.assert_equal(f(GoodPretty()), 'foo')
    # Just don't raise an exception for the following:
    f(BadPretty())

    f.pprint = False
    nt.assert_equal(f(A()), 'A()')
    nt.assert_equal(f(B()), 'B()')
    nt.assert_equal(f(GoodPretty()), 'GoodPretty()')


def test_deferred():
    f = PlainTextFormatter()

def test_precision():
    """test various values for float_precision."""
    f = PlainTextFormatter()
    nt.assert_equal(f(pi), repr(pi))
    f.float_precision = 0
    if numpy:
        po = numpy.get_printoptions()
        nt.assert_equal(po['precision'], 0)
    nt.assert_equal(f(pi), '3')
    f.float_precision = 2
    if numpy:
        po = numpy.get_printoptions()
        nt.assert_equal(po['precision'], 2)
    nt.assert_equal(f(pi), '3.14')
    f.float_precision = '%g'
    if numpy:
        po = numpy.get_printoptions()
        nt.assert_equal(po['precision'], 2)
    nt.assert_equal(f(pi), '3.14159')
    f.float_precision = '%e'
    nt.assert_equal(f(pi), '3.141593e+00')
    f.float_precision = ''
    if numpy:
        po = numpy.get_printoptions()
        nt.assert_equal(po['precision'], 8)
    nt.assert_equal(f(pi), repr(pi))

def test_bad_precision():
    """test various invalid values for float_precision."""
    f = PlainTextFormatter()
    def set_fp(p):
        f.float_precision=p
    nt.assert_raises(ValueError, set_fp, '%')
    nt.assert_raises(ValueError, set_fp, '%.3f%i')
    nt.assert_raises(ValueError, set_fp, 'foo')
    nt.assert_raises(ValueError, set_fp, -1)

def test_for_type():
    f = PlainTextFormatter()
    
    # initial return, None
    nt.assert_is(f.for_type(C, foo_printer), None)
    # no func queries
    nt.assert_is(f.for_type(C), foo_printer)
    # shouldn't change anything
    nt.assert_is(f.for_type(C), foo_printer)
    # None should do the same
    nt.assert_is(f.for_type(C, None), foo_printer)
    nt.assert_is(f.for_type(C, None), foo_printer)

def test_for_type_string():
    f = PlainTextFormatter()
    
    type_str = '%s.%s' % (C.__module__, 'C')
    
    # initial return, None
    nt.assert_is(f.for_type(type_str, foo_printer), None)
    # no func queries
    nt.assert_is(f.for_type(type_str), foo_printer)
    nt.assert_in(_mod_name_key(C), f.deferred_printers)
    nt.assert_is(f.for_type(C), foo_printer)
    nt.assert_not_in(_mod_name_key(C), f.deferred_printers)
    nt.assert_in(C, f.type_printers)

def test_for_type_by_name():
    f = PlainTextFormatter()
    
    mod = C.__module__
    
    # initial return, None
    nt.assert_is(f.for_type_by_name(mod, 'C', foo_printer), None)
    # no func queries
    nt.assert_is(f.for_type_by_name(mod, 'C'), foo_printer)
    # shouldn't change anything
    nt.assert_is(f.for_type_by_name(mod, 'C'), foo_printer)
    # None should do the same
    nt.assert_is(f.for_type_by_name(mod, 'C', None), foo_printer)
    nt.assert_is(f.for_type_by_name(mod, 'C', None), foo_printer)

def test_lookup():
    f = PlainTextFormatter()
    
    f.for_type(C, foo_printer)
    nt.assert_is(f.lookup(C()), foo_printer)
    with nt.assert_raises(KeyError):
        f.lookup(A())

def test_lookup_string():
    f = PlainTextFormatter()
    type_str = '%s.%s' % (C.__module__, 'C')
    
    f.for_type(type_str, foo_printer)
    nt.assert_is(f.lookup(C()), foo_printer)
    # should move from deferred to imported dict
    nt.assert_not_in(_mod_name_key(C), f.deferred_printers)
    nt.assert_in(C, f.type_printers)

def test_lookup_by_type():
    f = PlainTextFormatter()
    f.for_type(C, foo_printer)
    nt.assert_is(f.lookup_by_type(C), foo_printer)
    with nt.assert_raises(KeyError):
        f.lookup_by_type(A)

def test_lookup_by_type_string():
    f = PlainTextFormatter()
    type_str = '%s.%s' % (C.__module__, 'C')
    f.for_type(type_str, foo_printer)
    
    # verify insertion
    nt.assert_in(_mod_name_key(C), f.deferred_printers)
    nt.assert_not_in(C, f.type_printers)
    
    nt.assert_is(f.lookup_by_type(type_str), foo_printer)
    # lookup by string doesn't cause import
    nt.assert_in(_mod_name_key(C), f.deferred_printers)
    nt.assert_not_in(C, f.type_printers)
    
    nt.assert_is(f.lookup_by_type(C), foo_printer)
    # should move from deferred to imported dict
    nt.assert_not_in(_mod_name_key(C), f.deferred_printers)
    nt.assert_in(C, f.type_printers)

def test_in_formatter():
    f = PlainTextFormatter()
    f.for_type(C, foo_printer)
    type_str = '%s.%s' % (C.__module__, 'C')
    nt.assert_in(C, f)
    nt.assert_in(type_str, f)

def test_string_in_formatter():
    f = PlainTextFormatter()
    type_str = '%s.%s' % (C.__module__, 'C')
    f.for_type(type_str, foo_printer)
    nt.assert_in(type_str, f)
    nt.assert_in(C, f)

def test_pop():
    f = PlainTextFormatter()
    f.for_type(C, foo_printer)
    nt.assert_is(f.lookup_by_type(C), foo_printer)
    nt.assert_is(f.pop(C, None), foo_printer)
    f.for_type(C, foo_printer)
    nt.assert_is(f.pop(C), foo_printer)
    with nt.assert_raises(KeyError):
        f.lookup_by_type(C)
    with nt.assert_raises(KeyError):
        f.pop(C)
    with nt.assert_raises(KeyError):
        f.pop(A)
    nt.assert_is(f.pop(A, None), None)

def test_pop_string():
    f = PlainTextFormatter()
    type_str = '%s.%s' % (C.__module__, 'C')
    
    with nt.assert_raises(KeyError):
        f.pop(type_str)
    
    f.for_type(type_str, foo_printer)
    f.pop(type_str)
    with nt.assert_raises(KeyError):
        f.lookup_by_type(C)
    with nt.assert_raises(KeyError):
        f.pop(type_str)

    f.for_type(C, foo_printer)
    nt.assert_is(f.pop(type_str, None), foo_printer)
    with nt.assert_raises(KeyError):
        f.lookup_by_type(C)
    with nt.assert_raises(KeyError):
        f.pop(type_str)
    nt.assert_is(f.pop(type_str, None), None)
    

def test_error_method():
    f = HTMLFormatter()
    class BadHTML(object):
        def _repr_html_(self):
            raise ValueError("Bad HTML")
    bad = BadHTML()
    with capture_output() as captured:
        result = f(bad)
    nt.assert_is(result, None)
    nt.assert_in("Traceback", captured.stdout)
    nt.assert_in("Bad HTML", captured.stdout)
    nt.assert_in("_repr_html_", captured.stdout)

def test_nowarn_notimplemented():
    f = HTMLFormatter()
    class HTMLNotImplemented(object):
        def _repr_html_(self):
            raise NotImplementedError
    h = HTMLNotImplemented()
    with capture_output() as captured:
        result = f(h)
    nt.assert_is(result, None)
    nt.assert_equal("", captured.stderr)
    nt.assert_equal("", captured.stdout)

def test_warn_error_for_type():
    f = HTMLFormatter()
    f.for_type(int, lambda i: name_error)
    with capture_output() as captured:
        result = f(5)
    nt.assert_is(result, None)
    nt.assert_in("Traceback", captured.stdout)
    nt.assert_in("NameError", captured.stdout)
    nt.assert_in("name_error", captured.stdout)

def test_error_pretty_method():
    f = PlainTextFormatter()
    class BadPretty(object):
        def _repr_pretty_(self):
            return "hello"
    bad = BadPretty()
    with capture_output() as captured:
        result = f(bad)
    nt.assert_is(result, None)
    nt.assert_in("Traceback", captured.stdout)
    nt.assert_in("_repr_pretty_", captured.stdout)
    nt.assert_in("given", captured.stdout)
    nt.assert_in("argument", captured.stdout)


def test_bad_repr_traceback():
    f = PlainTextFormatter()
    bad = BadRepr()
    with capture_output() as captured:
        result = f(bad)
    # catches error, returns None
    nt.assert_is(result, None)
    nt.assert_in("Traceback", captured.stdout)
    nt.assert_in("__repr__", captured.stdout)
    nt.assert_in("ValueError", captured.stdout)


class MakePDF(object):
    def _repr_pdf_(self):
        return 'PDF'

def test_pdf_formatter():
    pdf = MakePDF()
    f = PDFFormatter()
    nt.assert_equal(f(pdf), 'PDF')

def test_print_method_bound():
    f = HTMLFormatter()
    class MyHTML(object):
        def _repr_html_(self):
            return "hello"
    with capture_output() as captured:
        result = f(MyHTML)
    nt.assert_is(result, None)
    nt.assert_not_in("FormatterWarning", captured.stderr)

    with capture_output() as captured:
        result = f(MyHTML())
    nt.assert_equal(result, "hello")
    nt.assert_equal(captured.stderr, "")

def test_print_method_weird():

    class TextMagicHat(object):
        def __getattr__(self, key):
            return key

    f = HTMLFormatter()
    
    text_hat = TextMagicHat()
    nt.assert_equal(text_hat._repr_html_, '_repr_html_')
    with capture_output() as captured:
        result = f(text_hat)
    
    nt.assert_is(result, None)
    nt.assert_not_in("FormatterWarning", captured.stderr)

    class CallableMagicHat(object):
        def __getattr__(self, key):
            return lambda : key
    
    call_hat = CallableMagicHat()
    with capture_output() as captured:
        result = f(call_hat)
    
    nt.assert_equal(result, None)

    class BadReprArgs(object):
        def _repr_html_(self, extra, args):
            return "html"
    
    bad = BadReprArgs()
    with capture_output() as captured:
        result = f(bad)
    
    nt.assert_is(result, None)
    nt.assert_not_in("FormatterWarning", captured.stderr)


def test_format_config():
    """config objects don't pretend to support fancy reprs with lazy attrs"""
    f = HTMLFormatter()
    cfg = Config()
    with capture_output() as captured:
        result = f(cfg)
    nt.assert_is(result, None)
    nt.assert_equal(captured.stderr, "")

    with capture_output() as captured:
        result = f(Config)
    nt.assert_is(result, None)
    nt.assert_equal(captured.stderr, "")

def test_pretty_max_seq_length():
    f = PlainTextFormatter(max_seq_length=1)
    lis = list(range(3))
    text = f(lis)
    nt.assert_equal(text, '[0, ...]')
    f.max_seq_length = 0
    text = f(lis)
    nt.assert_equal(text, '[0, 1, 2]')
    text = f(list(range(1024)))
    lines = text.splitlines()
    nt.assert_equal(len(lines), 1024)


def test_ipython_display_formatter():
    """Objects with _ipython_display_ defined bypass other formatters"""
    f = get_ipython().display_formatter
    catcher = []
    class SelfDisplaying(object):
        def _ipython_display_(self):
            catcher.append(self)

    class NotSelfDisplaying(object):
        def __repr__(self):
            return "NotSelfDisplaying"
        
        def _ipython_display_(self):
            raise NotImplementedError
    
    save_enabled = f.ipython_display_formatter.enabled
    f.ipython_display_formatter.enabled = True
    
    yes = SelfDisplaying()
    no = NotSelfDisplaying()
    
    d, md = f.format(no)
    nt.assert_equal(d, {'text/plain': repr(no)})
    nt.assert_equal(md, {})
    nt.assert_equal(catcher, [])
    
    d, md = f.format(yes)
    nt.assert_equal(d, {})
    nt.assert_equal(md, {})
    nt.assert_equal(catcher, [yes])

    f.ipython_display_formatter.enabled = save_enabled


def test_json_as_string_deprecated():
    class JSONString(object):
        def _repr_json_(self):
            return '{}'
    
    f = JSONFormatter()
    with warnings.catch_warnings(record=True) as w:
        d = f(JSONString())
    nt.assert_equal(d, {})
    nt.assert_equal(len(w), 1)


def test_repr_mime():
    class HasReprMime(object):
        def _repr_mimebundle_(self, include=None, exclude=None):
            return {
                'application/json+test.v2': {
                    'x': 'y'
                },
                'plain/text' : '<HasReprMime>',
                'image/png' : 'i-overwrite'
            }

        def _repr_png_(self):
            return 'should-be-overwritten'
        def _repr_html_(self):
            return '<b>hi!</b>'
    
    f = get_ipython().display_formatter
    html_f = f.formatters['text/html']
    save_enabled = html_f.enabled
    html_f.enabled = True
    obj = HasReprMime()
    d, md = f.format(obj)
    html_f.enabled = save_enabled
    
    nt.assert_equal(sorted(d), ['application/json+test.v2',
                                'image/png',
                                'plain/text',
                                'text/html',
                                'text/plain'])
    nt.assert_equal(md, {})

    d, md = f.format(obj, include={'image/png'})
    nt.assert_equal(list(d.keys()), ['image/png'],
                    'Include should filter out even things from repr_mimebundle')
    nt.assert_equal(d['image/png'], 'i-overwrite', '_repr_mimebundle_ take precedence')



def test_pass_correct_include_exclude():
    class Tester(object):

        def __init__(self, include=None, exclude=None):
            self.include = include
            self.exclude = exclude

        def _repr_mimebundle_(self, include, exclude, **kwargs):
            if include and (include != self.include):
                raise ValueError('include got modified: display() may be broken.')
            if exclude and (exclude != self.exclude):
                raise ValueError('exclude got modified: display() may be broken.')

            return None

    include = {'a', 'b', 'c'}
    exclude = {'c', 'e' , 'f'}

    f = get_ipython().display_formatter
    f.format(Tester(include=include, exclude=exclude), include=include, exclude=exclude)
    f.format(Tester(exclude=exclude), exclude=exclude)
    f.format(Tester(include=include), include=include)


def test_repr_mime_meta():
    class HasReprMimeMeta(object):
        def _repr_mimebundle_(self, include=None, exclude=None):
            data = {
                'image/png': 'base64-image-data',
            }
            metadata = {
                'image/png': {
                    'width': 5,
                    'height': 10,
                }
            }
            return (data, metadata)
    
    f = get_ipython().display_formatter
    obj = HasReprMimeMeta()
    d, md = f.format(obj)
    nt.assert_equal(sorted(d), ['image/png', 'text/plain'])
    nt.assert_equal(md, {
        'image/png': {
            'width': 5,
            'height': 10,
        }
    })

def test_repr_mime_failure():
    class BadReprMime(object):
        def _repr_mimebundle_(self, include=None, exclude=None):
            raise RuntimeError

    f = get_ipython().display_formatter
    obj = BadReprMime()
    d, md = f.format(obj)
    nt.assert_in('text/plain', d)

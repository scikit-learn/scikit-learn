"""
exception classes and constants handling test outcomes
as well as functions creating them
"""
from __future__ import absolute_import, division, print_function
import py
import sys


class OutcomeException(BaseException):
    """ OutcomeException and its subclass instances indicate and
        contain info about test and collection outcomes.
    """
    def __init__(self, msg=None, pytrace=True):
        BaseException.__init__(self, msg)
        self.msg = msg
        self.pytrace = pytrace

    def __repr__(self):
        if self.msg:
            val = self.msg
            if isinstance(val, bytes):
                val = py._builtin._totext(val, errors='replace')
            return val
        return "<%s instance>" % (self.__class__.__name__,)
    __str__ = __repr__


TEST_OUTCOME = (OutcomeException, Exception)


class Skipped(OutcomeException):
    # XXX hackish: on 3k we fake to live in the builtins
    # in order to have Skipped exception printing shorter/nicer
    __module__ = 'builtins'

    def __init__(self, msg=None, pytrace=True, allow_module_level=False):
        OutcomeException.__init__(self, msg=msg, pytrace=pytrace)
        self.allow_module_level = allow_module_level


class Failed(OutcomeException):
    """ raised from an explicit call to pytest.fail() """
    __module__ = 'builtins'


class Exit(KeyboardInterrupt):
    """ raised for immediate program exits (no tracebacks/summaries)"""
    def __init__(self, msg="unknown reason"):
        self.msg = msg
        KeyboardInterrupt.__init__(self, msg)

# exposed helper methods


def exit(msg):
    """ exit testing process as if KeyboardInterrupt was triggered. """
    __tracebackhide__ = True
    raise Exit(msg)


exit.Exception = Exit


def skip(msg="", **kwargs):
    """ skip an executing test with the given message.  Note: it's usually
    better to use the pytest.mark.skipif marker to declare a test to be
    skipped under certain conditions like mismatching platforms or
    dependencies.  See the pytest_skipping plugin for details.

    :kwarg bool allow_module_level: allows this function to be called at
        module level, skipping the rest of the module. Default to False.
    """
    __tracebackhide__ = True
    allow_module_level = kwargs.pop('allow_module_level', False)
    if kwargs:
        keys = [k for k in kwargs.keys()]
        raise TypeError('unexpected keyword arguments: {0}'.format(keys))
    raise Skipped(msg=msg, allow_module_level=allow_module_level)


skip.Exception = Skipped


def fail(msg="", pytrace=True):
    """ explicitly fail an currently-executing test with the given Message.

    :arg pytrace: if false the msg represents the full failure information
                  and no python traceback will be reported.
    """
    __tracebackhide__ = True
    raise Failed(msg=msg, pytrace=pytrace)


fail.Exception = Failed


class XFailed(fail.Exception):
    """ raised from an explicit call to pytest.xfail() """


def xfail(reason=""):
    """ xfail an executing test or setup functions with the given reason."""
    __tracebackhide__ = True
    raise XFailed(reason)


xfail.Exception = XFailed


def importorskip(modname, minversion=None):
    """ return imported module if it has at least "minversion" as its
    __version__ attribute.  If no minversion is specified the a skip
    is only triggered if the module can not be imported.
    """
    import warnings
    __tracebackhide__ = True
    compile(modname, '', 'eval')  # to catch syntaxerrors
    should_skip = False

    with warnings.catch_warnings():
        # make sure to ignore ImportWarnings that might happen because
        # of existing directories with the same name we're trying to
        # import but without a __init__.py file
        warnings.simplefilter('ignore')
        try:
            __import__(modname)
        except ImportError:
            # Do not raise chained exception here(#1485)
            should_skip = True
    if should_skip:
        raise Skipped("could not import %r" % (modname,), allow_module_level=True)
    mod = sys.modules[modname]
    if minversion is None:
        return mod
    verattr = getattr(mod, '__version__', None)
    if minversion is not None:
        try:
            from pkg_resources import parse_version as pv
        except ImportError:
            raise Skipped("we have a required version for %r but can not import "
                          "pkg_resources to parse version strings." % (modname,),
                          allow_module_level=True)
        if verattr is None or pv(verattr) < pv(minversion):
            raise Skipped("module %r has __version__ %r, required is: %r" % (
                          modname, verattr, minversion), allow_module_level=True)
    return mod

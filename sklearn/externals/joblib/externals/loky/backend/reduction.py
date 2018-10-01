###############################################################################
# Customizable Pickler with some basic reducers
#
# author: Thomas Moreau
#
# adapted from multiprocessing/reduction.py (17/02/2017)
#  * Replace the ForkingPickler with a similar _LokyPickler,
#  * Add CustomizableLokyPickler to allow customizing pickling process
#    on the fly.
#
import io
import sys
import functools
import warnings
from multiprocessing import util
try:
    # Python 2 compat
    from cPickle import loads
except ImportError:
    from pickle import loads
    import copyreg

if sys.platform == "win32":
    if sys.version_info[:2] > (3, 3):
        from multiprocessing.reduction import duplicate
    else:
        from multiprocessing.forking import duplicate

from pickle import HIGHEST_PROTOCOL
from . import LOKY_PICKLER

Pickler = None
try:
    if LOKY_PICKLER is None or LOKY_PICKLER == "":
        from pickle import Pickler
    elif LOKY_PICKLER == "cloudpickle":
        from cloudpickle import CloudPickler as Pickler
    elif LOKY_PICKLER == "dill":
        from dill import Pickler
    elif LOKY_PICKLER != "pickle":
        from importlib import import_module
        mpickle = import_module(LOKY_PICKLER)
        Pickler = mpickle.Pickler
    util.debug("Using default backend {} for pickling."
               .format(LOKY_PICKLER if LOKY_PICKLER is not None
                       else "pickle"))
except ImportError:
    warnings.warn("Failed to import {} as asked in LOKY_PICKLER. Make sure"
                  " it is correctly installed on your system. Falling back"
                  " to default builtin pickle.".format(LOKY_PICKLER))
except AttributeError:  # pragma: no cover
    warnings.warn("Failed to find Pickler object in module {}. The module "
                  "specified in LOKY_PICKLER should implement a Pickler "
                  "object. Falling back to default builtin pickle."
                  .format(LOKY_PICKLER))


if Pickler is None:
    from pickle import Pickler


###############################################################################
# Enable custom pickling in Loky.
# To allow instance customization of the pickling process, we use 2 classes.
# _LokyPickler gives module level customization and CustomizablePickler permits
# to use instance base custom reducers.  Only CustomizablePickler should be
# used.

class _LokyPickler(Pickler):
    """Pickler that uses custom reducers.

    HIGHEST_PROTOCOL is selected by default as this pickler is used
    to pickle ephemeral datastructures for interprocess communication
    hence no backward compatibility is required.

    """

    # We override the pure Python pickler as its the only way to be able to
    # customize the dispatch table without side effects in Python 2.6
    # to 3.2. For Python 3.3+ leverage the new dispatch_table
    # feature from http://bugs.python.org/issue14166 that makes it possible
    # to use the C implementation of the Pickler which is faster.

    if hasattr(Pickler, 'dispatch'):
        # Make the dispatch registry an instance level attribute instead of
        # a reference to the class dictionary under Python 2
        dispatch = Pickler.dispatch.copy()
    else:
        # Under Python 3 initialize the dispatch table with a copy of the
        # default registry
        dispatch_table = copyreg.dispatch_table.copy()

    @classmethod
    def register(cls, type, reduce_func):
        """Attach a reducer function to a given type in the dispatch table."""
        if hasattr(Pickler, 'dispatch'):
            # Python 2 pickler dispatching is not explicitly customizable.
            # Let us use a closure to workaround this limitation.
            def dispatcher(cls, obj):
                reduced = reduce_func(obj)
                cls.save_reduce(obj=obj, *reduced)
            cls.dispatch[type] = dispatcher
        else:
            cls.dispatch_table[type] = reduce_func


class CustomizableLokyPickler(Pickler):
    def __init__(self, writer, reducers=None, protocol=HIGHEST_PROTOCOL):
        Pickler.__init__(self, writer, protocol=protocol)
        if reducers is None:
            reducers = {}
        if hasattr(Pickler, 'dispatch'):
            # Make the dispatch registry an instance level attribute instead of
            # a reference to the class dictionary under Python 2
            self.dispatch = _LokyPickler.dispatch.copy()
        else:
            # Under Python 3 initialize the dispatch table with a copy of the
            # default registry
            self.dispatch_table = _LokyPickler.dispatch_table.copy()
        for type, reduce_func in reducers.items():
            self.register(type, reduce_func)

    def register(self, type, reduce_func):
        """Attach a reducer function to a given type in the dispatch table."""
        if hasattr(Pickler, 'dispatch'):
            # Python 2 pickler dispatching is not explicitly customizable.
            # Let us use a closure to workaround this limitation.
            def dispatcher(self, obj):
                reduced = reduce_func(obj)
                self.save_reduce(obj=obj, *reduced)
            self.dispatch[type] = dispatcher
        else:
            self.dispatch_table[type] = reduce_func

    @classmethod
    def loads(self, buf):
        if sys.version_info < (3, 3) and isinstance(buf, io.BytesIO):
            buf = buf.getvalue()
        return loads(buf)

    @classmethod
    def dumps(cls, obj, reducers=None, protocol=None):
        buf = io.BytesIO()
        p = cls(buf, reducers=reducers, protocol=protocol)
        p.dump(obj)
        if sys.version_info < (3, 3):
            return buf.getvalue()
        return buf.getbuffer()


def dump(obj, file, reducers=None, protocol=None):
    '''Replacement for pickle.dump() using LokyPickler.'''
    CustomizableLokyPickler(file, reducers=reducers,
                            protocol=protocol).dump(obj)


###############################################################################
# Registers extra pickling routines to improve picklization  for loky

register = _LokyPickler.register


# make methods picklable
def _reduce_method(m):
    if m.__self__ is None:
        return getattr, (m.__class__, m.__func__.__name__)
    else:
        return getattr, (m.__self__, m.__func__.__name__)


class _C:
    def f(self):
        pass

    @classmethod
    def h(cls):
        pass


register(type(_C().f), _reduce_method)
register(type(_C.h), _reduce_method)


if not hasattr(sys, "pypy_version_info"):
    # PyPy uses functions instead of method_descriptors and wrapper_descriptors
    def _reduce_method_descriptor(m):
        return getattr, (m.__objclass__, m.__name__)

    register(type(list.append), _reduce_method_descriptor)
    register(type(int.__add__), _reduce_method_descriptor)


# Make partial func pickable
def _reduce_partial(p):
    return _rebuild_partial, (p.func, p.args, p.keywords or {})


def _rebuild_partial(func, args, keywords):
    return functools.partial(func, *args, **keywords)


register(functools.partial, _reduce_partial)

if sys.platform != "win32":
    from ._posix_reduction import _mk_inheritable  # noqa: F401
else:
    from . import _win_reduction  # noqa: F401

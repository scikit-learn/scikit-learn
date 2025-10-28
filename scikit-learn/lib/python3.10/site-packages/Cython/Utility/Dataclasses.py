################### Dataclasses_fallback ###############################

# This is the fallback dataclass code if the stdlib module isn't available.
# It defines enough of the support types to be used with cdef classes
# and to fail if used on regular types.

# (Intended to be included as py code - not compiled)

from collections import namedtuple
try:
    from types import MappingProxyType
except ImportError:
    # mutable fallback if unavailable
    MappingProxyType = lambda x: x

class _MISSING_TYPE:
    pass
MISSING = _MISSING_TYPE()

_DataclassParams = namedtuple('_DataclassParams',
    ["init", "repr", "eq", "order", "unsafe_hash", "frozen",
     "match_args", "kw_only", "slots", "weakref_slot"])
class Field:
    __slots__ = ('name',
                 'type',
                 'default',
                 'default_factory',
                 'repr',
                 'hash',
                 'init',
                 'compare',
                 'metadata',
                 'kw_only',
                 '_field_type',  # Private: not to be used by user code.
                 )

    def __init__(self, default, default_factory, init, repr, hash, compare,
                 metadata, kw_only):
        self.name = None
        self.type = None
        self.default = default
        self.default_factory = default_factory
        self.init = init
        self.repr = repr
        self.hash = hash
        self.compare = compare
        # Be aware that if MappingProxyType is unavailable (i.e. py2?) then we
        # don't enforce non-mutability that the real module does
        self.metadata = (MappingProxyType({})
                         if metadata is None else
                         MappingProxyType(metadata))
        self.kw_only = kw_only
        self._field_type = None

    def __repr__(self):
        return ('Field('
                'name={!r},'
                'type={!r},'
                'default={!r},'
                'default_factory={!r},'
                'init={!r},'
                'repr={!r},'
                'hash={!r},'
                'compare={!r},'
                'metadata={!r},'
                'kwonly={!r},'
                ')'.format(self.name, self.type, self.default,
                           self.default_factory, self.init,
                           self.repr, self.hash, self.compare,
                           self.metadata, self.kw_only))

# A sentinel object for default values to signal that a default
# factory will be used.  This is given a nice repr() which will appear
# in the function signature of dataclasses' constructors.
class _HAS_DEFAULT_FACTORY_CLASS:
    def __repr__(self):
        return '<factory>'
_HAS_DEFAULT_FACTORY = _HAS_DEFAULT_FACTORY_CLASS()

def dataclass(*args, **kwds):
    raise NotImplementedError("Standard library 'dataclasses' module"
        "is unavailable, likely due to the version of Python you're using.")

# Markers for the various kinds of fields and pseudo-fields.
class _FIELD_BASE:
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name
_FIELD = _FIELD_BASE('_FIELD')
_FIELD_CLASSVAR = _FIELD_BASE('_FIELD_CLASSVAR')
_FIELD_INITVAR = _FIELD_BASE('_FIELD_INITVAR')

def field(*ignore, **kwds):
    default = kwds.pop("default", MISSING)
    default_factory = kwds.pop("default_factory", MISSING)
    init = kwds.pop("init", True)
    repr = kwds.pop("repr", True)
    hash = kwds.pop("hash", None)
    compare = kwds.pop("compare", True)
    metadata = kwds.pop("metadata", None)
    kw_only = kwds.pop("kw_only", None)

    if kwds:
        raise ValueError("field received unexpected keyword arguments: %s"
                         % list(kwds.keys()))
    if default is not MISSING and default_factory is not MISSING:
        raise ValueError('cannot specify both default and default_factory')
    if ignore:
        raise ValueError("'field' does not take any positional arguments")
    return Field(default, default_factory, init,
                 repr, hash, compare, metadata, kw_only)

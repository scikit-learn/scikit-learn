# encoding: utf-8
"""
A lightweight Traits like module.

This is designed to provide a lightweight, simple, pure Python version of
many of the capabilities of enthought.traits.  This includes:

* Validation
* Type specification with defaults
* Static and dynamic notification
* Basic predefined types
* An API that is similar to enthought.traits

We don't support:

* Delegation
* Automatic GUI generation
* A full set of trait types.  Most importantly, we don't provide container
  traits (list, dict, tuple) that can trigger notifications if their
  contents change.
* API compatibility with enthought.traits

There are also some important difference in our design:

* enthought.traits does not validate default values.  We do.

We choose to create this module because we need these capabilities, but
we need them to be pure Python so they work in all Python implementations,
including Jython and IronPython.

Inheritance diagram:

.. inheritance-diagram:: traitlets.traitlets
   :parts: 3
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
#
# Adapted from enthought.traits, Copyright (c) Enthought, Inc.,
# also under the terms of the Modified BSD License.

import contextlib
import inspect
import os
import re
import sys
import types
import enum
try:
    from types import ClassType, InstanceType
    ClassTypes = (ClassType, type)
except:
    ClassTypes = (type,)
from warnings import warn, warn_explicit

import six

from .utils.getargspec import getargspec
from .utils.importstring import import_item
from .utils.sentinel import Sentinel
from .utils.bunch import Bunch

SequenceTypes = (list, tuple, set, frozenset)

#-----------------------------------------------------------------------------
# Basic classes
#-----------------------------------------------------------------------------


Undefined = Sentinel('Undefined', 'traitlets',
'''
Used in Traitlets to specify that no defaults are set in kwargs
'''
)

All = Sentinel('All', 'traitlets',
'''
Used in Traitlets to listen to all types of notification or to notifications
from all trait attributes.
'''
)

# Deprecated alias
NoDefaultSpecified = Undefined

class TraitError(Exception):
    pass

#-----------------------------------------------------------------------------
# Utilities
#-----------------------------------------------------------------------------

from ipython_genutils.py3compat import cast_unicode_py2

_name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")

def isidentifier(s):
    if six.PY2:
        return bool(_name_re.match(s))
    else:
        return s.isidentifier()

_deprecations_shown = set()
def _should_warn(key):
    """Add our own checks for too many deprecation warnings.

    Limit to once per package.
    """
    env_flag = os.environ.get('TRAITLETS_ALL_DEPRECATIONS')
    if env_flag and env_flag != '0':
        return True

    if key not in _deprecations_shown:
        _deprecations_shown.add(key)
        return True
    else:
        return False

def _deprecated_method(method, cls, method_name, msg):
    """Show deprecation warning about a magic method definition.

    Uses warn_explicit to bind warning to method definition instead of triggering code,
    which isn't relevant.
    """
    warn_msg = "{classname}.{method_name} is deprecated in traitlets 4.1: {msg}".format(
        classname=cls.__name__, method_name=method_name, msg=msg
    )

    for parent in inspect.getmro(cls):
        if method_name in parent.__dict__:
            cls = parent
            break
    # limit deprecation messages to once per package
    package_name = cls.__module__.split('.', 1)[0]
    key = (package_name, msg)
    if not _should_warn(key):
        return
    try:
        fname = inspect.getsourcefile(method) or "<unknown>"
        lineno = inspect.getsourcelines(method)[1] or 0
    except (IOError, TypeError) as e:
        # Failed to inspect for some reason
        warn(warn_msg + ('\n(inspection failed) %s' % e), DeprecationWarning)
    else:
        warn_explicit(warn_msg, DeprecationWarning, fname, lineno)

def class_of(object):
    """ Returns a string containing the class name of an object with the
    correct indefinite article ('a' or 'an') preceding it (e.g., 'an Image',
    'a PlotValue').
    """
    if isinstance( object, six.string_types ):
        return add_article( object )

    return add_article( object.__class__.__name__ )


def add_article(name):
    """ Returns a string containing the correct indefinite article ('a' or 'an')
    prefixed to the specified string.
    """
    if name[:1].lower() in 'aeiou':
       return 'an ' + name

    return 'a ' + name


def repr_type(obj):
    """ Return a string representation of a value and its type for readable
    error messages.
    """
    the_type = type(obj)
    if six.PY2 and the_type is InstanceType:
        # Old-style class.
        the_type = obj.__class__
    msg = '%r %r' % (obj, the_type)
    return msg


def is_trait(t):
    """ Returns whether the given value is an instance or subclass of TraitType.
    """
    return (isinstance(t, TraitType) or
            (isinstance(t, type) and issubclass(t, TraitType)))


def parse_notifier_name(names):
    """Convert the name argument to a list of names.

    Examples
    --------

    >>> parse_notifier_name([])
    [All]
    >>> parse_notifier_name('a')
    ['a']
    >>> parse_notifier_name(['a', 'b'])
    ['a', 'b']
    >>> parse_notifier_name(All)
    [All]
    """
    if names is All or isinstance(names, six.string_types):
        return [names]
    else:
        if not names or All in names:
            return [All]
        for n in names:
            if not isinstance(n, six.string_types):
                raise TypeError("names must be strings, not %r" % n)
        return names


class _SimpleTest:
    def __init__ ( self, value ): self.value = value
    def __call__ ( self, test  ):
        return test == self.value
    def __repr__(self):
        return "<SimpleTest(%r)" % self.value
    def __str__(self):
        return self.__repr__()


def getmembers(object, predicate=None):
    """A safe version of inspect.getmembers that handles missing attributes.

    This is useful when there are descriptor based attributes that for
    some reason raise AttributeError even though they exist.  This happens
    in zope.inteface with the __provides__ attribute.
    """
    results = []
    for key in dir(object):
        try:
            value = getattr(object, key)
        except AttributeError:
            pass
        else:
            if not predicate or predicate(value):
                results.append((key, value))
    results.sort()
    return results

def _validate_link(*tuples):
    """Validate arguments for traitlet link functions"""
    for t in tuples:
        if not len(t) == 2:
            raise TypeError("Each linked traitlet must be specified as (HasTraits, 'trait_name'), not %r" % t)
        obj, trait_name = t
        if not isinstance(obj, HasTraits):
            raise TypeError("Each object must be HasTraits, not %r" % type(obj))
        if not trait_name in obj.traits():
            raise TypeError("%r has no trait %r" % (obj, trait_name))

class link(object):
    """Link traits from different objects together so they remain in sync.

    Parameters
    ----------
    source : (object / attribute name) pair
    target : (object / attribute name) pair

    Examples
    --------

    >>> c = link((src, 'value'), (tgt, 'value'))
    >>> src.value = 5  # updates other objects as well
    """
    updating = False

    def __init__(self, source, target):
        _validate_link(source, target)
        self.source, self.target = source, target
        try:
            setattr(target[0], target[1], getattr(source[0], source[1]))
        finally:
            source[0].observe(self._update_target, names=source[1])
            target[0].observe(self._update_source, names=target[1])

    @contextlib.contextmanager
    def _busy_updating(self):
        self.updating = True
        try:
            yield
        finally:
            self.updating = False

    def _update_target(self, change):
        if self.updating:
            return
        with self._busy_updating():
            setattr(self.target[0], self.target[1], change.new)

    def _update_source(self, change):
        if self.updating:
            return
        with self._busy_updating():
            setattr(self.source[0], self.source[1], change.new)

    def unlink(self):
        self.source[0].unobserve(self._update_target, names=self.source[1])
        self.target[0].unobserve(self._update_source, names=self.target[1])
        self.source, self.target = None, None


class directional_link(object):
    """Link the trait of a source object with traits of target objects.

    Parameters
    ----------
    source : (object, attribute name) pair
    target : (object, attribute name) pair
    transform: callable (optional)
        Data transformation between source and target.

    Examples
    --------

    >>> c = directional_link((src, 'value'), (tgt, 'value'))
    >>> src.value = 5  # updates target objects
    >>> tgt.value = 6  # does not update source object
    """
    updating = False

    def __init__(self, source, target, transform=None):
        self._transform = transform if transform else lambda x: x
        _validate_link(source, target)
        self.source, self.target = source, target
        try:
            setattr(target[0], target[1],
                    self._transform(getattr(source[0], source[1])))
        finally:
            self.source[0].observe(self._update, names=self.source[1])

    @contextlib.contextmanager
    def _busy_updating(self):
        self.updating = True
        try:
            yield
        finally:
            self.updating = False

    def _update(self, change):
        if self.updating:
            return
        with self._busy_updating():
            setattr(self.target[0], self.target[1],
                    self._transform(change.new))

    def unlink(self):
        self.source[0].unobserve(self._update, names=self.source[1])
        self.source, self.target = None, None

dlink = directional_link


#-----------------------------------------------------------------------------
# Base Descriptor Class
#-----------------------------------------------------------------------------


class BaseDescriptor(object):
    """Base descriptor class

    Notes
    -----
    This implements Python's descriptor prototol.

    This class is the base class for all such descriptors.  The
    only magic we use is a custom metaclass for the main :class:`HasTraits`
    class that does the following:

    1. Sets the :attr:`name` attribute of every :class:`BaseDescriptor`
       instance in the class dict to the name of the attribute.
    2. Sets the :attr:`this_class` attribute of every :class:`BaseDescriptor`
       instance in the class dict to the *class* that declared the trait.
       This is used by the :class:`This` trait to allow subclasses to
       accept superclasses for :class:`This` values.
    """

    name = None
    this_class = None

    def class_init(self, cls, name):
        """Part of the initialization which may depend on the underlying
        HasDescriptors class.

        It is typically overloaded for specific types.

        This method is called by :meth:`MetaHasDescriptors.__init__`
        passing the class (`cls`) and `name` under which the descriptor
        has been assigned.
        """
        self.this_class = cls
        self.name = name

    def instance_init(self, obj):
        """Part of the initialization which may depend on the underlying
        HasDescriptors instance.

        It is typically overloaded for specific types.

        This method is called by :meth:`HasTraits.__new__` and in the
        :meth:`BaseDescriptor.instance_init` method of descriptors holding
        other descriptors.
        """
        pass


class TraitType(BaseDescriptor):
    """A base class for all trait types.
    """

    metadata = {}
    default_value = Undefined
    allow_none = False
    read_only = False
    info_text = 'any value'

    def __init__(self, default_value=Undefined, allow_none=False, read_only=None, help=None,
        config=None, **kwargs):
        """Declare a traitlet.

        If *allow_none* is True, None is a valid value in addition to any
        values that are normally valid. The default is up to the subclass.
        For most trait types, the default value for ``allow_none`` is False.

        Extra metadata can be associated with the traitlet using the .tag() convenience method
        or by using the traitlet instance's .metadata dictionary.
        """
        if default_value is not Undefined:
            self.default_value = default_value
        if allow_none:
            self.allow_none = allow_none
        if read_only is not None:
            self.read_only = read_only
        self.help = help if help is not None else ''

        if len(kwargs) > 0:
            stacklevel = 1
            f = inspect.currentframe()
            # count supers to determine stacklevel for warning
            while f.f_code.co_name == '__init__':
                stacklevel += 1
                f = f.f_back
            mod = f.f_globals.get('__name__') or ''
            pkg = mod.split('.', 1)[0]
            key = tuple(['metadata-tag', pkg] + sorted(kwargs))
            if _should_warn(key):
                warn("metadata %s was set from the constructor. "
                     "With traitlets 4.1, metadata should be set using the .tag() method, "
                     "e.g., Int().tag(key1='value1', key2='value2')" % (kwargs,),
                     DeprecationWarning, stacklevel=stacklevel)
            if len(self.metadata) > 0:
                self.metadata = self.metadata.copy()
                self.metadata.update(kwargs)
            else:
                self.metadata = kwargs
        else:
            self.metadata = self.metadata.copy()
        if config is not None:
            self.metadata['config'] = config

        # We add help to the metadata during a deprecation period so that
        # code that looks for the help string there can find it.
        if help is not None:
            self.metadata['help'] = help

    def get_default_value(self):
        """DEPRECATED: Retrieve the static default value for this trait.

        Use self.default_value instead
        """
        warn("get_default_value is deprecated in traitlets 4.0: use the .default_value attribute", DeprecationWarning,
             stacklevel=2)
        return self.default_value

    def init_default_value(self, obj):
        """DEPRECATED: Set the static default value for the trait type.
        """
        warn("init_default_value is deprecated in traitlets 4.0, and may be removed in the future", DeprecationWarning,
             stacklevel=2)
        value = self._validate(obj, self.default_value)
        obj._trait_values[self.name] = value
        return value

    def _dynamic_default_callable(self, obj):
        """Retrieve a callable to calculate the default for this traitlet.

        This looks for:

        * default generators registered with the @default descriptor.
        * obj._{name}_default() on the class with the traitlet, or a subclass
          that obj belongs to.
        * trait.make_dynamic_default, which is defined by Instance

        If neither exist, it returns None
        """
        # Traitlets without a name are not on the instance, e.g. in List or Union
        if self.name:

            # Only look for default handlers in classes derived from self.this_class.
            mro = type(obj).mro()
            meth_name = '_%s_default' % self.name
            for cls in mro[:mro.index(self.this_class) + 1]:
                if hasattr(cls, '_trait_default_generators'):
                    default_handler = cls._trait_default_generators.get(self.name)
                    if default_handler is not None and default_handler.this_class == cls:
                        return types.MethodType(default_handler.func, obj)

                if meth_name in cls.__dict__:
                    method = getattr(obj, meth_name)
                    return method

        return getattr(self, 'make_dynamic_default', None)

    def instance_init(self, obj):
        # If no dynamic initialiser is present, and the trait implementation or
        # use provides a static default, transfer that to obj._trait_values.
        with obj.cross_validation_lock:
            if (self._dynamic_default_callable(obj) is None) \
                    and (self.default_value is not Undefined):
                v = self._validate(obj, self.default_value)
                if self.name is not None:
                    obj._trait_values[self.name] = v

    def get(self, obj, cls=None):
        try:
            value = obj._trait_values[self.name]
        except KeyError:
            # Check for a dynamic initializer.
            dynamic_default = self._dynamic_default_callable(obj)
            if dynamic_default is None:
                raise TraitError("No default value found for %s trait of %r"
                                 % (self.name, obj))
            value = self._validate(obj, dynamic_default())
            obj._trait_values[self.name] = value
            return value
        except Exception:
            # This should never be reached.
            raise TraitError('Unexpected error in TraitType: '
                             'default value not set properly')
        else:
            return value

    def __get__(self, obj, cls=None):
        """Get the value of the trait by self.name for the instance.

        Default values are instantiated when :meth:`HasTraits.__new__`
        is called.  Thus by the time this method gets called either the
        default value or a user defined value (they called :meth:`__set__`)
        is in the :class:`HasTraits` instance.
        """
        if obj is None:
            return self
        else:
            return self.get(obj, cls)

    def set(self, obj, value):
        new_value = self._validate(obj, value)
        try:
            old_value = obj._trait_values[self.name]
        except KeyError:
            old_value = self.default_value

        obj._trait_values[self.name] = new_value
        try:
            silent = bool(old_value == new_value)
        except:
            # if there is an error in comparing, default to notify
            silent = False
        if silent is not True:
            # we explicitly compare silent to True just in case the equality
            # comparison above returns something other than True/False
            obj._notify_trait(self.name, old_value, new_value)

    def __set__(self, obj, value):
        """Set the value of the trait by self.name for the instance.

        Values pass through a validation stage where errors are raised when
        impropper types, or types that cannot be coerced, are encountered.
        """
        if self.read_only:
            raise TraitError('The "%s" trait is read-only.' % self.name)
        else:
            self.set(obj, value)

    def _validate(self, obj, value):
        if value is None and self.allow_none:
            return value
        if hasattr(self, 'validate'):
            value = self.validate(obj, value)
        if obj._cross_validation_lock is False:
            value = self._cross_validate(obj, value)
        return value

    def _cross_validate(self, obj, value):
        if self.name in obj._trait_validators:
            proposal = Bunch({'trait': self, 'value': value, 'owner': obj})
            value = obj._trait_validators[self.name](obj, proposal)
        elif hasattr(obj, '_%s_validate' % self.name):
            meth_name = '_%s_validate' % self.name
            cross_validate = getattr(obj, meth_name)
            _deprecated_method(cross_validate, obj.__class__, meth_name,
                "use @validate decorator instead.")
            value = cross_validate(value, self)
        return value

    def __or__(self, other):
        if isinstance(other, Union):
            return Union([self] + other.trait_types)
        else:
            return Union([self, other])

    def info(self):
        return self.info_text

    def error(self, obj, value):
        if obj is not None:
            e = "The '%s' trait of %s instance must be %s, but a value of %s was specified." \
                % (self.name, class_of(obj),
                   self.info(), repr_type(value))
        else:
            e = "The '%s' trait must be %s, but a value of %r was specified." \
                % (self.name, self.info(), repr_type(value))
        raise TraitError(e)

    def get_metadata(self, key, default=None):
        """DEPRECATED: Get a metadata value.

        Use .metadata[key] or .metadata.get(key, default) instead.
        """
        if key == 'help':
            msg = "use the instance .help string directly, like x.help"
        else:
            msg = "use the instance .metadata dictionary directly, like x.metadata[key] or x.metadata.get(key, default)"
        warn("Deprecated in traitlets 4.1, " + msg, DeprecationWarning, stacklevel=2)
        return self.metadata.get(key, default)

    def set_metadata(self, key, value):
        """DEPRECATED: Set a metadata key/value.

        Use .metadata[key] = value instead.
        """
        if key == 'help':
            msg = "use the instance .help string directly, like x.help = value"
        else:
            msg = "use the instance .metadata dictionary directly, like x.metadata[key] = value"
        warn("Deprecated in traitlets 4.1, " + msg, DeprecationWarning, stacklevel=2)
        self.metadata[key] = value

    def tag(self, **metadata):
        """Sets metadata and returns self.

        This allows convenient metadata tagging when initializing the trait, such as:

        >>> Int(0).tag(config=True, sync=True)
        """
        maybe_constructor_keywords = set(metadata.keys()).intersection({'help','allow_none', 'read_only', 'default_value'})
        if maybe_constructor_keywords:
            warn('The following attributes are set in using `tag`, but seem to be constructor keywords arguments: %s '%
                    maybe_constructor_keywords, UserWarning, stacklevel=2)

        self.metadata.update(metadata)
        return self

    def default_value_repr(self):
        return repr(self.default_value)

#-----------------------------------------------------------------------------
# The HasTraits implementation
#-----------------------------------------------------------------------------

class _CallbackWrapper(object):
    """An object adapting a on_trait_change callback into an observe callback.

    The comparison operator __eq__ is implemented to enable removal of wrapped
    callbacks.
    """

    def __init__(self, cb):
        self.cb = cb
        # Bound methods have an additional 'self' argument.
        offset = -1 if isinstance(self.cb, types.MethodType) else 0
        self.nargs = len(getargspec(cb)[0]) + offset
        if (self.nargs > 4):
            raise TraitError('a trait changed callback must have 0-4 arguments.')

    def __eq__(self, other):
        # The wrapper is equal to the wrapped element
        if isinstance(other, _CallbackWrapper):
            return self.cb == other.cb
        else:
            return self.cb == other

    def __call__(self, change):
        # The wrapper is callable
        if self.nargs == 0:
            self.cb()
        elif self.nargs == 1:
            self.cb(change.name)
        elif self.nargs == 2:
            self.cb(change.name, change.new)
        elif self.nargs == 3:
            self.cb(change.name, change.old, change.new)
        elif self.nargs == 4:
            self.cb(change.name, change.old, change.new, change.owner)

def _callback_wrapper(cb):
    if isinstance(cb, _CallbackWrapper):
        return cb
    else:
        return _CallbackWrapper(cb)


class MetaHasDescriptors(type):
    """A metaclass for HasDescriptors.

    This metaclass makes sure that any TraitType class attributes are
    instantiated and sets their name attribute.
    """

    def __new__(mcls, name, bases, classdict):
        """Create the HasDescriptors class."""
        for k, v in classdict.items():
            # ----------------------------------------------------------------
            # Support of deprecated behavior allowing for TraitType types
            # to be used instead of TraitType instances.
            if inspect.isclass(v) and issubclass(v, TraitType):
                warn("Traits should be given as instances, not types (for example, `Int()`, not `Int`)."
                     " Passing types is deprecated in traitlets 4.1.",
                     DeprecationWarning, stacklevel=2)
                classdict[k] = v()
            # ----------------------------------------------------------------

        return super(MetaHasDescriptors, mcls).__new__(mcls, name, bases, classdict)

    def __init__(cls, name, bases, classdict):
        """Finish initializing the HasDescriptors class."""
        super(MetaHasDescriptors, cls).__init__(name, bases, classdict)
        cls.setup_class(classdict)

    def setup_class(cls, classdict):
        """Setup descriptor instance on the class

        This sets the :attr:`this_class` and :attr:`name` attributes of each
        BaseDescriptor in the class dict of the newly created ``cls`` before
        calling their :attr:`class_init` method.
        """
        for k, v in classdict.items():
            if isinstance(v, BaseDescriptor):
                v.class_init(cls, k)


class MetaHasTraits(MetaHasDescriptors):
    """A metaclass for HasTraits."""

    def setup_class(cls, classdict):
        cls._trait_default_generators = {}
        super(MetaHasTraits, cls).setup_class(classdict)


def observe(*names, **kwargs):
    """A decorator which can be used to observe Traits on a class.

    The handler passed to the decorator will be called with one ``change``
    dict argument. The change dictionary at least holds a 'type' key and a
    'name' key, corresponding respectively to the type of notification and the
    name of the attribute that triggered the notification.

    Other keys may be passed depending on the value of 'type'. In the case
    where type is 'change', we also have the following keys:
    * ``owner`` : the HasTraits instance
    * ``old`` : the old value of the modified trait attribute
    * ``new`` : the new value of the modified trait attribute
    * ``name`` : the name of the modified trait attribute.

    Parameters
    ----------
    *names
        The str names of the Traits to observe on the object.
    type: str, kwarg-only
        The type of event to observe (e.g. 'change')
    """
    if not names:
        raise TypeError("Please specify at least one trait name to observe.")
    for name in names:
        if name is not All and not isinstance(name, six.string_types):
            raise TypeError("trait names to observe must be strings or All, not %r" % name)
    return ObserveHandler(names, type=kwargs.get('type', 'change'))


def observe_compat(func):
    """Backward-compatibility shim decorator for observers

    Use with:

    @observe('name')
    @observe_compat
    def _foo_changed(self, change):
        ...

    With this, `super()._foo_changed(self, name, old, new)` in subclasses will still work.
    Allows adoption of new observer API without breaking subclasses that override and super.
    """
    def compatible_observer(self, change_or_name, old=Undefined, new=Undefined):
        if isinstance(change_or_name, dict):
            change = change_or_name
        else:
            clsname = self.__class__.__name__
            warn("A parent of %s._%s_changed has adopted the new (traitlets 4.1) @observe(change) API" % (
                clsname, change_or_name), DeprecationWarning)
            change = Bunch(
                type='change',
                old=old,
                new=new,
                name=change_or_name,
                owner=self,
            )
        return func(self, change)
    return compatible_observer


def validate(*names):
    """A decorator to register cross validator of HasTraits object's state
    when a Trait is set.

    The handler passed to the decorator must have one ``proposal`` dict argument.
    The proposal dictionary must hold the following keys:
    * ``owner`` : the HasTraits instance
    * ``value`` : the proposed value for the modified trait attribute
    * ``trait`` : the TraitType instance associated with the attribute

    Parameters
    ----------
    names
        The str names of the Traits to validate.

    Notes
    -----
    Since the owner has access to the ``HasTraits`` instance via the 'owner' key,
    the registered cross validator could potentially make changes to attributes
    of the ``HasTraits`` instance. However, we recommend not to do so. The reason
    is that the cross-validation of attributes may run in arbitrary order when
    exiting the ``hold_trait_notifications`` context, and such changes may not
    commute.
    """
    if not names:
        raise TypeError("Please specify at least one trait name to validate.")
    for name in names:
        if name is not All and not isinstance(name, six.string_types):
            raise TypeError("trait names to validate must be strings or All, not %r" % name)
    return ValidateHandler(names)


def default(name):
    """ A decorator which assigns a dynamic default for a Trait on a HasTraits object.

    Parameters
    ----------
    name
        The str name of the Trait on the object whose default should be generated.

    Notes
    -----
    Unlike observers and validators which are properties of the HasTraits
    instance, default value generators are class-level properties.

    Besides, default generators are only invoked if they are registered in
    subclasses of `this_type`.

    ::

        class A(HasTraits):
            bar = Int()

            @default('bar')
            def get_bar_default(self):
                return 11


        class B(A):
            bar = Float()  # This trait ignores the default generator defined in
                           # the base class A


        class C(B):

            @default('bar')
            def some_other_default(self):  # This default generator should not be
                return 3.0                 # ignored since it is defined in a
                                           # class derived from B.a.this_class.
    """
    if not isinstance(name, six.string_types):
        raise TypeError("Trait name must be a string or All, not %r" % name)
    return DefaultHandler(name)


class EventHandler(BaseDescriptor):

    def _init_call(self, func):
        self.func = func
        return self

    def __call__(self, *args, **kwargs):
        """Pass `*args` and `**kwargs` to the handler's function if it exists."""
        if hasattr(self, 'func'):
            return self.func(*args, **kwargs)
        else:
            return self._init_call(*args, **kwargs)

    def __get__(self, inst, cls=None):
        if inst is None:
            return self
        return types.MethodType(self.func, inst)


class ObserveHandler(EventHandler):

    def __init__(self, names, type):
        self.trait_names = names
        self.type = type

    def instance_init(self, inst):
        inst.observe(self, self.trait_names, type=self.type)


class ValidateHandler(EventHandler):

    def __init__(self, names):
        self.trait_names = names

    def instance_init(self, inst):
        inst._register_validator(self, self.trait_names)


class DefaultHandler(EventHandler):

    def __init__(self, name):
        self.trait_name = name

    def class_init(self, cls, name):
        super(DefaultHandler, self).class_init(cls, name)
        cls._trait_default_generators[self.trait_name] = self


class HasDescriptors(six.with_metaclass(MetaHasDescriptors, object)):
    """The base class for all classes that have descriptors.
    """

    def __new__(cls, *args, **kwargs):
        # This is needed because object.__new__ only accepts
        # the cls argument.
        new_meth = super(HasDescriptors, cls).__new__
        if new_meth is object.__new__:
            inst = new_meth(cls)
        else:
            inst = new_meth(cls, *args, **kwargs)
        inst.setup_instance(*args, **kwargs)
        return inst

    def setup_instance(self, *args, **kwargs):
        """
        This is called **before** self.__init__ is called.
        """
        self._cross_validation_lock = False
        cls = self.__class__
        for key in dir(cls):
            # Some descriptors raise AttributeError like zope.interface's
            # __provides__ attributes even though they exist.  This causes
            # AttributeErrors even though they are listed in dir(cls).
            try:
                value = getattr(cls, key)
            except AttributeError:
                pass
            else:
                if isinstance(value, BaseDescriptor):
                    value.instance_init(self)


class HasTraits(six.with_metaclass(MetaHasTraits, HasDescriptors)):

    def setup_instance(self, *args, **kwargs):
        self._trait_values = {}
        self._trait_notifiers = {}
        self._trait_validators = {}
        super(HasTraits, self).setup_instance(*args, **kwargs)

    def __init__(self, *args, **kwargs):
        # Allow trait values to be set using keyword arguments.
        # We need to use setattr for this to trigger validation and
        # notifications.
        super_args = args
        super_kwargs = {}
        with self.hold_trait_notifications():
            for key, value in kwargs.items():
                if self.has_trait(key):
                    setattr(self, key, value)
                else:
                    # passthrough args that don't set traits to super
                    super_kwargs[key] = value
        try:
            super(HasTraits, self).__init__(*super_args, **super_kwargs)
        except TypeError as e:
            arg_s_list = [ repr(arg) for arg in super_args ]
            for k, v in super_kwargs.items():
                arg_s_list.append("%s=%r" % (k, v))
            arg_s = ', '.join(arg_s_list)
            warn(
                "Passing unrecoginized arguments to super({classname}).__init__({arg_s}).\n"
                "{error}\n"
                "This is deprecated in traitlets 4.2."
                "This error will be raised in a future release of traitlets."
                .format(
                    arg_s=arg_s, classname=self.__class__.__name__,
                    error=e,
                ),
                DeprecationWarning,
                stacklevel=2,
            )

    def __getstate__(self):
        d = self.__dict__.copy()
        # event handlers stored on an instance are
        # expected to be reinstantiated during a
        # recall of instance_init during __setstate__
        d['_trait_notifiers'] = {}
        d['_trait_validators'] = {}
        return d

    def __setstate__(self, state):
        self.__dict__ = state.copy()

        # event handlers are reassigned to self
        cls = self.__class__
        for key in dir(cls):
            # Some descriptors raise AttributeError like zope.interface's
            # __provides__ attributes even though they exist.  This causes
            # AttributeErrors even though they are listed in dir(cls).
            try:
                value = getattr(cls, key)
            except AttributeError:
                pass
            else:
                if isinstance(value, EventHandler):
                    value.instance_init(self)

    @property
    @contextlib.contextmanager
    def cross_validation_lock(self):
        """
        A contextmanager for running a block with our cross validation lock set
        to True.

        At the end of the block, the lock's value is restored to its value
        prior to entering the block.
        """
        if self._cross_validation_lock:
            yield
            return
        else:
            try:
                self._cross_validation_lock = True
                yield
            finally:
                self._cross_validation_lock = False

    @contextlib.contextmanager
    def hold_trait_notifications(self):
        """Context manager for bundling trait change notifications and cross
        validation.

        Use this when doing multiple trait assignments (init, config), to avoid
        race conditions in trait notifiers requesting other trait values.
        All trait notifications will fire after all values have been assigned.
        """
        if self._cross_validation_lock:
            yield
            return
        else:
            cache = {}
            notify_change = self.notify_change

            def compress(past_changes, change):
                """Merges the provided change with the last if possible."""
                if past_changes is None:
                    return [change]
                else:
                    if past_changes[-1]['type'] == 'change' and change.type == 'change':
                        past_changes[-1]['new'] = change.new
                    else:
                        # In case of changes other than 'change', append the notification.
                        past_changes.append(change)
                    return past_changes

            def hold(change):
                name = change.name
                cache[name] = compress(cache.get(name), change)

            try:
                # Replace notify_change with `hold`, caching and compressing
                # notifications, disable cross validation and yield.
                self.notify_change = hold
                self._cross_validation_lock = True
                yield
                # Cross validate final values when context is released.
                for name in list(cache.keys()):
                    trait = getattr(self.__class__, name)
                    value = trait._cross_validate(self, getattr(self, name))
                    self.set_trait(name, value)
            except TraitError as e:
                # Roll back in case of TraitError during final cross validation.
                self.notify_change = lambda x: None
                for name, changes in cache.items():
                    for change in changes[::-1]:
                        # TODO: Separate in a rollback function per notification type.
                        if change.type == 'change':
                            if change.old is not Undefined:
                                self.set_trait(name, change.old)
                            else:
                                self._trait_values.pop(name)
                cache = {}
                raise e
            finally:
                self._cross_validation_lock = False
                # Restore method retrieval from class
                del self.notify_change

                # trigger delayed notifications
                for changes in cache.values():
                    for change in changes:
                        self.notify_change(change)

    def _notify_trait(self, name, old_value, new_value):
        self.notify_change(Bunch(
            name=name,
            old=old_value,
            new=new_value,
            owner=self,
            type='change',
        ))

    def notify_change(self, change):
        if not isinstance(change, Bunch):
            # cast to bunch if given a dict
            change = Bunch(change)
        name, type = change.name, change.type

        callables = []
        callables.extend(self._trait_notifiers.get(name, {}).get(type, []))
        callables.extend(self._trait_notifiers.get(name, {}).get(All, []))
        callables.extend(self._trait_notifiers.get(All, {}).get(type, []))
        callables.extend(self._trait_notifiers.get(All, {}).get(All, []))

        # Now static ones
        magic_name = '_%s_changed' % name
        if hasattr(self, magic_name):
            class_value = getattr(self.__class__, magic_name)
            if not isinstance(class_value, ObserveHandler):
                _deprecated_method(class_value, self.__class__, magic_name,
                    "use @observe and @unobserve instead.")
                cb = getattr(self, magic_name)
                # Only append the magic method if it was not manually registered
                if cb not in callables:
                    callables.append(_callback_wrapper(cb))

        # Call them all now
        # Traits catches and logs errors here.  I allow them to raise
        for c in callables:
            # Bound methods have an additional 'self' argument.

            if isinstance(c, _CallbackWrapper):
                c = c.__call__
            elif isinstance(c, EventHandler) and c.name is not None:
                c = getattr(self, c.name)

            c(change)

    def _add_notifiers(self, handler, name, type):
        if name not in self._trait_notifiers:
            nlist = []
            self._trait_notifiers[name] = {type: nlist}
        else:
            if type not in self._trait_notifiers[name]:
                nlist = []
                self._trait_notifiers[name][type] = nlist
            else:
                nlist = self._trait_notifiers[name][type]
        if handler not in nlist:
            nlist.append(handler)

    def _remove_notifiers(self, handler, name, type):
        try:
            if handler is None:
                del self._trait_notifiers[name][type]
            else:
                self._trait_notifiers[name][type].remove(handler)
        except KeyError:
            pass

    def on_trait_change(self, handler=None, name=None, remove=False):
        """DEPRECATED: Setup a handler to be called when a trait changes.

        This is used to setup dynamic notifications of trait changes.

        Static handlers can be created by creating methods on a HasTraits
        subclass with the naming convention '_[traitname]_changed'.  Thus,
        to create static handler for the trait 'a', create the method
        _a_changed(self, name, old, new) (fewer arguments can be used, see
        below).

        If `remove` is True and `handler` is not specified, all change
        handlers for the specified name are uninstalled.

        Parameters
        ----------
        handler : callable, None
            A callable that is called when a trait changes.  Its
            signature can be handler(), handler(name), handler(name, new),
            handler(name, old, new), or handler(name, old, new, self).
        name : list, str, None
            If None, the handler will apply to all traits.  If a list
            of str, handler will apply to all names in the list.  If a
            str, the handler will apply just to that name.
        remove : bool
            If False (the default), then install the handler.  If True
            then unintall it.
        """
        warn("on_trait_change is deprecated in traitlets 4.1: use observe instead",
             DeprecationWarning, stacklevel=2)
        if name is None:
            name = All
        if remove:
            self.unobserve(_callback_wrapper(handler), names=name)
        else:
            self.observe(_callback_wrapper(handler), names=name)

    def observe(self, handler, names=All, type='change'):
        """Setup a handler to be called when a trait changes.

        This is used to setup dynamic notifications of trait changes.

        Parameters
        ----------
        handler : callable
            A callable that is called when a trait changes. Its
            signature should be ``handler(change)``, where ``change`` is a
            dictionary. The change dictionary at least holds a 'type' key.
            * ``type``: the type of notification.
            Other keys may be passed depending on the value of 'type'. In the
            case where type is 'change', we also have the following keys:
            * ``owner`` : the HasTraits instance
            * ``old`` : the old value of the modified trait attribute
            * ``new`` : the new value of the modified trait attribute
            * ``name`` : the name of the modified trait attribute.
        names : list, str, All
            If names is All, the handler will apply to all traits.  If a list
            of str, handler will apply to all names in the list.  If a
            str, the handler will apply just to that name.
        type : str, All (default: 'change')
            The type of notification to filter by. If equal to All, then all
            notifications are passed to the observe handler.
        """
        names = parse_notifier_name(names)
        for n in names:
            self._add_notifiers(handler, n, type)

    def unobserve(self, handler, names=All, type='change'):
        """Remove a trait change handler.

        This is used to unregister handlers to trait change notifications.

        Parameters
        ----------
        handler : callable
            The callable called when a trait attribute changes.
        names : list, str, All (default: All)
            The names of the traits for which the specified handler should be
            uninstalled. If names is All, the specified handler is uninstalled
            from the list of notifiers corresponding to all changes.
        type : str or All (default: 'change')
            The type of notification to filter by. If All, the specified handler
            is uninstalled from the list of notifiers corresponding to all types.
        """
        names = parse_notifier_name(names)
        for n in names:
            self._remove_notifiers(handler, n, type)

    def unobserve_all(self, name=All):
        """Remove trait change handlers of any type for the specified name.
        If name is not specified, removes all trait notifiers."""
        if name is All:
            self._trait_notifiers = {}
        else:
            try:
                del self._trait_notifiers[name]
            except KeyError:
                pass

    def _register_validator(self, handler, names):
        """Setup a handler to be called when a trait should be cross validated.

        This is used to setup dynamic notifications for cross-validation.

        If a validator is already registered for any of the provided names, a
        TraitError is raised and no new validator is registered.

        Parameters
        ----------
        handler : callable
            A callable that is called when the given trait is cross-validated.
            Its signature is handler(proposal), where proposal is a Bunch (dictionary with attribute access)
            with the following attributes/keys:
                * ``owner`` : the HasTraits instance
                * ``value`` : the proposed value for the modified trait attribute
                * ``trait`` : the TraitType instance associated with the attribute
        names : List of strings
            The names of the traits that should be cross-validated
        """
        for name in names:
            magic_name = '_%s_validate' % name
            if hasattr(self, magic_name):
                class_value = getattr(self.__class__, magic_name)
                if not isinstance(class_value, ValidateHandler):
                    _deprecated_method(class_value, self.__class, magic_name,
                        "use @validate decorator instead.")
        for name in names:
            self._trait_validators[name] = handler

    def add_traits(self, **traits):
        """Dynamically add trait attributes to the HasTraits instance."""
        self.__class__ = type(self.__class__.__name__, (self.__class__,),
                              traits)
        for trait in traits.values():
            trait.instance_init(self)

    def set_trait(self, name, value):
        """Forcibly sets trait attribute, including read-only attributes."""
        cls = self.__class__
        if not self.has_trait(name):
            raise TraitError("Class %s does not have a trait named %s" %
                                (cls.__name__, name))
        else:
            getattr(cls, name).set(self, value)

    @classmethod
    def class_trait_names(cls, **metadata):
        """Get a list of all the names of this class' traits.

        This method is just like the :meth:`trait_names` method,
        but is unbound.
        """
        return list(cls.class_traits(**metadata))

    @classmethod
    def class_traits(cls, **metadata):
        """Get a ``dict`` of all the traits of this class.  The dictionary
        is keyed on the name and the values are the TraitType objects.

        This method is just like the :meth:`traits` method, but is unbound.

        The TraitTypes returned don't know anything about the values
        that the various HasTrait's instances are holding.

        The metadata kwargs allow functions to be passed in which
        filter traits based on metadata values.  The functions should
        take a single value as an argument and return a boolean.  If
        any function returns False, then the trait is not included in
        the output.  If a metadata key doesn't exist, None will be passed
        to the function.
        """
        traits = dict([memb for memb in getmembers(cls) if
                     isinstance(memb[1], TraitType)])

        if len(metadata) == 0:
            return traits

        result = {}
        for name, trait in traits.items():
            for meta_name, meta_eval in metadata.items():
                if type(meta_eval) is not types.FunctionType:
                    meta_eval = _SimpleTest(meta_eval)
                if not meta_eval(trait.metadata.get(meta_name, None)):
                    break
            else:
                result[name] = trait

        return result

    @classmethod
    def class_own_traits(cls, **metadata):
        """Get a dict of all the traitlets defined on this class, not a parent.

        Works like `class_traits`, except for excluding traits from parents.
        """
        sup = super(cls, cls)
        return {n: t for (n, t) in cls.class_traits(**metadata).items()
                if getattr(sup, n, None) is not t}

    def has_trait(self, name):
        """Returns True if the object has a trait with the specified name."""
        return isinstance(getattr(self.__class__, name, None), TraitType)

    def trait_names(self, **metadata):
        """Get a list of all the names of this class' traits."""
        return list(self.traits(**metadata))

    def traits(self, **metadata):
        """Get a ``dict`` of all the traits of this class.  The dictionary
        is keyed on the name and the values are the TraitType objects.

        The TraitTypes returned don't know anything about the values
        that the various HasTrait's instances are holding.

        The metadata kwargs allow functions to be passed in which
        filter traits based on metadata values.  The functions should
        take a single value as an argument and return a boolean.  If
        any function returns False, then the trait is not included in
        the output.  If a metadata key doesn't exist, None will be passed
        to the function.
        """
        traits = dict([memb for memb in getmembers(self.__class__) if
                     isinstance(memb[1], TraitType)])

        if len(metadata) == 0:
            return traits

        result = {}
        for name, trait in traits.items():
            for meta_name, meta_eval in metadata.items():
                if type(meta_eval) is not types.FunctionType:
                    meta_eval = _SimpleTest(meta_eval)
                if not meta_eval(trait.metadata.get(meta_name, None)):
                    break
            else:
                result[name] = trait

        return result

    def trait_metadata(self, traitname, key, default=None):
        """Get metadata values for trait by key."""
        try:
            trait = getattr(self.__class__, traitname)
        except AttributeError:
            raise TraitError("Class %s does not have a trait named %s" %
                                (self.__class__.__name__, traitname))
        metadata_name = '_' + traitname + '_metadata'
        if hasattr(self, metadata_name) and key in getattr(self, metadata_name):
            return getattr(self, metadata_name).get(key, default)
        else:
            return trait.metadata.get(key, default)

    @classmethod
    def class_own_trait_events(cls, name):
        """Get a dict of all event handlers defined on this class, not a parent.

        Works like ``event_handlers``, except for excluding traits from parents.
        """
        sup = super(cls, cls)
        return {n: e for (n, e) in cls.events(name).items()
                if getattr(sup, n, None) is not e}

    @classmethod
    def trait_events(cls, name=None):
        """Get a ``dict`` of all the event handlers of this class.

        Parameters
        ----------
        name: str (default: None)
            The name of a trait of this class. If name is ``None`` then all
            the event handlers of this class will be returned instead.

        Returns
        -------
        The event handlers associated with a trait name, or all event handlers.
        """
        events = {}
        for k, v in getmembers(cls):
            if isinstance(v, EventHandler):
                if name is None:
                    events[k] = v
                elif name in v.trait_names:
                    events[k] = v
                elif hasattr(v, 'tags'):
                    if cls.trait_names(**v.tags):
                        events[k] = v
        return events

#-----------------------------------------------------------------------------
# Actual TraitTypes implementations/subclasses
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# TraitTypes subclasses for handling classes and instances of classes
#-----------------------------------------------------------------------------


class ClassBasedTraitType(TraitType):
    """
    A trait with error reporting and string -> type resolution for Type,
    Instance and This.
    """

    def _resolve_string(self, string):
        """
        Resolve a string supplied for a type into an actual object.
        """
        return import_item(string)

    def error(self, obj, value):
        kind = type(value)
        if six.PY2 and kind is InstanceType:
            msg = 'class %s' % value.__class__.__name__
        else:
            msg = '%s (i.e. %s)' % ( str( kind )[1:-1], repr( value ) )

        if obj is not None:
            e = "The '%s' trait of %s instance must be %s, but a value of %s was specified." \
                % (self.name, class_of(obj),
                   self.info(), msg)
        else:
            e = "The '%s' trait must be %s, but a value of %r was specified." \
                % (self.name, self.info(), msg)

        raise TraitError(e)


class Type(ClassBasedTraitType):
    """A trait whose value must be a subclass of a specified class."""

    def __init__ (self, default_value=Undefined, klass=None, **kwargs):
        """Construct a Type trait

        A Type trait specifies that its values must be subclasses of
        a particular class.

        If only ``default_value`` is given, it is used for the ``klass`` as
        well. If neither are given, both default to ``object``.

        Parameters
        ----------
        default_value : class, str or None
            The default value must be a subclass of klass.  If an str,
            the str must be a fully specified class name, like 'foo.bar.Bah'.
            The string is resolved into real class, when the parent
            :class:`HasTraits` class is instantiated.
        klass : class, str [ default object ]
            Values of this trait must be a subclass of klass.  The klass
            may be specified in a string like: 'foo.bar.MyClass'.
            The string is resolved into real class, when the parent
            :class:`HasTraits` class is instantiated.
        allow_none : bool [ default False ]
            Indicates whether None is allowed as an assignable value.
        """
        if default_value is Undefined:
            new_default_value = object if (klass is None) else klass
        else:
            new_default_value = default_value

        if klass is None:
            if (default_value is None) or (default_value is Undefined):
                klass = object
            else:
                klass = default_value

        if not (inspect.isclass(klass) or isinstance(klass, six.string_types)):
            raise TraitError("A Type trait must specify a class.")

        self.klass = klass

        super(Type, self).__init__(new_default_value, **kwargs)

    def validate(self, obj, value):
        """Validates that the value is a valid object instance."""
        if isinstance(value, six.string_types):
            try:
                value = self._resolve_string(value)
            except ImportError:
                raise TraitError("The '%s' trait of %s instance must be a type, but "
                                 "%r could not be imported" % (self.name, obj, value))
        try:
            if issubclass(value, self.klass):
                return value
        except:
            pass

        self.error(obj, value)

    def info(self):
        """ Returns a description of the trait."""
        if isinstance(self.klass, six.string_types):
            klass = self.klass
        else:
            klass = self.klass.__module__ + '.' + self.klass.__name__
        result = "a subclass of '%s'" % klass
        if self.allow_none:
            return result + ' or None'
        return result

    def instance_init(self, obj):
        self._resolve_classes()
        super(Type, self).instance_init(obj)

    def _resolve_classes(self):
        if isinstance(self.klass, six.string_types):
            self.klass = self._resolve_string(self.klass)
        if isinstance(self.default_value, six.string_types):
            self.default_value = self._resolve_string(self.default_value)

    def default_value_repr(self):
        value = self.default_value
        if isinstance(value, six.string_types):
            return repr(value)
        else:
            return repr('{}.{}'.format(value.__module__, value.__name__))


class Instance(ClassBasedTraitType):
    """A trait whose value must be an instance of a specified class.

    The value can also be an instance of a subclass of the specified class.

    Subclasses can declare default classes by overriding the klass attribute
    """

    klass = None

    def __init__(self, klass=None, args=None, kw=None, **kwargs):
        """Construct an Instance trait.

        This trait allows values that are instances of a particular
        class or its subclasses.  Our implementation is quite different
        from that of enthough.traits as we don't allow instances to be used
        for klass and we handle the ``args`` and ``kw`` arguments differently.

        Parameters
        ----------
        klass : class, str
            The class that forms the basis for the trait.  Class names
            can also be specified as strings, like 'foo.bar.Bar'.
        args : tuple
            Positional arguments for generating the default value.
        kw : dict
            Keyword arguments for generating the default value.
        allow_none : bool [ default False ]
            Indicates whether None is allowed as a value.

        Notes
        -----
        If both ``args`` and ``kw`` are None, then the default value is None.
        If ``args`` is a tuple and ``kw`` is a dict, then the default is
        created as ``klass(*args, **kw)``.  If exactly one of ``args`` or ``kw`` is
        None, the None is replaced by ``()`` or ``{}``, respectively.
        """
        if klass is None:
            klass = self.klass

        if (klass is not None) and (inspect.isclass(klass) or isinstance(klass, six.string_types)):
            self.klass = klass
        else:
            raise TraitError('The klass attribute must be a class'
                                ' not: %r' % klass)

        if (kw is not None) and not isinstance(kw, dict):
            raise TraitError("The 'kw' argument must be a dict or None.")
        if (args is not None) and not isinstance(args, tuple):
            raise TraitError("The 'args' argument must be a tuple or None.")

        self.default_args = args
        self.default_kwargs = kw

        super(Instance, self).__init__(**kwargs)

    def validate(self, obj, value):
        if isinstance(value, self.klass):
            return value
        else:
            self.error(obj, value)

    def info(self):
        if isinstance(self.klass, six.string_types):
            klass = self.klass
        else:
            klass = self.klass.__name__
        result = class_of(klass)
        if self.allow_none:
            return result + ' or None'

        return result

    def instance_init(self, obj):
        self._resolve_classes()
        super(Instance, self).instance_init(obj)

    def _resolve_classes(self):
        if isinstance(self.klass, six.string_types):
            self.klass = self._resolve_string(self.klass)

    def make_dynamic_default(self):
        if (self.default_args is None) and (self.default_kwargs is None):
            return None
        return self.klass(*(self.default_args or ()),
                          **(self.default_kwargs or {}))

    def default_value_repr(self):
        return repr(self.make_dynamic_default())


class ForwardDeclaredMixin(object):
    """
    Mixin for forward-declared versions of Instance and Type.
    """
    def _resolve_string(self, string):
        """
        Find the specified class name by looking for it in the module in which
        our this_class attribute was defined.
        """
        modname = self.this_class.__module__
        return import_item('.'.join([modname, string]))


class ForwardDeclaredType(ForwardDeclaredMixin, Type):
    """
    Forward-declared version of Type.
    """
    pass


class ForwardDeclaredInstance(ForwardDeclaredMixin, Instance):
    """
    Forward-declared version of Instance.
    """
    pass


class This(ClassBasedTraitType):
    """A trait for instances of the class containing this trait.

    Because how how and when class bodies are executed, the ``This``
    trait can only have a default value of None.  This, and because we
    always validate default values, ``allow_none`` is *always* true.
    """

    info_text = 'an instance of the same type as the receiver or None'

    def __init__(self, **kwargs):
        super(This, self).__init__(None, **kwargs)

    def validate(self, obj, value):
        # What if value is a superclass of obj.__class__?  This is
        # complicated if it was the superclass that defined the This
        # trait.
        if isinstance(value, self.this_class) or (value is None):
            return value
        else:
            self.error(obj, value)


class Union(TraitType):
    """A trait type representing a Union type."""

    def __init__(self, trait_types, **kwargs):
        """Construct a Union  trait.

        This trait allows values that are allowed by at least one of the
        specified trait types. A Union traitlet cannot have metadata on
        its own, besides the metadata of the listed types.

        Parameters
        ----------
        trait_types: sequence
            The list of trait types of length at least 1.

        Notes
        -----
        Union([Float(), Bool(), Int()]) attempts to validate the provided values
        with the validation function of Float, then Bool, and finally Int.
        """
        self.trait_types = trait_types
        self.info_text = " or ".join([tt.info() for tt in self.trait_types])
        super(Union, self).__init__(**kwargs)

    def class_init(self, cls, name):
        for trait_type in self.trait_types:
            trait_type.class_init(cls, None)
        super(Union, self).class_init(cls, name)

    def instance_init(self, obj):
        for trait_type in self.trait_types:
            trait_type.instance_init(obj)
        super(Union, self).instance_init(obj)

    def validate(self, obj, value):
        with obj.cross_validation_lock:
            for trait_type in self.trait_types:
                try:
                    v = trait_type._validate(obj, value)
                    # In the case of an element trait, the name is None
                    if self.name is not None:
                        setattr(obj, '_' + self.name + '_metadata', trait_type.metadata)
                    return v
                except TraitError:
                    continue
        self.error(obj, value)

    def __or__(self, other):
        if isinstance(other, Union):
            return Union(self.trait_types + other.trait_types)
        else:
            return Union(self.trait_types + [other])

    def make_dynamic_default(self):
        if self.default_value is not Undefined:
            return self.default_value
        for trait_type in self.trait_types:
            if trait_type.default_value is not Undefined:
                return trait_type.default_value
            elif hasattr(trait_type, 'make_dynamic_default'):
                return trait_type.make_dynamic_default()


#-----------------------------------------------------------------------------
# Basic TraitTypes implementations/subclasses
#-----------------------------------------------------------------------------


class Any(TraitType):
    """A trait which allows any value."""
    default_value = None
    info_text = 'any value'


def _validate_bounds(trait, obj, value):
    """
    Validate that a number to be applied to a trait is between bounds.

    If value is not between min_bound and max_bound, this raises a
    TraitError with an error message appropriate for this trait.
    """
    if trait.min is not None and value < trait.min:
        raise TraitError(
            "The value of the '{name}' trait of {klass} instance should "
            "not be less than {min_bound}, but a value of {value} was "
            "specified".format(
                name=trait.name, klass=class_of(obj),
                value=value, min_bound=trait.min))
    if trait.max is not None and value > trait.max:
        raise TraitError(
            "The value of the '{name}' trait of {klass} instance should "
            "not be greater than {max_bound}, but a value of {value} was "
            "specified".format(
                name=trait.name, klass=class_of(obj),
                value=value, max_bound=trait.max))
    return value


class Int(TraitType):
    """An int trait."""

    default_value = 0
    info_text = 'an int'

    def __init__(self, default_value=Undefined, allow_none=False, **kwargs):
        self.min = kwargs.pop('min', None)
        self.max = kwargs.pop('max', None)
        super(Int, self).__init__(default_value=default_value,
                                  allow_none=allow_none, **kwargs)

    def validate(self, obj, value):
        if not isinstance(value, int):
            self.error(obj, value)
        return _validate_bounds(self, obj, value)


class CInt(Int):
    """A casting version of the int trait."""

    def validate(self, obj, value):
        try:
            value = int(value)
        except:
            self.error(obj, value)
        return _validate_bounds(self, obj, value)


if six.PY2:
    class Long(TraitType):
        """A long integer trait."""

        default_value = 0
        info_text = 'a long'

        def __init__(self, default_value=Undefined, allow_none=False, **kwargs):
            self.min = kwargs.pop('min', None)
            self.max = kwargs.pop('max', None)
            super(Long, self).__init__(
                default_value=default_value,
                allow_none=allow_none, **kwargs)

        def _validate_long(self, obj, value):
            if isinstance(value, long):
                return value
            if isinstance(value, int):
                return long(value)
            self.error(obj, value)

        def validate(self, obj, value):
            value = self._validate_long(obj, value)
            return _validate_bounds(self, obj, value)


    class CLong(Long):
        """A casting version of the long integer trait."""

        def validate(self, obj, value):
            try:
                value = long(value)
            except:
                self.error(obj, value)
            return _validate_bounds(self, obj, value)


    class Integer(TraitType):
        """An integer trait.

        Longs that are unnecessary (<= sys.maxint) are cast to ints."""

        default_value = 0
        info_text = 'an integer'

        def __init__(self, default_value=Undefined, allow_none=False, **kwargs):
            self.min = kwargs.pop('min', None)
            self.max = kwargs.pop('max', None)
            super(Integer, self).__init__(
                default_value=default_value,
                allow_none=allow_none, **kwargs)

        def _validate_int(self, obj, value):
            if isinstance(value, int):
                return value
            if isinstance(value, long):
                # downcast longs that fit in int:
                # note that int(n > sys.maxint) returns a long, so
                # we don't need a condition on this cast
                return int(value)
            if sys.platform == "cli":
                from System import Int64
                if isinstance(value, Int64):
                    return int(value)
            self.error(obj, value)

        def validate(self, obj, value):
            value = self._validate_int(obj, value)
            return _validate_bounds(self, obj, value)

else:
    Long, CLong = Int, CInt
    Integer = Int


class Float(TraitType):
    """A float trait."""

    default_value = 0.0
    info_text = 'a float'

    def __init__(self, default_value=Undefined, allow_none=False, **kwargs):
        self.min = kwargs.pop('min', -float('inf'))
        self.max = kwargs.pop('max', float('inf'))
        super(Float, self).__init__(default_value=default_value,
                                    allow_none=allow_none, **kwargs)

    def validate(self, obj, value):
        if isinstance(value, int):
            value = float(value)
        if not isinstance(value, float):
            self.error(obj, value)
        return _validate_bounds(self, obj, value)


class CFloat(Float):
    """A casting version of the float trait."""

    def validate(self, obj, value):
        try:
            value = float(value)
        except:
            self.error(obj, value)
        return _validate_bounds(self, obj, value)


class Complex(TraitType):
    """A trait for complex numbers."""

    default_value = 0.0 + 0.0j
    info_text = 'a complex number'

    def validate(self, obj, value):
        if isinstance(value, complex):
            return value
        if isinstance(value, (float, int)):
            return complex(value)
        self.error(obj, value)


class CComplex(Complex):
    """A casting version of the complex number trait."""

    def validate (self, obj, value):
        try:
            return complex(value)
        except:
            self.error(obj, value)

# We should always be explicit about whether we're using bytes or unicode, both
# for Python 3 conversion and for reliable unicode behaviour on Python 2. So
# we don't have a Str type.
class Bytes(TraitType):
    """A trait for byte strings."""

    default_value = b''
    info_text = 'a bytes object'

    def validate(self, obj, value):
        if isinstance(value, bytes):
            return value
        self.error(obj, value)


class CBytes(Bytes):
    """A casting version of the byte string trait."""

    def validate(self, obj, value):
        try:
            return bytes(value)
        except:
            self.error(obj, value)


class Unicode(TraitType):
    """A trait for unicode strings."""

    default_value = u''
    info_text = 'a unicode string'

    def validate(self, obj, value):
        if isinstance(value, six.text_type):
            return value
        if isinstance(value, bytes):
            try:
                return value.decode('ascii', 'strict')
            except UnicodeDecodeError:
                msg = "Could not decode {!r} for unicode trait '{}' of {} instance."
                raise TraitError(msg.format(value, self.name, class_of(obj)))
        self.error(obj, value)


class CUnicode(Unicode):
    """A casting version of the unicode trait."""

    def validate(self, obj, value):
        try:
            return six.text_type(value)
        except:
            self.error(obj, value)


class ObjectName(TraitType):
    """A string holding a valid object name in this version of Python.

    This does not check that the name exists in any scope."""
    info_text = "a valid object identifier in Python"

    if six.PY2:
        # Python 2:
        def coerce_str(self, obj, value):
            "In Python 2, coerce ascii-only unicode to str"
            if isinstance(value, unicode):
                try:
                    return str(value)
                except UnicodeEncodeError:
                    self.error(obj, value)
            return value
    else:
        coerce_str = staticmethod(lambda _,s: s)

    def validate(self, obj, value):
        value = self.coerce_str(obj, value)

        if isinstance(value, six.string_types) and isidentifier(value):
            return value
        self.error(obj, value)

class DottedObjectName(ObjectName):
    """A string holding a valid dotted object name in Python, such as A.b3._c"""
    def validate(self, obj, value):
        value = self.coerce_str(obj, value)

        if isinstance(value, six.string_types) and all(isidentifier(a)
          for a in value.split('.')):
            return value
        self.error(obj, value)


class Bool(TraitType):
    """A boolean (True, False) trait."""

    default_value = False
    info_text = 'a boolean'

    def validate(self, obj, value):
        if isinstance(value, bool):
            return value
        self.error(obj, value)


class CBool(Bool):
    """A casting version of the boolean trait."""

    def validate(self, obj, value):
        try:
            return bool(value)
        except:
            self.error(obj, value)


class Enum(TraitType):
    """An enum whose value must be in a given sequence."""

    def __init__(self, values, default_value=Undefined, **kwargs):
        self.values = values
        if kwargs.get('allow_none', False) and default_value is Undefined:
            default_value = None
        super(Enum, self).__init__(default_value, **kwargs)

    def validate(self, obj, value):
        if value in self.values:
                return value
        self.error(obj, value)

    def info(self):
        """ Returns a description of the trait."""
        result = 'any of ' + repr(self.values)
        if self.allow_none:
            return result + ' or None'
        return result

class CaselessStrEnum(Enum):
    """An enum of strings where the case should be ignored."""

    def __init__(self, values, default_value=Undefined, **kwargs):
        values = [cast_unicode_py2(value) for value in values]
        super(CaselessStrEnum, self).__init__(values, default_value=default_value, **kwargs)

    def validate(self, obj, value):
        if isinstance(value, str):
            value = cast_unicode_py2(value)
        if not isinstance(value, six.string_types):
            self.error(obj, value)

        for v in self.values:
            if v.lower() == value.lower():
                return v
        self.error(obj, value)

class Container(Instance):
    """An instance of a container (list, set, etc.)

    To be subclassed by overriding klass.
    """
    klass = None
    _cast_types = ()
    _valid_defaults = SequenceTypes
    _trait = None

    def __init__(self, trait=None, default_value=None, **kwargs):
        """Create a container trait type from a list, set, or tuple.

        The default value is created by doing ``List(default_value)``,
        which creates a copy of the ``default_value``.

        ``trait`` can be specified, which restricts the type of elements
        in the container to that TraitType.

        If only one arg is given and it is not a Trait, it is taken as
        ``default_value``:

        ``c = List([1, 2, 3])``

        Parameters
        ----------

        trait : TraitType [ optional ]
            the type for restricting the contents of the Container.  If unspecified,
            types are not checked.

        default_value : SequenceType [ optional ]
            The default value for the Trait.  Must be list/tuple/set, and
            will be cast to the container type.

        allow_none : bool [ default False ]
            Whether to allow the value to be None

        **kwargs : any
            further keys for extensions to the Trait (e.g. config)

        """
        # allow List([values]):
        if default_value is None and not is_trait(trait):
            default_value = trait
            trait = None

        if default_value is None:
            args = ()
        elif isinstance(default_value, self._valid_defaults):
            args = (default_value,)
        else:
            raise TypeError('default value of %s was %s' %(self.__class__.__name__, default_value))

        if is_trait(trait):
            if isinstance(trait, type):
                warn("Traits should be given as instances, not types (for example, `Int()`, not `Int`)."
                     " Passing types is deprecated in traitlets 4.1.",
                     DeprecationWarning, stacklevel=3)
            self._trait = trait() if isinstance(trait, type) else trait
        elif trait is not None:
            raise TypeError("`trait` must be a Trait or None, got %s" % repr_type(trait))

        super(Container,self).__init__(klass=self.klass, args=args, **kwargs)

    def element_error(self, obj, element, validator):
        e = "Element of the '%s' trait of %s instance must be %s, but a value of %s was specified." \
            % (self.name, class_of(obj), validator.info(), repr_type(element))
        raise TraitError(e)

    def validate(self, obj, value):
        if isinstance(value, self._cast_types):
            value = self.klass(value)
        value = super(Container, self).validate(obj, value)
        if value is None:
            return value

        value = self.validate_elements(obj, value)

        return value

    def validate_elements(self, obj, value):
        validated = []
        if self._trait is None or isinstance(self._trait, Any):
            return value
        for v in value:
            try:
                v = self._trait._validate(obj, v)
            except TraitError:
                self.element_error(obj, v, self._trait)
            else:
                validated.append(v)
        return self.klass(validated)

    def class_init(self, cls, name):
        if isinstance(self._trait, TraitType):
            self._trait.class_init(cls, None)
        super(Container, self).class_init(cls, name)

    def instance_init(self, obj):
        if isinstance(self._trait, TraitType):
            self._trait.instance_init(obj)
        super(Container, self).instance_init(obj)


class List(Container):
    """An instance of a Python list."""
    klass = list
    _cast_types = (tuple,)

    def __init__(self, trait=None, default_value=None, minlen=0, maxlen=sys.maxsize, **kwargs):
        """Create a List trait type from a list, set, or tuple.

        The default value is created by doing ``list(default_value)``,
        which creates a copy of the ``default_value``.

        ``trait`` can be specified, which restricts the type of elements
        in the container to that TraitType.

        If only one arg is given and it is not a Trait, it is taken as
        ``default_value``:

        ``c = List([1, 2, 3])``

        Parameters
        ----------

        trait : TraitType [ optional ]
            the type for restricting the contents of the Container.
            If unspecified, types are not checked.

        default_value : SequenceType [ optional ]
            The default value for the Trait.  Must be list/tuple/set, and
            will be cast to the container type.

        minlen : Int [ default 0 ]
            The minimum length of the input list

        maxlen : Int [ default sys.maxsize ]
            The maximum length of the input list
        """
        self._minlen = minlen
        self._maxlen = maxlen
        super(List, self).__init__(trait=trait, default_value=default_value,
                                   **kwargs)

    def length_error(self, obj, value):
        e = "The '%s' trait of %s instance must be of length %i <= L <= %i, but a value of %s was specified." \
            % (self.name, class_of(obj), self._minlen, self._maxlen, value)
        raise TraitError(e)

    def validate_elements(self, obj, value):
        length = len(value)
        if length < self._minlen or length > self._maxlen:
            self.length_error(obj, value)

        return super(List, self).validate_elements(obj, value)

    def validate(self, obj, value):
        value = super(List, self).validate(obj, value)
        value = self.validate_elements(obj, value)
        return value


class Set(List):
    """An instance of a Python set."""
    klass = set
    _cast_types = (tuple, list)

    # Redefine __init__ just to make the docstring more accurate.
    def __init__(self, trait=None, default_value=None, minlen=0, maxlen=sys.maxsize,
                 **kwargs):
        """Create a Set trait type from a list, set, or tuple.

        The default value is created by doing ``set(default_value)``,
        which creates a copy of the ``default_value``.

        ``trait`` can be specified, which restricts the type of elements
        in the container to that TraitType.

        If only one arg is given and it is not a Trait, it is taken as
        ``default_value``:

        ``c = Set({1, 2, 3})``

        Parameters
        ----------

        trait : TraitType [ optional ]
            the type for restricting the contents of the Container.
            If unspecified, types are not checked.

        default_value : SequenceType [ optional ]
            The default value for the Trait.  Must be list/tuple/set, and
            will be cast to the container type.

        minlen : Int [ default 0 ]
            The minimum length of the input list

        maxlen : Int [ default sys.maxsize ]
            The maximum length of the input list
        """
        super(Set, self).__init__(trait, default_value, minlen, maxlen, **kwargs)


class Tuple(Container):
    """An instance of a Python tuple."""
    klass = tuple
    _cast_types = (list,)

    def __init__(self, *traits, **kwargs):
        """Create a tuple from a list, set, or tuple.

        Create a fixed-type tuple with Traits:

        ``t = Tuple(Int(), Str(), CStr())``

        would be length 3, with Int,Str,CStr for each element.

        If only one arg is given and it is not a Trait, it is taken as
        default_value:

        ``t = Tuple((1, 2, 3))``

        Otherwise, ``default_value`` *must* be specified by keyword.

        Parameters
        ----------

        `*traits` : TraitTypes [ optional ]
            the types for restricting the contents of the Tuple.  If unspecified,
            types are not checked. If specified, then each positional argument
            corresponds to an element of the tuple.  Tuples defined with traits
            are of fixed length.

        default_value : SequenceType [ optional ]
            The default value for the Tuple.  Must be list/tuple/set, and
            will be cast to a tuple. If ``traits`` are specified,
            ``default_value`` must conform to the shape and type they specify.
        """
        default_value = kwargs.pop('default_value', Undefined)
        # allow Tuple((values,)):
        if len(traits) == 1 and default_value is Undefined and not is_trait(traits[0]):
            default_value = traits[0]
            traits = ()

        if default_value is Undefined:
            args = ()
        elif isinstance(default_value, self._valid_defaults):
            args = (default_value,)
        else:
            raise TypeError('default value of %s was %s' %(self.__class__.__name__, default_value))

        self._traits = []
        for trait in traits:
            if isinstance(trait, type):
                warn("Traits should be given as instances, not types (for example, `Int()`, not `Int`)"
                     " Passing types is deprecated in traitlets 4.1.",
                     DeprecationWarning, stacklevel=2)
            t = trait() if isinstance(trait, type) else trait
            self._traits.append(t)

        if self._traits and default_value is None:
            # don't allow default to be an empty container if length is specified
            args = None
        super(Container,self).__init__(klass=self.klass, args=args, **kwargs)

    def validate_elements(self, obj, value):
        if not self._traits:
            # nothing to validate
            return value
        if len(value) != len(self._traits):
            e = "The '%s' trait of %s instance requires %i elements, but a value of %s was specified." \
                % (self.name, class_of(obj), len(self._traits), repr_type(value))
            raise TraitError(e)

        validated = []
        for t, v in zip(self._traits, value):
            try:
                v = t._validate(obj, v)
            except TraitError:
                self.element_error(obj, v, t)
            else:
                validated.append(v)
        return tuple(validated)

    def class_init(self, cls, name):
        for trait in self._traits:
            if isinstance(trait, TraitType):
                trait.class_init(cls, None)
        super(Container, self).class_init(cls, name)

    def instance_init(self, obj):
        for trait in self._traits:
            if isinstance(trait, TraitType):
                trait.instance_init(obj)
        super(Container, self).instance_init(obj)


class Dict(Instance):
    """An instance of a Python dict."""
    _trait = None

    def __init__(self, trait=None, traits=None, default_value=Undefined,
                 **kwargs):
        """Create a dict trait type from a Python dict.

        The default value is created by doing ``dict(default_value)``,
        which creates a copy of the ``default_value``.

        Parameters
        ----------

        trait : TraitType [ optional ]
            The specified trait type to check and use to restrict contents of
            the Container. If unspecified, trait types are not checked.

        traits : Dictionary of trait types [ optional ]
            A Python dictionary containing the types that are valid for
            restricting the content of the Dict Container for certain keys.

        default_value : SequenceType [ optional ]
            The default value for the Dict.  Must be dict, tuple, or None, and
            will be cast to a dict if not None. If `trait` is specified, the
            `default_value` must conform to the constraints it specifies.
        """
        # Handling positional arguments
        if default_value is Undefined and trait is not None:
            if not is_trait(trait):
                default_value = trait
                trait = None

        # Handling default value
        if default_value is Undefined:
            default_value = {}
        if default_value is None:
            args = None
        elif isinstance(default_value, dict):
            args = (default_value,)
        elif isinstance(default_value, SequenceTypes):
            args = (default_value,)
        else:
            raise TypeError('default value of Dict was %s' % default_value)

        # Case where a type of TraitType is provided rather than an instance
        if is_trait(trait):
            if isinstance(trait, type):
                warn("Traits should be given as instances, not types (for example, `Int()`, not `Int`)"
                     " Passing types is deprecated in traitlets 4.1.",
                     DeprecationWarning, stacklevel=2)
            self._trait = trait() if isinstance(trait, type) else trait
        elif trait is not None:
            raise TypeError("`trait` must be a Trait or None, got %s" % repr_type(trait))

        self._traits = traits

        super(Dict, self).__init__(klass=dict, args=args, **kwargs)

    def element_error(self, obj, element, validator):
        e = "Element of the '%s' trait of %s instance must be %s, but a value of %s was specified." \
            % (self.name, class_of(obj), validator.info(), repr_type(element))
        raise TraitError(e)

    def validate(self, obj, value):
        value = super(Dict, self).validate(obj, value)
        if value is None:
            return value
        value = self.validate_elements(obj, value)
        return value

    def validate_elements(self, obj, value):
        use_dict = bool(self._traits)
        default_to = (self._trait or Any())
        if not use_dict and isinstance(default_to, Any):
            return value

        validated = {}
        for key in value:
            if use_dict and key in self._traits:
                validate_with = self._traits[key]
            else:
                validate_with = default_to
            try:
                v = value[key]
                if not isinstance(validate_with, Any):
                    v = validate_with._validate(obj, v)
            except TraitError:
                self.element_error(obj, v, validate_with)
            else:
                validated[key] = v

        return self.klass(validated)

    def class_init(self, cls, name):
        if isinstance(self._trait, TraitType):
            self._trait.class_init(cls, None)
        if self._traits is not None:
            for trait in self._traits.values():
                trait.class_init(cls, None)
        super(Dict, self).class_init(cls, name)

    def instance_init(self, obj):
        if isinstance(self._trait, TraitType):
            self._trait.instance_init(obj)
        if self._traits is not None:
            for trait in self._traits.values():
                trait.instance_init(obj)
        super(Dict, self).instance_init(obj)


class TCPAddress(TraitType):
    """A trait for an (ip, port) tuple.

    This allows for both IPv4 IP addresses as well as hostnames.
    """

    default_value = ('127.0.0.1', 0)
    info_text = 'an (ip, port) tuple'

    def validate(self, obj, value):
        if isinstance(value, tuple):
            if len(value) == 2:
                if isinstance(value[0], six.string_types) and isinstance(value[1], int):
                    port = value[1]
                    if port >= 0 and port <= 65535:
                        return value
        self.error(obj, value)

class CRegExp(TraitType):
    """A casting compiled regular expression trait.

    Accepts both strings and compiled regular expressions. The resulting
    attribute will be a compiled regular expression."""

    info_text = 'a regular expression'

    def validate(self, obj, value):
        try:
            return re.compile(value)
        except:
            self.error(obj, value)


class UseEnum(TraitType):
    """Use a Enum class as model for the data type description.
    Note that if no default-value is provided, the first enum-value is used
    as default-value.

    .. sourcecode:: python

        # -- SINCE: Python 3.4 (or install backport: pip install enum34)
        import enum
        from traitlets import HasTraits, UseEnum

        class Color(enum.Enum):
            red = 1         # -- IMPLICIT: default_value
            blue = 2
            green = 3

        class MyEntity(HasTraits):
            color = UseEnum(Color, default_value=Color.blue)

        entity = MyEntity(color=Color.red)
        entity.color = Color.green    # USE: Enum-value (preferred)
        entity.color = "green"        # USE: name (as string)
        entity.color = "Color.green"  # USE: scoped-name (as string)
        entity.color = 3              # USE: number (as int)
        assert entity.color is Color.green
    """
    default_value = None
    info_text = "Trait type adapter to a Enum class"

    def __init__(self, enum_class, default_value=None, **kwargs):
        assert issubclass(enum_class, enum.Enum), \
                          "REQUIRE: enum.Enum, but was: %r" % enum_class
        allow_none = kwargs.get("allow_none", False)
        if default_value is None and not allow_none:
            default_value = list(enum_class.__members__.values())[0]
        super(UseEnum, self).__init__(default_value=default_value, **kwargs)
        self.enum_class = enum_class
        self.name_prefix = enum_class.__name__ + "."

    def select_by_number(self, value, default=Undefined):
        """Selects enum-value by using its number-constant."""
        assert isinstance(value, int)
        enum_members = self.enum_class.__members__
        for enum_item in enum_members.values():
            if enum_item.value == value:
                return enum_item
        # -- NOT FOUND:
        return default

    def select_by_name(self, value, default=Undefined):
        """Selects enum-value by using its name or scoped-name."""
        assert isinstance(value, six.string_types)
        if value.startswith(self.name_prefix):
            # -- SUPPORT SCOPED-NAMES, like: "Color.red" => "red"
            value = value.replace(self.name_prefix, "", 1)
        return self.enum_class.__members__.get(value, default)

    def validate(self, obj, value):
        if isinstance(value, self.enum_class):
            return value
        elif isinstance(value, int):
            # -- CONVERT: number => enum_value (item)
            value2 = self.select_by_number(value)
            if value2 is not Undefined:
                return value2
        elif isinstance(value, six.string_types):
            # -- CONVERT: name or scoped_name (as string) => enum_value (item)
            value2 = self.select_by_name(value)
            if value2 is not Undefined:
                return value2
        elif value is None:
            if self.allow_none:
                return None
            else:
                return self.default_value
        self.error(obj, value)

    def info(self):
        """Returns a description of this Enum trait (in case of errors)."""
        result = "Any of: %s" % ", ".join(self.enum_class.__members__.keys())
        if self.allow_none:
            return result + " or None"
        return result

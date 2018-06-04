from __future__ import unicode_literals
from six import with_metaclass
from collections import defaultdict
import weakref

__all__ = (
    'CLIFilter',
    'SimpleFilter',
)

# Cache for _FilterTypeMeta. (Don't test the same __instancecheck__ twice as
# long as the object lives. -- We do this a lot and calling 'test_args' is
# expensive.)
_instance_check_cache = defaultdict(weakref.WeakKeyDictionary)


class _FilterTypeMeta(type):
    def __instancecheck__(cls, instance):
        cache = _instance_check_cache[tuple(cls.arguments_list)]

        def get():
            " The actual test. "
            if not hasattr(instance, 'test_args'):
                return False
            return instance.test_args(*cls.arguments_list)

        try:
            return cache[instance]
        except KeyError:
            result = get()
            cache[instance] = result
            return result


class _FilterType(with_metaclass(_FilterTypeMeta)):
    def __new__(cls):
        raise NotImplementedError('This class should not be initiated.')


class CLIFilter(_FilterType):
    """
    Abstract base class for filters that accept a
    :class:`~prompt_toolkit.interface.CommandLineInterface` argument. It cannot
    be instantiated, it's only to be used for instance assertions, e.g.::

        isinstance(my_filter, CliFilter)
    """
    arguments_list = ['cli']


class SimpleFilter(_FilterType):
    """
    Abstract base class for filters that don't accept any arguments.
    """
    arguments_list = []

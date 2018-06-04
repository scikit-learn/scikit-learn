"""
Prompt_toolkit is designed a way that the amount of changing state is reduced
to a minimum. Where possible, code is written in a pure functional way. In
general, this results in code where the flow is very easy to follow: the value
of a variable can be deducted from its first assignment.

However, often, practicality and performance beat purity and some classes still
have a changing state. In order to not having to care too much about
transferring states between several components we use some reactive
programming. Actually some kind of data binding.

We introduce two types:

- Filter: for binding a boolean state. They can be chained using & and |
  operators. Have a look in the ``filters`` module. Resolving the actual value
  of a filter happens by calling it.

- Integer: for binding integer values. Reactive operations (like addition and
  substraction) are not suppported. Resolving the actual value happens by
  casting it to int, like  ``int(integer)``. This way, it is possible to use
  normal integers as well for static values.
"""
from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from six import with_metaclass


class Integer(with_metaclass(ABCMeta, object)):
    """
    Reactive integer -- anything that can be resolved to an ``int``.
    """
    @abstractmethod
    def __int__(self):
        return 0

    @classmethod
    def from_callable(cls, func):
        """
        Create an Integer-like object that calls the given function when it is
        resolved to an int.
        """
        return _IntegerFromCallable(func)


Integer.register(int)


class _IntegerFromCallable(Integer):
    def __init__(self, func=0):
        self.func = func

    def __repr__(self):
        return 'Integer.from_callable(%r)' % self.func

    def __int__(self):
        return int(self.func())

"""
Backport of PEP 562.

https://pypi.org/search/?q=pep562

Licensed under MIT
Copyright (c) 2018 Isaac Muse <isaacmuse@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions
of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
"""
from __future__ import unicode_literals
import sys

__all__ = ('Pep562',)


class Pep562(object):
    """
    Backport of PEP 562 <https://pypi.org/search/?q=pep562>.

    Wraps the module in a class that exposes the mechanics to override `__dir__` and `__getattr__`.
    The given module will be searched for overrides of `__dir__` and `__getattr__` and use them when needed.
    """

    def __init__(self, name):
        """Acquire `__getattr__` and `__dir__`, but only replace module for versions less than Python 3.7."""

        self._module = sys.modules[name]
        self._get_attr = getattr(self._module, '__getattr__', None)
        self._get_dir = getattr(self._module, '__dir__', None)
        sys.modules[name] = self

    def __dir__(self):
        """Return the overridden `dir` if one was provided, else apply `dir` to the module."""

        return self._get_dir() if self._get_dir else dir(self._module)

    def __getattr__(self, name):
        """Attempt to retrieve the attribute from the module, and if missing, use the overridden function if present."""

        try:
            return getattr(self._module, name)
        except AttributeError:
            if self._get_attr:
                return self._get_attr(name)
            raise

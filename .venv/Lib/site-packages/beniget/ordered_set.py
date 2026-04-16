"""
Copied from https://github.com/bustawin/ordered-set-37
"""
# Unlicense
# This is free and unencumbered software released into the public domain.

# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.

# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

# For more information, please refer to <http://unlicense.org/>


import sys

from collections import OrderedDict
import itertools
from typing import TYPE_CHECKING, MutableSet

if TYPE_CHECKING:
    # trying to avoid polluting the global namespace with typing names.
    from typing import TypeVar, Iterator, Iterable, Optional
    T = TypeVar("T")

class ordered_set(MutableSet['T']):
    """
    A set that preserves insertion order by internally using a dict.
    """
    
    __slots__ = ('values',)

    def __init__(self, elements: 'Optional[Iterable[T]]' = None):
        self.values = OrderedDict.fromkeys(elements or [])

    def add(self, x: 'T') -> None:
        self.values[x] = None
    
    def update(self, values:'Iterable[T]') -> None:
        self.values.update((k, None) for k in values)

    def clear(self) -> None:
        self.values.clear()

    def discard(self, x: 'T') -> None:
        self.values.pop(x, None)

    def __getitem__(self, index:int) -> 'T':
        try:
            return next(itertools.islice(self.values, index, index + 1))
        except StopIteration:
            raise IndexError(f"index {index} out of range")

    def __contains__(self, x: object) -> bool:
        return self.values.__contains__(x)
    
    def __add__(self, other:'ordered_set[T]') -> 'ordered_set[T]':
        return ordered_set(itertools.chain(self, other))

    def __len__(self) -> int:
        return self.values.__len__()

    def __iter__(self) -> 'Iterator[T]':
        return self.values.__iter__()

    def __str__(self) -> str:
        return f"{{{', '.join(str(i) for i in self)}}}"

    def __repr__(self) -> str:
        return f"<ordered_set {self}>"
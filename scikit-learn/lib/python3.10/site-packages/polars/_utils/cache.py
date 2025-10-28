from __future__ import annotations

from collections import OrderedDict
from collections.abc import MutableMapping
from typing import TYPE_CHECKING, Any, TypeVar, overload

from polars._utils.various import no_default

if TYPE_CHECKING:
    import sys
    from collections.abc import ItemsView, Iterable, Iterator, KeysView, ValuesView

    from polars._utils.various import NoDefault

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

D = TypeVar("D")
K = TypeVar("K")
V = TypeVar("V")


class LRUCache(MutableMapping[K, V]):
    def __init__(self, maxsize: int) -> None:
        """
        Initialize an LRU (Least Recently Used) cache with a specified maximum size.

        Parameters
        ----------
        maxsize : int
            The maximum number of items the cache can hold.

        Examples
        --------
        >>> from polars._utils.cache import LRUCache
        >>> cache = LRUCache[str, int](maxsize=3)
        >>> cache["a"] = 1
        >>> cache["b"] = 2
        >>> cache["c"] = 3
        >>> cache["d"] = 4  # evicts the least recently used item ("a"), as maxsize=3
        >>> print(cache["b"])  # accessing "b" marks it as recently used
        2
        >>> print(list(cache.keys()))  # show the current keys in LRU order
        ['c', 'd', 'b']
        >>> cache.get("xyz", "not found")
        'not found'
        """
        self._items: OrderedDict[K, V] = OrderedDict()
        self.maxsize = maxsize

    def __bool__(self) -> bool:
        """Returns True if the cache is not empty, False otherwise."""
        return bool(self._items)

    def __contains__(self, key: Any) -> bool:
        """Check if the key is in the cache."""
        return key in self._items

    def __delitem__(self, key: K) -> None:
        """Remove the item with the specified key from the cache."""
        if key not in self._items:
            msg = f"{key!r} not found in cache"
            raise KeyError(msg)
        del self._items[key]

    def __getitem__(self, key: K) -> V:
        """Raises KeyError if the key is not found."""
        if key not in self._items:
            msg = f"{key!r} not found in cache"
            raise KeyError(msg)

        # moving accessed items to the end marks them as recently used
        self._items.move_to_end(key)
        return self._items[key]

    def __iter__(self) -> Iterator[K]:
        """Iterate over the keys in the cache."""
        yield from self._items

    def __len__(self) -> int:
        """Number of items in the cache."""
        return len(self._items)

    def __setitem__(self, key: K, value: V) -> None:
        """Insert a value into the cache."""
        if self._max_size == 0:
            return
        while len(self) >= self._max_size:
            self.popitem()
        if key in self:
            # moving accessed items to the end marks them as recently used
            self._items.move_to_end(key)
        self._items[key] = value

    def __repr__(self) -> str:
        """Return a string representation of the cache."""
        all_items = list(self._items.items())
        if len(self) > 4:
            items = (
                ", ".join(f"{k!r}: {v!r}" for k, v in all_items[:2])
                + " ..., "
                + ", ".join(f"{k!r}: {v!r}" for k, v in all_items[-2:])
            )
        else:
            items = ", ".join(f"{k!r}: {v!r}" for k, v in all_items)
        return f"{self.__class__.__name__}({{{items}}}, maxsize={self._max_size}, currsize={len(self)})"

    def clear(self) -> None:
        """Clear the cache, removing all items."""
        self._items.clear()

    @overload
    def get(self, key: K, default: None = None) -> V | None: ...

    @overload
    def get(self, key: K, default: D = ...) -> V | D: ...

    def get(self, key: K, default: D | V | None = None) -> V | D | None:
        """Return value associated with `key` if present, otherwise return `default`."""
        if key in self:
            # moving accessed items to the end marks them as recently used
            self._items.move_to_end(key)
            return self._items[key]
        return default

    @classmethod
    def fromkeys(cls, maxsize: int, *, keys: Iterable[K], value: V) -> Self:
        """Initialize cache with keys from an iterable, all set to the same value."""
        cache = cls(maxsize)
        for key in keys:
            cache[key] = value
        return cache

    def items(self) -> ItemsView[K, V]:
        """Return an iterable view of the cache's items (keys and values)."""
        return self._items.items()

    def keys(self) -> KeysView[K]:
        """Return an iterable view of the cache's keys."""
        return self._items.keys()

    @property
    def maxsize(self) -> int:
        return self._max_size

    @maxsize.setter
    def maxsize(self, n: int) -> None:
        """Set new maximum cache size; cache is trimmed if value is smaller."""
        if n < 0:
            msg = f"`maxsize` cannot be negative; found {n}"
            raise ValueError(msg)
        while len(self) > n:
            self.popitem()
        self._max_size = n

    def pop(self, key: K, default: D | NoDefault = no_default) -> V | D:
        """
        Remove specified key from the cache and return the associated value.

        If the key is not found, `default` is returned (if given).
        Otherwise, a KeyError is raised.
        """
        if (item := self._items.pop(key, default)) is no_default:
            msg = f"{key!r} not found in cache"
            raise KeyError(msg)
        return item

    def popitem(self) -> tuple[K, V]:
        """Remove the least recently used value; raises KeyError if cache is empty."""
        return self._items.popitem(last=False)

    def values(self) -> ValuesView[V]:
        """Return an iterable view of the cache's values."""
        return self._items.values()

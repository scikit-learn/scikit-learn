from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Generic, Hashable, TypeVar

try:
    from typing import Self  # pragma: ≥ 3.11 cover
except ImportError:
    from typing_extensions import Self  # pragma: < 3.11 cover

K = TypeVar("K", bound=Hashable)


class Cache(ABC, Generic[K]):
    """
    A generic cache interface.

    Add a close() method if the cache needs to perform any cleanup actions,
    and an __exit__ method to allow it to be used in a context manager.
    """

    @abstractmethod
    def get(self, key: K) -> Any | None:
        """
        Get a value from the cache.

        :param key: the key to retrieve
        :return: the cached value, or None if not found
        """
        raise NotImplementedError

    @abstractmethod
    def set(self, key: K, value: Any) -> None:
        """
        Set a value in the cache.

        :param key: the key to set
        :param value: the value to cache
        """
        raise NotImplementedError

    @abstractmethod
    def remove(self, key: K) -> None:
        """
        Remove a value from the cache.

        :param key: the key to remove
        """
        raise NotImplementedError

    @abstractmethod
    def clear(self) -> None:
        """Clear the entire cache."""
        raise NotImplementedError

    def __enter__(self) -> Self:
        return self


__all__ = [
    "Cache",
]

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Infrastructure for testing thread-safety."""

from collections import Counter
from collections.abc import Iterable
from hashlib import md5
from pickle import PicklingError, dumps
from types import MethodType
from typing import Any

from sklearn.utils.parallel import _FuncWrapper


class DetectChanges:
    """
    Heuristically detect if a set of objects have changed.

    This can be used to detect potential thread-safety issues.  If the
    arguments to tasks running in a thread pool have changed while running in
    the thread pool, that suggests unsafe mutation.
    """

    def __init__(self):
        self.counts = Counter()
        self.hashes = {}
        self.objects = {}

    def store_object(self, obj: Any) -> None:
        """Record an object's hash, if possible."""
        obj_id = id(obj)
        try:
            obj_hash = object_hash(obj)
        except (TypeError, PicklingError):
            return
        self.counts[obj_id] += 1
        self.hashes[obj_id] = obj_hash
        self.objects[obj_id] = obj

    def store_parallel_calls(
        self, iterable: Iterable[tuple[_FuncWrapper, tuple[Any], dict[str, Any]]]
    ) -> Iterable[tuple[_FuncWrapper, tuple[Any], dict[str, Any]]]:
        """
        Record hashes for inputs to ``Parallel.__call__``
        """
        for f, args, kwargs in iterable:
            if isinstance(f.function, MethodType):
                self.store_object(f.function.__self__)
            for obj in args:
                self.store_object(obj)
            for obj in kwargs.values():
                self.store_object(obj)
            yield (f, args, kwargs)

    def assert_objects_unchanged(self, min_count=1) -> None:
        """
        Raise an ``AssertionError`` if stored objects have changed.
        """
        for obj_id, count in self.counts.items():
            if count < min_count:
                continue
            obj = self.objects[obj_id]
            old_hash = self.hashes[obj_id]
            if object_hash(obj) != old_hash:
                raise AssertionError(
                    f"Hash for object of type {type(obj)} has changed, "
                    "suggesting potential thread unsafety"
                )


def object_hash(obj: Any) -> bytes:
    """
    Hash a Python object; mutating the object should change the result.

    Raises ``TypeError`` in some rare cases.
    """
    return md5(dumps(obj), usedforsecurity=False).digest()

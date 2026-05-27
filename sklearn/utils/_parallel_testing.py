# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Infrastructure for testing thread-safety."""

from collections import Counter
from collections.abc import Iterable
from copy import deepcopy
from dataclasses import dataclass
from inspect import getmro
from pickle import dumps, loads
from types import MethodType
from typing import Any, TypeVar

from deepdiff import DeepHash

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
        self.old_objects = {}
        self.objects = {}

    def store_object(self, obj: Any) -> None:
        """Record an object's hash, if possible."""
        if isinstance(obj, _ThreadSafe):
            return
        obj_id = id(obj)
        try:
            obj_hash = object_hash(obj)
        except TypeError:
            return
        self.counts[obj_id] += 1
        self.hashes[obj_id] = obj_hash
        self.objects[obj_id] = obj
        if hasattr(obj, "__deepcopy__"):
            # Fallback to pickle:
            copied = loads(dumps(obj))
        else:
            copied = deepcopy(obj)
        self.old_objects[obj_id] = obj

    def store_parallel_calls(
        self, iterable: Iterable[tuple[_FuncWrapper, tuple[Any], dict[str, Any]]]
    ) -> Iterable[tuple[_FuncWrapper, tuple[Any], dict[str, Any]]]:
        """
        Record hashes for inputs to ``Parallel.__call__``.

        Skip objects wrapped in ``_ThreadSafe``.
        """
        for f, args, kwargs in iterable:
            if isinstance(f.function, MethodType):
                self.store_object(f.function.__self__)
            new_args = []
            for obj in args:
                if isinstance(obj, _ThreadSafe):
                    obj = obj.unwrap()
                else:
                    self.store_object(obj)
                new_args.append(obj)
            new_kwargs = kwargs.copy()
            for key, obj in kwargs.items():
                if isinstance(obj, _ThreadSafe):
                    new_kwargs[key] = obj.unwrap()
                else:
                    self.store_object(obj)
            yield (f, tuple(new_args), new_kwargs)

    def assert_objects_unchanged(self, min_count=1) -> None:
        """
        Raise an ``AssertionError`` if stored objects have changed.
        """
        for obj_id, count in self.counts.items():
            if count < min_count:
                continue
            obj = self.objects[obj_id]
            old_hash = self.hashes[obj_id]
            orig_obj = self.old_objects[obj_id]
            if object_hash(obj) != old_hash:
                changed = set()
                if hasattr(orig_obj, "__dict__"):
                    attrs = orig_obj.__dict__.keys()
                else:
                    attrs = set()
                attrs |= set(getattr(orig_obj, "__slots__", ()))
                for attr in attrs:
                    if not attr.startswith("__") and object_hash(
                        getattr(orig_obj, attr)
                    ) != object_hash(getattr(obj, attr)):
                        changed.add(attr)

                expected_changed = get_thread_mutable_attributes(obj)
                if expected_changed and changed <= expected_changed:
                    continue

                raise AssertionError(
                    f"Hash for object of type {type(obj)} has changed, "
                    "suggesting potential thread unsafety. Attrs that "
                    f"changed unexpectedly: {changed - expected_changed}. "
                    f"Attrs whose change was expected: {changed & expected_changed}."
                )


def get_thread_mutable_attributes(obj: Any) -> set[str]:
    """
    Return attributes that are known to be mutated from multiple threads.

    These may be known thread-unsafe code that needs to be fixed, or known
    thread-safe code.
    """

    def _get(klass):
        return (
            getattr(klass, "__sklearn_thread_safe_attributes__", {})
            | getattr(klass, "__sklearn_thread_buggy_attributes__", {})
        ).keys()

    result = set()
    for klass in getmro(type(obj)):
        result |= _get(klass)
    return result


@dataclass
class _ThreadSafe:
    wrapped: object

    def unwrap(self):
        return self.wrapped


_T = TypeVar("_T")


def mark_thread_safe(obj: _T, reason: str = "") -> _T:
    """Mark object passed to ``delayed()`` as thread-safe."""
    from sklearn.utils.parallel import Parallel

    if Parallel._thread_safety_testing:
        # The return type doesn't match the signature, but this is only used in
        # debug mode, which has special handling for it in
        # DelayedChanges.store_parallel_calls().
        return _ThreadSafe(obj)  # type: ignore[return-value]

    return obj


def mark_thread_buggy(obj: _T, reason: str = "") -> _T:
    """
    Mark object passed to ``delayed()`` as thread-buggy.

    This is a separate function in order to aid in grepping.
    """
    return mark_thread_safe(obj)


def object_hash(obj: Any) -> str:
    """
    Hash a Python object; mutating the object should change the result.

    Raises ``TypeError`` in some rare cases.
    """
    return DeepHash(
        obj,
        exclude_paths=get_thread_mutable_attributes(obj),
        ignore_repetition=False,
        ignore_string_type_changes=False,
        ignore_iterable_order=False,
        # Apparently bytes are encoded into strings... this suggests long term
        # might want to recreate this library, or submit fixes.
        encodings=["utf-8", "latin-1"],
    )[obj]

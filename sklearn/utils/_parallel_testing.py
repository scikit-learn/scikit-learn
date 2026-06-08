# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Infrastructure for testing thread-safety."""

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import md5
from inspect import getmro
from io import BytesIO
from pickle import Pickler, PicklingError, dumps
from types import MethodType
from typing import Any, TypeVar

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
        self.old_attributes = {}
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
        self.old_attributes[obj_id] = attrs = {}
        for attr in _get_hopefully_immutable_attrs(obj):
            try:
                attrs[attr] = object_hash(getattr(obj, attr))
            except TypeError:
                continue

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
            if object_hash(obj) != old_hash:
                changed = set()
                for attr, old_attr_hash in self.old_attributes[obj_id].items():
                    if old_attr_hash != object_hash(getattr(obj, attr)):
                        changed.add(attr)

                raise AssertionError(
                    f"Hash for object of type {type(obj)} has changed, "
                    "suggesting potential thread unsafety. Attrs that "
                    f"changed unexpectedly: {changed}"
                )


def _get_thread_mutable_attributes(obj: Any) -> set[str]:
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


def mark_thread_safe(obj: _T, reason: str) -> _T:
    """Mark object passed to ``delayed()`` as thread-safe."""
    # If we import at top-level it's circular import.
    from sklearn.utils.parallel import Parallel

    del reason  # Only passed in for documentation purposes

    if Parallel._thread_safety_testing:
        # The return type doesn't match the signature, but this is only used in
        # free-threading debug mode, which hands it off for special handling
        # and unwrapping in DelayedChanges.store_parallel_calls().
        return _ThreadSafe(obj)  # type: ignore[return-value]

    return obj


def mark_thread_buggy(obj: _T, reason: str) -> _T:
    """
    Mark object passed to ``delayed()`` as thread-buggy.

    This is a separate function in order to aid in grepping.
    """
    return mark_thread_safe(obj, reason)


def _get_hopefully_immutable_attrs(obj: Any) -> set[str]:
    """
    Return attributes that shouldn't change.
    """
    return {
        attr
        for attr in (
            (getattr(obj, "__dict__", {}).keys() | getattr(obj, "__slots__", set()))
            - _get_thread_mutable_attributes(obj)
        )
    }


def _noop(_):
    """Meaningless placeholder for the pickler."""


class _HashingPickler(Pickler):
    """Write-only pickler designed for hashing."""

    # TODO Semantically identical but different objects might be hiding
    # thread-safety, should add id(obj) to the result, perhaps via
    # persistent_id().

    def reducer_override(self, obj: Any) -> Any:
        # If the object isn't pickleable, just hash a repr.
        try:
            _ = dumps(obj)
        except (PicklingError, TypeError):
            return _noop, (repr(obj),)

        # Check if there are any attributes to omit from the hash:
        changing_attrs = _get_thread_mutable_attributes(obj)
        if not changing_attrs:
            # Fall back to normal pickling:

            # TODO Attributes that don't get pickled might be hiding
            # thread-safety!
            return NotImplemented

        # TODO If you have a subclass of a built-in type that also has
        # __sklearn_thread_buggy/safe_attributes__, then this will give the
        # wrong result.
        return _noop, (
            {attr: getattr(obj, attr) for attr in _get_hopefully_immutable_attrs(obj)}
        )


def object_hash(obj: Any) -> bytes:
    """
    Hash a Python object; mutating the object should change the result.

    Raises ``TypeError`` in some rare cases.
    """
    f = BytesIO()
    hasher = _HashingPickler(f)
    try:
        hasher.dump(obj)
    except PicklingError:
        raise TypeError
    return md5(f.getvalue()).digest()

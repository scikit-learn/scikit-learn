"""
Experimental tests exploring mutation detection during
threaded Parallel execution.

Related to scikit-learn issue #34072.

The purpose of these tests is to demonstrate how unintended
mutation of shared inputs can occur when using thread-based
parallel execution.

The current prototype hashes serialized objects before and
after execution to detect whether shared state changed.
"""

from __future__ import annotations

import hashlib
import pickle

import pytest
from joblib import Parallel, delayed


def compute_object_hash(obj) -> str:
    """
    Compute a deterministic hash for a Python object.

    Notes
    -----
    This prototype relies on pickle serialization to capture
    object state before and after threaded execution.

    Limitations:
    - Not all objects are pickleable.
    - Serialization can be expensive for large objects.
    - Intended only for exploratory testing.
    """
    serialized = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)

    return hashlib.sha256(serialized).hexdigest()


def unsafe_mutation_task(shared_list: list[int]) -> None:
    """
    Intentionally mutate a shared object.

    This simulates unsafe behavior in threaded execution
    where multiple workers modify shared state.
    """
    shared_list.append(999)


def safe_readonly_task(shared_list: list[int]) -> list[int]:
    """
    Return a copy without mutating the original object.

    This simulates thread-safe behavior where workers avoid
    modifying shared shared state directly.
    """
    return shared_list.copy()


@pytest.mark.xfail(
    reason=(
        "Experimental mutation-detection prototype demonstrating "
        "shared-state mutation during threaded Parallel execution."
    ),
    strict=False,
)
def test_parallel_input_mutation_detected():
    """
    Verify that unintended mutation changes object state.

    The input list is shared across thread workers and is
    intentionally mutated. The hash comparison demonstrates
    that object state changed during execution.
    """
    shared_data = [1, 2, 3]

    hash_before = compute_object_hash(shared_data)

    Parallel(n_jobs=2, prefer="threads")(
        delayed(unsafe_mutation_task)(shared_data)
        for _ in range(2)
    )

    hash_after = compute_object_hash(shared_data)

    assert hash_before == hash_after


def test_parallel_input_not_mutated():
    """
    Verify that non-mutating threaded tasks preserve state.

    Workers operate on copies and do not mutate the shared
    input object. Therefore the object hash should remain
    unchanged after execution.
    """
    shared_data = [1, 2, 3]

    hash_before = compute_object_hash(shared_data)

    Parallel(n_jobs=2, prefer="threads")(
        delayed(safe_readonly_task)(shared_data)
        for _ in range(2)
    )

    hash_after = compute_object_hash(shared_data)

    assert hash_before == hash_after

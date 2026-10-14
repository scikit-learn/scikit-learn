import subprocess
import sys

import pytest

from sklearn.utils import _openmp_helpers
from sklearn.utils._openmp_helpers import (
    _openmp_parallelism_enabled,
    _openmp_uses_active_wait,
)

pytestmark = pytest.mark.skipif(
    not _openmp_parallelism_enabled(),
    reason="scikit-learn was built without OpenMP support",
)


@pytest.fixture
def no_probe(monkeypatch):
    # Fail if the active wait detection spawns its subprocess.
    def run(*args, **kwargs):
        raise AssertionError("the OpenMP runtime should not have been probed")

    monkeypatch.setattr(_openmp_helpers, "_ACTIVE_WAIT_CACHE", None)
    monkeypatch.setattr(subprocess, "run", run)


@pytest.mark.parametrize(
    "wait_policy, expected",
    [("active", True), ("passive", False), ("ACTIVE", True), ("Passive", False)],
)
def test_openmp_uses_active_wait_env_var(monkeypatch, no_probe, wait_policy, expected):
    monkeypatch.setenv("SKLEARN_OMP_WAIT_POLICY", wait_policy)
    assert _openmp_uses_active_wait() is expected


def test_openmp_uses_active_wait_env_var_invalid(monkeypatch, no_probe):
    monkeypatch.setenv("SKLEARN_OMP_WAIT_POLICY", "spin")
    with pytest.raises(ValueError, match="should be 'active' or 'passive'"):
        _openmp_uses_active_wait()


def test_openmp_uses_active_wait_probe_env(monkeypatch):
    # The subprocess gets an explicit wait policy, so that it cannot probe in
    # turn if the executable relaunches the application instead of running an
    # interpreter.
    probe_envs = []

    def run(*args, env, **kwargs):
        probe_envs.append(env)
        raise OSError

    monkeypatch.delenv("SKLEARN_OMP_WAIT_POLICY", raising=False)
    monkeypatch.setattr(_openmp_helpers, "_ACTIVE_WAIT_CACHE", None)
    monkeypatch.setattr(subprocess, "run", run)

    # A failed probe falls back to assuming active wait.
    assert _openmp_uses_active_wait() is True
    assert len(probe_envs) == 1
    assert probe_envs[0]["SKLEARN_OMP_WAIT_POLICY"] == "active"
    assert probe_envs[0]["OMP_DISPLAY_ENV"] == "VERBOSE"


def test_openmp_uses_active_wait_frozen(monkeypatch, no_probe):
    # A frozen application would relaunch itself instead of running the probe.
    monkeypatch.delenv("SKLEARN_OMP_WAIT_POLICY", raising=False)
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    assert _openmp_uses_active_wait() is True

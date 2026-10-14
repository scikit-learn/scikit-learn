import os
import subprocess
import sys
import warnings

import pytest
from threadpoolctl import threadpool_info

from sklearn.utils import _openmp_helpers
from sklearn.utils._openmp_helpers import (
    _openmp_parallelism_enabled,
    _openmp_uses_active_wait,
)

pytestmark = pytest.mark.skipif(
    not _openmp_parallelism_enabled(),
    reason="scikit-learn was built without OpenMP support",
)


def _unset_wait_policy_env(monkeypatch):
    # Setting the variable first makes monkeypatch restore its original state
    # at teardown, also undoing the value that the detection writes to
    # os.environ after probing.
    monkeypatch.setenv("SKLEARN_OMP_WAIT_POLICY", "")
    monkeypatch.delenv("SKLEARN_OMP_WAIT_POLICY")


@pytest.fixture
def raise_on_subprocess_probe(monkeypatch):
    # Fail if the active wait detection spawns its subprocess.
    def run(*args, **kwargs):
        raise AssertionError("the OpenMP runtime should not have been probed")

    monkeypatch.setattr(_openmp_helpers, "_ACTIVE_WAIT_CACHE", None)
    monkeypatch.setattr(subprocess, "run", run)


@pytest.mark.parametrize(
    "wait_policy, expected",
    [("active", True), ("passive", False), ("ACTIVE", True), ("Passive", False)],
)
def test_openmp_uses_active_wait_env_var(
    monkeypatch, raise_on_subprocess_probe, wait_policy, expected
):
    monkeypatch.setenv("SKLEARN_OMP_WAIT_POLICY", wait_policy)
    assert _openmp_uses_active_wait() is expected


def test_openmp_uses_active_wait_env_var_invalid(
    monkeypatch, raise_on_subprocess_probe
):
    monkeypatch.setenv("SKLEARN_OMP_WAIT_POLICY", "spin")
    with pytest.raises(ValueError, match="should be 'active' or 'passive'"):
        _openmp_uses_active_wait()


def _mock_probe(monkeypatch, result):
    # Replace the probe subprocess: `result` is either an exception to raise or
    # a (returncode, stderr) pair. Returns the list of environments it got.
    probe_envs = []

    def run(args, env, **kwargs):
        probe_envs.append(env)
        if isinstance(result, Exception):
            raise result
        returncode, stderr = result
        return subprocess.CompletedProcess(args, returncode, stdout="", stderr=stderr)

    _unset_wait_policy_env(monkeypatch)
    monkeypatch.setattr(_openmp_helpers, "_ACTIVE_WAIT_CACHE", None)
    monkeypatch.setattr(subprocess, "run", run)
    return probe_envs


def test_openmp_uses_active_wait_sets_env_var(monkeypatch):
    # The detected policy is passed on to subprocesses started later.
    _mock_probe(monkeypatch, (0, "  [host] KMP_BLOCKTIME='0ms'\n"))
    assert _openmp_uses_active_wait() is False
    assert os.environ["SKLEARN_OMP_WAIT_POLICY"] == "passive"


def test_openmp_uses_active_wait_probe_env(monkeypatch):
    # The subprocess gets an explicit wait policy, so that it cannot probe in
    # turn if the executable relaunches the application instead of running an
    # interpreter.
    probe_envs = _mock_probe(monkeypatch, (0, ""))
    _openmp_uses_active_wait()
    assert len(probe_envs) == 1
    assert probe_envs[0]["SKLEARN_OMP_WAIT_POLICY"] == "active"
    assert probe_envs[0]["OMP_DISPLAY_ENV"] == "VERBOSE"


@pytest.mark.parametrize(
    "stderr, expected",
    [
        ("  [host] KMP_BLOCKTIME='0ms'\n", False),
        ("  [host] KMP_BLOCKTIME='200ms'\n", True),
        ("  [host] GOMP_SPINCOUNT = '1'\n", False),
        ("  GOMP_SPINCOUNT = '300000'\n", True),
        # Runtimes that do not support OMP_DISPLAY_ENV, such as MSVC's vcomp.
        ("", True),
    ],
)
def test_openmp_uses_active_wait_probe_output(monkeypatch, stderr, expected):
    _mock_probe(monkeypatch, (0, stderr))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert _openmp_uses_active_wait() is expected


@pytest.mark.parametrize(
    "result, reason",
    [
        (OSError("boom"), "failed with OSError: boom"),
        (subprocess.TimeoutExpired("python", 5), "timed out after 5 s"),
        ((1, "Traceback ...\nImportError: boom\n"), "exited with code 1: ImportError"),
    ],
)
def test_openmp_uses_active_wait_probe_failure(monkeypatch, result, reason):
    # A failed probe falls back to assuming active wait, and says how to skip it.
    _mock_probe(monkeypatch, result)
    with pytest.warns(UserWarning, match="SKLEARN_OMP_WAIT_POLICY") as record:
        assert _openmp_uses_active_wait() is True
    assert reason in str(record[0].message)


def test_openmp_uses_active_wait_frozen(monkeypatch, raise_on_subprocess_probe):
    # A frozen application never attempt to run the subprocess probe to avoid
    # accidentally relaunching itself instead of running the probe.
    _unset_wait_policy_env(monkeypatch)
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    assert _openmp_uses_active_wait() is True


def _probe_openmp_runtime(monkeypatch):
    # Run the actual detection subprocess, failing on any warning except a
    # timeout, which can happen on overloaded CI runners.
    _unset_wait_policy_env(monkeypatch)
    monkeypatch.setattr(_openmp_helpers, "_ACTIVE_WAIT_CACHE", None)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        active_wait = _openmp_uses_active_wait()
    messages = [str(warning.message) for warning in record]
    if any("timed out" in message for message in messages):
        pytest.skip("The OpenMP runtime probe timed out.")
    assert not messages
    return active_wait


def test_openmp_uses_active_wait_probe_integration(monkeypatch):
    # The detection should never warn on our CI under nominal circumstances,
    # otherwise it would likely be too noisy for our users.
    assert isinstance(_probe_openmp_runtime(monkeypatch), bool)


@pytest.mark.parametrize(
    "runtime_env, expected",
    [
        ({"GOMP_SPINCOUNT": "0", "KMP_BLOCKTIME": "0"}, False),
        ({"GOMP_SPINCOUNT": "300000", "KMP_BLOCKTIME": "200"}, True),
    ],
)
def test_openmp_uses_active_wait_probe_integration_policy(
    monkeypatch, runtime_env, expected
):
    # The actual OpenMP runtime reports these settings in the format that the
    # detection parses, whichever of libgomp or libomp is used.
    openmp_prefixes = {
        info["prefix"] for info in threadpool_info() if info["internal_api"] == "openmp"
    }
    if "vcomp" in openmp_prefixes:
        pytest.skip("MSVC's vcomp does not support OMP_DISPLAY_ENV.")
    for name, value in runtime_env.items():
        monkeypatch.setenv(name, value)
    assert _probe_openmp_runtime(monkeypatch) is expected

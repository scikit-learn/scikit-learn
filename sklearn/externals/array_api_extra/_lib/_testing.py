"""
Private testing utilities.

See also ..testing for public testing utilities.
"""

from __future__ import annotations

import pytest


def xfail(
    request: pytest.FixtureRequest, *, reason: str, strict: bool | None = None
) -> None:
    """
    XFAIL the currently running test.

    Unlike ``pytest.xfail``, allow rest of test to execute instead of immediately
    halting it, so that it may result in a XPASS.
    xref https://github.com/pandas-dev/pandas/issues/38902

    Parameters
    ----------
    request : pytest.FixtureRequest
        ``request`` argument of the test function.
    reason : str
        Reason for the expected failure.
    strict : bool, optional
        If True, the test will be marked as failed if it passes.
        If False, the test will be marked as passed if it fails.
        Default: ``xfail_strict`` value in ``pyproject.toml``, or False if absent.
    """
    if strict is not None:
        marker = pytest.mark.xfail(reason=reason, strict=strict)
    else:
        marker = pytest.mark.xfail(reason=reason)
    request.node.add_marker(marker)

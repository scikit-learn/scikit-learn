import os

import pytest

from sklearn.utils._clibs import get_thread_limits, limit_threads_clibs
from sklearn.utils._clibs import _CLibsWrapper


@pytest.mark.parametrize("clib", _CLibsWrapper.SUPPORTED_CLIBS)
def test_limit_threads_clib_dict(clib):
    old_limits = get_thread_limits()

    if old_limits[clib] is not None:
        dynamic_scaling = limit_threads_clibs(limits={clib: 1})
        assert get_thread_limits()[clib] == 1
        assert dynamic_scaling[clib]

        limit_threads_clibs(limits={clib: 3})
        new_limits = get_thread_limits()
        assert new_limits[clib] in (3, os.cpu_count(), os.cpu_count() / 2)

        limit_threads_clibs(limits=old_limits)
        new_limits = get_thread_limits()
        assert new_limits[clib] == old_limits[clib]


@pytest.mark.parametrize("subset", ("all", "blas", "openmp"))
def test_limit_threads_clib_subset(subset):
    if subset == "all":
        clibs = _CLibsWrapper.SUPPORTED_CLIBS.keys()
    elif subset == "blas":
        clibs = ("openblas", "mkl", "mkl_win32")
    elif subset == "openmp":
        clibs = (c for c in _CLibsWrapper.SUPPORTED_CLIBS if "openmp" in c)

    old_limits = get_thread_limits()

    dynamic_scaling = limit_threads_clibs(limits=1, subset=subset)
    new_limits = get_thread_limits()
    for clib in clibs:
        if old_limits[clib] is not None:
            assert new_limits[clib] == 1
            assert dynamic_scaling[clib]

    limit_threads_clibs(limits=3, subset=subset)
    new_limits = get_thread_limits()
    for clib in clibs:
        if old_limits[clib] is not None:
            assert new_limits[clib] in (3, os.cpu_count(), os.cpu_count() / 2)

    limit_threads_clibs(limits=old_limits)
    new_limits = get_thread_limits()
    for clib in clibs:
        if old_limits[clib] is not None:
            assert new_limits[clib] == old_limits[clib]


def test_limit_threads_clib_bad_input():
    with pytest.raises(ValueError,
                       match="subset must be either 'all', 'blas' "
                             "or 'openmp'"):
        limit_threads_clibs(limits=1, subset="wrong")

    with pytest.raises(TypeError,
                       match="limits must either be an int or a dict"):
        limit_threads_clibs(limits=(1, 2, 3))

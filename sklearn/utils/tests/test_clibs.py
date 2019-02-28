import os

import pytest

from sklearn.utils.testing import SkipTest
from sklearn.utils._clibs import (get_thread_limits, _set_thread_limits,
                                  get_openblas_version, thread_pool_limits,
                                  _CLibsWrapper)
from sklearn.utils._clibs_helpers import check_num_threads


SKIP_OPENBLAS = get_openblas_version() is None


def test_openmp_enabled():
    # Check that an OpenMP library is loaded
    limits = get_thread_limits()

    assert not all([lib is None for lib in [limits['openmp_llvm'],
                                            limits['openmp_gnu'],
                                            limits['openmp_win32'],
                                            limits['openmp_intel']]])


@pytest.mark.parametrize("clib", _CLibsWrapper.SUPPORTED_CLIBS)
def test_set_thread_limits_dict(clib):
    # Check that the number of threads used by the multithreaded C-libs can be
    # modified dynamically.

    if clib == "openblas" and SKIP_OPENBLAS:
        raise SkipTest("Possible bug in getting maximum number of threads with"
                       " OpenBLAS < 0.2.16 and OpenBLAS does not expose it's "
                       "version before 0.3.4.")

    old_limits = get_thread_limits()

    if old_limits[clib] is None:
        raise SkipTest("The {} library is not present on this system."
                       .format(clib))

    dynamic_scaling = _set_thread_limits(limits={clib: 1})
    assert get_thread_limits()[clib] == 1
    assert dynamic_scaling[clib]

    thread_pool_limits(limits={clib: 3})
    new_limits = get_thread_limits()
    assert new_limits[clib] in (3, os.cpu_count(), os.cpu_count() / 2)

    thread_pool_limits(limits=old_limits)
    new_limits = get_thread_limits()
    assert new_limits[clib] == old_limits[clib]


@pytest.mark.parametrize("subset", ("all", "blas", "openmp"))
def test_set_thread_limits_subset(subset):
    # Check that the number of threads used by the multithreaded C-libs can be
    # modified dynamically.

    if subset == "all":
        clibs = list(_CLibsWrapper.SUPPORTED_CLIBS.keys())
    elif subset == "blas":
        clibs = ["openblas", "mkl", "mkl_win32"]
    elif subset == "openmp":
        clibs = list(c for c in _CLibsWrapper.SUPPORTED_CLIBS if "openmp" in c)

    if SKIP_OPENBLAS and "openblas" in clibs:
        clibs.remove("openblas")

    old_limits = get_thread_limits()

    dynamic_scaling = _set_thread_limits(limits=1, subset=subset)
    new_limits = get_thread_limits()
    for clib in clibs:
        if old_limits[clib] is not None:
            assert new_limits[clib] == 1
            assert dynamic_scaling[clib]

    thread_pool_limits(limits=3, subset=subset)
    new_limits = get_thread_limits()
    for clib in clibs:
        if old_limits[clib] is not None:
            assert new_limits[clib] in (3, os.cpu_count(), os.cpu_count() / 2)

    thread_pool_limits(limits=old_limits)
    new_limits = get_thread_limits()
    for clib in clibs:
        if old_limits[clib] is not None:
            assert new_limits[clib] == old_limits[clib]


def test_set_thread_limits_bad_input():
    # Check that appropriate errors are raised for invalid arguments

    with pytest.raises(ValueError,
                       match="subset must be either 'all', 'blas' "
                             "or 'openmp'"):
        thread_pool_limits(limits=1, subset="wrong")

    with pytest.raises(TypeError,
                       match="limits must either be an int, a dict"):
        thread_pool_limits(limits=(1, 2, 3))


@pytest.mark.parametrize("subset", (None, "all", "blas", "openmp"))
def test_thread_limit_context(subset):
    # Tests the thread limits context manager

    if subset in [None, "all"]:
        subset_clibs = list(_CLibsWrapper.SUPPORTED_CLIBS.keys())
    elif subset == "blas":
        subset_clibs = ["openblas", "mkl", "mkl_win32"]
    elif subset == "openmp":
        subset_clibs = list(c for c in _CLibsWrapper.SUPPORTED_CLIBS
                            if "openmp" in c)

    old_limits = get_thread_limits()

    with thread_pool_limits(limits=None, subset=subset):
        assert get_thread_limits() == old_limits

    with thread_pool_limits(limits=1, subset=subset):
        limits = get_thread_limits()
        if SKIP_OPENBLAS:
            del limits["openblas"]

        for clib in limits:
            if old_limits[clib] is None:
                assert limits[clib] is None
            elif clib in subset_clibs:
                assert limits[clib] == 1
            else:
                assert limits[clib] == old_limits[clib]

    assert get_thread_limits() == old_limits


def test_openmp_limit_num_threads():
    # checks that OpenMP effectively uses the number of threads requested by
    # the context manager

    old_num_threads = check_num_threads(100)

    for n_threads in [1, 2, 4]:
        with thread_pool_limits(limits=n_threads):
            assert check_num_threads(100) == n_threads

        assert check_num_threads(100) == old_num_threads

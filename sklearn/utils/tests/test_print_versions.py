import pytest

from sklearn.utils._print_versions import get_sys_info
from sklearn.utils._print_versions import get_deps_info
from sklearn.utils._print_versions import show_versions


def test_get_sys_info():
    sys_info = get_sys_info()

    assert 'python' in sys_info
    assert 'OS' in sys_info
    assert 'machine' in sys_info
    assert 'processor' in sys_info
    assert 'LANG' in sys_info
    assert 'LOCALE' in sys_info


def test_get_deps_info():
    deps_info = get_deps_info()

    assert 'pip' in deps_info
    assert 'setuptools' in deps_info
    assert 'numpy' in deps_info
    assert 'scipy' in deps_info
    assert 'Cython' in deps_info
    assert 'pandas' in deps_info
    assert 'matplotlib' in deps_info


def test_show_versions(capsys):
    show_versions()
    out, err = capsys.readouterr()
    assert 'python' in out
    assert 'numpy' in out


def test_show_versions_with_blas(capsys):
    show_versions(with_blas=True)
    out, err = capsys.readouterr()
    assert 'python' in out
    assert 'numpy' in out
    assert 'BLAS' in out

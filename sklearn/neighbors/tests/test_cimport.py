from sklearn.utils._testing import build_and_load_pyx_module


def test_cimport_works():
    pyx_module = "sklearn.neighbors.tests.cimport_test"
    assert build_and_load_pyx_module(pyx_module).all_ok()

from sklearn.utils._testing import build_and_load_pyx_module


def test_cimport_works():
    pyx_module = "sklearn.utils.tests.cimport_test"
    assert build_and_load_pyx_module(pyx_module).all_ok()


def test_cimport_cpp_works():
    pyx_module = "sklearn.utils.tests.cimport_cpp_test"
    assert build_and_load_pyx_module(pyx_module, cpp=True).all_ok()

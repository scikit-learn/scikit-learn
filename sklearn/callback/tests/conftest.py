import sysconfig

from sklearn.utils.fixes import _IS_WASM

_IS_FREE_THREADED = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))


def pytest_ignore_collect(collection_path, config):
    # Skip tests for the whole callack module on free-threaded Python and WASM/Pyodide.
    return _IS_WASM or _IS_FREE_THREADED

import pytest
import sysconfig

from sklearn.utils.fixes import _IS_WASM

_IS_FREE_THREADED = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))


def pytest_collection_modifyitems(config, items):
    if not (_IS_WASM or _IS_FREE_THREADED):
        return

    reason = "callback tests are skipped on free-threaded Python and WASM/Pyodide"
    skip_marker = pytest.mark.skip(reason=reason)
    for item in items:
        item.add_marker(skip_marker)

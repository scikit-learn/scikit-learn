from sklearn.callback.tests._utils import IS_WASM_OR_FREE_THREADED


def pytest_ignore_collect(collection_path, config):
    # Skip tests for the whole callback module on free-threaded Python and WASM/Pyodide.
    return IS_WASM_OR_FREE_THREADED

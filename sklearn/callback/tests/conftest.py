from sklearn.utils.fixes import _IS_WASM


def pytest_ignore_collect(collection_path, config):
    # Skip tests for the whole callback module on WASM/Pyodide.
    return _IS_WASM

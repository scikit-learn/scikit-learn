# Do not collect any tests in externals. This is more robust than using
# --ignore because --ignore needs a path and it is not convenient to pass in
# the externals path (very long install-dependent path in site-packages) when
# using --pyargs
import pytest

PYTEST_GTE_7 = hasattr(pytest, "version_tuple") and pytest.version_tuple >= (
    7,
    0,
)  # type: ignore[attr-defined]

if PYTEST_GTE_7:

    def pytest_ignore_collect(collection_path, config):
        return True

else:

    def pytest_ignore_collect(path, config):
        return True

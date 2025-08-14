import pytest


def test_utils_namespace():
    """Ensure objects are not unintentionally exposed in utils namespace."""
    with pytest.raises(ImportError):
        pass
    with pytest.raises(ImportError):
        pass
    with pytest.raises(ImportError):
        pass

"""
Test the Raetsch loader.
"""

from skdatasets.raetsch import fetch_raetsch


def check(data, shape):
    """Check dataset properties."""
    assert data.data.shape == shape
    assert data.target.shape[0] == shape[0]
    assert hasattr(data.inner_cv, '__iter__')
    assert hasattr(data.outer_cv, '__iter__')


def test_fetch_raetsch_banana():
    """Tests Gunnar Raetsch banana dataset."""
    data = fetch_raetsch('banana')
    check(data, (5300, 2))

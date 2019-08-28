# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause
import numpy as np
from numpy.testing import assert_allclose
import pytest

from sklearn.linear_model._glm.link import (
    IdentityLink,
    LogLink,
    LogitLink,
)


LINK_FUNCTIONS = [IdentityLink, LogLink, LogitLink]


@pytest.mark.parametrize('link', LINK_FUNCTIONS)
def test_link_properties(link):
    """Test link inverse and derivative."""
    rng = np.random.RandomState(42)
    x = rng.rand(100) * 100
    link = link()  # instantiate object
    if isinstance(link, LogitLink):
        # careful for large x, note expit(36) = 1
        # limit max eta to 15
        x = x / 100 * 15
    assert_allclose(link(link.inverse(x)), x)
    # if g(h(x)) = x, then g'(h(x)) = 1/h'(x)
    # g = link, h = link.inverse
    assert_allclose(link.derivative(link.inverse(x)),
                    1 / link.inverse_derivative(x))

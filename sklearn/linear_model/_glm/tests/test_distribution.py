# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#
# License: BSD 3 clause

from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
import pytest

from sklearn.linear_model._glm.distribution import (
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
)


@pytest.mark.parametrize(
    'family, expected',
    [(NormalDistribution(), [True, True, True]),
     (PoissonDistribution(), [False, True, True]),
     (TweedieDistribution(power=1.5), [False, True, True]),
     (GammaDistribution(), [False, False, True]),
     (InverseGaussianDistribution(), [False, False, True]),
     (TweedieDistribution(power=4.5), [False, False, True])])
def test_family_bounds(family, expected):
    """Test the valid range of distributions at -1, 0, 1."""
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, expected)


def test_tweedie_distribution_power():
    with pytest.raises(ValueError, match="no distribution exists"):
        TweedieDistribution(power=0.5)

    with pytest.raises(TypeError, match="must be a real number"):
        TweedieDistribution(power=1j)

    with pytest.raises(TypeError, match="must be a real number"):
        dist = TweedieDistribution()
        dist.power = 1j

    dist = TweedieDistribution()
    assert dist._include_lower_bound is False
    dist.power = 1
    assert dist._include_lower_bound is True


@pytest.mark.parametrize(
    'family, chk_values',
    [(NormalDistribution(), [-1.5, -0.1, 0.1, 2.5]),
     (PoissonDistribution(), [0.1, 1.5]),
     (GammaDistribution(), [0.1, 1.5]),
     (InverseGaussianDistribution(), [0.1, 1.5]),
     (TweedieDistribution(power=-2.5), [0.1, 1.5]),
     (TweedieDistribution(power=-1), [0.1, 1.5]),
     (TweedieDistribution(power=1.5), [0.1, 1.5]),
     (TweedieDistribution(power=2.5), [0.1, 1.5]),
     (TweedieDistribution(power=-4), [0.1, 1.5])])
def test_deviance_zero(family, chk_values):
    """Test deviance(y,y) = 0 for different families."""
    for x in chk_values:
        assert_allclose(family.deviance(x, x), 0, atol=1e-9)

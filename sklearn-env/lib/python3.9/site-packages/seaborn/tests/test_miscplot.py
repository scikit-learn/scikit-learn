import matplotlib.pyplot as plt

from .. import miscplot as misc
from ..palettes import color_palette
from .test_utils import _network


class TestPalPlot:
    """Test the function that visualizes a color palette."""
    def test_palplot_size(self):

        pal4 = color_palette("husl", 4)
        misc.palplot(pal4)
        size4 = plt.gcf().get_size_inches()
        assert tuple(size4) == (4, 1)

        pal5 = color_palette("husl", 5)
        misc.palplot(pal5)
        size5 = plt.gcf().get_size_inches()
        assert tuple(size5) == (5, 1)

        palbig = color_palette("husl", 3)
        misc.palplot(palbig, 2)
        sizebig = plt.gcf().get_size_inches()
        assert tuple(sizebig) == (6, 2)


class TestDogPlot:

    @_network(url="https://github.com/mwaskom/seaborn-data")
    def test_dogplot(self):
        misc.dogplot()
        ax = plt.gca()
        assert len(ax.images) == 1

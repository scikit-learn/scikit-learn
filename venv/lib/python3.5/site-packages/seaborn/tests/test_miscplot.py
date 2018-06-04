import nose.tools as nt
import numpy.testing as npt
import matplotlib.pyplot as plt

from . import PlotTestCase
from .. import miscplot as misc
from seaborn import color_palette


class TestPalPlot(PlotTestCase):
    """Test the function that visualizes a color palette."""
    def test_palplot_size(self):

        pal4 = color_palette("husl", 4)
        misc.palplot(pal4)
        size4 = plt.gcf().get_size_inches()
        nt.assert_equal(tuple(size4), (4, 1))

        pal5 = color_palette("husl", 5)
        misc.palplot(pal5)
        size5 = plt.gcf().get_size_inches()
        nt.assert_equal(tuple(size5), (5, 1))

        palbig = color_palette("husl", 3)
        misc.palplot(palbig, 2)
        sizebig = plt.gcf().get_size_inches()
        nt.assert_equal(tuple(sizebig), (6, 2))

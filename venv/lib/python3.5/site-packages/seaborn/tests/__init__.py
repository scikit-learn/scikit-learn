import numpy as np
from matplotlib.pyplot import close


class PlotTestCase(object):

    def setUp(self):
        np.random.seed(49)

    def tearDown(self):
        close('all')

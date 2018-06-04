from __future__ import absolute_import, division, print_function

from matplotlib.sankey import Sankey


def test_sankey():
    # lets just create a sankey instance and check the code runs
    sankey = Sankey()
    sankey.add()

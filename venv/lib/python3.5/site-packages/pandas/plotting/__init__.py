"""
Plotting api
"""

# flake8: noqa

from pandas.plotting._misc import (scatter_matrix, radviz,
                                   andrews_curves, bootstrap_plot,
                                   parallel_coordinates, lag_plot,
                                   autocorrelation_plot)
from pandas.plotting._core import boxplot
from pandas.plotting._style import plot_params
from pandas.plotting._tools import table
try:
    from pandas.plotting._converter import \
        register as register_matplotlib_converters
    from pandas.plotting._converter import \
        deregister as deregister_matplotlib_converters
except ImportError:
    pass

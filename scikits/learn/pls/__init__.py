"""
Projection to Latent Structures (aka Partial Least Squares) algorithm.

See http://scikit-learn.sourceforge.net/modules/pls.html for complete
documentation.

Author: Aman Thakral <aman.thakral@gmail.com> 

License: New BSD
"""

from .PLS import PLS, statistics_limits
from .PLS import plot_SPE_T2, plot_VIP, plot_WQ, plot_WQbar
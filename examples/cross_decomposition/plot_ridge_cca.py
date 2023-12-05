"""
==========================================
Ridge CCA
==========================================
This example illustrates
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cross_decomposition import RidgeCCA
from sklearn.model_selection import GridSearchCV

model = RidgeCCA(alpha_x=0.1, alpha_y=0.1)

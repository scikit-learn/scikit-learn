"""
Testing for the AverageRegressor module (sklearn.ensemble.AverageRegressor).
"""

# Author: Mohamed Ali Jamaoui m.ali.jamaoui@gmail.com 
# 
# License: BSD 3 clause

from sklearn import datasets
from sklearn.ensemble import AverageRegressor
from sklearn.utils.validation import check_random_state



rng = check_random_state(0)
# load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


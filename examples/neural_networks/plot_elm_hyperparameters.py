"""
===========================================================
Extreme Learning Machines: Effect of tuning hyperparameters
===========================================================

This example first demonstrates how sklearn.grid_search.GridSearchCV can
scan through the space of 2 significant hyperparameters C and weight_scale,
recording their cross-validation scores on the digits dataset.
Please note that n_hidden is also a significant parameter but was not
demonstrated here because of the involved computational costs.

The example then generates a color map to display the scores obtained from
GridSearchCV in a 2D plot. During this experiment, the number of hidden neurons
was set to 1024.

The resulting plot has a diagonal showing a set of high validation scores.
This is because, keeping the rest of the parameters fixed, the overfitting
caused by larger weight_scale is offset by smaller C, which regularizes
the trained output weights.
"""
print(__doc__)

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

from sklearn import grid_search
from sklearn.datasets import load_digits
from sklearn.neural_network import ELMClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

random_state = 0

# Generate sample data
digits = load_digits()
X, y = digits.data, digits.target

# Initialize parameters
parameters = {'elmclassifier__C': np.logspace(-2, 2, 5),
              'elmclassifier__weight_scale': np.logspace(-2, 2, 5)}

# Compute optimal parameters
elm = ELMClassifier(n_hidden=1024, random_state=random_state)

scaler = StandardScaler()
pipeline = make_pipeline(scaler, elm)

clf = grid_search.GridSearchCV(pipeline, parameters)
clf.fit(X, y)

# We extract just the scores
scores = [x[1] for x in clf.grid_scores_]
scores = np.array(scores).reshape(5, 5)

fig, ax = plt.subplots()

# Plot results functions
mat = ax.matshow(scores, origin="lower", cmap=plt.cm.Blues)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(mat, ax=ax, cax=cax)

ax.set_xticks(np.arange(5))
ax.set_yticks(np.arange(5))

ax.set_xlabel('C parameter')
ax.set_ylabel('weight_scale parameter')

ax.set_xticklabels(np.logspace(-2, 2, 5))
ax.xaxis.set_ticks_position('bottom')
ax.set_yticklabels(np.logspace(-2, 2, 5))
ax.set_title('Validation scores when n_hidden=1024', {'fontsize': 12})

plt.show()

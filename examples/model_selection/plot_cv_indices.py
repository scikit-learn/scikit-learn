"""
Visualizing cross-validation behavior in scikit-learn
=====================================================

Choosing the right cross-validation object is a crucial part of fitting a
model properly. There are many ways to split data into training and test
sets in order to avoid model overfitting, to standardize the number of
labels in test sets, etc.

This example visualizes the behavior of several common scikit-learn objects
for comparison.
"""

from sklearn.model_selection import (TimeSeriesSplit, KFold, ShuffleSplit,
                                     GroupShuffleSplit, GroupKFold)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

########################
# Visualize our data
# ------------------
#
# First, let's visualize the structure of our raw data. We have 100 data
# points in total, each one with a label attached to it, encoded by
# color. There are 10 labels in total.
#
# As we'll see, some cross-validation objects do specific things with
# labeled data, while others do not.

# Create our dataset
labels = np.hstack([[ii] * 10 for ii in range(5)])
labels = np.hstack([labels, labels])

# Visualize dataset labels
fig, ax = plt.subplots()
ax.scatter(range(len(labels)),  [.5] * len(labels), c=labels, marker='_',
           lw=50, cmap='rainbow')
ax.set(ylim=[-3, 4], title="Data labels (color = label)",
       xlabel="Sample index")
plt.setp(ax.get_yticklabels() + ax.yaxis.majorTicks, visible=False)


###############################################################################
# Define a function to visualize cross-validation behavior
# --------------------------------------------------------
#
# Next we'll define a function that lets us visualize the behavior of each
# cross-validation object. We'll perform 5 splits of the data. On each
# split, we'll visualize the indices chosen for the training set
# (in blue) and the test set (in red).

def plot_cv_indices(cv, data, ax, n_splits, lw=10):
    """Create a sample plot for indices of a cross-validation object."""

    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(cv.split(data, groups=data)):
        # Fill in indices with the training/test labels
        indices = np.array([np.nan] * len(data))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(range(len(indices)), [ii + .5] * len(indices),
                   c=indices, marker='_', lw=lw, cmap=plt.cm.coolwarm,
                   vmin=-.2, vmax=1.2)

        # Add white bars for group splits
        ixs_splits = np.arange(0, len(data), 10) - .5
        ax.scatter(ixs_splits, [ii + .5] * len(ixs_splits), marker='_',
                   lw=lw, c='w')


    # Plot the data at the end
    ax.scatter(range(len(data)), [ii + 1.5] * len(data),
               c=data, marker='_', lw=lw, cmap=plt.cm.rainbow)

    # Formatting
    yticklabels = list(range(n_splits)) + ['labels']
    ax.set(yticks=np.arange(n_splits+1) + .5, yticklabels=yticklabels,
           xlabel='Sample index', ylabel="CV iteration",
           ylim=[n_splits+1.2, -.2], xlim=[0, 100],
           title='{}'.format(type(cv).__name__))
    return ax


###############################################################################
# Let's see how it looks for the `KFold` cross-validation object:

fig, ax = plt.subplots()
n_splits = 5
cv = KFold(n_splits)
plot_cv_indices(cv, labels, ax, n_splits)

###############################################################################
# Visualize cross-validation indices for many CV objects
# ------------------------------------------------------
#
# Finally, let's visually-compare the cross validation behavior for many
# scikit-learn cross-validation objects. Below we will loop through several
# common cross-validation objects, visualizing the behavior of each.
#
# Note that some keep labels together, while others ignore
# label identity completely. Some have overlapping test sets between CV
# splits, while others do not.

cvs = [ShuffleSplit(n_splits=5), GroupShuffleSplit(n_splits=5),
       KFold(n_splits=5), GroupKFold(n_splits=5), TimeSeriesSplit(n_splits=5)]


fig, axs = plt.subplots(len(cvs), 1, figsize=(6, 3*len(cvs)), sharex=True)
for cv, ax in zip(cvs, axs.ravel()):
    plot_cv_indices(cv, labels, ax, n_splits)

cmap = plt.cm.coolwarm
axs[-1].legend([Patch(color=cmap(.8)), Patch(color=cmap(.2))],
               ['Testing set', 'Training set'], loc=(.7, .8))
plt.setp([ax for ax in axs[1:-1]], xlabel='')
plt.tight_layout()
plt.show()

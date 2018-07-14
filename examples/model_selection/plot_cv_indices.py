"""
Visualizing cross-validation behavior in scikit-learn
=====================================================

Choosing the right cross-validation object is a crucial part of fitting a
model properly. There are many ways to split data into training and test
sets in order to avoid model overfitting, to standardize the number of
groups in test sets, etc.

This example visualizes the behavior of several common scikit-learn objects
for comparison.
"""

from sklearn.model_selection import (TimeSeriesSplit, KFold, ShuffleSplit,
                                     StratifiedKFold, GroupShuffleSplit,
                                     GroupKFold)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

###############################################################################
# Visualize our data
# ------------------
#
# First, we must understand the structure of our data. We'll take a look at
# several datasets. They'll each have 100 randomly generated datapoints,
# but different "labels" that assign the data to groups.
#
# As we'll see, some cross-validation objects do specific things with
# labeled data, while others do not.
#
# To begin, we'll use a dataset in which each datapoint belongs to the same
# group (or, put another way, there are no groups).

n_points = 100
X = np.random.randn(100, 10)
y = np.ones(100)
y[50:] = 0

groups = np.ones(n_points)

def visualize_groups(groups, name):
    # Visualize dataset groups
    fig, ax = plt.subplots()
    ax.scatter(range(len(groups)),  [.5] * len(groups), c=groups, marker='_',
               lw=50, cmap='rainbow')
    ax.set(ylim=[-3, 4], title="Data groups (color = label)\n{}".format(name),
           xlabel="Sample index")
    plt.setp(ax.get_yticklabels() + ax.yaxis.majorTicks, visible=False)

visualize_groups(groups, 'no groups')

###############################################################################
# Define a function to visualize cross-validation behavior
# --------------------------------------------------------
#
# We'll define a function that lets us visualize the behavior of each
# cross-validation object. We'll perform 5 splits of the data. On each
# split, we'll visualize the indices chosen for the training set
# (in blue) and the test set (in red).

def plot_cv_indices(cv, X, y, group, ax, n_splits, lw=10):
    """Create a sample plot for indices of a cross-validation object."""

    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(cv.split(X=X, y=y, groups=group)):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(range(len(indices)), [ii + .5] * len(indices),
                   c=indices, marker='_', lw=lw, cmap=plt.cm.coolwarm,
                   vmin=-.2, vmax=1.2)


    # Plot the data at the end
    ax.scatter(range(len(X)), [ii + 1.5] * len(X),
               c=group, marker='_', lw=lw, cmap=plt.cm.rainbow)

    # Formatting
    yticklabels = list(range(n_splits)) + ['groups']
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
plot_cv_indices(cv, X, y, groups, ax, n_splits)

###############################################################################
# Visualize cross-validation indices for many CV objects
# ------------------------------------------------------
#
# Let's visually compare the cross validation behavior for many
# scikit-learn cross-validation objects. Below we will loop through several
# common cross-validation objects, visualizing the behavior of each.
#
# In this case, there is only a single group, so grouping behavior won't
# really apply.

cvs = [ShuffleSplit(n_splits=5), KFold(n_splits=5),
       TimeSeriesSplit(n_splits=5)]


fig, axs = plt.subplots(len(cvs), 1, figsize=(6, 3*len(cvs)), sharex=True)
for ax, cv in zip(axs, cvs):
    plot_cv_indices(cv, X, y, groups, ax, n_splits)

cmap = plt.cm.coolwarm
axs[-1].legend([Patch(color=cmap(.8)), Patch(color=cmap(.2))],
               ['Testing set', 'Training set'], loc=(.7, .8))
plt.setp([ax for ax in axs[1:-1]], xlabel='')
plt.tight_layout()

###############################################################################
# Using data with balanced groups
# -------------------------------
#
# Next we'll take a look at some data that has several groups, each with the
# same number of members. Here's what the data look like.

groups_even = np.hstack([[ii] * 10 for ii in range(5)])
groups_even = np.hstack([groups_even, groups_even])
y = groups_even
visualize_groups(groups_even, 'balanced groups')

###############################################################################
# We'll visualize these groups with several cross-validation objects that
# are relevant to grouped data.
#
# Note that some keep groups together, while others ignore
# label identity completely. Some have overlapping test sets between CV
# splits, while others do not.


cvs = [ShuffleSplit(n_splits=5), GroupShuffleSplit(n_splits=5),
       KFold(n_splits=5), GroupKFold(n_splits=5), StratifiedKFold(n_splits=5),
       TimeSeriesSplit(n_splits=5)]


fig, axs = plt.subplots(len(cvs), 1, figsize=(6, 3*len(cvs)), sharex=True)
for ax, cv in zip(axs, cvs):
    plot_cv_indices(cv, X, y, groups_even, ax, n_splits)

axs[-1].legend([Patch(color=cmap(.8)), Patch(color=cmap(.2))],
               ['Testing set', 'Training set'], loc=(.7, .8))
plt.setp([ax for ax in axs[1:-1]], xlabel='')
plt.tight_layout()

###############################################################################
# Using data in with imbalanced groups
# ------------------------------------
#
# Finally, let's see how these cross-validation objects behave with imbalanced
# groups.

percentiles = [.05, .1, .15, .2, .5]
groups_imbalanced = np.hstack([[ii] * int(100 * perc)
                               for ii, perc in enumerate(percentiles)])
y = groups_imbalanced
visualize_groups(groups_imbalanced, 'imbalanced groups')

###############################################################################
# We'll visualize these groups with several cross-validation objects that
# are relevant to grouped data with imbalanced groups.
#
# Several scikit-learn CV objects take special consideration to maintain
# the ratio of group membership present in the data.


cvs = [ShuffleSplit(n_splits=5), GroupShuffleSplit(n_splits=5),
       KFold(n_splits=5), GroupKFold(n_splits=5), StratifiedKFold(n_splits=5),
       TimeSeriesSplit(n_splits=5)]


fig, axs = plt.subplots(len(cvs), 1, figsize=(6, 3*len(cvs)), sharex=True)
for ax, cv in zip(axs, cvs):
    plot_cv_indices(cv, X, y, groups_imbalanced, ax, n_splits)

axs[-1].legend([Patch(color=cmap(.8)), Patch(color=cmap(.2))],
               ['Testing set', 'Training set'], loc=(.7, .8))
plt.setp([ax for ax in axs[1:-1]], xlabel='')
plt.tight_layout()
plt.show()

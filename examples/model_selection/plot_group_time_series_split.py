"""
==========================================================
A custom cross-validation splitter for grouped time series
==========================================================

This example shows how to write a custom :term:`CV splitter`. To be scikit-learn
compatible, a splitter must implement a `split` method and a `get_n_splits`
method.

As a motivating use case, we implement a group-aware time series split. However,
group structure can be handled in several ways, depending on the use case. For
instance, groups could be independent of one another while each is internally
time-ordered; multiple groups could instead be recorded simultaneously and share
the same timestamps; or a group's values could depend on another group's
historical data. No single built-in class covers every combination, which is why
it is often better to rely on a custom cross-validation splitter.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# The bike rides dataset
# ----------------------
#
# We use a cyclist's power meter recordings, sampled once per second during four
# separate rides recorded over the span of about a month.

import numpy as np
import pandas as pd

from sklearn.datasets import fetch_file

file_path = fetch_file(
    "https://raw.githubusercontent.com/INRIA/scikit-learn-mooc/"
    "refs/heads/main/datasets/bike_rides.csv"
)
cycling = pd.read_csv(file_path, index_col=0, parse_dates=True)
X, y = cycling.drop(columns="power"), cycling["power"]
cycling.head()

# %%
# Each ride is a natural group. Consecutive measurements inside a ride are
# highly correlated (speed, power output and heart-rate change smoothly from one
# second to the next), but there is no continuity between the end of one ride
# and the start of the next.

unique_ride_dates = np.unique(cycling.index.date)
print(f"There are {len(unique_ride_dates)} bike rides")

# %%
# Rides differ in duration, and each one is separated from the next by a
# different gap. This unequal, time-ordered structure is exactly what motivates
# our custom group- and time-aware splitter.

groups = pd.Series(cycling.index.date, index=cycling.index)
ride_bounds = cycling.index.to_series().groupby(groups).agg(["min", "max"])
ride_bounds["duration"] = ride_bounds["max"] - ride_bounds["min"]
ride_bounds["gap to next ride"] = ride_bounds["min"].shift(-1) - ride_bounds["max"]
ride_bounds[["duration", "gap to next ride"]]

# %%
# Plotting power output against elapsed time, rather than the actual timestamp,
# makes the four rides directly comparable even though they happened on
# different days and lasted different lengths of time.

import matplotlib.pyplot as plt

fig, axes = plt.subplots(
    len(unique_ride_dates), 1, figsize=(9, 8), sharex=True, sharey=True
)
for ax, (ride_date, ride) in zip(axes, cycling.groupby(groups)):
    power_per_minute = ride["power"].resample("1min").mean()
    ax.plot(np.arange(len(power_per_minute)), power_per_minute)
    ax.set_ylabel(str(ride_date))
axes[-1].set_xlabel("elapsed time (minutes)")
fig.supylabel("power (W)")
fig.suptitle("Power output for each of the four rides")
fig.tight_layout()

# %%
# Why not just use an existing splitter?
#
# * :class:`~sklearn.model_selection.TimeSeriesSplit` treats the whole dataset
#   as a single sequence, so it would happily cut a fold right through the
#   middle of a ride, letting samples from the same ride, only seconds apart,
#   end up on both sides of the train/test split.
# * :class:`~sklearn.model_selection.GroupKFold` (or
#   :class:`~sklearn.model_selection.LeaveOneGroupOut`) never splits a ride
#   between train and test, but it is also free to put a later ride in the
#   training set while an earlier ride is held out for testing, which does not
#   reflect how the model would actually be used, i.e. predicting a future ride
#   from past ones.
#
# What we need combines both constraints, that is, keep every ride whole, and
# only ever test on rides that come after the ones used for training. Since
# neither splitter does this alone, we write a small one that delegates to both.

# %%
# Writing a custom `GroupTimeSeriesSplit`
# ---------------------------------------
#
# Here we assume groups can be sorted in chronological order (they are not
# simultaneous).
#
# We reach for :class:`~sklearn.model_selection.LeaveOneGroupOut` to deal with
# the group structure, because it takes no `n_splits` argument. Instead, it
# yields exactly one fold per distinct group, which is the one-block-per-group
# decomposition our custom splitter needs. Then, we let
# :class:`~sklearn.model_selection.TimeSeriesSplit` deal with the chronological
# order.
#
# Notice that `groups` is passed to the constructor of our custom class, instead
# of being set by the method `split(X, y, groups)` like
# :class:`~sklearn.model_selection.GroupKFold` or
# :class:`~sklearn.model_selection.LeaveOneGroupOut` do. The reason is nested
# cross-validation. When a :class:`~sklearn.model_selection.GridSearchCV` is
# used as the inner loop of an outer evaluation loop, its splitter only sees the
# outer loop's training rows, with no group information attached unless it
# travels with the data itself.
#
# .. note::
#    If you are only composing built-in splitters, scikit-learn's metadata
#    routing solves this same nested-CV problem without a custom class

# %%
from sklearn.model_selection import LeaveOneGroupOut, TimeSeriesSplit


class GroupTimeSeriesSplit:
    """Time-ordered K-fold iterator that keeps each group intact.

    Parameters
    ----------
    groups : pandas.Series
        Group label for every sample of the full dataset, indexed
        exactly like the `X` that will later be passed to `split`,
        including any row subset of it, such as a training fold of an
        outer cross-validation loop.

    n_splits : int, default=5
        Number of folds. Must be strictly less than the number of
        groups seen in a given call to `split`.
    """

    def __init__(self, groups, n_splits=5):
        self.groups = groups
        self.n_splits = n_splits

    def split(self, X, y=None, groups=None):
        if groups is not None:
            raise ValueError(
                "GroupTimeSeriesSplit takes 'groups' once, at "
                "construction time, not through split()."
            )
        # Look up the group of each row of *this* X (the full dataset,
        # or any subset of it) by matching on its index.
        groups = self.groups.loc[X.index].to_numpy()

        # One block of sample indices per group, in chronological order.
        indices_per_group = [
            test_idx for _, test_idx in LeaveOneGroupOut().split(X, groups=groups)
        ]

        # We hand the list of per-group index blocks directly to
        # TimeSeriesSplit, one block standing in for one row/group.
        tss = TimeSeriesSplit(n_splits=self.n_splits)
        for train_groups, test_groups in tss.split(indices_per_group):
            train_idx = np.concatenate([indices_per_group[i] for i in train_groups])
            test_idx = np.concatenate([indices_per_group[i] for i in test_groups])
            yield np.sort(train_idx), np.sort(test_idx)

    def get_n_splits(self, X=None, y=None, groups=None):
        return self.n_splits


# %%
# Visualizing `GroupTimeSeriesSplit`
# ----------------------------------
#
# Because `GroupTimeSeriesSplit` always keeps a ride whole, a ride's status for
# a given fold (train, test, or unused) is the same for every one of its
# samples. The plot below reuses the actual power signal, now along the real
# timestamp axis so the multi-week gaps between rides show up as blank space,
# with one row per fold, colored blue where that fold trains, red where it
# tests, and light gray for a ride the fold does not use at all.

n_splits = 3  # at most n_rides - 1 = 3, since we only have 4 rides
cv = GroupTimeSeriesSplit(groups=groups, n_splits=n_splits)

power_per_minute = y.resample("1min").mean()
colors = {"train": "tab:blue", "test": "tab:red", "unused": "lightgray"}

fig, axes = plt.subplots(n_splits, 1, figsize=(9, 8), sharex=True, sharey=True)
for fold, (ax, (train_idx, test_idx)) in enumerate(zip(axes, cv.split(X))):
    status = pd.Series("unused", index=groups.index)
    status.iloc[train_idx] = "train"
    status.iloc[test_idx] = "test"
    status_per_minute = status.resample("1min").first().dropna()

    # One contiguous same-status run at a time, so a splitter that cuts
    # through the middle of a ride shows up as a color change partway.
    new_status = status_per_minute != status_per_minute.shift()
    gap_seconds = status_per_minute.index.to_series().diff().dt.total_seconds()
    run_id = (new_status | (gap_seconds > 60)).cumsum()
    for _, run in status_per_minute.groupby(run_id):
        segment = power_per_minute.loc[run.index]
        ax.plot(segment.index, segment, color=colors[run.iloc[0]])
    ax.set_ylabel(f"fold {fold}")
axes[-1].set_xlabel("date")
fig.supylabel("power (W)")
fig.suptitle("GroupTimeSeriesSplit: train (blue), test (red), unused (gray)")
fig.tight_layout()

# %%
# Each fold's colored segments stay entirely within their own ride's real
# timestamps: a red (test) segment never bleeds into the blank gap that follows
# it, confirming that `GroupTimeSeriesSplit` never splits a ride between train
# and test.
#
# Using the splitter to evaluate a model
# --------------------------------------
#
# `GroupTimeSeriesSplit` behaves like any other scikit-learn CV splitter and
# plugs directly into :func:`~sklearn.model_selection.cross_validate`. Notice
# that we no longer pass `groups=groups` to `cross_validate` either, as the
# splitter already carries it.

from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_validate

model = HistGradientBoostingRegressor(random_state=42)
cv_results = cross_validate(model, X, y, cv=cv, scoring="neg_mean_absolute_error")
print(
    f"MAE: {-cv_results['test_score'].mean():.1f} "
    f"+/- {cv_results['test_score'].std():.1f} Watts"
)

# %%
# Each fold trains on whole, chronologically earlier rides and evaluates on a
# whole, later ride, matching how this model would actually be used in
# production: predicting the power output of a future ride from past ones.

# %%
# Using the splitter in a hyperparameter search
# ---------------------------------------------
#
# Suppose we also want to tune `HistGradientBoostingRegressor`'s
# `min_samples_leaf`. The natural, leakage-free way to do this is to hold out
# the most recent ride as a final test set, and run
# :class:`~sklearn.model_selection.GridSearchCV` on the remaining, earlier
# rides, itself using a group- and time-aware split to choose the best
# hyperparameter.

from sklearn.model_selection import GridSearchCV

train_mask = groups.isin(unique_ride_dates[:-1])
test_mask = groups == unique_ride_dates[-1]
X_train, y_train = X.loc[train_mask], y.loc[train_mask]
X_test, y_test = X.loc[test_mask], y.loc[test_mask]

inner_cv = GroupTimeSeriesSplit(groups=groups, n_splits=2)
search = GridSearchCV(
    model,
    param_grid={"min_samples_leaf": [10, 20, 30, 50]},
    cv=inner_cv,
    scoring="neg_mean_absolute_error",
)
search.fit(X_train, y_train)
print(f"Best min_samples_leaf: {search.best_params_['min_samples_leaf']}")
print(f"Held-out ride MAE: {-search.score(X_test, y_test):.1f} Watts")

# %%
# Here, no `groups=` argument is passed to the grid search either, as `inner_cv`
# already knows how to recover the group of every row of `X_train`, a subset of
# `X`, on its own, exactly as it would if nested inside an outer cross_validate
# loop instead of this single final holdout.
#
# With only four rides total, three is barely enough data to tune a
# hyperparameter on. A larger dataset would let both the outer evaluation and
# the inner search use more folds, and would let a further level of nesting
# evaluate time windows of a new ride.
#
# Wrapping existing splitters like we did in this example already covers many
# custom cross-validation needs. All a class needs is a `split(X, y, groups)`
# generator of `(train_idx, test_idx)` pairs and a matching `get_n_splits`
# method to satisfy scikit-learn's cross-validator API.

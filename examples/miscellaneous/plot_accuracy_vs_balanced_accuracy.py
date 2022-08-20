#!/usr/bin/env python
# coding: utf-8

# Comparison accuracy_score - balanced_accuracy_score
#
# We will look at a binary classification scenario.
# We simulate ground truth data and predictions with a specific degree of alignment.
# In other words, we plant a level of accuracy.
# The goal is to show how the metrics `accuracy` and `balanced_accuracy`
# vary according to such degree of alignmen between ground truth (i.e., y_true)
# and predicted data (i.e., y_pred).
# Also we show how that relationship (metric vs. proportion of alignment)
# varies with different degrees of imbalance in the data, that is, it depends
# on the underlying ground truth data distribution (i.e., proportion of positive/negative labels).


from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import accuracy_score, balanced_accuracy_score

sns.set_context(context="notebook", font_scale=1.2)
sns.set_style(style="whitegrid")


def shuffle_vector(vector):
    """Return a shuffled version of vector"""
    shuffled_idx = np.random.permutation(len(vector))
    return vector[shuffled_idx]


def compute_metrics(y_true, y_pred):
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "balanced_accuracy_adjusted": balanced_accuracy_score(
            y_true, y_pred, adjusted=True
        ),
    }


# We can experiment with these metrics by varying the degree in which`y_pred` aligns with `y_true`.
# In other words, creating data for with planted accuracy, so that we can observe how `accuracy_score` relates to `balanced_accuracy_score` (and its adjusted version).


def create_aligned_y_pred(y_true, alignment_prop=0):
    n_aligned = int(len(y_true) * alignment_prop)
    aligned = y_true[:n_aligned]
    misaligned = y_true[n_aligned:]
    misaligned = shuffle_vector(misaligned)
    return np.concatenate((aligned, misaligned))


def evaluate_preds_for_ground_truth(y_true):
    metrics = defaultdict(list)
    # Pick a range of alignment proportions to evaluate
    alignments = np.arange(0, 1.1, 0.1).round(3)
    for a in alignments:
        metrics["alignment"].append(a)
        y_pred = create_aligned_y_pred(y_true, alignment_prop=a)
        for metric, value in compute_metrics(y_true, y_pred).items():
            metrics[metric].append(value)
    return dict(metrics)


def make_samples(prob_pos, n_samples=100_000):
    return np.random.binomial(1, prob_pos, n_samples)


# Let's fake a target binary variable `y_true` (our ground truth) and predictions `y_pred`.
# In order to make sure we have the exact same distrubution in the fake predictions, we just shuffle the y_true values.
# That way we know for sure that the maximum accuracy possible is 1.
y_true = make_samples(prob_pos=0.5)
y_pred = shuffle_vector(y_true)
assert y_pred.sum() == y_true.sum()

# Now we can repeat that more systematically for different proportions of positive labels in our ground truth:
prob_pos_values = np.arange(0.05, 1.0, 0.05).round(3)
dfs = []
for prob_pos in prob_pos_values:
    y_true = make_samples(prob_pos=prob_pos)
    metrics = evaluate_preds_for_ground_truth(y_true)
    df = pd.DataFrame(metrics)
    df["prob_pos"] = prob_pos
    dfs.append(df)
results = pd.concat(dfs)  # Collect results in a dataframe
results.head(3)

melted = results.melt(
    value_vars=["accuracy", "balanced_accuracy", "balanced_accuracy_adjusted"],
    id_vars=["alignment", "prob_pos"],
    var_name="metric",
)

sns.catplot(
    data=melted,
    y="value",
    x="alignment",
    col="metric",
    hue="prob_pos",
    palette="viridis",
    height=6,
    aspect=1.2,
)
sns.catplot(
    data=melted,
    y="value",
    x="prob_pos",
    col="metric",
    hue="alignment",
    palette="viridis",
    height=6,
    aspect=1.2,
)

# We can observe that both `balanced_accuracy` and `balanced_accuracy_adjusted` have a
# consistent behaviour which only depends on the level of alignment and not on the underlying labels distribution.

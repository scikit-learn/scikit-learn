"""
==================================
Using ensemble of Classifier Chain
==================================
This example shows how using ensemble of :class:`~sklearn.multioutput.ClassifierChain`
can help to achieve higher quality on a multilabel dataset,
compared to using a single instance of a chain.

:class:`~sklearn.multioutput.ClassifierChain` is the metamodel used for multilabel
classification tasks. It fits a base model for each label. The models in chain are
ordered randomly. In addition to input features in the dataset, each model gets the
predictions of the preceding models in the chain as features. These additional features
allow each chain to exploit correlations among the classes.
The :ref:`Jaccard similarity <jaccard_similarity_score>` score for chain tends to be
greater than that of the set independent base models.
"""

# Author: Adam Kleczewski
# License: BSD 3 clause

# %%
# Loading a dataset
# -----------------
# For this example we will use the `yeast
# <https://www.openml.org/d/40597>`_ dataset which contains
# 2417 datapoints each with 103 features and 14 possible labels. Each
# data point has at least one label. As a baseline we first train a logistic
# regression classifier for each of the 14 labels. To evaluate the performance of
# these classifiers we predict on a held-out test set and calculate the
# Jaccard similarity for each sample.

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split

# Load a multi-label dataset from https://www.openml.org/d/40597
X, Y = fetch_openml("yeast", version=4, return_X_y=True, parser="pandas")
Y = Y == "TRUE"
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=0)

# %%
# Fit models
# ----------
# We will fit :class:`~sklearn.linear_model.LogisticRegression` wrapped by
# :class:`~sklearn.multiclass.OneVsRestClassifier` and ensemble of multiple
# :class:`~sklearn.multioutput.ClassifierChain`.
#
# LogisticRegression wrapped by OneVsRestClassifier
# **************************************************
# Since by default :class:`~sklearn.linear_model.LogisticRegression` can't
# handle data with multiple targets, we need to use
# :class:`~sklearn.multiclass.OneVsRestClassifier`.
# After fitting the model we calculate Jaccard similarity.

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import jaccard_score
from sklearn.multiclass import OneVsRestClassifier

base_lr = LogisticRegression(random_state=0)
ovr = OneVsRestClassifier(base_lr)
ovr.fit(X_train, Y_train)
Y_pred_ovr = ovr.predict(X_test)
ovr_jaccard_score = jaccard_score(Y_test, Y_pred_ovr, average="samples")

# %%
# Ensemble of ClassifierChain
# ***************************
# Because the models in each chain are arranged randomly there is significant
# variation in performance among the chains. Presumably there is an optimal
# ordering of the classes in a chain that will yield the best performance.
# However, we do not know that ordering a priori. Instead, we can construct a
# voting ensemble of classifier chains by averaging the binary predictions of
# the chains and apply a threshold of 0.5. The Jaccard similarity score of the
# ensemble is greater than that of the independent models and tends to exceed
# the score of each chain in the ensemble (although this is not guaranteed
# with randomly ordered chains).

from sklearn.multioutput import ClassifierChain

chains = [ClassifierChain(base_lr, order="random", random_state=i) for i in range(10)]
for chain in chains:
    chain.fit(X_train, Y_train)

Y_pred_chains = np.array([chain.predict(X_test) for chain in chains])
chain_jaccard_scores = [
    jaccard_score(Y_test, Y_pred_chain >= 0.5, average="samples")
    for Y_pred_chain in Y_pred_chains
]

Y_pred_ensemble = Y_pred_chains.mean(axis=0)
ensemble_jaccard_score = jaccard_score(
    Y_test, Y_pred_ensemble >= 0.5, average="samples"
)

# %%
# Plot results
# ------------
# Plot the Jaccard similarity scores for the independent model, each of the
# chains, and the ensemble (note that the vertical axis on this plot does
# not begin at 0).

model_scores = [ovr_jaccard_score] + chain_jaccard_scores + [ensemble_jaccard_score]

model_names = (
    "Independent",
    "Chain 1",
    "Chain 2",
    "Chain 3",
    "Chain 4",
    "Chain 5",
    "Chain 6",
    "Chain 7",
    "Chain 8",
    "Chain 9",
    "Chain 10",
    "Ensemble",
)

x_pos = np.arange(len(model_names))

fig, ax = plt.subplots(figsize=(7, 4))
ax.grid(True)
ax.set_title("Classifier Chain Ensemble Performance Comparison")
ax.set_xticks(x_pos)
ax.set_xticklabels(model_names, rotation="vertical")
ax.set_ylabel("Jaccard Similarity Score")
ax.set_ylim([min(model_scores) * 0.9, max(model_scores) * 1.1])
colors = ["r"] + ["b"] * len(chain_jaccard_scores) + ["g"]
ax.bar(x_pos, model_scores, alpha=0.5, color=colors)
plt.tight_layout()
plt.show()

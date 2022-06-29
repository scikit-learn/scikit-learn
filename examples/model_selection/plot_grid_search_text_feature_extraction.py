"""
==========================================================
Sample pipeline for text feature extraction and evaluation
==========================================================

The dataset used in this example is the 20 newsgroups dataset which will be
automatically downloaded, cached and reused for the document
classification example.

In this example we tune the hyperparameters of a particular classifier using a
:class:`~sklearn.model_selection.GridSearchCV`. For a demo on the performance of
some other classifiers, see the
:ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
notebook.

"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Mathieu Blondel <mathieu@mblondel.org>
#         Arturo Amor <david-arturo.amor-quiroz@inria.fr>
# License: BSD 3 clause

# %%
# Data loading
# ------------
# We load two categories from the training set. You can adjust the number of
# categories by giving their names to the dataset loader or setting them to None
# to get the 20 of them.

from sklearn.datasets import fetch_20newsgroups

categories = [
    "alt.atheism",
    "talk.religion.misc",
]

data_train = fetch_20newsgroups(
    subset="train",
    categories=categories,
    shuffle=True,
    random_state=42,
    remove=("headers", "footers", "quotes"),
)

data_test = fetch_20newsgroups(
    subset="test",
    categories=categories,
    shuffle=True,
    random_state=42,
    remove=("headers", "footers", "quotes"),
)

print(f"Loading 20 newsgroups dataset for {len(data_train.target_names)} categories:")
print(data_train.target_names)
print(f"{len(data_train.data)} documents")

# %%
# Pipeline with hyperparameter tuning
# -----------------------------------
# We define a pipeline combining a text feature vectorizer with a simple
# classifier. We also define a set of parameters to use for grid search.
# Uncommenting more parameters will give better exploring power but will
# increase processing time in a combinatorial way.

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.naive_bayes import ComplementNB
from sklearn.pipeline import Pipeline

pipeline = Pipeline(
    [
        ("vect", TfidfVectorizer()),
        ("clf", ComplementNB()),
    ]
)
pipeline

# %%
import numpy as np

parameters = {
    "vect__max_df": (0.5, 0.75, 1.0),
    "vect__min_df": (1, 3, 5),
    # 'vect__max_features': (None, 5000, 10000, 50000),
    "vect__ngram_range": ((1, 1), (1, 2)),  # unigrams or bigrams
    # 'vect__norm': ('l1', 'l2'),
    "clf__alpha": np.logspace(-2.5, -1, 20),
    # (0.01, 0.1),
}

# %%
# We search for the best parameters for both the feature extraction and the
# classifier.

from pprint import pprint
from sklearn.model_selection import RandomizedSearchCV

grid_search = RandomizedSearchCV(
    estimator=pipeline,
    param_distributions=parameters,
    n_iter=40,
    random_state=0,
    n_jobs=2,
    verbose=1,
)

print("Performing grid search...")
print("parameters:")
pprint(parameters)

# %%
from time import time

t0 = time()
grid_search.fit(data_train.data, data_train.target)
print(f"done in {time() - t0:.3f}s")

# %%
print("Best parameters set:")
best_parameters = grid_search.best_estimator_.get_params()
for param_name in sorted(parameters.keys()):
    print(f"{param_name}: {best_parameters[param_name]}")

# %%
test_accuracy = grid_search.score(data_test.data, data_test.target)
print(
    "Accuracy of the best parameters using the inner CV of "
    f"the grid search: {grid_search.best_score_:.3f}"
)
print(f"Accuracy on test set: {test_accuracy:.3f}")

# %%
# The prefixes `vect` and `clf` are required to avoid possible ambiguities in
# the pipeline, but are not necessary for visualizing the results. Because of
# this, we define a function that will rename the tuned hyperparameters and
# improve the readability.

import pandas as pd


def shorten_param(param_name):
    if "__" in param_name:
        return param_name.rsplit("__", 1)[1]
    return param_name


cv_results = pd.DataFrame(grid_search.cv_results_)
cv_results = cv_results.rename(shorten_param, axis=1)
# unigrams are mapped to index 1 and bigrams to index 2
cv_results["ngram_range"] = cv_results["ngram_range"].apply(lambda x: x[1])

# %%
# We can use a `plotly.express.scatter
# <https://plotly.com/python-api-reference/generated/plotly.express.scatter.html>`_
# to visualize the trade-off between scoring time and mean test score. Passing
# the cursor over a given point displays the corresponding parameters.

import plotly.express as px

param_names = [shorten_param(name) for name in parameters.keys()]
fig = px.scatter(
    cv_results, x="mean_score_time", y="mean_test_score", hover_data=param_names
)
fig

# %%
# Notice that the cluster of models in the upper-left corner of the plot are the
# most optimal in terms of accuracy and scoring time. In this case, using
# bigrams increases the required scoring time without improving considerably the
# accuracy of the pipeline.
#
# For more information on how to customize an automated tuning to maximize score
# and minimize scoring time, see the example notebook
# :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_digits.py`.
#
# We can also use a `plotly.express.parallel_coordinates
# <https://plotly.com/python-api-reference/generated/plotly.express.parallel_coordinates.html>`_
# to further visualize the mean test score as a function of the tuned
# hyperparameters. This helps finding interactions between (more than two)
# hyperparameters and provide an intuition on the relevance they have for
# maximizing the performance of a pipeline.

column_results = param_names + ["mean_test_score"]

fig = px.parallel_coordinates(
    cv_results[column_results].apply(dict.fromkeys(column_results, lambda x: x)),
    color="mean_test_score",
    color_continuous_scale=px.colors.sequential.Jet,
)
fig

# %%
# The parallel coordinates plot displays the values of the hyperparameters on
# different columns while the performance metric is color coded. It is possible
# to select a range of results by clicking and holding on any axis of the
# parallel coordinate plot. You can then slide (move) the range selection and
# cross two selections to see the intersections. You can undo a selection by
# clicking once again on the same axis.
#
# In particular for this hyperparameter search, it is interesting to notice that
# the top performing models (mean test score > 0.82) are reached when `min_df=1`
# and `alpha` is close to 0.01, regardless of the values of `max_df` and
# `ngram_range`.
#
# The parallel coordinates plot does not provide information on the variability
# of the score accross the different folds of the cross-validation. For such
# purpose we can further inspect the `cv_results` as follows:

column_results += ["std_test_score"]
cv_results = (
    cv_results[column_results]
    .rename(shorten_param, axis=1)
    .sort_values("mean_test_score", ascending=False)
)
cv_results.reset_index(drop=True)

# %%
# By a manual inspection of the results, one can notice that the top performing
# models overlap within one standard deviation of their test score, showing that
# `max_df` and `ngram_range` are indeed not meaningful parameters in this
# particular case.

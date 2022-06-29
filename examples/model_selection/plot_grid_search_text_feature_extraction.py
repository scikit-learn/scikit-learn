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
parameters = {
    "vect__max_df": (0.5, 0.75, 1.0),
    "vect__min_df": (1, 3, 5),
    # 'vect__max_features': (None, 5000, 10000, 50000),
    "vect__ngram_range": ((1, 1), (1, 2)),  # unigrams or bigrams
    # 'vect__norm': ('l1', 'l2'),
    "clf__alpha": (0.01, 0.1),
}

# %%
# Find the best parameters for both the feature extraction and the classifier.

from pprint import pprint
from time import time
from sklearn.model_selection import GridSearchCV

grid_search = GridSearchCV(
    estimator=pipeline, param_grid=parameters, n_jobs=2, verbose=1
)

print("Performing grid search...")
print("parameters:")
pprint(parameters)
t0 = time()
grid_search.fit(data_train.data, data_train.target)
print(f"done in {time() - t0:.3f}s")
print()
print("Best parameters set:")
best_parameters = grid_search.best_estimator_.get_params()
for param_name in sorted(parameters.keys()):
    print(f"{param_name}: {best_parameters[param_name]}")

accuracy = grid_search.score(data_test.data, data_test.target)
print(
    "Accuracy of the best parameters using the inner CV of "
    f"the grid search: {grid_search.best_score_:.3f}"
)
print(f"Accuracy on test set: {accuracy:.3f}")

# %%
# Finally, we use a `plotly.express.parallel_coordinates
# <https://plotly.com/python-api-reference/generated/plotly.express.parallel_coordinates.html>`_
# to visualize the results from the
# :class:`~sklearn.model_selection.GridSearchCV`.

import pandas as pd
import plotly.express as px


def shorten_param(param_name):
    if "__" in param_name:
        return param_name.rsplit("__", 1)[1]
    return param_name


cv_results = pd.DataFrame(grid_search.cv_results_)
# unigrams are mapped to index 1 and bigrams to index 2
cv_results["param_vect__ngram_range"] = cv_results["param_vect__ngram_range"].apply(
    lambda x: x[1]
)

column_results = [f"param_{name}" for name in parameters.keys()]
column_results += ["mean_test_score"]

fig = px.parallel_coordinates(
    cv_results[column_results]
    .rename(shorten_param, axis=1)
    .apply(dict.fromkeys(map(shorten_param, column_results), lambda x: x)),
    color="mean_test_score",
    color_continuous_scale=px.colors.sequential.Jet,
)
fig

# %%
# The parallel coordinates plot displays the values of the hyperparameters on
# different columns while the performance metric is color coded. It allows us to
# quickly inspect the combinations of hyperparameters that maximize the
# performance.
#
# It is possible to select a range of results by clicking and holding on any
# axis of the parallel coordinate plot. You can then slide (move) the range
# selection and cross two selections to see the intersections. You can undo a
# selection by clicking once again on the same axis.
#
# In particular for this hyperparameter search, it is interesting to notice that
# the top performing models (dark red lines with mean test score > 0.82) are
# reached when `min_df=1` and `alpha=0.01`, regardless of the value of
# `max_df` and the `ngram_range`.

# %%
column_results += ["std_test_score"]
cv_results = (
    cv_results[column_results]
    .rename(shorten_param, axis=1)
    .sort_values("mean_test_score", ascending=False)
)
cv_results

# %%
# By a manual inspection of the results, one can notice that the top performing
# models overlap within one standard deviation of their test score, showing that
# `max_df` and `ngram_range` are indeed not meaningful parameters in this
# particular case. For more information on how to customize a
# :class:`~sklearn.model_selection.GridSearchCV`, see the example notebook
# :ref:`sphx_glr_auto_examples_model_selection_plot_grid_search_digits.py`.

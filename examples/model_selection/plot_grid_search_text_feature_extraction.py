"""
==========================================================
Sample pipeline for text feature extraction and evaluation
==========================================================

The dataset used in this example is the 20 newsgroups dataset which will be
automatically downloaded and then cached and reused for the document
classification example.

You can adjust the number of categories by giving their names to the dataset
loader or setting them to None to get the 20 of them.

"""

# Author: Olivier Grisel <olivier.grisel@ensta.org>
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Mathieu Blondel <mathieu@mblondel.org>
# License: BSD 3 clause

# %%
# Data loading
# ------------
# We load two categories from the training set. Set `categories=None` to do the
# analysis on all the categories.

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
    estimator=pipeline, param_grid=parameters, n_jobs=-1, verbose=1
)

print("Performing grid search...")
print("pipeline:", [name for name, _ in pipeline.steps])
print("parameters:")
pprint(parameters)
t0 = time()
grid_search.fit(data_train.data, data_train.target)
print(f"done in {time() - t0:.3f}s")
print()

print("Best score: {grid_search.best_score_:.3f}")
print("Best parameters set:")
best_parameters = grid_search.best_estimator_.get_params()
for param_name in sorted(parameters.keys()):
    print(f"{param_name}: {best_parameters[param_name]}")

# %%

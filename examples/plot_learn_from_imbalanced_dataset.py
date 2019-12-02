"""
==============================
Learn from imbalanced datasets
==============================

This example illustrates the problem induced by learning on datasets having
imbalanced classes. Subsequently, we compare different approaches alleviating
these negative effects using scikit-learn.

"""

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: MIT

print(__doc__)

###############################################################################
# Problem definition
###############################################################################

from sklearn.datasets import fetch_openml

df, y = fetch_openml('adult', version=2, as_frame=True, return_X_y=True)
# we are dropping the following features:
# - "fnlwgt": this feature was created while studying the "adult" dataset.
#   Thus, we will not use this feature which is not acquired during the survey.
# - "education-num": it is encoding the same information than "education".
#   Thus, we are removing one of these 2 features.
df = df.drop(columns=['fnlwgt', 'education-num'])

###############################################################################
# The "adult" dataset as a class ratio of about 3:1

classes_count = y.value_counts()
classes_count

###############################################################################
# This dataset is only slightly imbalanced. To better highlight the effect of
# learning from an imbalanced dataset, we will increase its ratio to 100:1

import numpy as np
import pandas as pd

rng = np.random.RandomState(0)

# we define a ratio 100:1
n_samples_minority_class = classes_count.max() // 100

mask_minority_class = y == classes_count.idxmin()
indices_minority_class = np.flatnonzero(mask_minority_class)
indices_minority_class_subsampled = rng.choice(
    indices_minority_class, size=n_samples_minority_class, replace=False
)

# sample the dataframe
df_res = pd.concat([df.loc[~mask_minority_class, :],
                    df.loc[indices_minority_class_subsampled, :]])
# sample the target
y_res = pd.concat([y.loc[~mask_minority_class],
                   y.loc[indices_minority_class_subsampled]])
y_res.value_counts()

###############################################################################
# For the rest of the notebook, we will make a single split to get training
# and testing data. Note that you should use cross-validation to have an
# estimate of the performance variation in practice.

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    df_res, y_res, stratify=y_res, random_state=42
)

###############################################################################
# As a baseline, we use a classifier which will always predict the majority
# class independently of the features provided.

from sklearn.dummy import DummyClassifier

dummy_clf = DummyClassifier(strategy="most_frequent")
score = dummy_clf.fit(X_train, y_train).score(X_test, y_test)
print("Accuracy score of a dummy classifier: {:.3f}".format(score))

##############################################################################
# Instead of using the accuracy, we can use the balanced accuracy which will
# take into account the balancing issue.

from sklearn.metrics import balanced_accuracy_score

y_pred = dummy_clf.predict(X_test)
score = balanced_accuracy_score(y_test, y_pred)
print("Balanced accuracy score of a dummy classifier: {:.3f}".format(score))

###############################################################################
# Strategies to learn from an imbalanced dataset
###############################################################################

###############################################################################
# We will first define a helper function which will train a given model
# and compute both accuracy and balanced accuracy. The results will be stored
# in a dataframe


def evaluate_classifier(clf, clf_name=None):
    from sklearn.pipeline import Pipeline
    if clf_name is None:
        if isinstance(clf, Pipeline):
            clf_name = clf[-1].__class__.__name__
        else:
            clf_name = clf.__class__.__name__
    acc = clf.fit(X_train, y_train).score(X_test, y_test)
    y_pred = clf.predict(X_test)
    bal_acc = balanced_accuracy_score(y_test, y_pred)
    clf_score = pd.DataFrame(
        {"Accuracy": acc, "Balanced accuracy": bal_acc},
        index=[clf_name]
    )
    # to avoid passing df_scores and returning it, we make it a global variable
    global df_scores
    df_scores = pd.concat([df_scores, clf_score], axis=0).round(decimals=3)


# Let's define an empty dataframe to store the results
df_scores = pd.DataFrame()

###############################################################################
# Dummy baseline
# ..............
#
# Before to train a real machine learning model, we can store the results
# obtained with our :class:`sklearn.dummy.DummyClassifier`.

evaluate_classifier(dummy_clf)
df_scores

###############################################################################
# Linear classifier baseline
# ..........................
#
# We will create a machine learning pipeline using a
# :class:`sklearn.linear_model.LogisticRegression` classifier. In this regard,
# we will need to one-hot encode the categorical columns and standardized the
# numerical columns before to inject the data into the
# :class:`sklearn.linear_model.LogisticRegression` classifier.
#
# First, we define our numerical and categorical pipelines.

from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import make_pipeline

num_pipe = make_pipeline(
    StandardScaler(), SimpleImputer(strategy="mean", add_indicator=True)
)
cat_pipe = make_pipeline(
    SimpleImputer(strategy="constant", fill_value="missing"),
    OneHotEncoder(handle_unknown="ignore")
)

###############################################################################
# Then, we can create a preprocessor which will dispatch the categorical
# columns to the categorical pipeline and the numerical columns to the
# numerical pipeline

from pandas.api.types import CategoricalDtype
from sklearn.compose import ColumnTransformer
from sklearn.compose import make_column_selector as selector

preprocessor_linear = ColumnTransformer(
    [("num-pipe", num_pipe, selector(dtype_include=np.number)),
     ("cat-pipe", cat_pipe, selector(dtype_include=CategoricalDtype))],
    n_jobs=2
)

###############################################################################
# Finally, we connect our preprocessor with our `LogisticRegression`. We can
# then evaluate our model.

from sklearn.linear_model import LogisticRegression

lr_clf = make_pipeline(
    preprocessor_linear, LogisticRegression(max_iter=1000)
)
evaluate_classifier(lr_clf)
df_scores

###############################################################################
# We can see that our linear model is learning slightly better than our dummy
# baseline. However, it is impacted by the class imbalance.
#
# We can verify that something similar is happening with a tree-based model
# such as :class:`sklearn.ensemble.RandomForestClassifier`. With this type of
# classifier, we will not need to scale the numerical data, and we will only
# need to ordinal encode the categorical data.

from sklearn.preprocessing import OrdinalEncoder
from sklearn.ensemble import RandomForestClassifier

cat_pipe = make_pipeline(
    SimpleImputer(strategy="constant", fill_value="missing"),
    OrdinalEncoder()
)

preprocessor_tree = ColumnTransformer(
    [("num-pipe", num_pipe, selector(dtype_include=np.number)),
     ("cat-pipe", cat_pipe, selector(dtype_include=CategoricalDtype))],
    n_jobs=2
)

rf_clf = make_pipeline(
    preprocessor_tree, RandomForestClassifier(random_state=42, n_jobs=2)
)

evaluate_classifier(rf_clf)
df_scores

###############################################################################
# The :class:`sklearn.ensemble.RandomForestClassifier` is as well affected by
# the class imbalanced, slightly less than the linear model. Now, we will
# present different approach to improve the performance of these 2 models.
#
# Use `class_weight`
# ..................
#
# Most of the models in `scikit-learn` have a parameter `class_weight`. This
# parameter will affect the computation of the loss in linear model or the
# criterion in the tree-based model to penalize differently a false
# classification from the minority and majority class. We can set
# `class_weight="balanced"` such that the weight applied is inversely
# proportional to the class frequency. We test this parametrization in both
# linear model and tree-based model.

lr_clf.set_params(logisticregression__class_weight="balanced")
evaluate_classifier(
    lr_clf, "LogisticRegression with class weight='balanced'"
)
df_scores

###############################################################################
# This weighting strategy is particularly efficient for the logistic
# regression. The balanced accuracy increased significantly.

rf_clf.set_params(randomforestclassifier__class_weight="balanced")
evaluate_classifier(
    rf_clf, "RandomForestClassifier with class weight='balanced'"
)
df_scores

###############################################################################
# However, the same weighting strategy is not efficient with random forest.
# Indeed, the chosen criteria (e.g. entropy) is known to be sensitive to class
# imbalanced.

###############################################################################
# From a random-forest toward a balanced random-forest
# ....................................................
#
# One way to improve the accuracy of the tree-based method is to perform
# some under-sampling such that each tree in the ensemble is learning from a
# balanced set.
#
# The :class:`sklearn.ensemble.RandomForestClassifier` provide an option
# `class_weight="balanced_bootstrap"` such that each tree will learn from a
# bootstrap sample with equal number of samples for each classes. This
# algorithm is also known as a balanced random-forest.

rf_clf.set_params(randomforestclassifier__class_weight="balanced_bootstrap")
evaluate_classifier(
    rf_clf, "RandomForestClassifier with class_weight='balanced_bootstrap'"
)
df_scores

###############################################################################
# We can observe by taking a balanced bootstrap for each tree alleviate the
# overfitting in the random-forest.

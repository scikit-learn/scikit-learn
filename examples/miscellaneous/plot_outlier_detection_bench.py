"""
==========================================
Evaluation of outlier detection estimators
==========================================

This example benchmarks two outlier detection algorithms, namely
:ref:`local_outlier_factor` (LOF) and :ref:`isolation_forest` (IForest), using
ROC curves on classical anomaly detection datasets. The goal is to show that
different algorithms perform well on different datasets.

The algorithm performance is assessed in an outlier detection context:

1. The algorithms are trained on the whole dataset which is assumed to
contain outliers.

2. The ROC curve from :class:`~sklearn.metrics.RocCurveDisplay` is computed
on the same dataset using the knowledge of the labels.

"""

# Author: Pharuj Rajborirug <pharuj.ra@kmitl.ac.th>
#         Arturo Amor <david-arturo.amor-quiroz@inria.fr>
# License: BSD 3 clause

# %%
# Datasets preprocessing and model training
# =========================================
#
# This example uses real-world datasets available in :class:`sklearn.datasets`.
# Due to computational constraints of the scikit-learn documentation, the sample
# size of some datasets is reduced using a stratified
# :class:`~sklearn.model_selection.train_test_split`.
# After the data preprocessing, the datasets' targets will have two classes, 0
# representing inliers and 1 representing outliers.
#
# Different datasets require different preprocessing. One has to take into
# account that tree based models such as
# :class:`~sklearn.ensemble.IsolationForest` can deal with categorical variables
# encoded using an :class:`~sklearn.preprocessing.OrdinalEncoder`, whereas
# neighbors based models such as :class:`~sklearn.neighbors.LocalOutlierFactor`
# require a :class:`~sklearn.preprocessing.OneHotEncoder`.
#
# Neighbors based models may also require scaling of the numerical features (see
# for instance :ref:`neighbors_scaling`). In the presence of outliers, a good
# option is to use a :class:`~sklearn.preprocessing.RobustScaler`.
#
# The following `fit_predict` function returns average outlier score of X.

from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import (
    OneHotEncoder,
    OrdinalEncoder,
    RobustScaler,
)
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
from time import perf_counter


def fit_predict(X, model_name, categorical_columns=(), n_neighbors=20):

    tic = perf_counter()
    if model_name == "LOF":
        preprocessor = ColumnTransformer(
            [("categorical", OneHotEncoder(), categorical_columns)],
            remainder=RobustScaler(),
        )
        clf = LocalOutlierFactor(n_neighbors=n_neighbors, n_jobs=2)
        model = make_pipeline(preprocessor, clf)
        model.fit(X)
        y_pred = model[-1].negative_outlier_factor_

    if model_name == "IForest":
        ordinal_encoder = OrdinalEncoder(
            handle_unknown="use_encoded_value", unknown_value=-1
        )
        clf = IsolationForest(random_state=rng)
        preprocessor = ColumnTransformer(
            [("categorical", ordinal_encoder, categorical_columns)],
            remainder="passthrough",
        )
        model = make_pipeline(preprocessor, clf)
        y_pred = model.fit(X).decision_function(X)
    toc = perf_counter()
    print(f"Duration for {model_name}: {toc - tic:.2f} s")
    return y_pred


# %%
# On the rest of the example we process one dataset per section and summarize
# the results in a final plotting section.
#
# KDDCup99 - SA dataset
# ---------------------
#
# The :ref:`kddcup99_dataset` was generated using a closed network and
# hand-injected attacks. The SA dataset is a subset of it obtained by simply
# selecting all the normal data and an anomaly proportion of around 3%.

# %%
import numpy as np
from sklearn.datasets import fetch_kddcup99
from sklearn.model_selection import train_test_split

rng = np.random.RandomState(42)

X, y = fetch_kddcup99(
    subset="SA", percent10=True, random_state=rng, return_X_y=True, as_frame=True
)
y = (y != b"normal.").astype(np.int32)
X, _, y, _ = train_test_split(X, y, train_size=0.1, stratify=y, random_state=rng)

n_samples, anomaly_frac = X.shape[0], y.mean()
print(f"{n_samples} datapoints with anomaly propotion of {anomaly_frac:.02%}")

# %%
# The SA dataset contains 41 features out of which 3 are categorical:
# "protocol_type", "service" and "flag".

# %%
y_true = {}
y_pred = {"LOF": {}, "IForest": {}}
model_names = ["LOF", "IForest"]
cat_columns = ["protocol_type", "service", "flag"]

y_true["KDDCup99 - SA"] = y
for model_name in model_names:
    y_pred[model_name]["KDDCup99 - SA"] = fit_predict(
        X,
        model_name=model_name,
        categorical_columns=cat_columns,
        n_neighbors=int(n_samples * anomaly_frac),
    )

# %%
# In this example we set `n_neighbors` to match the number of anomalies on the
# dataset. This is a good heuristic as long as the proportion of outliers is not
# very low. Notice that this means the optimal number of neighbors scales with
# the number of samples. Therefore the fit time of LOF models increases linearly
# with the size of the dataset. If one had access to ground truth labels and was
# to additionally tune the number of neighbors, the whole computation would
# result quadratic on `n_samples`.
#
# Forest covertypes dataset
# -------------------------
#
# The :ref:`covtype_dataset` is a multiclass dataset where the target is the
# dominant species of tree in a given patch of forest. It contains 54 features,
# some of wich ("Wilderness_Area" and "Soil_Type") are already binary encoded.
# Though originally meant as a classification task, one can regard inliers as
# samples encoded with label 2 and outliers as those with label 4.

# %%
from sklearn.datasets import fetch_covtype

X, y = fetch_covtype(return_X_y=True, as_frame=True)
s = (y == 2) + (y == 4)
X = X.loc[s]
y = y.loc[s]
y = (y != 2).astype(np.int32)

X, _, y, _ = train_test_split(X, y, train_size=0.05, stratify=y, random_state=rng)
X_forestcover = X  # save X for later use

n_samples, anomaly_frac = X.shape[0], y.mean()
print(f"{n_samples} datapoints with anomaly propotion of {anomaly_frac:.02%}")

# %%
y_true["forestcover"] = y
for model_name in model_names:
    y_pred[model_name]["forestcover"] = fit_predict(
        X,
        model_name=model_name,
        n_neighbors=int(n_samples * anomaly_frac),
    )

# %%
# Breast Cancer dataset
# ---------------------
#
# The :ref:`breast_cancer_dataset` (WDBC) is a binary classification dataset
# where the class is whether a tumor is malignant or benign. It contains 212
# malignant samples, and 357 benign. All the 30 features are continuous and
# positive.

# %%
import pandas as pd
from sklearn.datasets import load_breast_cancer

X, y = load_breast_cancer(return_X_y=True, as_frame=True)
y = np.logical_not(y).astype(np.int32)  # make label 1 to be the minority class

idx_benign = y.index[y == 0]
idx_malign = rng.choice(y.index[y == 1], size=85)  # downsample to 85 points
X = pd.concat([X.iloc[idx_benign], X.iloc[idx_malign]])
y = pd.concat([y.iloc[idx_benign], y.iloc[idx_malign]])

n_samples, anomaly_frac = X.shape[0], y.mean()
print(f"{n_samples} datapoints with anomaly propotion of {anomaly_frac:.02%}")

# %%
y_true["breast cancer"] = y
for model_name in model_names:
    y_pred[model_name]["breast cancer"] = fit_predict(
        X,
        model_name=model_name,
        n_neighbors=int(n_samples * anomaly_frac),
    )

# %%
# Cardiotocography dataset
# ------------------------
#
# The `Cardiotocography dataset <http://www.openml.org/d/1466>`_ is a multiclass
# dataset of fetal cardiotocograms, the classes being the fetal heart rate (FHR)
# pattern encoded with labels from 1 to 10. Here we set class 3 (the minority
# class) to represent the outliers. It contains 30 numerical features.

# %%
from sklearn.datasets import fetch_openml

X, y = fetch_openml(
    name="cardiotocography", version=1, return_X_y=True, as_frame=False, parser="pandas"
)
s = y == "3"
y = s.astype(np.int32)

n_samples, anomaly_frac = X.shape[0], y.mean()
print(f"{n_samples} datapoints with anomaly propotion of {anomaly_frac:.02%}")

# %%
y_true["cardiotocography"] = y
for model_name in model_names:
    y_pred[model_name]["cardiotocography"] = fit_predict(
        X,
        model_name=model_name,
        n_neighbors=int(n_samples * anomaly_frac),
    )

# %%
# Plot and interpret results
# ==========================
#
# The algorithm performance relates to how good the true positive rate (TPR)
# is at low value of the false positive rate (FPR). The best algorithms
# have the curve on the top-left of the plot and the area under curve (AUC)
# close to 1. The diagonal dashed line represents a random classification
# of outliers and inliers.

# %%
import math
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay

# plotting parameters
cols = 2
linewidth = 1
pos_label = 0  # mean 0 belongs to positive class
datasets_names = y_true.keys()
rows = math.ceil(len(datasets_names) / cols)

fig, axs = plt.subplots(nrows=rows, ncols=cols, squeeze=False, figsize=(10, rows * 4))

for i, dataset_name in enumerate(datasets_names):
    for model_idx, model_name in enumerate(model_names):
        display = RocCurveDisplay.from_predictions(
            y_true[dataset_name],
            y_pred[model_name][dataset_name],
            pos_label=pos_label,
            name=model_name,
            linewidth=linewidth,
            ax=axs[i // cols, i % cols],
            plot_chance_level=(model_idx == len(model_names) - 1),
            chance_level_kw={
                "linewidth": linewidth,
                "linestyle": ":",
            },
        )
    axs[i // cols, i % cols].set_title(dataset_name)
_ = plt.tight_layout(pad=2.0)  # spacing between subplots

# %%
# We observe that once the number of neighbors is tuned, LOF and IForest perform
# similarly in terms of ROC AUC for the forestcover and cardiotocography
# datasets. The score for IForest is slightly better for the SA dataset and LOF
# performs considerably better on WDBC than IForest.
#
# Ablation study
# ==============
#
# In this section we explore the impact of the hyperparameter `n_neighbors` and
# the choice of scaling the numerical variables on the LOF model. Here we use
# the :ref:`covtype_dataset` dataset as the binary encoded categories introduce
# a natural scale of euclidean distances between 0 and 1. We then want a scaling
# method to avoid granting a privilege to non-binary features and that is robust
# enough to outliers so that the task of finding them does not become too
# difficult.

# %%
X = X_forestcover
y = y_true["forestcover"]

n_samples = X.shape[0]
n_neighbors_list = (n_samples * np.array([0.02, 0.01, 0.001])).astype(np.int32)
model = make_pipeline(RobustScaler(), LocalOutlierFactor())

fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1], linewidth=linewidth, linestyle=":")
for n_neighbors in n_neighbors_list:
    model.set_params(localoutlierfactor__n_neighbors=n_neighbors)
    model.fit(X)
    y_pred = model[-1].negative_outlier_factor_
    display = RocCurveDisplay.from_predictions(
        y,
        y_pred,
        pos_label=pos_label,
        name=f"n_neighbors = {n_neighbors}",
        linewidth=linewidth,
        ax=ax,
    )
    ax.plot([0, 1], [0, 1], linewidth=linewidth, linestyle=":")
_ = ax.set_title("RobustScaler with varying n_neighbors")

# %%
# We observe that the number of neighbors has a big impact on the performance of
# the model. If one has access to ground truth labels it is then important to
# tune `n_neighbors` accordingly. A convenient way to do so is to explore a
# number of samples proportional to the size of the dataset and possibly over a
# range of values of the order of magnitud of the expected contamination.

# %%
from sklearn.preprocessing import StandardScaler, MinMaxScaler, FunctionTransformer

scaler_list = [RobustScaler(), FunctionTransformer(), StandardScaler(), MinMaxScaler()]
clf = LocalOutlierFactor(n_neighbors=int(0.02 * n_samples))

fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1], linewidth=linewidth, linestyle=":")
for scaler in scaler_list:
    model = make_pipeline(scaler, clf)
    model.fit(X)
    y_pred = model[-1].negative_outlier_factor_
    display = RocCurveDisplay.from_predictions(
        y,
        y_pred,
        pos_label=pos_label,
        name=str(scaler).split("(")[0],
        linewidth=linewidth,
        ax=ax,
    )
    ax.plot([0, 1], [0, 1], linewidth=linewidth, linestyle=":")
ax.set_title("Fixed n_neighbors with varying scaler")
plt.show()

# %%
# On the one hand, :class:`~sklearn.preprocessing.RobustScaler` scales each
# feature independently by using the interquartile range (IQR) by default, which
# is the range between the 25th and 75th percentiles of the data. It centers the
# data by subtracting the median and then scale it by dividing by the IQR. The
# IQR is robust to outliers: the median and interquartile range are less
# affected by extreme values than the range, the mean and the standard
# deviation.
#
# On the other hand, the :class:`~sklearn.preprocessing.MinMaxScaler` scales
# each feature individually such that its range maps into the range between zero
# and one. If there are outliers in the data, they can skew it towards either
# the minimum or maximum values, leading to suboptimal performance. In the case
# of this example, it makes the extreme values present in the continuous
# variables indistinct from the binary encoded categories, which themselves
# remain unchanged by this scaling method.
#
# Because of the above, an IQR scaling more effective for outlier detection than
# a min-max range scaling. We can additionally observe that scaling by the
# standard deviation (as done by :class:`~sklearn.preprocessing.StandardScaler`)
# and not scaling at all (using the
# :class:`~sklearn.preprocessing.FunctionTransformer` with default argument)
# give a relatively good result once the number of neighbors is tuned, but still
# does not perform as well as the :class:`~sklearn.preprocessing.RobustScaler`.

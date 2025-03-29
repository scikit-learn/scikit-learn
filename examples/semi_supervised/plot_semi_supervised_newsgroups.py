"""
================================================
Semi-supervised Classification on a Text Dataset
================================================

This example demonstrates the effectiveness of semi-supervised learning
in text classification when labeled data is scarce.
We compare four different approaches:

1. Supervised learning using 100% of labeled data (baseline)
   - Uses SGDClassifier with TF-IDF features
   - Represents the best possible performance with full supervision

2. Supervised learning using only 20% of labeled data
   - Same model as baseline but with limited training data
   - Shows the performance degradation due to limited labeled data

3. SelfTrainingClassifier (semi-supervised)
   - Uses 20% labeled data + 80% unlabeled data
   - Iteratively predicts labels for unlabeled data
   - Demonstrates how self-training can improve performance

4. LabelSpreading (semi-supervised)
   - Uses 20% labeled data + 80% unlabeled data
   - Propagates labels through the data manifold
   - Shows how graph-based methods can leverage unlabeled data

The example uses the 20 newsgroups dataset, focusing on five computer-related
categories. The results demonstrate how semi-supervised methods can achieve
better performance than supervised learning with limited labeled data by
effectively utilizing unlabeled samples.

You can adjust the number of categories by giving their names to the dataset
loader or setting them to `None` to get all 20 of them.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer
from sklearn.semi_supervised import LabelSpreading, SelfTrainingClassifier

# Loading dataset containing first five categories
data = fetch_20newsgroups(
    subset="train",
    categories=[
        "alt.atheism",
        "comp.graphics",
        "comp.os.ms-windows.misc",
        "comp.sys.ibm.pc.hardware",
        "comp.sys.mac.hardware",
    ],
)

print(
    f"Dataset contains {len(data.filenames)} documents from "
    f"{len(data.target_names)} categories\n"
)

# Parameters
sdg_params = dict(alpha=1e-5, penalty="l2", loss="log_loss")
vectorizer_params = dict(ngram_range=(1, 2), min_df=5, max_df=0.8)

# Supervised Pipeline
pipeline = Pipeline(
    [
        ("vect", CountVectorizer(**vectorizer_params)),
        ("tfidf", TfidfTransformer()),
        ("clf", SGDClassifier(**sdg_params)),
    ]
)
# SelfTraining Pipeline
st_pipeline = Pipeline(
    [
        ("vect", CountVectorizer(**vectorizer_params)),
        ("tfidf", TfidfTransformer()),
        ("clf", SelfTrainingClassifier(SGDClassifier(**sdg_params), verbose=False)),
    ]
)
# LabelSpreading Pipeline
ls_pipeline = Pipeline(
    [
        ("vect", CountVectorizer(**vectorizer_params)),
        ("tfidf", TfidfTransformer()),
        # LabelSpreading does not support dense matrices
        ("toarray", FunctionTransformer(lambda x: x.toarray())),
        ("clf", LabelSpreading()),
    ]
)


def eval_and_get_f1(clf, X_train, y_train, X_test, y_test):
    """Evaluate model performance and return F1 score"""
    print(f"   Number of training samples: {len(X_train)}")
    print(f"   Unlabeled samples in training set: {sum(1 for x in y_train if x == -1)}")
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    f1 = f1_score(y_test, y_pred, average="micro")
    print(f"   Micro-averaged F1 score on test set: {f1:.3f}")
    print("\n")
    return f1


X, y = data.data, data.target
X_train, X_test, y_train, y_test = train_test_split(X, y)

f1_scores = {}

# Evaluate supervised model with 100% of training data
print("1. Supervised SGDClassifier on 100% of the data:")
f1_scores["Supervised (100%)"] = eval_and_get_f1(
    pipeline, X_train, y_train, X_test, y_test
)

# Evaluate supervised model with 20% of training data
print("2. Supervised SGDClassifier on 20% of the training data:")
y_mask = np.random.rand(len(y_train)) < 0.2
# X_20 and y_20 are the subset of the train dataset indicated by the mask
X_20, y_20 = map(list, zip(*((x, y) for x, y, m in zip(X_train, y_train, y_mask) if m)))
f1_scores["Supervised (20%)"] = eval_and_get_f1(pipeline, X_20, y_20, X_test, y_test)

# Evaluate semi-supervised approaches
print(
    "3. SelfTrainingClassifier (semi-supervised) using 20% labeled "
    "+ 80% unlabeled data):"
)
y_train_semi = y_train.copy()
y_train_semi[~y_mask] = -1  # Mark unlabeled data with -1
f1_scores["SelfTraining"] = eval_and_get_f1(
    st_pipeline, X_train, y_train_semi, X_test, y_test
)
print("4. LabelSpreading (semi-supervised) using 20% labeled + 80% unlabeled data:")
f1_scores["LabelSpreading"] = eval_and_get_f1(
    ls_pipeline, X_train, y_train_semi, X_test, y_test
)

# Create the plot
plt.figure(figsize=(10, 6))

models = list(f1_scores.keys())
scores = list(f1_scores.values())

# Define colors for each bar
colors = ["royalblue", "royalblue", "forestgreen", "royalblue"]
bars = plt.bar(models, scores, color=colors)

# Customize plot
plt.title("Comparison of Classification Approaches")
plt.ylabel("Micro-averaged F1 Score")
plt.xticks()

# Add value labels on bars
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2.0,
        height,
        f"{height:.2f}",
        ha="center",
        va="bottom",
    )

# Add note below plot
plt.figtext(
    0.5,
    0.02,
    "SelfTraining classifier shows improved performance over "
    "supervised learning with limited data",
    ha="center",
    va="bottom",
    fontsize=10,
    style="italic",
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.show()

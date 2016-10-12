"""
============================
Classifier Chain
============================

Demonstrate that a single classifier chain
out performs 10 independent classifiers
and that an ensemble of 10 randomly ordered
classifier chains out performs a single classifier chain

"""

import numpy as np
from sklearn.multi_label import ClassifierChain
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import jaccard_similarity_score
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import fetch_rcv1

# Get the Reuters Corpus Volume I dataset
rcv1 = fetch_rcv1()

# RCV1 contains 804414 rows. By convention the first 23149 rows
# are designated as training data
split = 23149
X = rcv1.data[:split, :].toarray()

# RCV1 is a multi-label with 103 labels. To keep this example simple
# we will only try to predict 10 of those labels.
selected_labels = np.random.choice(list(range(rcv1.target.shape[1])), 10)
Y = rcv1.target[:split, selected_labels].toarray()

# But since we are only using a subset of the labels many
# of the rows in our training data will not contain any of
# the 10 labels we are trying to predict.
# So let's only use the rows that contain at least one of the
# labels we are trying to predict.
mask = Y.sum(axis=1) != 0

# Now split the data into training and test datasets.
# There needs to be at least one positive example of
# each label in the training dataset or the train method
# will throw an error. So let's ensure that we get
# a good train/test split.
valid_split_found = False
count = 0
while not valid_split_found:
    X_train, X_test, Y_train, Y_test = \
        train_test_split(X[mask, :], Y[mask, :])
    valid_split_found = all(Y_train.sum(axis=0) > 0)
    if valid_split_found:
        print("valid split found")
    if count >= 10:
        print("valid split not found")
        break
    count += 1


# 1. Fit an independent LogisticRegression model for each class
# using the OneVsRestClassifier wrapper
ovr = OneVsRestClassifier(LogisticRegression())
ovr.fit(X_train, Y_train)
Y_pred_ovr = ovr.predict(X_test)
print("one vs rest jaccard similarity",
      jaccard_similarity_score(Y_test, Y_pred_ovr))

# 2. fit a single chain of LogisticRegression models
classifier_chain = ClassifierChain(LogisticRegression())
classifier_chain.fit(X_train, Y_train)
Y_chain = classifier_chain.predict(X_test)
print("classifier chain jaccard similarity",
      jaccard_similarity_score(Y_test, Y_chain))

# 3. fit an ensemble of classifier chains
chains = [ClassifierChain(LogisticRegression())
          for i in range(10)]
for chain in chains:
    chain.fit(X_train, Y_train)
Y_pred_ensemble = \
    np.array([chain.predict(X_test) for chain in chains]).mean(axis=0)
print("classifier chain ensemble jaccard similarity",
      jaccard_similarity_score(Y_test, Y_pred_ensemble >= .5))


# The three cases should have increasing Jaccard similarity scores.

"""
Name gender classification using symbolic features

This example demonstrates how to use symbolic, rather than numeric, features
for a classification task. It is inspired by, and uses the same dataset as,
the gender classification in Steven Bird, Ewan Klein, Edward Loper (2009),
Natural Language Processing with Python,
http://nltk.googlecode.com/svn/trunk/doc/book/ch06.html.
"""

import numpy as np

from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import zero_one_score
from sklearn.naive_bayes import BernoulliNB
from sklearn.svm import LinearSVC

with open("female.txt") as input:
    female_names = [ln for ln in input if ln[0] != '#']
with open("male.txt") as input:
    male_names = [ln for ln in input if ln[0] != '#']
names = female_names + male_names


def gender_features(name):
    """Extract features from a name for gender identification.

    Returns a dict mapping feature names to booleans.
    """
    x = {"first": name[0],
         "first2": name[:2],
         "first3": name[:3],
         "last": name[-1],
         "last2": name[-2:],
         "last3": name[-3:]}
    for c in "abcdefghijklmnopqrstuvwzyx":
        x["count(%s)" % c] = name.count(c)
    return x


dv = DictVectorizer()
X = dv.fit_transform(gender_features(n) for n in names)
# TODO scale/center X
X = X.tocsr()
print("%d samples, %d features\n" % X.shape)

y = np.array([0] * len(female_names) + [1] * len(male_names))

# Instead of splitting our data into training and test sets,
# we perform 10-fold cross validation.

for clf in (BernoulliNB(), LinearSVC()):
    print("Training and testing %r" % clf)
    for i, (train, test) in enumerate(StratifiedKFold(y, k=10)):
        clf.fit(X[train], y[train])

        y_pred = clf.predict(X[test])
        acc = zero_one_score(y[test], y_pred)
        print("  Fold: %d  Accuracy: %.2f%%" % (i, acc * 100))

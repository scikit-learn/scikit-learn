from collections import defaultdict
from itertools import chain
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import SGDClassifier


def gender_features(name):
    """Extract features from name for gender identification.

    Returns a dict mapping feature names to values, intended as input to
    DictVectorizer.
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


def read_census_file(gender):
    """Parse US census file. Yields (name, gender, frequency) pairs.

    The names will be used for feature extraction; gender becomes part of the
    target vector y; and frequency will be used to generate a sample_weight
    vector.
    """
    filename = "dist.%s.first" % gender

    with open(filename) as f:
        freq = defaultdict(float)
        for ln in f:
            name, freq, _, _ = ln.split()
            yield name, gender, float(freq)


census_data = chain(read_census_file("female"), read_census_file("male"))
names, y, sample_weight = zip(*census_data)

vect = DictVectorizer()
X = vect.fit_transform(gender_features(n) for n in names)

clf = SGDClassifier(verbose=True).fit(X, y)
clf.fit(X, y, sample_weight=sample_weight)

import numpy as np
print(np.mean(clf.predict(X) == np.asarray(y)))

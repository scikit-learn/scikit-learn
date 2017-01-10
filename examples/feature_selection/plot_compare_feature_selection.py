"""
=========================================
Comparison of feature selection functions
=========================================

This example illustrates performance of different univariate feature selection
functions on the text classification task (the 20 newsgroups dataset).

The plot shows the accuracy of a multinomial Naive Bayes classifier as a
function of the amount of the best features selected for training it using five
methods: chi-square, information gain, information gain ratio, F-test and
Kraskov et al's mutual information based on k-nearest neighbor distances.

The script prints to stdout the best 10 features selected using each feature
selection method.
"""

from __future__ import print_function

import time
import numpy as np
import matplotlib.pyplot as plt
from logging import basicConfig, INFO

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.feature_selection import info_gain, info_gain_ratio
from sklearn.feature_selection import mutual_info_classif, f_classif
from sklearn.naive_bayes import MultinomialNB
from sklearn import metrics


print(__doc__)
basicConfig(level=INFO, format='%(asctime)s %(levelname)s %(message)s')


# Load data

print("Loading 20 newsgroups dataset for categories: all")
remove = ('headers', 'footers', 'quotes')
data_train = fetch_20newsgroups(subset='train', categories=None,
                                shuffle=True, random_state=42,
                                remove=remove)
data_test = fetch_20newsgroups(subset='test', categories=None,
                               shuffle=True, random_state=42,
                               remove=remove)
print('data loaded')

# Train-test split

y_train, y_test = data_train.target, data_test.target
categories = data_train.target_names    # for case categories == None

vectorizer = CountVectorizer(max_df=0.5, stop_words='english')
X_train = vectorizer.fit_transform(data_train.data)
X_test = vectorizer.transform(data_test.data)
feature_names = vectorizer.get_feature_names()
cutoffs = [int(x) for x in np.logspace(np.log10(1000.0),
           np.log10(X_train.shape[1]), num=10)]

results = {}

clf = MultinomialNB(alpha=.01)

for func in [chi2, info_gain, info_gain_ratio, f_classif, mutual_info_classif]:

    results[func.__name__] = []

    for k in cutoffs:

        print('_' * 80)
        print("%s, %s, %s features" % (clf, func.__name__, k))

        # apply feature selection

        t0 = time.time()
        selector = SelectKBest(func, k)
        X_train2 = selector.fit_transform(X_train, y_train)
        X_test2 = selector.transform(X_test)
        duration = time.time() - t0

        # keep selected feature names
        feature_names2 = [feature_names[i] for i
                          in selector.get_support(indices=True)]
        feature_names2 = np.asarray(feature_names2)

        # train and evaluate a classifier

        clf.fit(X_train2, y_train)
        pred = clf.predict(X_test2)
        score = metrics.accuracy_score(y_test, pred)
        print("Accuracy:   %0.3f" % score)

        if hasattr(clf, 'coef_') and feature_names is not None:
            print("Top 10 keywords per class:")
            for i, category in enumerate(categories):
                top10 = np.argsort(clf.coef_[i])[-10:]
                print("%s: %s" % (category, " ".join(feature_names2[top10])))

        results[func.__name__].append((score, duration))

# Plot results

f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(12, 8))
ax1.set_title('20 newsgroups dataset')

ax1.set_xlabel('#Features')
ax1.set_ylabel('Accuracy')
ax2.set_ylabel('Time, secs')
colors = "bgrcmyk"
plt.ticklabel_format(useOffset=False)

for i, (name, results) in enumerate(results.items()):
    scores, durations = zip(*results)
    ax1.plot(cutoffs, scores, color=colors[i], label=name)
    ax2.plot(cutoffs, durations, color=colors[i], label=name)

ax1.grid(True)
ax2.grid(True)
ax1.legend(loc='best')
ax2.legend(loc='best')

plt.show()

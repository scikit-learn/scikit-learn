"""
=========================================
Comparison of feature selection functions
=========================================

This example illustrates performance of different feature selection functions
on the text classification task (the 20 newsgroups dataset).

The plot shows the accuracy of a multinomial Naive Bayes classifier as a
function of the amount of the best features selected for training it using five
methods: chi-square (CHI2), information gain (IG), information gain ratio
(IGR), mutual information (MI) and F-test (F).

The script prints to stdout the best 10 features selected using each feature
selection method.
"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from logging import basicConfig, INFO

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import TfidfVectorizer
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

vectorizer = TfidfVectorizer(sublinear_tf=True, max_df=0.5,
                             stop_words='english')
X_train = vectorizer.fit_transform(data_train.data)
X_test = vectorizer.transform(data_test.data)
feature_names = vectorizer.get_feature_names()
k_cutoffs = [int(x) for x in np.logspace(np.log10(1000.0),
             np.log10(X_train.shape[1]), num=10)]

results = {}

clf = MultinomialNB(alpha=.01)

for func, name in [(chi2, "CHI2"), (info_gain, "IG"), (info_gain_ratio, "IGR"),
                   (mutual_info_classif, "MI"), (f_classif, "F")]:

    results[name] = []

    for k in k_cutoffs:

        print('_' * 80)
        print("%s, %s, %s features" % (clf, name, k))

        # apply feature selection

        selector = SelectKBest(func, k)
        X_train2 = selector.fit_transform(X_train, y_train)
        X_test2 = selector.transform(X_test)

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

        results[name].append(score)

# Plot results

plt.figure(figsize=(12, 8))
plt.title("20 newsgroups dataset")
plt.xlabel('#Features')
plt.ylabel('Accuracy')
colors = "bgrcmyk"
plt.ticklabel_format(useOffset=False)

for i, (name, scores) in enumerate(results.items()):
    plt.plot(k_cutoffs, scores, color=colors[i], label=name)

plt.grid(True)
plt.legend(loc='best')

plt.show()

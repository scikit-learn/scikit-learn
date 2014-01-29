#!/usr/bin/env python
"""
================================
Learning to Rank with LambdaMART
================================

LambdaMART is a state-of-the-art algorithm for learn to rank problems.
This example uses the MQ2008 dataset and compares LambdaMART with
linear regression and least squares gradient boosting.

The MQ2008 dataset is available at the following links,
download MQ2008.rar and extract the files to 'data' folder.

https://research.microsoft.com/en-us/um/beijing/projects/letor/letor4download.aspx
https://research.microsoft.com/en-us/um/beijing/projects/letor/LETOR4.0/Data/MQ2008.rar
"""

# Author: Jacques Kvam <jwkvam@gmail.com>
# License: BSD 3 clause

from sklearn.datasets import load_svmlight_file
from sklearn.ensemble import LambdaMART
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from time import time

X_train, y_train, group_train = load_svmlight_file('data/train.txt',
                                                   query_id=True)
X_valid, y_valid, group_valid = load_svmlight_file('data/vali.txt',
                                                   query_id=True)

X_train = X_train.todense()
X_valid = X_valid.todense()

# transform labels to NDCG relevance values
y_train = 2**y_train - 1
y_valid = 2**y_valid - 1

t0 = time()
lm = LambdaMART().fit(X_train, y_train, group=group_train)
print("LambdaMART fit in %fs" % (time() - t0))
gb = GradientBoostingRegressor().fit(X_train, y_train)
lr = LinearRegression().fit(X_train, y_train)
for reg in [lm, gb, lr]:
    print("%s training score is %f" % (reg.__class__.__name__,
          lm.loss_(y_train, reg.predict(X_train)[:, None], group_train)))
    print("%s validation score is %f" % (reg.__class__.__name__,
          lm.loss_(y_valid, reg.predict(X_valid)[:, None], group_valid)))

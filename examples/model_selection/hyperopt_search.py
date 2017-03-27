# Test HyperoptSearchCV

import numpy as np
from sklearn.metrics import make_scorer
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

from sklearn.model_selection import HyperoptSearchCV
from sklearn.model_selection import train_test_split

from hyperopt import hp


###################################
## Regression
print ("\nRegression task...")
from sklearn.datasets import load_boston
from sklearn.linear_model import LinearRegression

x, y = load_boston(return_X_y=True)

train, test = train_test_split(range(x.shape[0]), random_state=1234, train_size=0.7)

x_train = x[train]
y_train = y[train]
x_test = x[test]
y_test = y[test]

pipeline = Pipeline([
    ('pca', PCA()),
    ('model', LinearRegression())
])
space = {
    'pca__n_components': hp.choice('pca__n_components', [1, 2, 3, 4]),
    'model__fit_intercept': hp.choice('model__fit_intercept', [True, False])
}

### Minimize
from sklearn.metrics import median_absolute_error
scorer = make_scorer(median_absolute_error, greater_is_better=False, needs_proba=False, needs_threshold=False)

hs = HyperoptSearchCV(estimator=pipeline,
                      space=space,
                      scoring=scorer,
                      fit_params=None)
hs.fit(x_train, y_train)
for k in hs.report.keys():
    print (k, np.median(hs.report[k]['test_scores']), hs.report[k]['parameters'])
print (hs.best_index_)


print (hs.best_params_)
y_pred = hs.predict(x_test)
print ("MAE (min): {}\n".format(median_absolute_error(y_true=y_test, y_pred=y_pred)))

### Maximize
from sklearn.metrics import explained_variance_score
scorer = make_scorer(explained_variance_score, greater_is_better=True, needs_proba=False, needs_threshold=False)

hs = HyperoptSearchCV(estimator=pipeline,
                      space=space,
                      scoring=scorer,
                      fit_params=None)
hs.fit(x_train, y_train)
for k in hs.report.keys():
    print (k, np.median(hs.report[k]['test_scores']), hs.report[k]['parameters'])
print (hs.best_index_)


print (hs.best_params_)
y_pred = hs.predict(x_test)
print ("Explained variance (max): {}\n".format(explained_variance_score(y_true=y_test, y_pred=y_pred)))


###################################
## Classification
print ("\nClassification task...")
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelBinarizer

x, y = load_iris(return_X_y=True)
lb = LabelBinarizer()
y = lb.fit_transform(y)[:, 1]

train, test = train_test_split(range(x.shape[0]), random_state=1234, train_size=0.7)

x_train = x[train]
y_train = y[train]
x_test = x[test]
y_test = y[test]

pipeline = Pipeline([
    ('pca', PCA()),
    ('model', LogisticRegression())
])
space = {
    'pca__n_components': hp.choice('pca__n_components', [1, 2, 3, 4]),
    'model__fit_intercept': hp.choice('model__fit_intercept', [True, False]),
    'model__penalty': hp.choice('model__penalty',['l1', 'l2']),
    'model__C': hp.uniform('model__C', 0.00001, 100000)
}

### Minimize

#### With probability
from sklearn.metrics import brier_score_loss
scorer = make_scorer(lambda y_true, y_pred: brier_score_loss(y_true, y_pred[:, 1]),
                     greater_is_better=False, needs_proba=True, needs_threshold=False)

hs = HyperoptSearchCV(estimator=pipeline,
                      space=space,
                      scoring=scorer,
                      fit_params=None)
hs.fit(x_train, y_train)
for k in hs.report.keys():
    print (k, np.median(hs.report[k]['test_scores']), hs.report[k]['parameters'])
print (hs.best_index_)


print (hs.best_params_)
y_prob = hs.predict_proba(x_test)
print ("Brier score loss (min): {}\n".format(brier_score_loss(y_true=y_test, y_prob=y_prob[:, 1])))

#### Without probability
from sklearn.metrics import hamming_loss
scorer = make_scorer(hamming_loss, greater_is_better=False, needs_proba=False, needs_threshold=False)

hs = HyperoptSearchCV(estimator=pipeline,
                      space=space,
                      scoring=scorer,
                      fit_params=None)
hs.fit(x_train, y_train)
for k in hs.report.keys():
    print(k, np.median(hs.report[k]['test_scores']), hs.report[k]['parameters'])
print (hs.best_index_)


print (hs.best_params_)
y_pred = hs.predict(x_test)
print ("Hamming loss (min): {}\n".format(hamming_loss(y_true=y_test, y_pred=y_pred)))

### Maximize

#### With probability
from sklearn.metrics import roc_auc_score
scorer = make_scorer(lambda y_true, y_score: roc_auc_score(y_true, y_score),
                     greater_is_better=True, needs_proba=False, needs_threshold=True)

hs = HyperoptSearchCV(estimator=pipeline,
                      space=space,
                      scoring=scorer,
                      fit_params=None)
hs.fit(x_train, y_train)
for k in hs.report.keys():
    print (k, np.median(hs.report[k]['test_scores']), hs.report[k]['parameters'])
print (hs.best_index_)


print (hs.best_params_)
y_prob = hs.predict_proba(x_test)
print ("Roc Auc (max): {}\n".format(roc_auc_score(y_true=y_test, y_score=y_prob[:, 1])))

#### Without probability
from sklearn.metrics import f1_score
scorer = make_scorer(f1_score, greater_is_better=True, needs_proba=False, needs_threshold=False)

hs = HyperoptSearchCV(estimator=pipeline,
                      space=space,
                      scoring=scorer,
                      fit_params=None)
hs.fit(x_train, y_train)
for k in hs.report.keys():
    print (k, np.median(hs.report[k]['test_scores']), hs.report[k]['parameters'])
print (hs.best_index_)


print (hs.best_params_)
y_pred = hs.predict(x_test)
print ("F1_score (max): {}\n".format(f1_score(y_true=y_test, y_pred=y_pred)))


###################################

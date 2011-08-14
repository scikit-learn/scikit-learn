"""
==================================================
Plot Label Propagation versus SVM with iris dataset
==================================================

Performance comparison between Label Propagation in the semi-supervised setting
to SVM in the supervised setting in the iris dataset.

First 9 experiments: SVM (SVM), Label Propagation (LP), Label Spreading (LS)
operate in the "inductive setting". That is, the system is trained with some
percentage of data and then queried against unseen datapoints to infer a label.

The final 10th experiment is in the transductive setting. Using a label
spreading algorithm, the system is trained with approximately 24 percent of the
data labeled and during training, unlabeled points are transductively assigned
values. The test precision, recall, and F1 scores are based on these
transductively assigned labels.
"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import datasets
from scikits.learn import svm
from scikits.learn import label_propagation
from scikits.learn.metrics.metrics import precision_score
from scikits.learn.metrics.metrics import recall_score
from scikits.learn.metrics.metrics import f1_score

iris = datasets.load_iris()

X = iris.data
Y = iris.target

# title for the plots
titles = ['Label Propagation 20% labeled data',
          'Support Vector Machine 80% labeled data']

# 80% data to keep
hold_80 = np.random.rand(len(Y)) < 0.8
train, = np.where(hold_80)

# 20% test data
test, = np.where(hold_80 == False)

X_all = X[train]
Y_all = Y[train]

svc = svm.SVC(kernel='rbf')
svc.fit(X_all, Y_all)
print "Limited Label data example"
print "Test name\tprecision\trecall   \tf1"
print "SVM 80.0pct\t%0.6f\t%0.6f\t%0.6f" %\
        (precision_score(svc.predict(X[test]), Y[test]),\
         recall_score(svc.predict(X[test]), Y[test]),\
         f1_score(svc.predict(X[test]), Y[test]))

print "-------"

for num in [0.2, 0.3, 0.4, 1.0]:
    lp = label_propagation.LabelPropagation()
    hold_new = np.random.rand(len(train)) > num
    train_new, = np.where(hold_new)
    Y_dup = np.copy(Y_all)
    Y_dup[train_new] = -1
    lp.fit(X_all, Y_dup)
    print "LP %0.1fpct\t%0.6f\t%0.6f\t%0.6f" % \
            (80 * num, precision_score(lp.predict(X[test]), Y[test]),\
             recall_score(lp.predict(X[test]), Y[test]),\
             f1_score(lp.predict(X[test]), Y[test]))

# label spreading
for num in [0.2, 0.3, 0.4, 1.0]:
    lspread = label_propagation.LabelSpreading()
    hold_new = np.random.rand(len(train)) > num
    train_new, = np.where(hold_new)
    Y_dup = np.copy(Y_all)
    Y_dup[train_new] = -1
    lspread.fit(X_all, Y_dup)
    print "LS %0.1fpct\t%0.6f\t%0.6f\t%0.6f" % \
            (80 * num, precision_score(lspread.predict(X[test]), Y[test]),\
             recall_score(lspread.predict(X[test]), Y[test]), \
             f1_score(lspread.predict(X[test]), Y[test]))

print "-------"
lspread = label_propagation.LabelSpreading(alpha=0.8)
Y_dup = np.copy(Y)
hold_new = np.random.rand(len(train)) > 0.3
train_new, = np.where(hold_new)
Y_dup = np.copy(Y)
Y_dup[train_new] = -1
lspread.fit(X, Y, suppress_warning=True)
trans_result = np.asarray(lspread.transduction_)
print "LS 20tran\t%0.6f\t%0.6f\t%0.6f" % \
        (precision_score(trans_result[test], Y[test]),
         recall_score(trans_result[test], Y[test]), \
         f1_score(trans_result[test], Y[test]))

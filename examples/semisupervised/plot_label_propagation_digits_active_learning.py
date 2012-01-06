"""
==============================================
Label Propagation digits active learning demo
==============================================

Demonstrates an active learning technique to learn handwritten digits
using label propagation.

We start by training a label propagation model with only 10 labeled points,
then we select the top five most uncertain points to label. Next, we train
with 15 labeled points (original 10 + 5 new ones). We repeat this process
four times to have a model trained with 30 labeled examples.

A plot will appear showing the top 5 most uncertain digits for each iteration
of training. These may or may not contain mistakes, but we will train the next
model with their true labels.

It may take several minutes to finish.
"""
print __doc__

import numpy as np
import pylab as pl

from scipy import stats

from sklearn import datasets
from sklearn import label_propagation

from sklearn.metrics import metrics
from sklearn.metrics.metrics import confusion_matrix

np.random.RandomState(42)

digits = datasets.load_digits()
X = digits.data
Y = digits.target

n_total_samples = len(Y)
n_labeled_points = 10

unlabeled_indicies = np.arange(n_total_samples)[n_labeled_points:]
f = pl.figure()

for i in range(5):
    Y_train = np.copy(Y)
    Y_train[unlabeled_indicies] = -1

    lp_model = label_propagation.LabelSpreading()
    lp_model.fit(X, Y_train, gamma=0.25, max_iters=5)

    y_pred = lp_model.transduction_[unlabeled_indicies]
    y_true = Y[unlabeled_indicies]

    cm = confusion_matrix(y_true, y_pred, labels=lp_model.unq_labels)

    print "Label Spreading model: %d labeled & %d unlabeled (%d total)" %\
        (n_labeled_points, n_total_samples - n_labeled_points, n_total_samples)

    print metrics.classification_report(y_true, y_pred)

    print "Confusion matrix"
    print cm

    pred_entropies = stats.distributions.entropy(lp_model.y_.T)

    arg = zip(pred_entropies, np.arange(n_total_samples))
    arg.sort()
    uncertain_idx = [x[1] for x in arg[-5:]]
    del_inds = np.array([])

    f.text(.05, (1 - (i + 1) * .15), "model %d" % (i + 1))
    for index, im_ind in enumerate(uncertain_idx):
        image = digits.images[im_ind]

        sub = f.add_subplot(5, 5, index + 1 + (5 * i))
        sub.imshow(image, cmap=pl.cm.gray_r)
        sub.set_title('%i' % lp_model.transduction_[im_ind])
        sub.axis('off')

        # labeling 5 points
        del_ind, = np.where(unlabeled_indicies == im_ind)
        del_inds = np.concatenate((del_inds, del_ind))

    unlabeled_indicies = np.delete(unlabeled_indicies, del_inds)
    n_labeled_points += 5

pl.show()

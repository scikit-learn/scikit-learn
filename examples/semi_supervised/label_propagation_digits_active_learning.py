"""
========================================
Label Propagation digits active learning
========================================

Demonstrates an active learning technique to learn handwritten digits
using label propagation.

We start by training a label propagation model with only 10 labeled points,
then we select the top five most uncertain points to label. Next, we train
with 15 labeled points (original 10 + 5 new ones). We repeat this process
four times to have a model trained with 30 labeled examples.

A plot will appear showing the top 5 most uncertain digits for each iteration
of training. These may or may not contain mistakes, but we will train the next
model with their true labels.
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
X = digits.data[:330]
Y = digits.target[:330]

n_total_samples = len(Y)
n_labeled_points = 10

unlabeled_indices = np.arange(n_total_samples)[n_labeled_points:]
f = pl.figure()

for i in range(5):
    Y_train = np.copy(Y)
    Y_train[unlabeled_indices] = -1

    lp_model = label_propagation.LabelSpreading(gamma=0.25, max_iter=5)
    lp_model.fit(X, Y_train)

    y_pred = lp_model.transduction_[unlabeled_indices]
    y_true = Y[unlabeled_indices]

    cm = confusion_matrix(y_true, y_pred, labels=lp_model.unique_labels)

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

    f.text(.05, (1 - (i + 1) * .18), "model %d\n\nfit with\n%d labels" % ((i + 1), 
        i * 5 + 10), size=10)
    for index, im_ind in enumerate(uncertain_idx):
        image = digits.images[im_ind]

        sub = f.add_subplot(5, 5, index + 1 + (5 * i))
        sub.imshow(image, cmap=pl.cm.gray_r)
        sub.set_title('pred: %i\ntrue: %i' % (lp_model.transduction_[im_ind],
            Y[im_ind]), size=10)
        sub.axis('off')

        # labeling 5 points
        del_ind, = np.where(unlabeled_indices == im_ind)
        del_inds = np.concatenate((del_inds, del_ind))

    unlabeled_indices = np.delete(unlabeled_indices, del_inds)
    n_labeled_points += 5

pl.subplots_adjust(0.12, 0.03, 0.9, 0.9, 0.2, 0.45)
pl.show()

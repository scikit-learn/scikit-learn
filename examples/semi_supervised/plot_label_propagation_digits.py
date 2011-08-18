"""
===================================================
Label Propagation digits: Demonstrating performance
===================================================

This example demonstrates the power of semisupervised learning by
training a Label Spreading model to classify handwritten digits
with sets of very few labels.

The handwritten digit dataset has 1797 total points. The model will
be trained using all points, but only 30 will be labeled. Results 
in the form of a confusion matrix and a series of metrics over each 
class will be very good.

At the end, the top 10 most uncertain predictions will be shown.
"""
print __doc__

import numpy as np
import pylab as pl

from scipy import stats

from scikits.learn import datasets
from scikits.learn import label_propagation

from scikits.learn.metrics import metrics
from scikits.learn.metrics.metrics import confusion_matrix

np.random.RandomState(42)

digits = datasets.load_digits()
X = digits.data[:330]
Y = digits.target[:330]

n_total_samples = len(Y)
n_labeled_points = 30

indicies = np.arange(n_total_samples)

unlabeled_set = indicies[n_labeled_points:]

# shuffle everything around
Y_train = np.copy(Y)
Y_train[unlabeled_set] = -1

lp_model = label_propagation.LabelSpreading()
lp_model.fit(X, Y_train, gamma=0.25, max_iters=5)
y_pred = lp_model.transduction_[unlabeled_set]
y_true = Y[unlabeled_set]

cm = confusion_matrix(y_true, y_pred, labels=lp_model.unq_labels)

print "Label Spreading model: %d labeled & %d unlabeled points (%d total)" %\
        (n_labeled_points, n_total_samples - n_labeled_points, n_total_samples)

print metrics.classification_report(y_true, y_pred)

print "Confusion matrix"
print cm

pred_entropies = stats.distributions.entropy(lp_model.y_.T)

arg = zip(pred_entropies, np.arange(n_total_samples))
arg.sort()
uncertain_idx = [x[1] for x in arg[-10:]]

for index, im_ind in enumerate(uncertain_idx):
    image = digits.images[im_ind]

    pl.subplot(2, 5, index + 1)
    pl.imshow(image, cmap=pl.cm.gray_r)
    pl.title('%i' % lp_model.transduction_[im_ind])

pl.show()

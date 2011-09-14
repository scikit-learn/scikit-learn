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

from sklearn import datasets
from sklearn import label_propagation

from sklearn.metrics import metrics
from sklearn.metrics.metrics import confusion_matrix

np.random.RandomState(42)

digits = datasets.load_digits()
X = digits.data[:330]
Y = digits.target[:330]

n_total_samples = len(Y)
n_labeled_points = 30

indices = np.arange(n_total_samples)

unlabeled_set = indices[n_labeled_points:]
f = pl.figure()

# shuffle everything around
Y_train = np.copy(Y)
Y_train[unlabeled_set] = -1

lp_model = label_propagation.LabelSpreading(gamma=0.25, max_iter=5)
lp_model.fit(X, Y_train)
predicted_labels = lp_model.transduction_[unlabeled_set]
true_labels = Y[unlabeled_set]

cm = confusion_matrix(true_labels, predicted_labels,
        labels=lp_model.unique_labels_)

print "Label Spreading model: %d labeled & %d unlabeled points (%d total)" %\
        (n_labeled_points, n_total_samples - n_labeled_points, n_total_samples)

print metrics.classification_report(true_labels, predicted_labels)

print "Confusion matrix"
print cm

pred_entropies = stats.distributions.entropy(lp_model.label_distributions_.T)

arg = zip(pred_entropies, np.arange(n_total_samples))
arg.sort()
uncertainty_index = [x[1] for x in arg[-10:]]

for index, image_index in enumerate(uncertainty_index):
    image = digits.images[image_index]

    sub = f.add_subplot(2, 5, index + 1)
    sub.imshow(image, cmap=pl.cm.gray_r)
    sub.set_title('predict: %i\ntrue: %i' % (
        lp_model.transduction_[image_index], Y[image_index]))

f.suptitle('Learning with small amount of labeled data')
pl.show()

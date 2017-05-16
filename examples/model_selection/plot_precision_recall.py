"""
================
Precision-Recall
================

Example of Precision-Recall metric to evaluate classifier output quality.

In information retrieval, precision is a measure of result relevancy, while
recall is a measure of how many truly relevant results are returned. A high
area under the curve represents both high recall and high precision, where high
precision relates to a low false positive rate, and high recall relates to a
low false negative rate. High scores for both show that the classifier is
returning accurate results (high precision), as well as returning a majority of
all positive results (high recall).

A system with high recall but low precision returns many results, but most of
its predicted labels are incorrect when compared to the training labels. A
system with high precision but low recall is just the opposite, returning very
few results, but most of its predicted labels are correct when compared to the
training labels. An ideal system with high precision and high recall will
return many results, with all results labeled correctly.

Precision (:math:`P`) is defined as the number of true positives (:math:`T_p`)
over the number of true positives plus the number of false positives
(:math:`F_p`).

:math:`P = \\frac{T_p}{T_p+F_p}`

Recall (:math:`R`) is defined as the number of true positives (:math:`T_p`)
over the number of true positives plus the number of false negatives
(:math:`F_n`).

:math:`R = \\frac{T_p}{T_p + F_n}`

These quantities are also related to the (:math:`F_1`) score, which is defined
as the harmonic mean of precision and recall.

:math:`F1 = 2\\frac{P \\times R}{P+R}`

It is important to note that the precision may not decrease with recall. The
definition of precision (:math:`\\frac{T_p}{T_p + F_p}`) shows that lowering
the threshold of a classifier may increase the denominator, by increasing the
number of results returned. If the threshold was previously set too high, the
new results may all be true positives, which will increase precision. If the
previous threshold was about right or too low, further lowering the threshold
will introduce false positives, decreasing precision.

Recall is defined as :math:`\\frac{T_p}{T_p+F_n}`, where :math:`T_p+F_n` does
not depend on the classifier threshold. This means that lowering the classifier
threshold may increase recall, by increasing the number of true positive
results. It is also possible that lowering the threshold may leave recall
unchanged, while the precision fluctuates.

The relationship between recall and precision can be observed in the
stairstep area of the plot - at the edges of these steps a small change
in the threshold considerably reduces precision, with only a minor gain in
recall.

**Average precision** summarizes such a plot as the weighted mean of precisions
achieved at each threshold, with the increase in recall from the previous
threshold used as the weight:

:math:`\\text{AP} = \\sum_n (R_n - R_{n-1}) P_n`

where :math:`P_n` and :math:`R_n` are the precision and recall at the
:math:`n`th threshold. A pair :math:`(R_k, P_k)` is referred to as an
*operating point*.

In *interpolated* average precision, a set of desired recall values is
specified and for each desired value, we average the best precision scores
possible with a recall value at least equal to the target value.
The most common choice is 'eleven point' interpolated precision, where the
desired recall values are [0, 0.1, 0.2, ..., 1.0]. This is the metric used in
`The PASCAL Visual Object Classes (VOC) Challenge <http://citeseerx.ist.psu.edu
/viewdoc/download?doi=10.1.1.157.5766&rep=rep1&type=pdf>`_. In the example
below, the eleven precision values are indicated with an arrow to pointing to
the best precision possible while meeting or exceeding the desired recall.
Note that it's possible that the same operating point might correspond to
multiple desired recall values.

Precision-recall curves are typically used in binary classification to study
the output of a classifier. In order to extend the Precision-recall curve and
average precision to multi-class or multi-label classification, it is necessary
to binarize the output. One curve can be drawn per label, but one can also draw
a precision-recall curve by considering each element of the label indicator
matrix as a binary prediction (micro-averaging).

.. note::

    See also :func:`sklearn.metrics.average_precision_score`,
             :func:`sklearn.metrics.recall_score`,
             :func:`sklearn.metrics.precision_score`,
             :func:`sklearn.metrics.f1_score`
"""
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

from sklearn import svm, datasets
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier

# import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target

# Binarize the output
y = label_binarize(y, classes=[0, 1, 2])
n_classes = y.shape[1]

# Add noisy features
random_state = np.random.RandomState(0)
n_samples, n_features = X.shape
X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]

# Split into training and test
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                    random_state=random_state)

# Run classifier
classifier = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True,
                                 random_state=random_state))
y_score = classifier.fit(X_train, y_train).decision_function(X_test)

# Compute Precision-Recall and plot curve
precision = dict()
recall = dict()
average_precision = dict()
for i in range(n_classes):
    precision[i], recall[i], _ = precision_recall_curve(y_test[:, i],
                                                        y_score[:, i])
    average_precision[i] = average_precision_score(y_test[:, i], y_score[:, i])

# Compute micro-average ROC curve and ROC area
precision["micro"], recall["micro"], _ = precision_recall_curve(y_test.ravel(),
    y_score.ravel())
average_precision["micro"] = average_precision_score(y_test, y_score,
                                                     average="micro")


###############################################################################
# Plot micro-averaged Precision-Recall curve
# ------------------------------------------
#

plt.clf()
plt.step(recall['micro'], precision['micro'], color='b', alpha=0.2, where='post',
         label='Precision-recall curve (area = {:0.3f})'
               ''.format(average_precision['micro']))
plt.fill_between(recall["micro"], precision["micro"], step='post', alpha=0.2,
                 color='b')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('Precision-Recall example: AUC={0:0.2f}'.format(
        average_precision["micro"]))
plt.legend(loc="lower left", prop={'size': 8})
plt.show()


###############################################################################
# Plot Precision-Recall curve for each class
# ------------------------------------------
#
plt.clf()
colors = ['r', 'b', 'g']
for i in range(n_classes):
    plt.step(recall[i], precision[i], color=colors[i], where='post',
             label='Precision-recall curve of class {0} (area = {1:0.2f})'
                   ''.format(i, average_precision[i]))
    plt.fill_between(recall[i], precision[i], step='post', alpha=0.2,
                     color=colors[i])
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Extension of Precision-Recall curve to multi-class')
plt.legend(loc="lower right", prop={'size': 8})
plt.show()


###############################################################################
# Compute eleven-point average precision
# --------------------------------------
#
# Illustrative example of selecting the eleven precision values to use in
# computing eleven-point interpolated average precision for class 2 above.

def pick_eleven_points(recall_, precision_):
    """Choose the eleven operating points that correspond 
    to the best precision for any ``recall >= r`` for r in 
    [0, 0.1, 0.2, ..., 1.0]
    """
    operating_points = list()
    for target_recall in np.arange(0, 1.1, 0.1):
        operating_points_to_consider = [pair for pair in zip(recall_, precision_)
                                        if pair[0] >= target_recall]
        operating_points.append(max(operating_points_to_consider, key=itemgetter(1)))
    return operating_points

iris_cls = 2
plt.clf()
plt.step(recall[iris_cls], precision[iris_cls], color='g', where='post', alpha=0.5,
         linewidth=2,
         label='Precision-recall curve of class {0} (area = {1:0.2f})'
               ''.format(iris_cls, average_precision[iris_cls]))
plt.fill_between(recall[iris_cls], precision[iris_cls], step='post', alpha=0.1,
                 color='g')

eleven_points = pick_eleven_points(recall[iris_cls], precision[iris_cls])
interpolated_average_precision = np.mean([e[1] for e in eleven_points])

print("Target recall    Selected recall   Precision")
for i in range(11):
    plt.annotate('',
                 xy=(eleven_points[i][0], eleven_points[i][1]), xycoords='data',
                 xytext=(i / 10, 0), textcoords='data',
                 arrowprops=dict(arrowstyle="->", alpha=0.7,
                                 connectionstyle="angle3,angleA=90,angleB=45"))
    print("  >= {}           {:3.3f}             {:3.3f}".format(i/10, *eleven_points[i]))

print("  Average:", " "*9, "-", " "*13, "{:3.3f}".format(interpolated_average_precision))

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xticks(np.arange(0, 1.1, 0.1))
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Eleven point Precision Recall for class {}'.format(iris_cls))
plt.legend(loc="upper right", prop={'size': 8})

plt.show()

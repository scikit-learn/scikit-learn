"""
=======================================
PAUC - Partial Area Under Receiver Operating Characteristic (ROC) curve
=======================================

Example of setting False Positive Rate (FPR) or
True Positive Rate (TPR) limitations,
and calculate partial AUC to evaluate binary classifier output quality.

According to the definition of ROC, and as explained in
:ref:`sphx_glr_auto_examples_model_selection_plot_roc.py`,
a larger area under the curve (AUC) is usually better.
However, in real life cases, we sometimes have strict requirements
on keeping low false positive rate,
or keeping high true positive rate, or both of them.
In such scenarios, we are interested in a portion of the ROC curve
which satisfy certain restrictions.

For example, in fraud detection domain,
there are requirements on building a classifier
to distinguish fraud behavior from normal behavior.
To reduce the workload of manually checking the model output on
false positive cases (normal behavior),
we need to keep a low level of FPR.
At the same time, to mitigate the risk of missing fraud cases,
we need to keep a low level of false negative rate (FNR).
Please note that since `TPR = TP/P = (TP)/(TP + FN) = 1 - FNR`
, the limitation on low FNR is
equivalent to the limitation on high TPR.


In :func:`sklearn.metrics.roc_auc_score`,
with the options ``max_fpr`` and ``min_tpr``,
you could define the restrictions to ROC and
get a standardized partial AUC value returned,
which could be used as your model selection criteria.

Let's visualize a ROC example to clarify the definitions.
You could follow below examples and manually calculate
standardized pAUC for a better understanding.


.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :ref:`sphx_glr_auto_examples_model_selection_plot_roc.py`

"""
# %%
# Plot partial ROC area on one example

print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

# Make a fake y_true and y_score
y_true = np.array([0, 0, 1, 1])
y_score = np.array([0.1,  0,  0.1, 0.01])
max_fpr = 0.9
min_tpr = 0.6

# Calculate the fpr, tpr list
fpr, tpr, thresholds = metrics.roc_curve(y_true=y_true,
                                         y_score=y_score,
                                         pos_label=1)
spauc = metrics.roc_auc_score(y_true=y_true,
                              y_score=y_score,
                              max_fpr=max_fpr,
                              min_tpr=min_tpr)
print(f'The fpr list is: {fpr}')
print(f'The tpr list is: {tpr}')
print(f'With max_fpr {max_fpr} and min_tpr '
      f'{min_tpr}, the standardized partial AUC is {spauc}.')
plt.rcParams["figure.figsize"] = (10, 10)

# Plot the min_tpr horizontal line
if min_tpr is not None:
    plt.axhline(y=min_tpr, color='g', linestyle='--',
                alpha=0.5, label='min_tpr')

# Plot the max_fpr vertical line
if max_fpr is not None:
    plt.axvline(x=max_fpr, color='r', linestyle='--',
                alpha=0.5, label='max_fpr')

plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve')

# Plot the non-discriminant ROC line for reference:
plt.plot([0, 1], [0, 1], color='navy', lw=2,
         linestyle='--', label='non-discriminant ROC')

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.grid()
plt.gca().set_aspect('equal', adjustable='box')

if min_tpr is None:
    min_tpr = 0.0
if max_fpr is None:
    max_fpr = 1.0
# Fill the region which conforms to max_fpr with red:
plt.fill_betweenx(tpr, x1=max_fpr, x2=fpr, where=(fpr <= max_fpr),
                  alpha=0.1, interpolate=True, color='red')
# Fill the region which conforms to min_tpr with green:
plt.fill_between(fpr, y1=min_tpr, y2=tpr, where=(tpr >= min_tpr),
                 alpha=0.1, interpolate=True, color='green')
plt.show()

# %%
# Manual calculation on standardized pAUC for the above example
# ..........................................
# The partial AUC is the area of the region formed by "max_fpr",
# "min_tpr", and "ROC curve" lines. The pAUC value ``A = 0.4 * 0.4 = 0.16``.
#
# According to the Formula 7 in [1]_,
# the standardized pAUC is
#
# .. math::
#    \frac{1}{2}\left ( 1 + \frac{A - min}{max - min} \right )
#
# where the max/min means the maximum/minimum value of AUC
# under the `max_fpr` and `min_tpr` restrictions.
#
# Following the ideas proposed by [1]_ (McClish in 1989),
# under ``max_fpr`` and ``min_tpr`` restrictions,
# we know the maximum area could attain is the rectangle with 4 endpoints:
# ``(0, min_tpr), (0, 1), (max_fpr, min_tpr), (max_fpr, 1)`` .
#
# Thus, the max area is ``(1 - min_tpr) * max_fpr = 0.4 * 0.9 = 0.36``.
#
# The minimum area could attain is the triangle with 3 endpoints:
# ``(min_tpr, min_tpr), (max_fpr, min_tpr), and (max_fpr, max_fpr)``.
# Thus, the min area is
# ``0.5 * (max_fpr - min_tpr) * (max_fpr - min_tpr) = 0.045``.
#
# Finally, we get the standardized pAUC as
#
# .. math::
#    spauc =
#    \frac{1}{2}\left ( 1 + \frac{0.16 - 0.045}{0.36 - 0.045} \right )
#    \approx  0.6825
#

# %%
# .. topic:: References
#
#    .. [1] `Analyzing a portion of the ROC curve. McClish, 1989
#             <https://www.ncbi.nlm.nih.gov/pubmed/2668680>`_

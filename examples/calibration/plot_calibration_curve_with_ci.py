"""
========================================================
Probability calibration curves with confidence intervals
========================================================

In practice, it's important to emphasize where along the range of predicted
probabilities that a model's calibration cannot be confidently assessed.
Plotting a calibration curve with confidence intervals is one way of
quantifying how uncertain one should be in assessing a model's calibration.

This example fits a logistic regression to a small and slightly noisy
dataset. The first figure (without confidence intervals) indicates reasonably
good calibration. While its histogram helpfully suggests that predicted
probabilities in the middle range are scarce, it doesn't directly display how
much this scarcity affects our guess as to where a model is well calibrated.
Where is there enough data to evaluate calibration? Where is there not enough?

Confidence intervals display how varied the fraction of positive cases could be
for each bin of predicted probabilities. The interpretation of a 95% confidence
interval for each bin is that out of 100 future test bins, 95 estimated
intervals are expected to capture the true fraction of positives. For example,
consider the last bin in the second figure. The 95% confidence interval,
(0.84, 0.98), was estimated from test samples in the bin corresponding to
predicted probabilities around 0.90. This interval may be too wide in some
sensitive applications. Analogous interpretations for the rest of the bins may
suggest that calibration can't be confidently assessed for predicted
probabilities in any range except for those around 0.06 (the first bin).

The conclusion we got from confidence intervals may not have been gotten from
only looking at the histogram. It's tempting to be confident that a model is
well calibrated for predicted probabilities around 0.90 (the last bin) by only
looking at the histogram, as predictions are peaked at the ends. The confidence
interval at the last bin is less ambiguous. In some sensitive applications,
the interval that the second figure provides may be too wide. It suggests that
more data needs to be collected to evaluate this model.
"""
print(__doc__)

# Author: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD Style.

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import calibration_curve
from sklearn.model_selection import train_test_split


# Create dataset of classification task with seemingly good calibration but few
# predicted probabilities in middle of 0-1 line
X, y = datasets.make_classification(n_samples=1000, flip_y=0.1,
                                    random_state=9139)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
                                                    random_state=9245)

# Logistic regression with no calibration
lr = LogisticRegression(C=1.)
lr.fit(X_train, y_train)
prob_pos = lr.predict_proba(X_test)[:, 1]
fraction_of_positives, mean_predicted_value, fraction_of_positives_std = \
    calibration_curve(y_test, prob_pos, n_bins=5, return_std=True)


def plot_calibration_curve(fraction_of_positives, mean_predicted_value,
                           fig_index, with_ci=False):
    _ = plt.figure(fig_index, figsize=(10, 10))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))

    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")
    if with_ci:
        ax1.errorbar(mean_predicted_value, fraction_of_positives,
                     yerr=2*fraction_of_positives_std, fmt="s-",
                     label="Model with 95% CI", capsize=5)
    else:
        ax1.plot(mean_predicted_value, fraction_of_positives, "s-",
                 label="Model")

    ax2.hist(prob_pos, range=(0, 1), bins=10, label="Model",
             histtype="step", lw=2)

    ax1.set_ylabel("Fraction of positives")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="lower right")
    ax1.set_title("Calibration plot (reliability curve)")

    ax2.set_xlabel("Mean predicted value")
    ax2.set_ylabel("Count")
    ax2.legend(loc="upper center", ncol=2)

    plt.tight_layout()

# Plot calibration curve without 95% CIs
plot_calibration_curve(fraction_of_positives, mean_predicted_value, 1)

# Plot calibration curve with 95% CIs
plot_calibration_curve(fraction_of_positives, mean_predicted_value, 2,
                       with_ci=True)

plt.show()

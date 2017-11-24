# -*- coding: utf-8 -*-

"""
================================================================
Using KBinsDiscretizer to discretize continuous features
================================================================

The example compares prediction result of linear regression (linear model)
and decision tree (tree based model) with and without discretization of
real-valued features.

As is shown in the result before discretization, linear model is fast to
build and relatively straightforward to interpret, but can only model
linear relationships, while decision tree can build a much more complex model
of the data. One way to make linear model more powerful on continuous data
is to use discretization (also known as binning). In the example, we
discretize the feature and one-hot encode the transformed data. Note that if
the bins are not reasonably wide, there would appear to be a substantially
increased risk of overfitting, so the discretiser parameters need to be tuned
under cv.

After discretization, linear regression and decision tree make exactly the
same prediction. As features are constant within each bin, any model must
predict the same value for all points within a bin. Compared with the result
before discretization, linear model become much more flexible while decision
tree gets much less flexible. Note that binning features generally has no
beneficial effect for tree-based models, as these models can learn to split
up the data anywhere.

"""

# Author: Andreas Müller
#         Hanmin Qin <qinhanmin2005@sina.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.tree import DecisionTreeRegressor

print(__doc__)

# construct the dataset
rnd = np.random.RandomState(42)
X = rnd.uniform(-3, 3, size=100)
y = np.sin(X) + rnd.normal(size=len(X)) / 3
X = X.reshape(-1, 1)

# transform the dataset with KBinsDiscretizer
enc = KBinsDiscretizer(n_bins=10, encode='onehot')
X_binned = enc.fit_transform(X)

# predict with original dataset
fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 4))
plt.sca(axes[0])
line = np.linspace(-3, 3, 1000, endpoint=False).reshape(-1, 1)
reg = LinearRegression().fit(X, y)
plt.plot(line, reg.predict(line), linewidth=2, color='green',
         label="linear regression")
reg = DecisionTreeRegressor(min_samples_split=3, random_state=0).fit(X, y)
plt.plot(line, reg.predict(line), linewidth=2, color='red',
         label="decision tree")
plt.plot(X[:, 0], y, 'o', c='k')
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.legend(loc="best")
plt.title("Result before discretization")

# predict with transformed dataset
plt.sca(axes[1])
line_binned = enc.transform(line)
reg = LinearRegression().fit(X_binned, y)
plt.plot(line, reg.predict(line_binned), linewidth=2, color='green',
         linestyle='-', label='linear regression')
reg = DecisionTreeRegressor(min_samples_split=3,
                            random_state=0).fit(X_binned, y)
plt.plot(line, reg.predict(line_binned), linewidth=2, color='red',
         linestyle=':', label='decision tree')
plt.plot(X[:, 0], y, 'o', c='k')
bins = enc.offset_[0] + enc.bin_width_[0] * np.arange(1, enc.n_bins_[0])
plt.vlines(bins, *plt.gca().get_ylim(), linewidth=1, alpha=.2)
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.title("Result after discretization")
plt.tight_layout()
plt.show()

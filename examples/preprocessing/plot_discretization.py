# -*- coding: utf-8 -*-

"""

================================================================
Using KBinsDiscretizer to discretize continuous features
================================================================

The example compares prediction result of linear regression (linear model)
and decision tree (tree based model) before and after discretization.

As is shown in the result before discretization, linear model can only model
linear relationships, while decision tree can build a much more complex model
of the data. One way to make linear model more powerful on continuous data
is to use discretization (also known as binning).

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
y_no_noise = (np.sin(4 * X) + X)
y = (y_no_noise + rnd.normal(size=len(X))) / 2
X = X.reshape(-1, 1)

# transform the dataset with KBinsDiscretizer
enc = KBinsDiscretizer(n_bins=10, encode='onehot-dense')
X_binned = enc.fit_transform(X)

# predict with original dataset
plt.figure(figsize=(10, 4))
plt.subplot(121)
line = np.linspace(-3, 3, 1000, endpoint=False).reshape(-1, 1)
reg = DecisionTreeRegressor(min_samples_split=3).fit(X, y)
plt.plot(line, reg.predict(line), linewidth=2, color='green',
         label="decision tree")
reg = LinearRegression().fit(X, y)
plt.plot(line, reg.predict(line), linewidth=2, color='red',
         label="linear regression")
plt.plot(X[:, 0], y, 'o', c='k')
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.legend(loc="best")
plt.title("Result before discretization")

# predict with transformed dataset
plt.subplot(122)
line_binned = enc.transform(line)
reg = LinearRegression().fit(X_binned, y)
plt.plot(line, reg.predict(line_binned), linewidth=2, color='green',
         linestyle='-', label='linear regression')
reg = DecisionTreeRegressor(min_samples_split=3).fit(X_binned, y)
plt.plot(line, reg.predict(line_binned), linewidth=2, color='red',
         linestyle=':', label='decision tree')
plt.plot(X[:, 0], y, 'o', c='k')
bins = enc.offset_[0] + enc.bin_width_[0] * np.arange(1, enc.n_bins_[0])
plt.vlines(bins, -3, 3, linewidth=1, alpha=.2)
plt.legend(loc="best")
plt.ylabel("Regression output")
plt.xlabel("Input feature")
plt.title("Result after discretization")
plt.tight_layout()
plt.show()

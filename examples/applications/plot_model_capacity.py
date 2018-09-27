print(__doc__)

# Author: Danilo Bzdok
# License: BSD 3 clause

import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline

rs = np.random.RandomState(123)
X = rs.randn(1500, 2) * 100
Y = rs.randn(1500) * 100

target_func = make_pipeline(PolynomialFeatures(5), Ridge(random_state=0))

Y = target_func.fit(X, Y).predict(X)
X += rs.randn(X.shape[0], X.shape[1]) * 25  # add some noise

X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=500, random_state=0)

results = []
pins = np.arange(1000) + 1
for n_train in pins:
    print(n_train)
    perm_results = []
    for n_perm in range(50):
        # draw samples of always bigger shares of training data
        subset_inds = rs.randint(0, len(X_train), n_train)

        cur_X = X_train[subset_inds]
        cur_Y = Y_train[subset_inds]

        lowcap_func = make_pipeline(PolynomialFeatures(2),
                                    Ridge(random_state=0))
        lowcap_func.fit(cur_X, cur_Y)
        is_low_MSE = mean_squared_error(lowcap_func.predict(cur_X), cur_Y)
        oos_low_MSE = mean_squared_error(lowcap_func.predict(X_test), Y_test)

        optcap_func = make_pipeline(PolynomialFeatures(5),
                                    Ridge(random_state=0))
        optcap_func.fit(cur_X, cur_Y)
        is_opt_MSE = mean_squared_error(optcap_func.predict(cur_X), cur_Y)
        oos_opt_MSE = mean_squared_error(optcap_func.predict(X_test), Y_test)
        
        perm_results.append([is_low_MSE, oos_low_MSE, is_opt_MSE, oos_opt_MSE])
    
    results.append(perm_results)

results = np.array(results)
means = results.mean(axis=1)
std = results.std(axis=1)

# compute the optimal accuracy based on the target function (Bayes error rate)
BER = mean_squared_error(target_func.predict(X), Y)

from matplotlib import pylab as plt
plt.figure(figsize=(6, 5))
plt.plot(pins[::3], means[:, 0][::3], label='Train (quadratic)',
         color='slateblue')
plt.plot(pins[::3], means[:, 1][::3], label='Test (quadratic)',
         color='darkblue')
plt.plot(pins[::3], means[:, 2][::3], label='Train (optimal)', color='salmon')
plt.plot(pins[::3], means[:, 3][::3], label='Test (optimal)', color='red')
plt.ylim(0, 125)
# plt.xticks(np.arange(pins), pins)
plt.xlabel('Number of training examples')
plt.ylabel('Error (MSE)')
plt.xticks(
    np.arange(len(pins))[::100],
    np.linspace(100, 1000, 10, dtype=np.int))
plt.hlines(BER, 0 - 85, len(pins) + 85, linestyles='dashed', color='grey',
           label='Bayes error rate', linewidth=3)
plt.legend(loc='middle right')
plt.savefig('polynomial_fit.png', DPI=500)
plt.show()


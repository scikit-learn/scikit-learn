import numpy as np
import pylab as pl
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
from scipy.stats.mstats import mquantiles

from sklearn.cross_validation import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import gradient_boosting
from sklearn.datasets.cal_housing import fetch_cal_housing

# fetch California housing dataset
cal_housing = fetch_cal_housing()

# split 80/20 train-test
X_train, X_test, y_train, y_test = train_test_split(cal_housing.data,
                                                    cal_housing.target,
                                                    test_size=0.2,
                                                    random_state=1)
names = cal_housing.feature_names

clf = GradientBoostingRegressor(n_estimators=800, max_depth=4,
                                learn_rate=0.1, loss='huber',
                                random_state=1)
print('_' * 80)
print("Training GBRT...")
clf.fit(X_train, y_train)
print("fin.")

print('_' * 80)
print('Plot train-test deviance')
print


def aae(y_true, y_pred):
    return np.mean(np.abs(y_true.ravel() - y_pred.ravel()))


test_score = np.zeros((clf.n_estimators,), dtype=np.float64)
for i, y_pred in enumerate(clf.staged_decision_function(X_test)):
    test_score[i] = aae(y_test, y_pred)

train_score = np.zeros((clf.n_estimators,), dtype=np.float64)
for i, y_pred in enumerate(clf.staged_decision_function(X_train)):
    train_score[i] = aae(y_train, y_pred)

pl.figure()
pl.title('Training and Test Absolute Error')
pl.plot(np.arange(clf.n_estimators) + 1, train_score, 'b-',
        label='Train error')
pl.plot(np.arange(clf.n_estimators) + 1, test_score, 'r-',
        label='Test error')
pl.legend(loc='upper right')
pl.xlabel('Boosting Iterations')
pl.ylabel('Deviance')

print('_' * 80)
print('One-way partial dependence plots')
print

pl.figure()
sub_plots = []

for i, fx in enumerate([0, 5, 1, 2]):
    name = names[fx]
    target_feature = np.array([fx], dtype=np.int32)

    ax = pl.subplot(2, 2, i + 1)

    # plot partial dependency
    pdp, (axis,) = gradient_boosting.partial_dependency(clf, target_feature,
                                                        X=X_train)
    ax.plot(axis, pdp.ravel(), 'g-')

    # plot data deciles
    deciles = mquantiles(X_train[:, fx], prob=np.arange(0.1, 1.0, 0.1))
    trans = matplotlib.transforms.blended_transform_factory(ax.transData,
                                                            ax.transAxes)
    ax.vlines(deciles, 0.0, 0.05, transform=trans, color='red')

    pl.xlabel(name)
    pl.ylabel('Partial Dependency')

    sub_plots.append(ax)

# set common ylim
y_min = min((ax.get_ylim()[0] for ax in sub_plots))
y_max = max((ax.get_ylim()[1] for ax in sub_plots))
for ax in sub_plots:
    ax.set_ylim((y_min, y_max))

pl.tight_layout()

print('_' * 80)
print('Two-way partial dependence plot')
print

fig = pl.figure()

target_feature = np.array([1, 5], dtype=np.int32)
pdp, (x_axis, y_axis) = gradient_boosting.partial_dependency(clf,
                                                             target_feature,
                                                             X=X_train,
                                                             grid_resolution=50)
XX, YY = np.meshgrid(x_axis, y_axis)

Z = pdp.T.reshape(XX.shape).T
ax = Axes3D(fig)
surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=pl.cm.BuPu)

ax.set_xlabel(names[target_feature[0]])
ax.set_ylabel(names[target_feature[1]])
ax.set_zlabel('Partial dependency')

#  pretty init view
ax.view_init(elev=22, azim=122)

pl.colorbar(surf)

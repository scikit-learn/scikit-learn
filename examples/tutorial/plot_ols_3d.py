import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from scikits.learn import datasets, linear_model

diabetes = datasets.load_diabetes()
indices = (0, 1)

X_train = diabetes.data[:-20, indices]
X_test  = diabetes.data[-20:, indices]
y_train = diabetes.target[:-20]
y_test  = diabetes.target[-20:]

ols = linear_model.LinearRegression()
ols.fit(X_train, y_train)

fig = pl.figure(1, figsize=(4, 3))
pl.clf()
ax = Axes3D(fig, elev=43.5, azim=-110)
n = 100
ax.scatter(X_train[:, 0], X_train[:, 1], y_train, c='k', marker='+')

ax.plot_surface(np.array([[-.1, -.1], [.15, .15]]), 
                np.array([[-.1, .15], [-.1, .15]]),
                ols.predict(np.array([[-.1, -.1, .15, .15],
                                      [-.1, .15, -.1, .15]]).T
                            ).reshape((2, 2)),
                alpha=.5)
ax.set_xlabel('X_1')
ax.set_ylabel('X_2')
ax.set_zlabel('Y')
ax.set_xticks(())
ax.set_yticks(())
ax.set_zticks(())
#pl.savefig('diabetes_ols_diag.png')
ax.elev = -.5
ax.azim = 0
#pl.savefig('diabetes_ols_x1.png')
ax.elev = -.5
ax.azim = 90
#pl.savefig('diabetes_ols_x2.png')



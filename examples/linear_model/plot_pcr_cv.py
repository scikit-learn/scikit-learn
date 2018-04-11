from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import scale
from sklearn import datasets
from sklearn.linear_model import PCR
from sklearn.cross_validation import train_test_split
import matplotlib.pyplot as plt
import numpy as np

diabetes = datasets.load_diabetes()

X = diabetes.data
y = diabetes.target
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.33, random_state=42)

pcr = PCR()
kf_10 = KFold(n_splits=10, shuffle=True, random_state=1)

pcr.fit(scale(X_train), y_train)

X_reduced_train = pcr.get_transformed_data()  # get the reduced dataset
regr = pcr.get_regression_model()  # get the undelying regression model
n = len(X_reduced_train)

mse = list()

score = -1 * cross_val_score(regr, np.ones((n, 1)), y_train.ravel(),
                             cv=kf_10, scoring='neg_mean_squared_error').mean()
mse.append(score)

for i in np.arange(1, X.shape[1]):
    score = -1 * cross_val_score(regr, X_reduced_train[:, :i], y_train.ravel(
    ), cv=kf_10, scoring='neg_mean_squared_error').mean()
    mse.append(score)

fig, ax = plt.subplots()
ax.plot(mse, '-v')
ax.xaxis.set_ticks(np.arange(1, X.shape[1], 1))
plt.show()

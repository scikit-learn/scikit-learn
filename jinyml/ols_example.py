

import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

# load the diabets dataset

diabetes = datasets.load_diabetes()

# Use only one feature

diabetes_X = diabetes.data[:, np.newaxis, 2]

# Split the data into training/testing sets
diabetes_x_train = diabetes_X[:-20]
diabetes_x_test = diabetes_X[-20:]

diabetes_y_train = diabetes.target[:-20]
diabetes_y_test = diabetes.target[-20:]

# create linear regression object
regr = linear_model.LinearRegression()
# train the model using the training sets
regr.fit(diabetes_x_train, diabetes_y_train)

# make predictions using the testing set
diabetes_y_pred = regr.predict(diabetes_x_test)

# the coefficients
print("Coefficients: \n'", regr.coef_)
# The mean squared error
print("Variance score: %.2f " % r2_score(diabetes_y_test, diabetes_y_pred))

# Plot outputs

plt.scatter(diabetes_x_test, diabetes_y_test, color='black')
plt.plot(diabetes_x_test, diabetes_y_pred, color='blue', linewidth=3)

plt.xticks(())
plt.yticks(())
plt.show()
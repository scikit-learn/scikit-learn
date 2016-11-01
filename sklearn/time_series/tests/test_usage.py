from numpy.testing import assert_equal
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.time_series.time_series_estimator import *
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import ElasticNet

X =np.vstack((np.linspace(1,1000,1000),np.linspace(1,1000,1000)*10)).transpose()
y = np.linspace(1,1000,1000)*np.linspace(1,1000,1000)
X_train, X_test = time_series_split(X,test_size=.5)
y_train, y_test = time_series_split(y,test_size=.5)
train_size = X_train.shape[0]
test_size = X_test.shape[0]

def test_fit_predict():
    elastic = ElasticNet()
    elastic_tsr = TimeSeriesRegressor(elastic, alpha=1, max_iter=10000)
    # Note that ElasticNet parameters can be input in the ElasticNet() or the TimeSeriesRegressor()
    # This allows grid search and other meta estimation techniques to access the parameters of the base_estimator
    elastic_tsr.fit(X_train, y_train)
    pred_train_2 = elastic_tsr.predict(X_train)  # outputs a numpy array of length: len(X_train)-n_prev
    pred_test_2 = elastic_tsr.predict(X_test)
    assert_equal(pred_train_2.shape,(train_size-1,))
    assert_equal(pred_test_2.shape, (test_size-1,))

def test_n_prev():
    n_prev=2
    linear_model = LinearRegression()
    linear_tsr = TimeSeriesRegressor(linear_model, n_prev=n_prev)
    linear_tsr.fit(X_train, y_train)
    pred_train = linear_tsr.predict(X_train) #outputs a numpy array of length: len(X_train)-n_prev
    pred_test = linear_tsr.predict(X_test)
    assert_equal(pred_train.shape,(train_size-n_prev,))
    assert_equal(pred_test.shape, (test_size-n_prev,))

def test_n_ahead():
    n_ahead = 3
    linear_model = LinearRegression()
    linear_tsr = TimeSeriesRegressor(linear_model, n_ahead=n_ahead)
    linear_tsr.fit(X_train, y_train)
    pred_train_2 = linear_tsr.predict(X_train)
    pred_test_2 = linear_tsr.predict(X_test)
    assert_equal(pred_train_2.shape,(train_size-n_ahead,))
    assert_equal(pred_test_2.shape, (test_size-n_ahead,))

def test_forecast():
    n_prev = 3
    tsr = TimeSeriesRegressor(LinearRegression(), n_prev=n_prev)
    tsr.fit(X_train)
    fc = tsr.forecast(X_train, len(X_test))
    assert_equal(fc.shape,(test_size,2L))

def test_pipeline():
    n_prev = 3
    tsr = TimeSeriesRegressor(ElasticNet(max_iter=1000), n_prev=n_prev)

    param_grid = [{'alpha': [.01, .05],
                   'l1_ratio': [0, .25, 1]}]
    cv = cascade_cv(len(X_train), n_folds=3)
    grid = GridSearchCV(tsr, param_grid, cv=cv)
    grid.fit(X_train, y_train)
    pred_train_3 = grid.predict(X_train)  # outputs a numpy array of length: len(X_train)-n_prev
    pred_test_3 = grid.predict(X_test)
    assert_equal(pred_train_3.shape,(train_size-n_prev,))
    assert_equal(pred_test_3.shape, (test_size - n_prev,))


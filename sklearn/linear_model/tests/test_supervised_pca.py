# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD 3 clause


from SupervisedPCA import SupervisedPCARegressor
from SupervisedPCA import SupervisedPCAClassifier
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn import datasets
import numpy as np

def test_supervisedpcaRegressor_fit():
    # Test LinearRegression on a simple dataset.
    # a simple dataset
    diabetes=datasets.load_diabetes()
    X = diabetes.data
    Y = diabetes.target

    spca = SupervisedPCARegressor()
    spca.fit(X, Y,threshold=300,n_components=2)

    assert_array_equal(spca._leavouts,[1,5])
    assert_almost_equal(spca._model.coef_[0], [-537.7584256])

def test_supervisedpcaRegressor_predict():
    diabetes=datasets.load_diabetes()
    X = diabetes.data
    Y = diabetes.target

    spca = SupervisedPCARegressor()
    spca.fit(X, Y)
    
    predictions = spca.predict(X)
    mae=np.mean(abs(predictions-Y))
    assert_almost_equal(mae,51.570682097)
    
def test_supervisedpcaClassifier():
    iris=datasets.load_iris()
    X = iris.data
    Y = iris.target

    spca = SupervisedPCAClassifier()
    spca.fit(X, Y,threshold=1,n_components=2)

    assert_array_equal(spca._leavouts,[0,1])
    assert_almost_equal(spca._model.coef_[0][0], -2.43973048)
    
def test_supervisedpcaClassifier_predict():
    iris=datasets.load_iris()
    X = iris.data
    Y = iris.target

    spca = SupervisedPCAClassifier()
    spca.fit(X, Y)
    
    predictions = spca.predict(X)
    error=np.mean(sum(abs(predictions-Y))/float(len(predictions)))
    assert_almost_equal(error,0.08666666)
    
    
    
test_supervisedpcaRegressor_fit()
test_supervisedpcaRegressor_predict()
test_supervisedpcaClassifier()
test_supervisedpcaClassifier_predict()
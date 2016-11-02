import numpy as np

from sklearn import svm, multiclass

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raise_message

def test_ovr_decision_function():
    train_base_points = np.array([[1, 2], [2, 1]])
    
    # For all the quadrants (classes)
    x = np.vstack((
        train_base_points * [1, 1],    # Q1
        train_base_points * [-1, 1],   # Q2
        train_base_points * [-1, -1],  # Q3
        train_base_points * [1, -1]    # Q4
        ))
    
    y = [0] * 2 + [1] * 2 + [2] * 2 + [3] * 2
    
    # First point is closer to the decision boundaries than the second point
    test_base_points = np.array([[5, 5], [10, 10]])
    
    # For all the quadrants (classes)
    x_test = np.vstack((
        test_base_points * [1, 1],    # Q1
        test_base_points * [-1, 1],   # Q2
        test_base_points * [-1, -1],  # Q3
        test_base_points * [1, -1],   # Q4
        ))
    
    y_test = [0] * 2 + [1] * 2 + [2] * 2 + [3] * 2
    
    svc = svm.SVC(kernel='linear', decision_function_shape='ovr')
    svc.fit(x, y)
    
    y_pred = svc.predict(x_test)
    
    # Test if the prediction is the same as y
    assert_array_equal(y_pred, y_test)
    
    deci_vals = svc.decision_function(x_test)
    
    # Assert that the predicted class has the maximum value
    assert_array_equal(np.argmax(deci_vals, axis=1), y_pred)
    
    # Get decision value at test points for the predicted class
    pred_class_deci_vals = deci_vals[range(8), y_pred].reshape((4, 2))
    
    # Assert pred_class_deci_vals > 0 here
    assert_greater(np.min(pred_class_deci_vals), 0.0)
    
    # Test if the first point has lower decision value on every quadrant
    # compared to the second point
    assert_true(np.all(pred_class_deci_vals[:, 0] < pred_class_deci_vals[:, 1]))

from sklearn.datasets import make_classification
from sklearn.naive_bayes import MultinomialNB
from sklearn.mimic_calibration import _MimicCalibration
import numpy as np


def test_mimic_example():
    n_samples = 1000
    X, y = make_classification(n_samples=3 * n_samples, n_features=6,
                               random_state=42)
    X -= X.min()  # MultinomialNB only allows positive X
    # split train and test
    X_train, y_train = X[:n_samples], y[:n_samples]
    X, y = X[n_samples:2 * n_samples], y[n_samples:2 * n_samples]
    clf = MultinomialNB().fit(X_train, y_train)
    y_prob = clf.predict_proba(X)
    y_pre_calib_prob = np.array([y[1] for y in y_prob])
    mimic_obj = _MimicCalibration(threshold_pos=5, boundary_choice=2, record_history=False)

    # X: probability prediction from MultinomialNB
    # y: binary target, [0, 1]
    X = y_pre_calib_prob
    # use mimi calibration to calibrate X
    mimic_obj.fit(X, y)
    # y_calib_prob: the mimic-calibrated probaility
    y_calib_prob = mimic_obj.predict(X)
    assert (y_calib_prob.shape[0] == X.shape[0]), "The length of calibrated prob must be the same as pre-calibrated prob."


test_mimic_example()

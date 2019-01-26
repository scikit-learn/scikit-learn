from sklearn.datasets import make_classification
from sklearn.calibration import CalibratedClassifierCV
from sklearn.naive_bayes import MultinomialNB
from sklearn.mimic_calibration import _MimicCalibration
from copy import copy


def test_mimic_calibration():
    n_samples = 50
    X, y = make_classification(n_samples=3 * n_samples, n_features=6,
                               random_state=42)
    X -= X.min()  # MultinomialNB only allows positive X
    # split train and test
    X_train, y_train = X[:n_samples], y[:n_samples]
    X_calib, y_calib = X[n_samples:2 * n_samples], y[n_samples:2 * n_samples]
    clf = MultinomialNB().fit(X_train, y_train)
    clf_prob = CalibratedClassifierCV(clf, method="mimic", cv="prefit")
    clf_prob.fit(X_calib, y_calib)
    y_mimic_calibrated_score = clf_prob.predict(y_calib)

def test_mimic_calibration_2():
    n_samples = 50
    X, y = make_classification(n_samples=3 * n_samples, n_features=6,
                               random_state=42)
    # MultinomialNB only allows positive X
    X -= X.min()
    # split train and test
    X_train, y_train = X[:n_samples], y[:n_samples]
    X_calib, y_calib = X[n_samples:2 * n_samples], y[n_samples:2 * n_samples]
    clf = MultinomialNB().fit(X_train, y_train)
    y_calib_score = clf.predict(X_calib)
    mimicObject = _MimicCalibration()
    sorted_index = y_calib_score.argsort()
    y_score = y_calib_score[sorted_index]
    y_target = y_calib[sorted_index]
    bin_info, total_number_pos = mimicObject.construct_initial_bin(y_score, y_target, 5)
    assert(total_number_pos == y_target.sum()), "The number of positive is inconsistent.{x} vs {y}"\
        .format(x = total_number_pos, y = y_target.sum())

    current_binning = bin_info
    # test 2.
    new_bin_temp, increasing_flag = mimicObject.merge_bins(copy(current_binning), True)

    # test 3.
    result = mimicObject.run_merge_function(current_binning, record_history=False)

    # test 4.
    history_result = mimicObject.run_merge_function(current_binning, record_history=True)
    import matplotlib.pylab as plt
    for hist_trail in history_result:
        print(hist_trail)
        plt.plot(hist_trail)

# test_mimic_calibration()
test_mimic_calibration_2()

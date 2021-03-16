"""
Testing for the tree_m5 module (sklearn.tree_m5).
"""
import numpy as np
import pytest

from sklearn import ensemble
from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.metrics import mean_squared_error

from sklearn.tree import export_text_m5
from sklearn.tree._classes_m5 import M5Prime, predict_from_leaves


@pytest.mark.parametrize("smoothing_constant", [75, 100], ids="smoothing_constant={}".format)
def test_m5p_smoothing(smoothing_constant):
    """ tests that the M5P smoothing feature works correctly """

    import matplotlib.pyplot as plt
    plt.ion()  # comment to debug
    debug_prints = False  # change to true to debug the smoothing procedure

    # use a given random seed so as to enforce deterministic output
    # (and 100 gives a particularly interesting result with constant leaves, thats why we use it)
    np.random.seed(100)

    x = np.random.uniform(-20, -2, (500, 1))
    piece_1_mask = x < -12
    piece_2_mask = ~piece_1_mask
    piece_2_train_indices = np.argsort(x, axis=0)[-5:]  # last 5 samples
    # training_mask = np.ma.mask_or(piece_1_mask, piece_2_train_mask)
    training_mask = piece_1_mask.copy()
    training_mask[piece_2_train_indices] = True
    test_mask = ~training_mask
    X_train = x[training_mask].reshape(-1, 1)
    X_test = x[test_mask].reshape(-1, 1)

    # create a y that is piecewise linear
    def piece1(x):
        return 12.3 * x + 52.1
    y = piece1(x) + 5 * np.random.randn(*x.shape)

    z = 0
    def piece2(x):
        return 5.3 * x + z
    z = piece1(-12) - piece2(-12)

    y[piece_2_mask] -= piece1(x[piece_2_mask])
    y[piece_2_mask] += piece2(x[piece_2_mask])
    y_train = y[training_mask].reshape(-1, 1)
    y_test = y[test_mask].reshape(-1, 1)
    plt.plot(X_train, y_train, 'xb')
    plt.plot(X_test, y_test, 'xg')
    plt.draw()
    plt.pause(0.001)

    # Fit M5P without smoothing
    print("Without smoothing: ")
    model_no_smoothing = M5Prime(use_smoothing=False, debug_prints=False)
    model_no_smoothing.fit(X_train, y_train)
    print(export_text_m5(model_no_smoothing, out_file=None, node_ids=True))

    # Fit M5P with smoothing
    print("With smoothing: ")
    model_smoothing = M5Prime(use_smoothing=True, debug_prints=debug_prints, smoothing_constant=smoothing_constant)
    model_smoothing.fit(X_train, y_train)
    print(export_text_m5(model_smoothing, out_file=None, node_ids=True))

    # Predict
    # --no smoothing: compare the fast and slow methods
    y_pred_no_smoothing = model_no_smoothing.predict(x).reshape(-1, 1)
    y_pred_no_smoothing_slow = predict_from_leaves(model_no_smoothing, x, smoothing=False).reshape(-1, 1)

    # --smoothing: compare pre-computed and late methods
    y_pred_installed_smoothing = model_smoothing.predict(x).reshape(-1, 1)
    y_pred_late_smoothing = model_no_smoothing.predict(x, smooth_predictions=True,
                                                       smoothing_constant=smoothing_constant).reshape(-1, 1)

    plt.plot(x, y_pred_no_smoothing, '.r', label="no smoothing")
    # plt.plot(x, y_pred_no_smoothing_slow, '.r', label="no smoothing slow")
    plt.plot(x, y_pred_late_smoothing, '.m', label="smoothing (on prediction)")
    plt.plot(x, y_pred_installed_smoothing, '.g', label="smoothing (installed models)")
    plt.legend()
    plt.draw()
    plt.pause(0.001)

    # make sure both ways to smooth give the same output
    np.testing.assert_array_equal(y_pred_no_smoothing, y_pred_no_smoothing_slow)
    np.testing.assert_array_almost_equal(y_pred_late_smoothing,
                                         y_pred_installed_smoothing, decimal=4)

    # compare performances without/with smoothing
    mse_no_smoothing = mean_squared_error(y_test, y_pred_no_smoothing[test_mask])
    mse_smoothing = mean_squared_error(y_test, y_pred_installed_smoothing[test_mask])
    print("M5P MSE: %.4f (no smoothing) %.4f (smoothing)" % (mse_no_smoothing, mse_smoothing))

    # simple assert: smoothing improves performance (in most cases but not always - thats why we fixed our random seed)
    assert mse_smoothing < mse_no_smoothing

    plt.close('all')


def test_boston_housing():
    """
    A copy of the Scikit Learn Gradient Boosting regression example, with the M5P in addition

    http://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_regression.html
    :return:
    """
    import matplotlib.pyplot as plt
    plt.ion()  # comment to debug

    # #############################################################################
    # Load data
    boston = datasets.load_boston()
    X, y = shuffle(boston.data, boston.target, random_state=13)
    X = X.astype(np.float32)
    offset = int(X.shape[0] * 0.9)
    X_train, y_train = X[:offset], y[:offset]
    X_test, y_test = X[offset:], y[offset:]

    # #############################################################################
    # Fit regression model
    params = {'n_estimators': 500, 'max_depth': 4, 'min_samples_split': 2,
              'learning_rate': 0.01, 'loss': 'ls'}
    clf = ensemble.GradientBoostingRegressor(**params)
    clf.fit(X_train, y_train)
    y_predicted = clf.predict(X_test)
    mse = mean_squared_error(y_test, y_predicted)
    print("XGBoost MSE: %.4f" % mse)

    # #############################################################################
    # Plot predictions
    plt.figure(figsize=(18, 12))
    plt.subplot(2, 3, 1)
    plt.title('Predictions on test set (RMSE = {:2f})'.format(np.sqrt(mse)))
    plt.plot(y_test, y_predicted, '.')
    plt.xlabel('true y')
    plt.ylabel('predicted_y')

    # #############################################################################
    # Plot training deviance

    # compute test set deviance
    test_score = np.zeros((params['n_estimators'],), dtype=np.float64)

    for i, y_pred in enumerate(clf.staged_predict(X_test)):
        test_score[i] = clf.loss_(y_test, y_pred)

    plt.subplot(2, 3, 2)
    plt.title('Deviance')
    plt.plot(np.arange(params['n_estimators']) + 1, clf.train_score_, 'b-',
             label='Training Set Deviance')
    plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
             label='Test Set Deviance')
    plt.legend(loc='upper right')
    plt.xlabel('Boosting Iterations')
    plt.ylabel('Deviance')

    # #############################################################################
    # Plot feature importance
    feature_importance = clf.feature_importances_
    # make importances relative to max importance
    feature_importance = 100.0 * (feature_importance / feature_importance.max())
    sorted_idx = np.argsort(feature_importance)
    pos = np.arange(sorted_idx.shape[0]) + .5
    plt.subplot(2, 3, 3)
    plt.barh(pos, feature_importance[sorted_idx], align='center')
    plt.yticks(pos, boston.feature_names[sorted_idx])
    plt.xlabel('Relative Importance')
    plt.title('Variable Importance')

    # ------- M5P
    # #############################################################################
    # Fit regression model
    params = {}
    clf = M5Prime(**params)
    clf.fit(X_train, y_train)

    # Print the tree
    print(export_text_m5(clf, out_file=None))
    print(export_text_m5(clf, out_file=None, feature_names=boston.feature_names))

    # Predict
    y_predicted = clf.predict(X_test)
    mse = mean_squared_error(y_test, y_predicted)
    print("M5P MSE: %.4f" % mse)

    # #############################################################################
    # Plot predictions
    plt.subplot(2, 3, 4)
    plt.title('Predictions on test set (RMSE = {:2f})'.format(np.sqrt(mse)))
    plt.plot(y_test, y_predicted, '.')
    plt.xlabel('true y')
    plt.ylabel('predicted_y')

    # Compress the tree (features-wise)
    idx = clf.compress_features()
    # Print the tree
    print(export_text_m5(clf, out_file=None))
    print(export_text_m5(clf, out_file=None, feature_names=boston.feature_names))

    # Predict
    y_predicted2 = clf.predict(X_test[:, idx])
    mse2 = mean_squared_error(y_test, y_predicted2)
    print("M5P2 MSE: %.4f" % mse)

    # #############################################################################
    # Plot predictions
    plt.subplot(2, 3, 5)
    plt.title('Predictions on test set (RMSE = {:2f})'.format(np.sqrt(mse2)))
    plt.plot(y_test, y_predicted2, '.')
    plt.xlabel('true y')
    plt.ylabel('predicted_y')

    # #############################################################################
    # Plot feature importance
    feature_importance = clf.feature_importances_
    # make importances relative to max importance
    feature_importance = 100.0 * (feature_importance / feature_importance.max())
    sorted_idx = np.argsort(feature_importance)
    pos = np.arange(sorted_idx.shape[0]) + .5
    plt.subplot(2, 3, 6)
    plt.barh(pos, feature_importance[sorted_idx], align='center')
    # do not forget that we now work on reindexed features
    plt.yticks(pos, boston.feature_names[idx][sorted_idx])
    plt.xlabel('Relative Importance')
    plt.title('Variable Importance')

    plt.close('all')


@pytest.mark.parametrize("use_smoothing", [None,  # None (default) = True = 'installed'
                                           False, 'on_prediction'])
def test_default_smoothing_modes(use_smoothing):
    """ Tests with default constant/ratio, depending on smoothing mode """

    half_x = np.random.random(25)
    X = np.r_[half_x, -half_x]
    Y = np.r_[2.8 * half_x + 5, -0.2 * half_x + 5]

    X = X.reshape(-1, 1)
    Y = Y.reshape(-1, 1)

    regr = M5Prime(use_smoothing=use_smoothing)
    regr.fit(X, Y)
    print(export_text_m5(regr, out_file=None, node_ids=True))
    preds = regr.predict(X)

    import matplotlib.pyplot as plt
    plt.ion()  # comment to debug
    plt.plot(X, Y, 'x', label='true')
    plt.plot(X, preds, 'x', label='prediction')
    plt.close('all')

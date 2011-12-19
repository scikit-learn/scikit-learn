from unittest import TestCase
import numpy as np
from sklearn.boosting import FunctionalGradientBoosting
from sklearn.tree import DecisionTreeRegressor
from sklearn.datasets.base import load_boston

class TestFunctionalGradientBoosting(TestCase):
    def setUp(self):
        self.task = load_boston()
        self.base_est = DecisionTreeRegressor(max_depth=2, min_split=4)
        self.boosting = FunctionalGradientBoosting(
                base_estimator=DecisionTreeRegressor(
                    max_depth=2,
                    min_split=4),
                n_estimators=5)

    def test_fit_returns_self(self):
        r = self.boosting.fit(self.task['data'], self.task['target'])
        assert r is self.boosting

    def test_1_estimator_matches_base(self):
        self.boosting = FunctionalGradientBoosting(
                base_estimator=DecisionTreeRegressor(
                    max_depth=2,
                    min_split=4),
                n_estimators=1)
        self.base_est.fit(self.task['data'], self.task['target'])
        self.boosting.fit(self.task['data'], self.task['target'])
        pred1 = self.base_est.predict(self.task['data'])
        pred2 = self.boosting.predict(self.task['data'])
        self.assert_(np.allclose(pred1, pred2))

    def test_n_estimators(self):
        assert len(self.boosting.estimators_) == 0
        self.boosting.fit(self.task['data'], self.task['target'])
        assert len(self.boosting.estimators_) == self.boosting.n_estimators

    def test_int_y_not_implemented(self):
        self.assertRaises(NotImplementedError,
                self.boosting.fit,
                np.ones((4, 5)), np.arange(4).astype('int'))

    def test_mse_always_goes_down(self):
        model = self.boosting
        task = self.task
        mse_list = []
        for fit_iter in model.fit_iter(task['data'], task['target']):
            mse_list.append(np.mean(fit_iter.residual ** 2))
            if len(mse_list) > 1:
                self.assert_(mse_list[-1] < mse_list[-2])

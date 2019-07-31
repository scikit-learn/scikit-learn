import pytest
from scipy.stats import norm

from sklearn.datasets import make_classification
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import HalvingGridSearchCV
from sklearn.model_selection import HalvingRandomSearchCV


class FastClassifier(DummyClassifier):
    """Dummy classifier that accepts parameters a, b, ... z.

    These parameter don't affect the predictions and are useful for fast
    grid searching."""

    def __init__(self, strategy='stratified', random_state=None,
                 constant=None, **kwargs):
        super().__init__(strategy=strategy, random_state=random_state,
                         constant=constant)

    def get_params(self, deep=False):
        params = super().get_params(deep=deep)
        for char in range(ord('a'), ord('z') + 1):
            params[chr(char)] = 'whatever'
        return params


def test_aggressive_elimination():
    # Test the aggressive_elimination parameter.

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': ('l1', 'l2'), 'b': list(range(30))}
    base_estimator = FastClassifier()
    ratio = 3

    # aggressive_elimination is only really relevant when there is not enough
    # budget.
    max_resources = 180

    # aggressive_elimination=True
    # In this case, the first iterations only use min_resources_ resources
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             aggressive_elimination=True,
                             max_resources=max_resources, ratio=ratio)
    sh.fit(X, y)
    assert sh.n_iterations_ == 4
    assert sh.n_required_iterations_ == 4
    assert sh.n_possible_iterations_ == 3
    assert sh._r_i_list == [20, 20, 60, 180]  # see how it loops at the start
    assert sh.n_remaining_candidates_ == 1

    # Make sure we get the same results with randomized search
    sh = HalvingRandomSearchCV(base_estimator, parameters,
                               n_candidates=60, cv=5,
                               aggressive_elimination=True,
                               max_resources=max_resources, ratio=ratio)
    sh.fit(X, y)
    assert sh.n_iterations_ == 4
    assert sh.n_required_iterations_ == 4
    assert sh.n_possible_iterations_ == 3
    assert sh._r_i_list == [20, 20, 60, 180]  # see how it loops at the start
    assert sh.n_remaining_candidates_ == 1

    # aggressive_elimination=False
    # In this case we don't loop at the start, and might end up with a lot of
    # candidates at the last iteration
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             aggressive_elimination=False,
                             max_resources=max_resources, ratio=ratio)
    sh.fit(X, y)

    assert sh.n_iterations_ == 3
    assert sh.n_required_iterations_ == 4
    assert sh.n_possible_iterations_ == 3
    assert sh._r_i_list == [20, 60, 180]
    assert sh.n_remaining_candidates_ == 3

    max_resources = n_samples
    # with enough budget, aggressive_elimination has no effect since it is not
    # needed

    # aggressive_elimination=True
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             aggressive_elimination=True,
                             max_resources=max_resources, ratio=ratio)
    sh.fit(X, y)

    assert sh.n_iterations_ == 4
    assert sh.n_required_iterations_ == 4
    assert sh.n_possible_iterations_ == 4
    assert sh._r_i_list == [20, 60, 180, 540]
    assert sh.n_remaining_candidates_ == 1

    # aggressive_elimination=False
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             aggressive_elimination=False,
                             max_resources=max_resources, ratio=ratio)
    sh.fit(X, y)

    assert sh.n_iterations_ == 4
    assert sh.n_required_iterations_ == 4
    assert sh.n_possible_iterations_ == 4
    assert sh._r_i_list == [20, 60, 180, 540]
    assert sh.n_remaining_candidates_ == 1


def test_force_exhaust_resources_false():
    # Test the force_exhaust_resources parameter when it's false or ignored.
    # This is the default case: we start at the beginning no matter what since
    # we do not overwrite min_resources_

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': [1, 2], 'b': [1, 2, 3]}
    base_estimator = FastClassifier()
    ratio = 3

    # with enough budget
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             force_exhaust_resources=False, ratio=ratio)
    sh.fit(X, y)
    assert sh.n_iterations_ == 2
    assert sh.n_required_iterations_ == 2
    assert sh.n_possible_iterations_ == 4
    assert sh._r_i_list == [20, 60]

    # with enough budget but min_resources!='auto': ignored
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             force_exhaust_resources=False, ratio=ratio,
                             min_resources=50)
    sh.fit(X, y)
    assert sh.n_iterations_ == 2
    assert sh.n_required_iterations_ == 2
    assert sh.n_possible_iterations_ == 3
    assert sh._r_i_list == [50, 150]

    # without enough budget (budget is exhausted anyway)
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             force_exhaust_resources=False, ratio=ratio,
                             max_resources=30)
    sh.fit(X, y)
    assert sh.n_iterations_ == 1
    assert sh.n_required_iterations_ == 2
    assert sh.n_possible_iterations_ == 1
    assert sh._r_i_list == [20]


@pytest.mark.parametrize('max_resources, r_i_list', [
    ('auto', [333, 999]),
    (1000, [333, 999]),
    (999, [333, 999]),
    (600, [200, 600]),
    (599, [199, 597]),
    (300, [100, 300]),
    (60, [20, 60]),
    (50, [20]),
    (20, [20]),
])
def test_force_exhaust_budget_true(max_resources, r_i_list):
    # Test the force_exhaust_resources parameter when it's true
    # in this case we need to change min_resources so that the last iteration uses as
    # much budget as possible

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': [1, 2], 'b': [1, 2, 3]}
    base_estimator = FastClassifier()
    ratio = 3
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=5,
                             force_exhaust_resources=True, ratio=ratio,
                             max_resources=max_resources)
    sh.fit(X, y)

    assert sh.n_possible_iterations_ == sh.n_iterations_ == len(sh._r_i_list)
    assert sh._r_i_list == r_i_list

    # Test same for randomized search
    sh = HalvingRandomSearchCV(base_estimator, parameters, n_candidates=6,
                               cv=5, force_exhaust_resources=True,
                               ratio=ratio, max_resources=max_resources)
    sh.fit(X, y)

    assert sh.n_possible_iterations_ == sh.n_iterations_ == len(sh._r_i_list)
    assert sh._r_i_list == r_i_list


@pytest.mark.parametrize(
    'max_resources, n_iterations, n_possible_iterations', [
        ('auto', 5, 9),  # whole budget is used
        (1024, 5, 9),
        (700, 5, 8),
        (512, 5, 8),
        (511, 5, 7),
        (32, 4, 4),
        (31, 3, 3),
        (16, 3, 3),
        (4, 1, 1),   # max_resources == min_resources, only one iteration is possible
    ])
def test_n_iterations(max_resources, n_iterations, n_possible_iterations):
    # test the number of actual iterations that were run depending on
    # max_resources

    n_samples = 1024
    X, y = make_classification(n_samples=n_samples, random_state=1)
    parameters = {'a': [1, 2], 'b': list(range(10))}
    base_estimator = FastClassifier()
    ratio = 2

    sh = HalvingGridSearchCV(base_estimator, parameters, cv=2, ratio=ratio,
                             max_resources=max_resources, min_resources=4)
    sh.fit(X, y)
    assert sh.n_required_iterations_ == 5
    assert sh.n_iterations_ == n_iterations
    assert sh.n_possible_iterations_ == n_possible_iterations


def test_resource_parameter():
    # Test the resource parameter

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': [1, 2], 'b': list(range(10))}
    base_estimator = FastClassifier()
    sh = HalvingGridSearchCV(base_estimator, parameters, cv=2,
                             resource='c', max_resources=10, ratio=3)
    sh.fit(X, y)
    assert set(sh._r_i_list) == set([1, 3, 9])
    for r_i, params, param_c in zip(sh.cv_results_['r_i'],
                                    sh.cv_results_['params'],
                                    sh.cv_results_['param_c']):
        assert r_i == params['c'] == param_c

    with pytest.raises(
            ValueError,
            match='Cannot use resource=1234 which is not supported '):
        sh = HalvingGridSearchCV(base_estimator, parameters, cv=2,
                                 resource='1234', max_resources=10)
        sh.fit(X, y)

    with pytest.raises(
            ValueError,
            match='Cannot budget on parameter c since it is part of the '
                  'searched parameters.'):
        parameters = {'a': [1, 2], 'b': [1, 2], 'c': [1, 3]}
        sh = HalvingGridSearchCV(base_estimator, parameters, cv=2,
                                 resource='c', max_resources=10)
        sh.fit(X, y)


@pytest.mark.parametrize(
    'max_resources, n_candidates, expected_n_candidates_', [
        (512, 'auto', 128),  # generate exactly as much as needed
        (32, 'auto', 8),
        (32, 8, 8),
        (32, 7, 7),  # ask for less than what we could
        (32, 9, 9),  # ask for more than 'reasonable'
    ])
def test_random_search(max_resources, n_candidates, expected_n_candidates_):
    # Test random search and make sure the number of generated candidates is as
    # expected

    n_samples = 1024
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': norm, 'b': norm}
    base_estimator = FastClassifier()
    sh = HalvingRandomSearchCV(base_estimator, parameters,
                               n_candidates=n_candidates,
                               cv=2,
                               max_resources=max_resources, ratio=2, min_resources=4)
    sh.fit(X, y)
    assert sh.n_candidates_ == expected_n_candidates_
    if n_candidates == 'auto':
        # Make sure 'auto' makes the last iteration use as much budget as we
        # can
        assert sh._r_i_list[-1] == max_resources


def test_groups_not_supported():
    base_estimator = FastClassifier()
    param_grid = {'a': [1]}
    sh = HalvingRandomSearchCV(base_estimator, param_grid)
    X, y = make_classification(n_samples=10)
    groups = [0] * 10

    with pytest.raises(ValueError, match="groups are not supported"):
        sh.fit(X, y, groups)

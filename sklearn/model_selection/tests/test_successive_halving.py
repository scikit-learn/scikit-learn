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


@pytest.mark.parametrize('klass', (HalvingGridSearchCV, HalvingRandomSearchCV))
@pytest.mark.parametrize(
    ('aggressive_elimination,'
     'max_resources,'
     'expected_n_iterations,'
     'expected_n_required_iterations,'
     'expected_n_possible_iterations,'
     'expected_n_remaining_candidates,'
     'expected_r_i_list,'), [
         # notice how it loops at the beginning
         (True, 'small', 4, 4, 3, 1, [20, 20, 60, 180]),
         # no aggressive elimination: we end up with less iterations and more
         # candidates at the end
         (False, 'small', 3, 4, 3, 3, [20, 60, 180]),
         # When the amount of resource isn't limited, aggressive_elimination
         # doesn't matter.
         (True, 'high', 4, 4, 4, 1, [20, 60, 180, 540]),
         (False, 'high', 4, 4, 4, 1, [20, 60, 180, 540]),
     ]
)
def test_aggressive_elimination(
        klass, aggressive_elimination, max_resources, expected_n_iterations,
        expected_n_required_iterations, expected_n_possible_iterations,
        expected_n_remaining_candidates, expected_r_i_list):
    # Test the aggressive_elimination parameter.

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': ('l1', 'l2'), 'b': list(range(30))}
    base_estimator = FastClassifier()

    if max_resources == 'small':
        max_resources = 180
    else:
        max_resources = n_samples

    sh = klass(base_estimator, parameters,
               aggressive_elimination=aggressive_elimination,
               max_resources=max_resources, ratio=3)

    if klass is HalvingRandomSearchCV:
        sh.set_params(n_candidates=2 * 30)  # same number as with the grid

    sh.fit(X, y)

    assert sh.n_iterations_ == expected_n_iterations
    assert sh.n_required_iterations_ == expected_n_required_iterations
    assert sh.n_possible_iterations_ == expected_n_possible_iterations
    assert sh._r_i_list == expected_r_i_list
    assert sh.n_remaining_candidates_ == expected_n_remaining_candidates


@pytest.mark.parametrize('klass', (HalvingGridSearchCV, HalvingRandomSearchCV))
@pytest.mark.parametrize(
    ('min_resources,'
     'max_resources,'
     'expected_n_iterations,'
     'expected_n_required_iterations,'
     'expected_n_possible_iterations,'
     'expected_r_i_list,'), [
         # with enough resources
         ('auto', 'auto', 2, 2, 4, [20, 60]),
         # with enough resources but min_resources!='auto': ignored
         (50, 'auto', 2, 2, 3, [50, 150]),
         # without enough resources (resources are exhausted anyway)
         ('auto', 30, 1, 2, 1, [20]),
     ]
)
def test_force_exhaust_resources_false(
        klass, min_resources, max_resources, expected_n_iterations,
        expected_n_required_iterations, expected_n_possible_iterations,
        expected_r_i_list):
    # Test the force_exhaust_resources parameter when it's false or ignored.
    # This is the default case: we start at the beginning no matter what since
    # we do not overwrite min_resources_
    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': [1, 2], 'b': [1, 2, 3]}
    base_estimator = FastClassifier()

    sh = klass(base_estimator, parameters, force_exhaust_resources=False,
               ratio=3, min_resources=min_resources,
               max_resources=max_resources)
    if klass is HalvingRandomSearchCV:
        sh.set_params(n_candidates=6)  # same number as with the grid

    sh.fit(X, y)
    assert sh.n_iterations_ == expected_n_iterations
    assert sh.n_required_iterations_ == expected_n_required_iterations
    assert sh.n_possible_iterations_ == expected_n_possible_iterations
    assert sh._r_i_list == expected_r_i_list


@pytest.mark.parametrize('klass', (HalvingRandomSearchCV, HalvingGridSearchCV))
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
def test_force_exhaust_resources_true(klass, max_resources, r_i_list):
    # Test the force_exhaust_resources parameter when it's true
    # in this case we need to change min_resources so that the last iteration
    # uses as much resources as possible

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': [1, 2], 'b': [1, 2, 3]}
    base_estimator = FastClassifier()

    sh = klass(base_estimator, parameters, force_exhaust_resources=True,
               ratio=3, max_resources=max_resources)
    if klass is HalvingRandomSearchCV:
        sh.set_params(n_candidates=6)  # same as for HalvingGridSearchCV
    sh.fit(X, y)

    assert sh.n_possible_iterations_ == sh.n_iterations_ == len(sh._r_i_list)
    assert sh._r_i_list == r_i_list


@pytest.mark.parametrize('klass', (HalvingRandomSearchCV, HalvingGridSearchCV))
@pytest.mark.parametrize(
    'max_resources, n_iterations, n_possible_iterations', [
        ('auto', 5, 9),  # all resources are used
        (1024, 5, 9),
        (700, 5, 8),
        (512, 5, 8),
        (511, 5, 7),
        (32, 4, 4),
        (31, 3, 3),
        (16, 3, 3),
        (4, 1, 1),  # max_resources == min_resources, only one iteration is
                    # possible
    ])
def test_n_iterations(klass, max_resources, n_iterations,
                      n_possible_iterations):
    # test the number of actual iterations that were run depending on
    # max_resources

    n_samples = 1024
    X, y = make_classification(n_samples=n_samples, random_state=1)
    parameters = {'a': [1, 2], 'b': list(range(10))}
    base_estimator = FastClassifier()
    ratio = 2

    sh = klass(base_estimator, parameters, cv=2, ratio=ratio,
               max_resources=max_resources, min_resources=4)
    if klass is HalvingRandomSearchCV:
        sh.set_params(n_candidates=20)  # same as for HalvingGridSearchCV
    sh.fit(X, y)
    assert sh.n_required_iterations_ == 5
    assert sh.n_iterations_ == n_iterations
    assert sh.n_possible_iterations_ == n_possible_iterations


@pytest.mark.parametrize('klass', (HalvingRandomSearchCV, HalvingGridSearchCV))
def test_resource_parameter(klass):
    # Test the resource parameter

    n_samples = 1000
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': [1, 2], 'b': list(range(10))}
    base_estimator = FastClassifier()
    sh = klass(base_estimator, parameters, cv=2, resource='c',
               max_resources=10, ratio=3)
    sh.fit(X, y)
    assert set(sh._r_i_list) == set([1, 3, 9])
    for r_i, params, param_c in zip(sh.cv_results_['resource_iter'],
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
            match='Cannot use parameter c as the resource since it is part '
                  'of the searched parameters.'):
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
    # Test random search and make sure the number of generated candidates is
    # as expected

    n_samples = 1024
    X, y = make_classification(n_samples=n_samples, random_state=0)
    parameters = {'a': norm, 'b': norm}
    base_estimator = FastClassifier()
    sh = HalvingRandomSearchCV(base_estimator, parameters,
                               n_candidates=n_candidates, cv=2,
                               max_resources=max_resources, ratio=2,
                               min_resources=4)
    sh.fit(X, y)
    assert sh.n_candidates_[0] == expected_n_candidates_
    if n_candidates == 'auto':
        # Make sure 'auto' makes the last iteration use as much resources as
        # we can
        assert sh._r_i_list[-1] == max_resources


@pytest.mark.parametrize('klass', (HalvingRandomSearchCV, HalvingGridSearchCV))
def test_groups_not_supported(klass):
    base_estimator = FastClassifier()
    param_grid = {'a': [1]}
    sh = klass(base_estimator, param_grid)
    X, y = make_classification(n_samples=10)
    groups = [0] * 10

    with pytest.raises(ValueError, match="groups are not supported"):
        sh.fit(X, y, groups)


@pytest.mark.parametrize('klass', (HalvingGridSearchCV, HalvingRandomSearchCV))
@pytest.mark.parametrize('params, expected_error_message', [
    ({'scoring': {'accuracy', 'accuracy'}},
     'Multimetric scoring is not supported'),
    ({'resource': 'not_a_parameter'},
     'Cannot use resource=not_a_parameter which is not supported'),
    ({'resource': 'a', 'max_resources': 100},
     'Cannot use parameter a as the resource since it is part of'),
    ({'max_resources': 'not_auto'},
     'max_resources must be either'),
    ({'max_resources': 100.5},
     'max_resources must be either'),
    ({'max_resources': -10},
     'max_resources must be either'),
    ({'min_resources': 'not_auto'},
     'min_resources must be either'),
    ({'min_resources': 0.5},
     'min_resources must be either'),
    ({'min_resources': -10},
     'min_resources must be either'),
    ({'force_exhaust_resources': True, 'min_resources': 5},
     'min_resources must be set to auto if '),
    ({'max_resources': 'auto', 'resource': 'b'},
     "max_resources can only be 'auto' if resource='n_samples'"),
    ({'min_resources': 15, 'max_resources': 14},
     "min_resources_=15 is greater than max_resources_=14"),
])
def test_input_errors(klass, params, expected_error_message):
    base_estimator = FastClassifier()
    param_grid = {'a': [1]}
    X, y = make_classification(100)

    sh = klass(base_estimator, param_grid, **params)

    with pytest.raises(ValueError, match=expected_error_message):
        sh.fit(X, y)

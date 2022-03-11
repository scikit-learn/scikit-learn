import pytest
from os import environ
from random import Random
from datetime import datetime

# Definition of the random_seed fixture. See the docstring of the fixture for
# more information.

RANDOM_SEED_RANGE = list(range(100))  # All seeds in [0, 99] should be valid.
random_seed_var = environ.get("SKLEARN_TESTS_GLOBAL_RANDOM_SEED")
if random_seed_var is None:
    # If the environment variable is not defined, pick-up one seed at random in
    # the range of admissible random seeds. Note, to make sure that all
    # pytest-xdist workers see the same seed, we seed the meta-random number
    # generator with a value derived from the year and the day.
    rng = Random(int(datetime.now().strftime("%Y%j")))
    random_seeds = [rng.choice(RANDOM_SEED_RANGE)]
elif random_seed_var == "all":
    random_seeds = RANDOM_SEED_RANGE
else:
    if "-" in random_seed_var:
        start, stop = random_seed_var.split("-")
        random_seeds = list(range(int(start), int(stop) + 1))
    else:
        random_seeds = [int(random_seed_var)]

    if min(random_seeds) < 0 or max(random_seeds) > 99:
        raise ValueError(
            "The value(s) of the environment variable "
            "SKLEARN_TESTS_GLOBAL_RANDOM_SEED must be in the range [0, 99] "
            f"(or 'all'), got: {random_seed_var}"
        )


def pytest_report_header(config):
    if len(random_seeds) == 1:
        repro_value = random_seeds[0]
    else:
        repro_value = random_seed_var
    return [
        "To reproduce this test run, set the following environment variable:",
        f'    SKLEARN_TESTS_GLOBAL_RANDOM_SEED="{repro_value}"',
    ]


@pytest.fixture(params=random_seeds)
def global_random_seed(request):
    """Fixture to ask for a random yet controllable random seed.

    All tests that use this fixture accept the contract that they should
    deterministically pass for any seed value from 0 to 99 included.

    If the SKLEARN_TESTS_GLOBAL_RANDOM_SEED environment variable is not set
    (which should be the default, in particular on the CI), the fixture will
    choose an arbitrary seed in the above range and all fixtured tests will run
    for that specific seed. This ensures that over time, our CI will run all
    tests with different seeds while keeping the test duration of a single run
    of the full test suite limited. This will enforce that the tests assertions
    of tests written to use this fixture are not dependent on a specific seed
    value.

    The range of admissible seed values is limited to [0, 99] because it is
    often not possible to write a test that can work for any possible seed and
    we want to avoid having tests that randomly fail on the CI.

    Valid values for SKLEARN_TESTS_GLOBAL_RANDOM_SEED:

    - SKLEARN_TESTS_GLOBAL_RANDOM_SEED="42": run tests with a fixed seed of 42
    - SKLEARN_TESTS_GLOBAL_RANDOM_SEED="40-42": run the tests with all seeds
      between 40 and 42 included
    - SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all": run the tests with all seeds
      between 0 and 99 included

    When writing a new test function that uses this fixture, please use the
    following command to make sure that it passes deterministically for all
    admissible seeds on your local machine:

        SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all" pytest -v -k test_your_test_name
    """
    yield request.param

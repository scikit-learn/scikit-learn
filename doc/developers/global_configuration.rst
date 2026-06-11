Global configuration and environment variables
==============================================

.. _environment_variable:

:func:`sklearn.set_config` and :func:`sklearn.config_context` can be used to
change global scikit-learn configuration at runtime.

The environment variables documented below should be set before importing
scikit-learn.

Performance and memory
----------------------

.. _envvar_SKLEARN_ASSUME_FINITE:

`SKLEARN_ASSUME_FINITE`
~~~~~~~~~~~~~~~~~~~~~~~

Sets the default value for the `assume_finite` argument of
:func:`sklearn.set_config`.

.. _envvar_SKLEARN_WORKING_MEMORY:

`SKLEARN_WORKING_MEMORY`
~~~~~~~~~~~~~~~~~~~~~~~~

Sets the default value for the `working_memory` argument of
:func:`sklearn.set_config`. The global default value is 1024.

.. _envvar_SKLEARN_PAIRWISE_DIST_CHUNK_SIZE:

`SKLEARN_PAIRWISE_DIST_CHUNK_SIZE`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This sets the size of chunk to be used by the underlying `PairwiseDistancesReductions`
implementations. The default value is `256` which has been showed to be adequate on
most machines.

Users looking for the best performance might want to tune this variable using
powers of 2 so as to get the best parallelism behavior for their hardware,
especially with respect to their caches' sizes.

Testing and CI
--------------

`SKLEARN_SEED`
~~~~~~~~~~~~~~

Sets the seed of the global random generator when running the tests, for
reproducibility.

Note that scikit-learn tests are expected to run deterministically with
explicit seeding of their own independent RNG instances instead of relying on
the numpy or Python standard library RNG singletons to make sure that test
results are independent of the test execution order. However some tests might
forget to use explicit seeding and this variable is a way to control the initial
state of the aforementioned singletons.

`SKLEARN_TESTS_GLOBAL_RANDOM_SEED`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Controls the seeding of the random number generator used in tests that rely on
the `global_random_seed` fixture.

All tests that use this fixture accept the contract that they should
deterministically pass for any seed value from 0 to 99 included.

In nightly CI builds, the `SKLEARN_TESTS_GLOBAL_RANDOM_SEED` environment
variable is drawn randomly in the above range and all fixtured tests will run
for that specific seed. The goal is to ensure that, over time, our CI will run
all tests with different seeds while keeping the test duration of a single run
of the full test suite limited. This will check that the assertions of tests
written to use this fixture are not dependent on a specific seed value.

The range of admissible seed values is limited to [0, 99] because it is often
not possible to write a test that can work for any possible seed and we want to
avoid having tests that randomly fail on the CI.

Valid values for `SKLEARN_TESTS_GLOBAL_RANDOM_SEED`:

- `SKLEARN_TESTS_GLOBAL_RANDOM_SEED="42"`: run tests with a fixed seed of 42
- `SKLEARN_TESTS_GLOBAL_RANDOM_SEED="40-42"`: run the tests with all seeds
  between 40 and 42 included
- `SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all"`: run the tests with all seeds
  between 0 and 99 included. This can take a long time: only use for individual
  tests, not the full test suite!

If the variable is not set, then 42 is used as the global seed in a
deterministic manner. This ensures that, by default, the scikit-learn test
suite is as deterministic as possible to avoid disrupting our friendly
third-party package maintainers. Similarly, this variable should not be set in
the CI config of pull-requests to make sure that our friendly contributors are
not the first people to encounter a seed-sensitivity regression in a test
unrelated to the changes of their own PR. Only the scikit-learn maintainers who
watch the results of the nightly builds are expected to be annoyed by this.

When writing a new test function that uses this fixture, please use the
following command to make sure that it passes deterministically for all
admissible seeds on your local machine:

.. prompt:: bash $

    SKLEARN_TESTS_GLOBAL_RANDOM_SEED="all" pytest -v -k test_your_test_name

`SKLEARN_SKIP_NETWORK_TESTS`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When this environment variable is set to ``0``, tests that need network access
are enabled. When unset or set to a non-zero value, network tests are skipped.

`SKLEARN_RUN_FLOAT32_TESTS`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When this environment variable is set to '1', the tests using the
`global_dtype` fixture are also run on float32 data.
When this environment variable is not set, the tests are only run on
float64 data.

`SKLEARN_WARNINGS_AS_ERRORS`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This environment variable is used to turn warnings into errors in tests and
documentation build.

Some CI (Continuous Integration) builds set `SKLEARN_WARNINGS_AS_ERRORS=1`, for
example to make sure that we catch deprecation warnings from our dependencies
and that we adapt our code.

To locally run with the same "warnings as errors" setting as in these CI builds
you can set `SKLEARN_WARNINGS_AS_ERRORS=1`.

By default, warnings are not turned into errors. This is the case if
`SKLEARN_WARNINGS_AS_ERRORS` is unset, or `SKLEARN_WARNINGS_AS_ERRORS=0`.

This environment variable uses specific warning filters to ignore some warnings,
since sometimes warnings originate from third-party libraries and there is not
much we can do about it. You can see the warning filters in the
`_get_warnings_filters_info_list` function in `sklearn/utils/_testing.py`.

Note that for documentation build, `SKLEARN_WARNINGS_AS_ERRORS=1` is checking
that the documentation build, in particular running examples, does not produce
any warnings. This is different from the `-W` `sphinx-build` argument that
catches syntax warnings in the rst files.

Build and debug
---------------

`SKLEARN_ENABLE_DEBUG_CYTHON_DIRECTIVES`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When this environment variable is set to a non zero value, the `Cython`
derivative, `boundscheck` is set to `True`. This is useful for finding
segfaults.

`SKLEARN_BUILD_ENABLE_DEBUG_SYMBOLS`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When this environment variable is set to a non zero value, the debug symbols
will be included in the compiled C extensions. Only debug symbols for POSIX
systems are configured.

`SKLEARN_FAIL_NO_OPENMP`
~~~~~~~~~~~~~~~~~~~~~~~~

If OpenMP is not supported by the compiler, the build will be done with OpenMP
functionalities disabled. This is not recommended since it will force some
estimators to run in sequential mode instead of leveraging thread-based
parallelism. Setting this environment variable (before cythonization) will
force the build to fail if OpenMP is not supported.

`SKLEARN_SKIP_OPENMP_TEST`
~~~~~~~~~~~~~~~~~~~~~~~~~~

When set to any value, skips the test that checks whether scikit-learn was built
with OpenMP-based parallelism enabled.

`SKLEARN_TEST_NO_OPENMP`
~~~~~~~~~~~~~~~~~~~~~~~~

When set to ``true``, used in CI to build and test scikit-learn on macOS without
OpenMP support.

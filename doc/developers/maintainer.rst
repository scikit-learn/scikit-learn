    Maintainer / core-developer information
========================================

Merging Pull Requests
---------------------

Individual commits are squashed when a Pull Request (PR) is merged on Github.
Before merging,

- the resulting commit title can be edited if necessary. Note
  that this will rename the PR title by default.
- the detailed description, containing the titles of all the commits, can
  be edited or deleted.
- for PRs with multiple code contributors care must be taken to keep
  the `Co-authored-by: name <name@example.com>` tags in the detailed
  description. This will mark the PR as having `multiple co-authors
  <https://help.github.com/en/github/committing-changes-to-your-project/creating-a-commit-with-multiple-authors>`_.
  Whether code contributions are significanly enough to merit co-authorship is
  left to the maintainer's discretion, same as for the "what's new" entry.


The scikit-learn.org web site
-----------------------------

The scikit-learn web site (https://scikit-learn.org) is hosted at GitHub,
but should rarely be updated manually by pushing to the
https://github.com/scikit-learn/scikit-learn.github.io repository. Most
updates can be made by pushing to master (for /dev) or a release branch
like 0.99.X, from which Circle CI builds and uploads the documentation
automatically.

Travis Cron jobs
----------------

From `<https://docs.travis-ci.com/user/cron-jobs>`_: Travis CI cron jobs work
similarly to the cron utility, they run builds at regular scheduled intervals
independently of whether any commits were pushed to the repository. Cron jobs
always fetch the most recent commit on a particular branch and build the project
at that state. Cron jobs can run daily, weekly or monthly, which in practice
means up to an hour after the selected time span, and you cannot set them to run
at a specific time.

For scikit-learn, Cron jobs are used for builds that we do not want to run in
each PR. As an example the build with the dev versions of numpy and scipy is
run as a Cron job. Most of the time when this numpy-dev build fail, it is
related to a numpy change and not a scikit-learn one, so it would not make sense
to blame the PR author for the Travis failure.

The definition of what gets run in the Cron job is done in the .travis.yml
config file, exactly the same way as the other Travis jobs. We use a ``if: type
= cron`` filter in order for the build to be run only in Cron jobs.

The branch targeted by the Cron job and the frequency of the Cron job is set
via the web UI at https://www.travis-ci.org/scikit-learn/scikit-learn/settings.

Experimental features
---------------------

The :mod:`sklearn.experimental` module was introduced in 0.21 and contains
experimental features / estimators that are subject to change without
deprecation cycle.

To create an experimental module, you can just copy and modify the content of
`enable_hist_gradient_boosting.py
<https://github.com/scikit-learn/scikit-learn/blob/c9c89cfc85dd8dfefd7921c16c87327d03140a06/sklearn/experimental/enable_hist_gradient_boosting.py>`__,
or
`enable_iterative_imputer.py
<https://github.com/scikit-learn/scikit-learn/blob/c9c89cfc85dd8dfefd7921c16c87327d03140a06/sklearn/experimental/enable_iterative_imputer.py>`_.

.. note::

  These are permalink as in 0.24, where these estimators are still
  experimental. They might be stable at the time of reading - hence the
  permalink. See below for instructions on the transition from experimental
  to stable.

Note that the public import path must be to a public subpackage (like
``sklearn/ensemble`` or ``sklearn/impute``), not just a ``.py`` module.
Also, the (private) experimental features that are imported must be in a
submodule/subpackage of the public subpackage, e.g.
``sklearn/ensemble/_hist_gradient_boosting/`` or
``sklearn/impute/_iterative.py``. This is needed so that pickles still work
in the future when the features aren't experimental anymore.

To avoid type checker (e.g. mypy) errors a direct import of experimental
estimators should be done in the parent module, protected by the
``if typing.TYPE_CHECKING`` check. See `sklearn/ensemble/__init__.py
<https://github.com/scikit-learn/scikit-learn/blob/c9c89cfc85dd8dfefd7921c16c87327d03140a06/sklearn/ensemble/__init__.py>`_,
or `sklearn/impute/__init__.py
<https://github.com/scikit-learn/scikit-learn/blob/c9c89cfc85dd8dfefd7921c16c87327d03140a06/sklearn/impute/__init__.py>`_
for an example.

Please also write basic tests following those in
`test_enable_hist_gradient_boosting.py
<https://github.com/scikit-learn/scikit-learn/blob/c9c89cfc85dd8dfefd7921c16c87327d03140a06/sklearn/experimental/tests/test_enable_hist_gradient_boosting.py>`__.


Make sure every user-facing code you write explicitly mentions that the feature
is experimental, and add a ``# noqa`` comment to avoid pep8-related warnings::

    # To use this experimental feature, we need to explicitly ask for it:
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    from sklearn.ensemble import HistGradientBoostingRegressor

For the docs to render properly, please also import
``enable_my_experimental_feature`` in ``doc/conf.py``, else sphinx won't be
able to import the corresponding modules. Note that using ``from
sklearn.experimental import *`` **does not work**.

Note that some experimental classes / functions are not included in the
:mod:`sklearn.experimental` module: ``sklearn.datasets.fetch_openml``.

Once the feature become stable, remove all `enable_my_experimental_feature`
in the scikit-learn code (even feature highlights etc.) and make the
`enable_my_experimental_feature` a no-op that just raises a warning:
`enable_hist_gradient_boosting.py
<https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/experimental/enable_hist_gradient_boosting.py>`__.
The file should stay there indefinitely as we don't want to break users code:
we just incentivize them to remove that import with the warning.

Also update the tests accordingly: `test_enable_hist_gradient_boosting.py
<https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/experimental/tests/test_enable_hist_gradient_boosting.py>`__.

Maintainer / core-developer information
========================================

A major release, incrementing the minor version, or a bug fix release
incrementing the patch version. Our convention is that we release one or more
release candidates (0.RRrcN) before releasing the final distributions. Creating
the first RC involves creating the release branch (0.99.X for instance)
directly on the main repo. This is considered the *feature freeze*. The rest of
the changes (required to fix some of the issues detected during the release
process are a good example) should be done as a PR to the master branch, if
it's possible, or done through a PR to the release branch, for the other
maintainers to review before the release. Note that all development for the
minor releases and the major release of for instance 0.99, happens under the
branch called ``0.99.X``, where ``X`` is a literal, and independent of minor
bug fix releases. Each release candidate, major, or minor release is a tag
under the same branch. This is further explained under :ref:
`preparing-a-release-pr`. We follow the `PEP101
<https://www.python.org/dev/peps/pep-0101/>`_ to indicate release candidates,
post, and minor releases.

Before a release
----------------

1. Update authors table::

    $ cd build_tools; make authors; cd ..

   and commit. This is only needed if the authors have changed since the last
   release. This step is sometimes done independent of the release. This
   updates the maintainer list and is not the contributor list for the release.

2. Confirm any blockers tagged for the milestone are resolved, and that other
   issues tagged for the milestone can be postponed.

3. Ensure the change log and commits correspond (within reason!), and that the
   change log is reasonably well curated. Some tools for these tasks include:

   - ``maint_tools/sort_whats_new.py`` can put what's new entries into
     sections. It's not perfect, and requires manual checking of the changes.
     If the whats new list if well curated, it may not be necessary.

   - The ``maint_tools/whats_missing.sh`` script may be used to identify pull
     requests that were merged but likely missing from What's New.

Permissions
...........

The release manager requires a set of permissions on top of the usual
permissions given to maintainers, which includes:

- *maintainer* role on ``scikit-learn`` projects on ``pypi.org`` and
  ``test.pypi.org``, separately.
- become a member of the *scikit-learn* team on conda-forge.
- *maintainer* on ``https://github.com/MacPython/scikit-learn-wheels``


.. _preparing-a-release-pr:

Preparing a release PR
......................

Since any commits to an existing branch (e.g. 0.999.X) will automatically
update the web site documentation, it is best to develop a release with a pull
request in which 0.999.X is the base. It also allows you to keep track of any
tasks towards release with a TO DO list. Since the documentation for the branch
is updated live, it should ideally match the version available to users. Having
all the changes for a release in one PR minimizes the gap between the time when
the website points to a release, and the release being available to users.

Most development of the release, and its documentation, should happen in master
to avoid asynchrony. To select commits from master for use in the bug fix
(version 0.999.3), you can use::

    $ # # assuming master and upstream/master are the same
    $ git checkout -b release-0.999.3 master
    $ git rebase -i upstream/0.999.X

It's nice to have a copy of the ``git rebase -i`` log in the PR to help others
understand what's included.

In terms of including changes, the first RC is a pure copy of the ``master``
with a release version number, and ideally counts as a *feature freeze*. Each
coming release candidate and the final release afterwards will include minor
documentation changes and bug fixes. Any major enhancement or feature should be
excluded.

The minor releases should include bug fixes and some relevant documentation
changes only. Any PR resulting in a behavior change which is not a bug fix
should be excluded.

Then pick the commits for release and resolve any issues, and create a pull
request with 0.999.X as base. Add a commit updating ``sklearn.__version__``.
Additional commits can be cherry-picked into the ``release-0.999.3`` branch
while preparing the release.

Making a release
----------------

1. Update docs. Note that this is for the final release, not necessarily for
   the RC releases.

   - Edit the doc/whats_new.rst file to add release title and commit
     statistics. You can retrieve commit statistics with::

        $ git shortlog -s 0.99.33.. | cut -f2- | sort --ignore-case | tr '\n' ';' | sed 's/;/, /g;s/, $//'

   - Update the release date in whats_new.rst

   - Edit the doc/templates/index.html to change the 'News' entry of the front page.

   - Note that these changes should be made in master and cherry-picked into
     the release branch, only before the final release.

2. On the branch for releasing, update the version number in
   sklearn/__init__.py, the ``__version__`` variable by removing ``dev*`` only
   when ready to release. On master, increment the version in the same place
   (when branching for release). This means while we're in the release
   candidate period, the latest stable is two versions behind the master
   branch, instead of one.

3. Create the source tarball:

   - Wipe clean your repo::

       $ git clean -xfd

   - Generate the tarball::

       $ python setup.py sdist

   - You can also test a binary dist build using::

       $ python setup.py bdist_wheel

   - You can test if PyPi is going to accept the package using::

       $ twine check dist/*

   You can run ``twine check`` after step 5 (fetching artifacts) as well.

   The result should be in the `dist/` folder. We will upload it later
   with the wheels. Check that you can install it in a new virtualenv and
   that the tests pass.

4. Proceed with caution. Ideally, tags should be created when you're almost
   certain that the release is ready to avoid premature tagging. You can test
   upload the sdist to test.pypi.org, and test the next step by setting
   ``BUILD_COMMIT`` to the branch name (``0.22.X`` for instance) in a PR to the
   wheel building repo. Once all works, you can proceed with tagging. Create
   the tag and push it (if it's an RC, it can be 0.xxrc1 for instance)::

    $ git tag -a 0.999

    $ git push git@github.com:scikit-learn/scikit-learn.git --tags

5. Update the dependency versions and set ``BUILD_COMMIT`` variable to the
   release tag at:

   https://github.com/MacPython/scikit-learn-wheels

   Once the CI has completed successfully, collect the generated binary wheel
   packages and upload them to PyPI by running the following commands in the
   scikit-learn source folder (checked out at the release tag)::

       $ rm -r dist # only if there's anything other than the sdist tar.gz there
       $ pip install -U wheelhouse_uploader twine
       $ python setup.py fetch_artifacts

6. Check the content of the `dist/` folder: it should contain all the wheels
   along with the source tarball ("scikit-learn-RRR.tar.gz").

   Make sure that you do not have developer versions or older versions of
   the scikit-learn package in that folder.

   Before uploading to pypi, you can test upload to test.pypi.org::

       $ twine upload --verbose --repository-url https://test.pypi.org/legacy/ dist/*

   Upload everything at once to https://pypi.org::

       $ twine upload dist/*

7. For major/minor (not bug-fix release), update the symlink for ``stable``
   and the ``latestStable`` variable in
   https://github.com/scikit-learn/scikit-learn.github.io::

       $ cd /tmp
       $ git clone --depth 1 --no-checkout git@github.com:scikit-learn/scikit-learn.github.io.git
       $ cd scikit-learn.github.io
       $ echo stable > .git/info/sparse-checkout
       $ git checkout master
       $ rm stable
       $ ln -s 0.999 stable
       $ sed -i "s/latestStable = '.*/latestStable = '0.999';/g" versionwarning.js
       $ git commit -am "Update stable to point to 0.999"
       $ git push origin master

The following GitHub checklist might be helpful in a release PR::

    * [ ] update news and what's new date in master and release branch
    * [ ] create tag
    * [ ] update dependencies and release tag at https://github.com/MacPython/scikit-learn-wheels
    * [ ] twine the wheels to PyPI when that's green
    * [ ] https://github.com/scikit-learn/scikit-learn/releases draft
    * [ ] confirm bot detected at https://github.com/conda-forge/scikit-learn-feedstock and wait for merge
    * [ ] https://github.com/scikit-learn/scikit-learn/releases publish
    * [ ] announce on mailing list and on twitter
    * [ ] (regenerate Dash docs: https://github.com/Kapeli/Dash-User-Contributions/tree/master/docsets/Scikit)

The scikit-learn.org web site
-----------------------------

The scikit-learn web site (http://scikit-learn.org) is hosted at GitHub,
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
<https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/experimental/enable_hist_gradient_boosting.py>`_,
or
`enable_iterative_imputer.py
<https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/experimental/enable_iterative_imputer.py>`_.

Note that the public import path must be to a public subpackage (like
``sklearn/ensemble`` or ``sklearn/impute``), not just a ``.py`` module.
Also, the (private) experimental features that are imported must be in a
submodule/subpackage of the public subpackage, e.g.
``sklearn/ensemble/_hist_gradient_boosting/`` or
``sklearn/impute/_iterative.py``. This is needed so that pickles still work
in the future when the features aren't experimental anymore

Please also write basic tests following those in
`test_enable_hist_gradient_boosting.py
<https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/experimental/tests/test_enable_hist_gradient_boosting.py>`_.

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

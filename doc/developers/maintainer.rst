Maintainer / core-developer information
========================================

Before a release
----------------

1. Update authors table::

    $ cd build_tools; make authors; cd ..

   and commit.

2. Confirm any blockers tagged for the milestone are resolved, and that other
   issues tagged for the milestone can be postponed.

3. Ensure the change log and commits correspond (within reason!), and that the
   change log is reasonably well curated. Some tools for these tasks include:

   - ``maint_tools/sort_whats_new.py`` can put what's new entries into
     sections.

   - The ``maint_tools/whats_missing.sh`` script may be used to identify pull
     requests that were merged but likely missing from What's New.

Preparing a bug-fix-release
...........................

Since any commits to a released branch (e.g. 0.999.X) will automatically update
the web site documentation, it is best to develop a bug-fix release with a pull
request in which 0.999.X is the base. It also allows you to keep track of any
tasks towards release with a TO DO list.

Most development of the bug fix release, and its documentation, should
happen in master to avoid asynchrony. To select commits from master for use in
the bug fix (version 0.999.3), you can use::

    $ git checkout -b release-0.999.3 master
    $ git rebase -i 0.999.X

Then pick the commits for release and resolve any issues, and create a pull
request with 0.999.X as base. Add a commit updating ``sklearn.__version__``.
Additional commits can be cherry-picked into the ``release-0.999.3`` branch
while preparing the release.

Making a release
----------------

1. Update docs:

   - Edit the doc/whats_new.rst file to add release title and commit
     statistics. You can retrieve commit statistics with::

        $ git shortlog -ns 0.998..

   - Update the release date in whats_new.rst

   - Edit the doc/index.rst to change the 'News' entry of the front page.

   - Note that these changes should be made in master and cherry-picked into
     the release branch.

2. On the branch for releasing, update the version number in
   sklearn/__init__.py, the ``__version__`` variable by removing ``dev*`` only
   when ready to release.
   On master, increment the verson in the same place (when branching for
   release).

3. Create the tag and push it::

    $ git tag -a 0.999

    $ git push git@github.com:scikit-learn/scikit-learn.git --tags

4. Create the source tarball:

   - Wipe clean your repo::

       $ git clean -xfd

   - Generate the tarball::

       $ python setup.py sdist

   The result should be in the `dist/` folder. We will upload it later
   with the wheels. Check that you can install it in a new virtualenv and
   that the tests pass.

5. Update the dependency versions and set ``BUILD_COMMIT`` variable to the
   release tag at:

   https://github.com/MacPython/scikit-learn-wheels

   Once the CI has completed successfully, collect the generated binary wheel
   packages and upload them to PyPI by running the following commands in the
   scikit-learn source folder (checked out at the release tag)::

       $ pip install -U wheelhouse_uploader twine
       $ python setup.py fetch_artifacts

6. Check the content of the `dist/` folder: it should contain all the wheels
   along with the source tarball ("scikit-learn-XXX.tar.gz").

   Make sure that you do not have developer versions or older versions of
   the scikit-learn package in that folder.

   Upload everything at once to https://pypi.org::

       $ twine upload dist/

7. For major/minor (not bug-fix release), update the symlink for ``stable``
   in https://github.com/scikit-learn/scikit-learn.github.io::

       $ cd /tmp
       $ git clone --depth 1 --no-checkout git@github.com:scikit-learn/scikit-learn.github.io.git
       $ cd scikit-learn.github.io
       $ echo stable > .git/info/sparse-checkout
       $ git checkout master
       $ ln -sf 0.999 stable
       $ git push origin master

The following GitHub checklist might be helpful in a release PR::

    * [ ] update news and what's new date in master and release branch
    * [ ] create tag
    * [ ] update dependencies and release tag at https://github.com/MacPython/scikit-learn-wheels
    * [ ] twine the wheels to PyPI when that's green
    * [ ] https://github.com/scikit-learn/scikit-learn/releases draft
    * [ ] confirm bot detected at https://github.com/conda-forge/scikit-learn-feedstock and wait for merge
    * [ ] https://github.com/scikit-learn/scikit-learn/releases publish
    * [ ] announce on mailing list
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

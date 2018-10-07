Maintainer / core-developer information
========================================

Before a release
----------------

1. Update authors table::

    $ cd build_tools; make authors; cd ..

   and commit.

Making a release
----------------
For more information see https://github.com/scikit-learn/scikit-learn/wiki/How-to-make-a-release


1. Update docs:

   - Edit the doc/whats_new.rst file to add release title and commit
     statistics. You can retrieve commit statistics with::

        $ git shortlog -ns 0.998..

   - Edit the doc/index.rst to change the 'News' entry of the front page.

2. Update the version number in sklearn/__init__.py, the __version__
   variable

3. Create the tag and push it::

    $ git tag 0.999

    $ git push origin --tags

4. create tarballs:

   - Wipe clean your repo::

       $ git clean -xfd

   - Register and upload on PyPI::

       $ python setup.py sdist register upload


5. Push the documentation to the website. Circle CI should do this
   automatically for master and <N>.<N>.X branches.

6. Build binaries using dedicated CI servers by updating the git submodule
   reference to the new scikit-learn tag of the release at:

   https://github.com/MacPython/scikit-learn-wheels

   Once the CI has completed successfully, collect the generated binary wheel
   packages and upload them to PyPI by running the following commands in the
   scikit-learn source folder (checked out at the release tag)::

       $ pip install -U wheelhouse_uploader
       $ python setup.py sdist fetch_artifacts upload_all


7. FOR FINAL RELEASE: Update the release date in What's New

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

.. _contributing:

============
Contributing
============

This project is a community effort, and everyone is welcome to
contribute.

The project is hosted on https://github.com/scikit-learn/scikit-learn

The decision making process and governance structure of scikit-learn is laid
out in the governance document: :ref:`governance`.

Scikit-learn is somewhat :ref:`selective <selectiveness>` when it comes to
adding new algorithms, and the best way to contribute and to help the project
is to start working on known issues.
See :ref:`new_contributors` to get started.

.. topic:: **Our community, our values**

    We are a community based on openness and friendly, didactic,
    discussions.

    We aspire to treat everybody equally, and value their contributions.

    Decisions are made based on technical merit and consensus.

    Code is not the only way to help the project. Reviewing pull
    requests, answering questions to help others on mailing lists or
    issues, organizing and teaching tutorials, working on the website,
    improving the documentation, are all priceless contributions.

    We abide by the principles of openness, respect, and consideration of
    others of the Python Software Foundation:
    https://www.python.org/psf/codeofconduct/


In case you experience issues using this package, do not hesitate to submit a
ticket to the
`GitHub issue tracker
<https://github.com/scikit-learn/scikit-learn/issues>`_. You are also
welcome to post feature requests or pull requests.

Ways to contribute
==================

There are many ways to contribute to scikit-learn, with the most common ones
being contribution of code or documentation to the project. Improving the
documentation is no less important than improving the library itself.  If you
find a typo in the documentation, or have made improvements, do not hesitate to
send an email to the mailing list or preferably submit a GitHub pull request.
Full documentation can be found under the doc/ directory.

But there are many other ways to help. In particular answering queries on the
`issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_,
investigating bugs, and :ref:`reviewing other developers' pull requests
<code_review>` are very valuable contributions that decrease the burden on the
project maintainers.

Another way to contribute is to report issues you're facing, and give a "thumbs
up" on issues that others reported and that are relevant to you.  It also helps
us if you spread the word: reference the project from your blog and articles,
link to it from your website, or simply star to say "I use it":

In case a contribution/issue involves changes to the API principles
or changes to dependencies or supported versions, it must be backed by a
:ref:`slep`, where a SLEP must be submitted as a pull-request to
`enhancement proposals <https://scikit-learn-enhancement-proposals.readthedocs.io>`_
using the `SLEP template <https://scikit-learn-enhancement-proposals.readthedocs.io/en/latest/slep_template.html>`_
and follows the decision-making process outlined in :ref:`governance`.

.. raw:: html

   <a class="github-button" href="https://github.com/scikit-learn/scikit-learn"
   data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star
   scikit-learn/scikit-learn on GitHub">Star</a>
   <script async defer src="https://buttons.github.io/buttons.js"></script>

.. topic:: Contributing to related projects

   Scikit-learn thrives in an ecosystem of several related projects, which also
   may have relevant issues to work on, including smaller projects such as:

   * `scikit-learn-contrib <https://github.com/search?q=org%3Ascikit-learn-contrib+is%3Aissue+is%3Aopen+sort%3Aupdated-desc&type=Issues>`__
   * `joblib <https://github.com/joblib/joblib/issues>`__
   * `sphinx-gallery <https://github.com/sphinx-gallery/sphinx-gallery/issues>`__
   * `numpydoc <https://github.com/numpy/numpydoc/issues>`__
   * `liac-arff <https://github.com/renatopp/liac-arff>`__

   and larger projects:

   * `numpy <https://github.com/numpy/numpy/issues>`__
   * `scipy <https://github.com/scipy/scipy/issues>`__
   * `matplotlib <https://github.com/matplotlib/matplotlib/issues>`__
   * and so on.

   Look for issues marked "help wanted" or similar.
   Helping these projects may help Scikit-learn too.
   See also :ref:`related_projects`.


Submitting a bug report or a feature request
============================================

We use GitHub issues to track all bugs and feature requests; feel free to open
an issue if you have found a bug or wish to see a feature implemented.

In case you experience issues using this package, do not hesitate to submit a
ticket to the
`Bug Tracker <https://github.com/scikit-learn/scikit-learn/issues>`_. You are
also welcome to post feature requests or pull requests.

It is recommended to check that your issue complies with the
following rules before submitting:

-  Verify that your issue is not being currently addressed by other
   `issues <https://github.com/scikit-learn/scikit-learn/issues?q=>`_
   or `pull requests <https://github.com/scikit-learn/scikit-learn/pulls?q=>`_.

-  If you are submitting an algorithm or feature request, please verify that
   the algorithm fulfills our
   `new algorithm requirements
   <http://scikit-learn.org/stable/faq.html#what-are-the-inclusion-criteria-for-new-algorithms>`_.

-  If you are submitting a bug report, we strongly encourage you to follow the guidelines in
   :ref:`filing_bugs`.

.. _filing_bugs:

How to make a good bug report
-----------------------------

When you submit an issue to `Github
<https://github.com/scikit-learn/scikit-learn/issues>`__, please do your best to
follow these guidelines! This will make it a lot easier to provide you with good
feedback:

- The ideal bug report contains a **short reproducible code snippet**, this way
  anyone can try to reproduce the bug easily (see `this
  <https://stackoverflow.com/help/mcve>`_ for more details). If your snippet is
  longer than around 50 lines, please link to a `gist
  <https://gist.github.com>`_ or a github repo.

- If not feasible to include a reproducible snippet, please be specific about
  what **estimators and/or functions are involved and the shape of the data**.

- If an exception is raised, please **provide the full traceback**.

- Please include your **operating system type and version number**, as well as
  your **Python, scikit-learn, numpy, and scipy versions**. This information
  can be found by running the following code snippet::

    >>> import sklearn
    >>> sklearn.show_versions()  # doctest: +SKIP

  .. note::

    This utility function is only available in scikit-learn v0.20+.
    For previous versions, one has to explicitly run::

     import platform; print(platform.platform())
     import sys; print("Python", sys.version)
     import numpy; print("NumPy", numpy.__version__)
     import scipy; print("SciPy", scipy.__version__)
     import sklearn; print("Scikit-Learn", sklearn.__version__)

- Please ensure all **code snippets and error messages are formatted in
  appropriate code blocks**.  See `Creating and highlighting code blocks
  <https://help.github.com/articles/creating-and-highlighting-code-blocks>`_
  for more details.

Contributing code
=================

.. note::

  To avoid duplicating work, it is highly advised that you search through the
  `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_ and
  the `PR list <https://github.com/scikit-learn/scikit-learn/pulls>`_.
  If in doubt about duplicated work, or if you want to work on a non-trivial
  feature, it's recommended to first open an issue in
  the `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_
  to get some feedbacks from core developers.

How to contribute
-----------------

The preferred way to contribute to scikit-learn is to fork the `main
repository <https://github.com/scikit-learn/scikit-learn/>`__ on GitHub,
then submit a "pull request" (PR).

In the first few steps, we explain how to locally install scikit-learn, and
how to set up your git repository:

1. `Create an account <https://github.com/join>`_ on
   GitHub if you do not already have one.

2. Fork the `project repository
   <https://github.com/scikit-learn/scikit-learn>`__: click on the 'Fork'
   button near the top of the page. This creates a copy of the code under your
   account on the GitHub user account. For more details on how to fork a
   repository see `this guide <https://help.github.com/articles/fork-a-repo/>`_.

3. Clone your fork of the scikit-learn repo from your GitHub account to your
   local disk::

       $ git clone git@github.com:YourLogin/scikit-learn.git
       $ cd scikit-learn

4. Install the development dependencies::

       $ pip install cython pytest flake8

5. Install scikit-learn in editable mode::

       $ pip install --editable .

   for more details about advanced installation, see the
   :ref:`install_bleeding_edge` section.

6. Add the ``upstream`` remote. This saves a reference to the main
   scikit-learn repository, which you can use to keep your repository
   synchronized with the latest changes::

    $ git remote add upstream https://github.com/scikit-learn/scikit-learn.git

You should now have a working installation of scikit-learn, and your git
repository properly configured. The next steps now describe the process of
modifying code and submitting a PR:

7. Synchronize your master branch with the upstream master branch::

        $ git checkout master
        $ git pull upstream master

8. Create a feature branch to hold your development changes::

        $ git checkout -b my_feature

   and start making changes. Always use a feature branch. It's good
   practice to never work on the ``master`` branch!

9. Develop the feature on your feature branch on your computer, using Git to
   do the version control. When you're done editing, add changed files using
   ``git add`` and then ``git commit``::

       $ git add modified_files
       $ git commit

   to record your changes in Git, then push the changes to your GitHub
   account with::

       $ git push -u origin my-feature

10. Follow `these
    <https://help.github.com/articles/creating-a-pull-request-from-a-fork>`_
    instructions to create a pull request from your fork. This will send an
    email to the committers. You may want to consider sending an email to the
    mailing list for more visibility.

.. note::

  If you are modifying a Cython module, you have to re-run step 5 after modifications
  and before testing them.

It is often helpful to keep your local feature branch synchronized with the
latest changes of the main scikit-learn repository::

    $ git fetch upstream
    $ git merge upstream/master

Subsequently, you might need to solve the conflicts. You can refer to the
`Git documentation related to resolving merge conflict using the command
line
<https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/>`_.

.. topic:: Learning git:

    The `Git documentation <https://git-scm.com/documentation>`_ and
    http://try.github.io are excellent resources to get started with git,
    and understanding all of the commands shown here.

Pull request checklist
----------------------

Before a PR can be merged, it needs to be approved by two core developers.
Please prefix the title of your pull request with ``[MRG]`` if the
contribution is complete and should be subjected to a detailed review. An
incomplete contribution -- where you expect to do more work before receiving
a full review -- should be prefixed ``[WIP]`` (to indicate a work in
progress) and changed to ``[MRG]`` when it matures. WIPs may be useful to:
indicate you are working on something to avoid duplicated work, request
broad review of functionality or API, or seek collaborators. WIPs often
benefit from the inclusion of a `task list
<https://github.com/blog/1375-task-lists-in-gfm-issues-pulls-comments>`_ in
the PR description.

In order to ease the reviewing process, we recommend that your contribution
complies with the following rules before marking a PR as ``[MRG]``. The
**bolded** ones are especially important:

1. **Give your pull request a helpful title** that summarises what your
   contribution does. This title will often become the commit message once
   merged so it should summarise your contribution for posterity. In some
   cases "Fix <ISSUE TITLE>" is enough. "Fix #<ISSUE NUMBER>" is never a
   good title.

2. **Make sure your code passes the tests**. The whole test suite can be run
   with `pytest`, but it is usually not recommended since it takes a long
   time. It is often enough to only run the test related to your changes:
   for example, if you changed something in
   `sklearn/linear_model/logistic.py`, running the following commands will
   usually be enough:

   - `pytest sklearn/linear_model/logistic.py` to make sure the doctest
     examples are correct
   - `pytest sklearn/linear_model/tests/test_logistic.py` to run the tests
     specific to the file
   - `pytest sklearn/linear_model` to test the whole `linear_model` module
   - `pytest sklearn/doc/linear_model.rst` to make sure the user guide
     examples are correct.
   - `pytest sklearn/tests/test_common.py -k LogisticRegression` to run all our
     estimator checks (specifically for `LogisticRegression`, if that's the
     estimator you changed).

   There may be other failing tests, but they will be caught by the CI so
   you don't need to run the whole test suite locally. You can read more in
   :ref:`testing_coverage`.

3. **Make sure your code is properly commented and documented**, and **make
   sure the documentation renders properly**. To build the documentation, please
   refer to our :ref:`contribute_documentation` guidelines. The CI will also
   build the docs: please refer to :ref:`generated_doc_CI`.

4. **Tests are necessary for enhancements to be
   accepted**. Bug-fixes or new features should be provided with
   `non-regression tests
   <https://en.wikipedia.org/wiki/Non-regression_testing>`_. These tests
   verify the correct behavior of the fix or feature. In this manner, further
   modifications on the code base are granted to be consistent with the
   desired behavior. In the case of bug fixes, at the time of the PR, the
   non-regression tests should fail for the code base in the master branch
   and pass for the PR code.

5. **Make sure that your PR does not add PEP8 violations**. On a Unix-like
   system, you can run `make flake8-diff`. `flake8 path_to_file`, would work
   for any system, but please avoid reformatting parts of the file that your
   pull request doesn't change, as it distracts from code review.

6. Follow the `coding-guidelines`_ (see below).

7. When applicable, use the validation tools and scripts in the
   ``sklearn.utils`` submodule.  A list of utility routines available
   for developers can be found in the :ref:`developers-utils` page.

8. Often pull requests resolve one or more other issues (or pull requests).
   If merging your pull request means that some other issues/PRs should
   be closed, you should `use keywords to create link to them
   <https://github.com/blog/1506-closing-issues-via-pull-requests/>`_
   (e.g., ``Fixes #1234``; multiple issues/PRs are allowed as long as each
   one is preceded by a keyword). Upon merging, those issues/PRs will
   automatically be closed by GitHub. If your pull request is simply
   related to some other issues/PRs, create a link to them without using
   the keywords (e.g., ``See also #1234``).

9. PRs should often substantiate the change, through benchmarks of
   performance and efficiency or through examples of usage. Examples also
   illustrate the features and intricacies of the library to users. Have a
   look at other examples in the `examples/
   <https://github.com/scikit-learn/scikit-learn/tree/master/examples>`_
   directory for reference. Examples should demonstrate why the new
   functionality is useful in practice and, if possible, compare it to other
   methods available in scikit-learn.

10. New features often need to be illustrated with narrative documentation in
    the user guide, with small code snipets. If relevant, please also add
    references in the literature, with PDF links when possible.

11. The user guide should also include expected time and space complexity
    of the algorithm and scalability, e.g. "this algorithm can scale to a
    large number of samples > 100000, but does not scale in dimensionality:
    n_features is expected to be lower than 100".

You can also check our :ref:`code_review` to get an idea of what reviewers
will expect.

You can check for common programming errors with the following tools:

* Code with a good unittest coverage (at least 80%, better 100%), check
  with::

    $ pip install pytest pytest-cov
    $ pytest --cov sklearn path/to/tests_for_package

  see also :ref:`testing_coverage`

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the GitHub issue).

Also check out the :ref:`performance-howto` guide for more details on
profiling and Cython optimizations.

.. note::

  The current state of the scikit-learn code base is not compliant with
  all of those guidelines, but we expect that enforcing those constraints
  on all new contributions will get the overall code base quality in the
  right direction.

.. note::

   For two very well documented and more detailed guides on development
   workflow, please pay a visit to the `Scipy Development Workflow
   <https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html>`_ -
   and the `Astropy Workflow for Developers
   <https://astropy.readthedocs.io/en/latest/development/workflow/development_workflow.html>`_
   sections.

Continuous Integration (CI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Azure pipelines are used for testing scikit-learn on Linux, Mac and Windows,
  with different dependencies and settings.
* CircleCI is used to build the docs for viewing, for linting with flake8, and
  for testing with PyPy on Linux

Please note that if one of the following markers appear in the latest commit
message, the following actions are taken.

    ====================== ===================
    Commit Message Marker  Action Taken by CI
    ---------------------- -------------------
    [scipy-dev]            Add a Travis build with our dependencies (numpy, scipy, etc ...) development builds
    [ci skip]              CI is skipped completely
    [doc skip]             Docs are not built
    [doc quick]            Docs built, but excludes example gallery plots
    [doc build]            Docs built including example gallery plots
    ====================== ===================

Stalled pull requests
^^^^^^^^^^^^^^^^^^^^^

As contributing a feature can be a lengthy process, some
pull requests appear inactive but unfinished. In such a case, taking
them over is a great service for the project.

A good etiquette to take over is:

* **Determine if a PR is stalled**

  * A pull request may have the label "stalled" or "help wanted" if we
    have already identified it as a candidate for other contributors.

  * To decide whether an inactive PR is stalled, ask the contributor if
    she/he plans to continue working on the PR in the near future.
    Failure to respond within 2 weeks with an activity that moves the PR
    forward suggests that the PR is stalled and will result in tagging
    that PR with "help wanted".

    Note that if a PR has received earlier comments on the contribution
    that have had no reply in a month, it is safe to assume that the PR
    is stalled and to shorten the wait time to one day.

    After a sprint, follow-up for un-merged PRs opened during sprint will
    be communicated to participants at the sprint, and those PRs will be
    tagged "sprint". PRs tagged with "sprint" can be reassigned or
    declared stalled by sprint leaders.

* **Taking over a stalled PR**: To take over a PR, it is important to
  comment on the stalled PR that you are taking over and to link from the
  new PR to the old one. The new PR should be created by pulling from the
  old one.

.. _new_contributors:

Issues for New Contributors
---------------------------

New contributors should look for the following tags when looking for issues.  We
strongly recommend that new contributors tackle "easy" issues first: this helps
the contributor become familiar with the contribution workflow, and for the core
devs to become acquainted with the contributor; besides which, we frequently
underestimate how easy an issue is to solve!

.. topic:: good first issue tag

    A great way to start contributing to scikit-learn is to pick an item from
    the list of `good first issues
    <https://github.com/scikit-learn/scikit-learn/labels/good%20first%20issue>`_
    in the issue tracker. Resolving these issues allow you to start contributing
    to the project without much prior knowledge. If you have already contributed
    to scikit-learn, you should look at Easy issues instead.

.. topic:: Easy tag

    If you have already contributed to scikit-learn, another great way to contribute
    to scikit-learn is to pick an item from the list of `Easy issues
    <https://github.com/scikit-learn/scikit-learn/labels/Easy>`_ in the issue
    tracker. Your assistance in this area will be greatly appreciated by the
    more experienced developers as it helps free up their time to concentrate on
    other issues.

.. topic:: help wanted tag

    We often use the help wanted tag to mark issues regardless of difficulty. Additionally,
    we use the help wanted tag to mark Pull Requests which have been abandoned
    by their original contributor and are available for someone to pick up where the original
    contributor left off. The list of issues with the help wanted tag can be found
    `here <https://github.com/scikit-learn/scikit-learn/labels/help%20wanted>`__ .

    Note that not all issues which need contributors will have this tag.

.. _contribute_documentation:

Documentation
=============

We are glad to accept any sort of documentation: function docstrings,
reStructuredText documents (like this one), tutorials, etc. reStructuredText
documents live in the source code repository under the ``doc/`` directory.

You can edit the documentation using any text editor, and then generate the
HTML output by typing ``make`` from the ``doc/`` directory. Alternatively,
``make html`` may be used to generate the documentation **with** the example
gallery (which takes quite some time). The resulting HTML files will be
placed in ``_build/html/stable`` and are viewable in a web browser.


Building the documentation
--------------------------

First, make sure you have :ref:`properly installed <install_bleeding_edge>`
the development version.

Building the documentation requires installing some additional packages::

    pip install sphinx sphinx-gallery numpydoc matplotlib Pillow pandas scikit-image

To build the documentation, you need to be in the ``doc`` folder::

    cd doc

In the vast majority of cases, you only need to generate the full web site,
without the example gallery::

    make

The documentation will be generated in the ``_build/html/stable`` directory.
To also generate the example gallery you can use::

    make html

This will run all the examples, which takes a while. If you only want to
generate a few examples, you can use::

    EXAMPLES_PATTERN=your_regex_goes_here make html

This is particularly useful if you are modifying a few examples.

Set the environment variable `NO_MATHJAX=1` if you intend to view
the documentation in an offline setting.

To build the PDF manual, run::

    make latexpdf

.. warning:: **Sphinx version**

   While we do our best to have the documentation build under as many
   versions of Sphinx as possible, the different versions tend to
   behave slightly differently. To get the best results, you should
   use the same version as the one we used on CircleCI. Look at this
   `github search <https://github.com/search?utf8=%E2%9C%93&q=sphinx+repo%3Ascikit-learn%2Fscikit-learn+extension%3Ash+path%3Abuild_tools%2Fcircle&type=Code>`_
   to know the exact version.

Guidelines for writing documentation
------------------------------------

It is important to keep a good compromise between mathematical and algorithmic
details, and give intuition to the reader on what the algorithm does.

Basically, to elaborate on the above, it is best to always
start with a small paragraph with a hand-waving explanation of what the
method does to the data. Then, it is very helpful to point out why the feature is
useful and when it should be used - the latter also including "big O"
(:math:`O\left(g\left(n\right)\right)`) complexities of the algorithm, as opposed
to just *rules of thumb*, as the latter can be very machine-dependent. If those
complexities are not available, then rules of thumb may be provided instead.

Secondly, a generated figure from an example (as mentioned in the previous
paragraph) should then be included to further provide some intuition.

Next, one or two small code examples to show its use can be added.

Next, any math and equations, followed by references,
can be added to further the documentation. Not starting the
documentation with the maths makes it more friendly towards
users that are just interested in what the feature will do, as
opposed to how it works "under the hood".

Finally, follow the formatting rules below to make it consistently good:

* Add "See also" in docstrings for related classes/functions.

* "See also" in docstrings should be one line per reference,
  with a colon and an explanation, for example::

    See also
    --------
    SelectKBest : Select features based on the k highest scores.
    SelectFpr : Select features based on a false positive rate test.

* When documenting the parameters and attributes, here is a list of some
  well-formatted examples::

    n_clusters : int, default=3
        The number of clusters detected by the algorithm.

    some_param : {'hello', 'goodbye'}, bool or int, default=True
        The parameter description goes here, which can be either a string
        literal (either `hello` or `goodbye`), a bool, or an int. The default
        value is True.

    array_parameter : {array-like, sparse matrix, dataframe} of shape (n_samples, n_features) or (n_samples,)
        This parameter accepts data in either of the mentioned forms, with one
        of the mentioned shapes. The default value is
        `np.ones(shape=(n_samples,))`.

In general have the following in mind:

    1. Use Python basic types. (``bool`` instead of ``boolean``)
    2. Use parenthesis for defining shapes: ``array-like of shape (n_samples,)``
       or ``array-like of shape (n_samples, n_features)``
    3. For strings with multiple options, use brackets:
       ``input: {'log', 'squared', 'multinomial'}``
    4. 1D or 2D data can be a subset of
       ``{array-like, ndarray, sparse matrix, dataframe}``. Note that ``array-like``
       can also be a ``list``, while ``ndarray`` is explicitly only a ``numpy.ndarray``.

* For unwritten formatting rules, try to follow existing good works:

    * For "References" in docstrings, see the Silhouette Coefficient
      (:func:`sklearn.metrics.silhouette_score`).

* When editing reStructuredText (``.rst``) files, try to keep line length under
  80 characters when possible (exceptions include links and tables).

.. _generated_doc_CI:

Generated documentation on CircleCI
-----------------------------------

When you change the documentation in a pull request, CircleCI automatically
builds it. To view the documentation generated by CircleCI:

* navigate to the bottom of your pull request page to see the CI
  statuses. You may need to click on "Show all checks" to see all the CI
  statuses.
* click on the CircleCI status with "doc" in the title.
* add ``#artifacts`` at the end of the URL. Note: you need to wait for the
  CircleCI build to finish before being able to look at the artifacts.
* once the artifacts are visible, navigate to ``doc/_changed.html`` to see a
  list of documentation pages that are likely to be affected by your pull
  request. Navigate to ``doc/index.html`` to see the full generated html
  documentation.

If you often need to look at the documentation generated by CircleCI, e.g. when
reviewing pull requests, you may find :ref:`this tip
<viewing_rendered_html_documentation>` very handy.

.. _testing_coverage:

Testing and improving test coverage
===================================

High-quality `unit testing <https://en.wikipedia.org/wiki/Unit_testing>`_
is a corner-stone of the scikit-learn development process. For this
purpose, we use the `pytest <https://docs.pytest.org>`_
package. The tests are functions appropriately named, located in `tests`
subdirectories, that check the validity of the algorithms and the
different options of the code.

The full scikit-learn tests can be run using 'make' in the root folder.
Alternatively, running 'pytest' in a folder will run all the tests of
the corresponding subpackages.

We expect code coverage of new features to be at least around 90%.

For guidelines on how to use ``pytest`` efficiently, see the
:ref:`pytest_tips`.

Writing matplotlib related tests
--------------------------------

Test fixtures ensure that a set of tests will be executing with the appropriate
initialization and cleanup. The scikit-learn test suite implements a fixture
which can be used with ``matplotlib``.

``pyplot``
    The ``pyplot`` fixture should be used when a test function is dealing with
    ``matplotlib``. ``matplotlib`` is a soft dependency and is not required.
    This fixture is in charge of skipping the tests if ``matplotlib`` is not
    installed. In addition, figures created during the tests will be
    automatically closed once the test function has been executed.

To use this fixture in a test function, one needs to pass it as an
argument::

    def test_requiring_mpl_fixture(pyplot):
        # you can now safely use matplotlib

Workflow to improve test coverage
---------------------------------

To test code coverage, you need to install the `coverage
<https://pypi.org/project/coverage/>`_ package in addition to pytest.

1. Run 'make test-coverage'. The output lists for each file the line
    numbers that are not tested.

2. Find a low hanging fruit, looking at which lines are not tested,
    write or adapt a test specifically for these lines.

3. Loop.

Issue Tracker Tags
==================

All issues and pull requests on the
`GitHub issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_
should have (at least) one of the following tags:

:Bug / Crash:
    Something is happening that clearly shouldn't happen.
    Wrong results as well as unexpected errors from estimators go here.

:Cleanup / Enhancement:
    Improving performance, usability, consistency.

:Documentation:
    Missing, incorrect or sub-standard documentations and examples.

:New Feature:
    Feature requests and pull requests implementing a new feature.

There are four other tags to help new contributors:

:good first issue:
    This issue is ideal for a first contribution to scikit-learn. Ask for help
    if the formulation is unclear. If you have already contributed to
    scikit-learn, look at Easy issues instead.

:Easy:
    This issue can be tackled without much prior experience.

:Moderate:
    Might need some knowledge of machine learning or the package,
    but is still approachable for someone new to the project.

:help wanted:
    This tag marks an issue which currently lacks a contributor or a
    PR that needs another contributor to take over the work. These
    issues can range in difficulty, and may not be approachable
    for new contributors. Note that not all issues which need
    contributors will have this tag.


.. _coding-guidelines:

Coding guidelines
=================

The following are some guidelines on how new code should be written. Of
course, there are special cases and there will be exceptions to these
rules. However, following these rules when submitting new code makes
the review easier so new code can be integrated in less time.

Uniformly formatted code makes it easier to share code ownership. The
scikit-learn project tries to closely follow the official Python guidelines
detailed in `PEP8 <https://www.python.org/dev/peps/pep-0008>`_ that
detail how code should be formatted and indented. Please read it and
follow it.

In addition, we add the following guidelines:

* Use underscores to separate words in non class names: ``n_samples``
  rather than ``nsamples``.

* Avoid multiple statements on one line. Prefer a line return after
  a control flow statement (``if``/``for``).

* Use relative imports for references inside scikit-learn.

* Unit tests are an exception to the previous rule;
  they should use absolute imports, exactly as client code would.
  A corollary is that, if ``sklearn.foo`` exports a class or function
  that is implemented in ``sklearn.foo.bar.baz``,
  the test should import it from ``sklearn.foo``.

* **Please don't use** ``import *`` **in any case**. It is considered harmful
  by the `official Python recommendations
  <https://docs.python.org/2/howto/doanddont.html#from-module-import>`_.
  It makes the code harder to read as the origin of symbols is no
  longer explicitly referenced, but most important, it prevents
  using a static analysis tool like `pyflakes
  <https://divmod.readthedocs.io/en/latest/products/pyflakes.html>`_ to automatically
  find bugs in scikit-learn.

* Use the `numpy docstring standard
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  in all your docstrings.


A good example of code that we like can be found `here
<https://gist.github.com/nateGeorge/5455d2c57fb33c1ae04706f2dc4fee01>`_.

.. _backwards-compatibility:

Maintaining backwards compatibility
===================================

.. _contributing_deprecation:

Deprecation
-----------

If any publicly accessible method, function, attribute or parameter
is renamed, we still support the old one for two releases and issue
a deprecation warning when it is called/passed/accessed.
E.g., if the function ``zero_one`` is renamed to ``zero_one_loss``,
we add the decorator ``deprecated`` (from ``sklearn.utils``)
to ``zero_one`` and call ``zero_one_loss`` from that function::

    from ..utils import deprecated

    def zero_one_loss(y_true, y_pred, normalize=True):
        # actual implementation
        pass

    @deprecated("Function 'zero_one' was renamed to 'zero_one_loss' "
                "in version 0.13 and will be removed in release 0.15. "
                "Default behavior is changed from 'normalize=False' to "
                "'normalize=True'")
    def zero_one(y_true, y_pred, normalize=False):
        return zero_one_loss(y_true, y_pred, normalize)

If an attribute is to be deprecated,
use the decorator ``deprecated`` on a property. Please note that the
``property`` decorator should be placed before the ``deprecated``
decorator for the docstrings to be rendered properly.
E.g., renaming an attribute ``labels_`` to ``classes_`` can be done as::

    @deprecated("Attribute labels_ was deprecated in version 0.13 and "
                "will be removed in 0.15. Use 'classes_' instead")
    @property
    def labels_(self):
        return self.classes_

If a parameter has to be deprecated, use ``DeprecationWarning`` appropriately.
In the following example, k is deprecated and renamed to n_clusters::

    import warnings

    def example_function(n_clusters=8, k='deprecated'):
        if k != 'deprecated':
            warnings.warn("'k' was renamed to n_clusters in version 0.13 and "
                          "will be removed in 0.15.", DeprecationWarning)
            n_clusters = k

When the change is in a class, we validate and raise warning in ``fit``::

  import warnings

  class ExampleEstimator(BaseEstimator):
      def __init__(self, n_clusters=8, k='deprecated'):
          self.n_clusters = n_clusters
          self.k = k

      def fit(self, X, y):
          if self.k != 'deprecated':
              warnings.warn("'k' was renamed to n_clusters in version 0.13 and "
                            "will be removed in 0.15.", DeprecationWarning)
              self._n_clusters = self.k
          else:
              self._n_clusters = self.n_clusters

As in these examples, the warning message should always give both the
version in which the deprecation happened and the version in which the
old behavior will be removed. If the deprecation happened in version
0.x-dev, the message should say deprecation occurred in version 0.x and
the removal will be in 0.(x+2), so that users will have enough time to
adapt their code to the new behaviour. For example, if the deprecation happened
in version 0.18-dev, the message should say it happened in version 0.18
and the old behavior will be removed in version 0.20.

In addition, a deprecation note should be added in the docstring, recalling the
same information as the deprecation warning as explained above. Use the
``.. deprecated::`` directive::

  .. deprecated:: 0.13
     ``k`` was renamed to ``n_clusters`` in version 0.13 and will be removed
     in 0.15.

What's more, a deprecation requires a test which ensures that the warning is
raised in relevant cases but not in other cases. The warning should be caught
in all other tests (using e.g., ``@pytest.mark.filterwarnings``),
and there should be no warning in the examples.


Change the default value of a parameter
---------------------------------------

If the default value of a parameter needs to be changed, please replace the
default value with a specific value (e.g., ``warn``) and raise
``FutureWarning`` when users are using the default value. In the following
example, we change the default value of ``n_clusters`` from 5 to 10
(current version is 0.20)::

    import warnings

    def example_function(n_clusters='warn'):
        if n_clusters == 'warn':
            warnings.warn("The default value of n_clusters will change from "
                          "5 to 10 in 0.22.", FutureWarning)
            n_clusters = 5

When the change is in a class, we validate and raise warning in ``fit``::

  import warnings

  class ExampleEstimator:
      def __init__(self, n_clusters='warn'):
          self.n_clusters = n_clusters

      def fit(self, X, y):
          if self.n_clusters == 'warn':
            warnings.warn("The default value of n_clusters will change from "
                          "5 to 10 in 0.22.", FutureWarning)
            self._n_clusters = 5

Similar to deprecations, the warning message should always give both the
version in which the change happened and the version in which the old behavior
will be removed. The docstring needs to be updated accordingly. We need a test
which ensures that the warning is raised in relevant cases but not in other
cases. The warning should be caught in all other tests
(using e.g., ``@pytest.mark.filterwarnings``), and there should be no warning
in the examples.

.. currentmodule:: sklearn

.. _code_review:

Code Review Guidelines
======================
Reviewing code contributed to the project as PRs is a crucial component of
scikit-learn development. We encourage anyone to start reviewing code of other
developers. The code review process is often highly educational for everybody
involved. This is particularly appropriate if it is a feature you would like to
use, and so can respond critically about whether the PR meets your needs. While
each pull request needs to be signed off by two core developers, you can speed
up this process by providing your feedback.

Here are a few important aspects that need to be covered in any code review,
from high-level questions to a more detailed check-list.

- Do we want this in the library? Is it likely to be used? Do you, as
  a scikit-learn user, like the change and intend to use it? Is it in
  the scope of scikit-learn? Will the cost of maintaining a new
  feature be worth its benefits?

- Is the code consistent with the API of scikit-learn? Are public
  functions/classes/parameters well named and intuitively designed?

- Are all public functions/classes and their parameters, return types, and
  stored attributes named according to scikit-learn conventions and documented clearly?

- Is any new functionality described in the user-guide and illustrated with examples?

- Is every public function/class tested? Are a reasonable set of
  parameters, their values, value types, and combinations tested? Do
  the tests validate that the code is correct, i.e. doing what the
  documentation says it does? If the change is a bug-fix, is a
  non-regression test included? Look at `this
  <https://jeffknupp.com/blog/2013/12/09/improve-your-python-understanding-unit-testing>`__
  to get started with testing in Python.

- Do the tests pass in the continuous integration build? If
  appropriate, help the contributor understand why tests failed.

- Do the tests cover every line of code (see the coverage report in the build
  log)? If not, are the lines missing coverage good exceptions?

- Is the code easy to read and low on redundancy? Should variable names be
  improved for clarity or consistency? Should comments be added? Should comments
  be removed as unhelpful or extraneous?

- Could the code easily be rewritten to run much more efficiently for
  relevant settings?

- Is the code backwards compatible with previous versions? (or is a
  deprecation cycle necessary?)

- Will the new code add any dependencies on other libraries? (this is
  unlikely to be accepted)

- Does the documentation render properly (see the
  :ref:`contribute_documentation` section for more details), and are the plots
  instructive?

:ref:`saved_replies` includes some frequent comments that reviewers may make.

Reading the existing code base
==============================

Reading and digesting an existing code base is always a difficult exercise
that takes time and experience to master. Even though we try to write simple
code in general, understanding the code can seem overwhelming at first,
given the sheer size of the project. Here is a list of tips that may help
make this task easier and faster (in no particular order).

- Get acquainted with the :ref:`api_overview`: understand what :term:`fit`,
  :term:`predict`, :term:`transform`, etc. are used for.
- Before diving into reading the code of a function / class, go through the
  docstrings first and try to get an idea of what each parameter / attribute
  is doing. It may also help to stop a minute and think *how would I do this
  myself if I had to?*
- The trickiest thing is often to identify which portions of the code are
  relevant, and which are not. In scikit-learn **a lot** of input checking
  is performed, especially at the beginning of the :term:`fit` methods.
  Sometimes, only a very small portion of the code is doing the actual job.
  For example looking at the ``fit()`` method of
  :class:`sklearn.linear_model.LinearRegression`, what you're looking for
  might just be the call the ``scipy.linalg.lstsq``, but it is buried into
  multiple lines of input checking and the handling of different kinds of
  parameters.
- Due to the use of `Inheritance
  <https://en.wikipedia.org/wiki/Inheritance_(object-oriented_programming)>`_,
  some methods may be implemented in parent classes. All estimators inherit
  at least from :class:`BaseEstimator <sklearn.base.BaseEstimator>`, and
  from a ``Mixin`` class (e.g. :class:`ClassifierMixin
  <sklearn.base.ClassifierMixin>`) that enables default behaviour depending
  on the nature of the estimator (classifier, regressor, transformer, etc.).
- Sometimes, reading the tests for a given function will give you an idea of
  what its intended purpose is. You can use ``git grep`` (see below) to find
  all the tests written for a function. Most tests for a specific
  function/class are placed under the ``tests/`` folder of the module
- You'll often see code looking like this:
  ``out = Parallel(...)(delayed(some_function)(param) for param in
  some_iterable)``. This runs ``some_function`` in parallel using `Joblib
  <https://joblib.readthedocs.io/>`_. ``out`` is then an iterable containing
  the values returned by ``some_function`` for each call.
- We use `Cython <https://cython.org/>`_ to write fast code. Cython code is
  located in ``.pyx`` and ``.pxd`` files. Cython code has a more C-like
  flavor: we use pointers, perform manual memory allocation, etc. Having
  some minimal experience in C / C++ is pretty much mandatory here.
- Master your tools.

  - With such a big project, being efficient with your favorite editor or
    IDE goes a long way towards digesting the code base. Being able to quickly
    jump (or *peek*) to a function/class/attribute definition helps a lot.
    So does being able to quickly see where a given name is used in a file.
  - `git <https://git-scm.com/book/en>`_ also has some built-in killer
    features. It is often useful to understand how a file changed over time,
    using e.g. ``git blame`` (`manual
    <https://git-scm.com/docs/git-blame>`_). This can also be done directly
    on GitHub. ``git grep`` (`examples
    <https://git-scm.com/docs/git-grep#_examples>`_) is also extremely
    useful to see every occurrence of a pattern (e.g. a function call or a
    variable) in the code base.


.. _plotting_api:

Plotting API
============

Scikit-learn defines a simple API for creating visualizations for machine
learning. The key features of this API is to run calculations once and to have
the flexibility to adjust the visualizations after the fact. This logic is
encapsulated into a display object where the computed data is stored and
the plotting is done in a `plot` method. The display object's `__init__`
method contains only the data needed to create the visualization. The `plot`
method takes in parameters that only have to do with visualization, such as a
matplotlib axes. The `plot` method will store the matplotlib artists as
attributes allowing for style adjustments through the display object. A
`plot_*` helper function accepts parameters to do the computation and the
parameters used for plotting. After the helper function creates the display
object with the computed values, it calls the display's plot method. Note
that the `plot` method defines attributes related to matplotlib, such as the
line artist. This allows for customizations after calling the `plot` method.

For example, the `RocCurveDisplay` defines the following methods and
attributes:

.. code-block:: python

   class RocCurveDisplay:
       def __init__(self, fpr, tpr, roc_auc, estimator_name):
           ...
           self.fpr = fpr
           self.tpr = tpr
           self.roc_auc = roc_auc
           self.estimator_name = estimator_name

       def plot(self, ax=None, name=None, **kwargs):
           ...
           self.line_ = ...
           self.ax_ = ax
           self.figure_ = ax.figure_

   def plot_roc_curve(estimator, X, y, pos_label=None, sample_weight=None,
                      drop_intermediate=True, response_method="auto",
                      name=None, ax=None, **kwargs):
       # do computation
       viz = RocCurveDisplay(fpr, tpr, roc_auc, 
                                estimator.__class__.__name__)
       return viz.plot(ax=ax, name=name, **kwargs)

Read more in the :ref:`User Guide <visualizations>`.

.. _contributing:

============
Contributing
============

This project is a community effort, and everyone is welcome to
contribute.

The project is hosted on https://github.com/scikit-learn/scikit-learn

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

|


In case you experience issues using this package, do not hesitate to submit a
ticket to the
`GitHub issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_. You are
also welcome to post feature requests or pull requests.


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
then submit a "pull request" (PR):

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

 4. Install library in editable mode::

        $ pip install --editable .

    for more details about advanced installation, see the
    :ref:`install_bleeding_edge` section.

 5. Create a branch to hold your development changes::

        $ git checkout -b my-feature

    and start making changes. Always use a ``feature`` branch. It's good practice to
    never work on the ``master`` branch!

 6. Develop the feature on your feature branch on your computer, using Git to do the
    version control. When you're done editing, add changed files using ``git add``
    and then ``git commit`` files::

        $ git add modified_files
        $ git commit

    to record your changes in Git, then push the changes to your GitHub account with::

        $ git push -u origin my-feature

 7. Follow `these
    <https://help.github.com/articles/creating-a-pull-request-from-a-fork>`_
    instructions to create a pull request from your fork. This will send an
    email to the committers. You may want to consider sending an email to the
    mailing list for more visibility.

.. note::

  If you are modifying a Cython module, you have to re-run step 4 after modifications
  and before testing them.

.. note::

  In the above setup, your ``origin`` remote repository points to
  YourLogin/scikit-learn.git. If you wish to fetch/merge from the main
  repository instead of your forked one, you will need to add another remote
  to use instead of ``origin``. If we choose the name ``upstream`` for it, the
  command will be::

        $ git remote add upstream https://github.com/scikit-learn/scikit-learn.git

If any of the above seems like magic to you, then look up the `Git documentation
<https://git-scm.com/documentation>`_ and the `Git development workflow
<http://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html>`_ on the
web, or ask a friend or another contributor for help.

If some conflicts arise between your branch and the ``master`` branch, you need
to merge ``master``. The command will be::

  $ git merge master

with ``master`` being synchronized with the ``upstream``.

Subsequently, you need to solve the conflicts. You can refer to the `Git
documentation related to resolving merge conflict using the command line
<https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/>`_.

.. note::

   In the past, the policy to resolve conflicts was to rebase your branch on
   ``master``. GitHub interface deals with merging ``master`` better than in
   the past.


Contributing pull requests
--------------------------

We recommend that that your contribution complies with the following
rules before submitting a pull request:

* Follow the `coding-guidelines`_ (see below). To make sure that
  your PR does not add PEP8 violations you can run
  `./build_tools/circle/flake8_diff.sh` or `make flake8-diff` on a
  Unix-like system.

* When applicable, use the validation tools and scripts in the
  ``sklearn.utils`` submodule.  A list of utility routines available
  for developers can be found in the :ref:`developers-utils` page.

* Give your pull request a helpful title that summarises what your
  contribution does. In some cases "Fix <ISSUE TITLE>" is enough.
  "Fix #<ISSUE NUMBER>" is not enough.

* Often pull requests resolve one or more other issues (or pull requests).
  If merging your pull request means that some other issues/PRs should
  be closed, you should `use keywords to create link to them
  <https://github.com/blog/1506-closing-issues-via-pull-requests/>`_
  (e.g., ``Fixes #1234``; multiple issues/PRs are allowed as long as each
  one is preceded by a keyword). Upon merging, those issues/PRs will
  automatically be closed by GitHub. If your pull request is simply
  related to some other issues/PRs, create a link to them without using
  the keywords (e.g., ``See also #1234``).

* All public methods should have informative docstrings with sample
  usage presented as doctests when appropriate.

* Please prefix the title of your pull request with ``[MRG]`` if the
  contribution is complete and should be subjected to a detailed review.
  Two core developers will review your code and change the prefix of the pull
  request to ``[MRG + 1]`` and ``[MRG + 2]`` on approval, making it eligible
  for merging. An incomplete contribution -- where you expect to do more
  work before receiving a full review -- should be prefixed ``[WIP]`` (to
  indicate a work in progress) and changed to ``[MRG]`` when it matures.
  WIPs may be useful to: indicate you are working on something to avoid
  duplicated work, request broad review of functionality or API, or seek
  collaborators. WIPs often benefit from the inclusion of a
  `task list
  <https://github.com/blog/1375-task-lists-in-gfm-issues-pulls-comments>`_
  in the PR description.

* All other tests pass when everything is rebuilt from scratch. On
  Unix-like systems, check with (from the toplevel source folder)::

    $ make

* When adding additional functionality, provide at least one example script
  in the ``examples/`` folder. Have a look at other examples for reference.
  Examples should demonstrate why the new functionality is useful in
  practice and, if possible, compare it to other methods available in
  scikit-learn.

* Documentation and high-coverage tests are necessary for enhancements to be
  accepted. Bug-fixes or new features should be provided with
  `non-regression tests
  <https://en.wikipedia.org/wiki/Non-regression_testing>`_. These tests
  verify the correct behavior of the fix or feature. In this manner, further
  modifications on the code base are granted to be consistent with the
  desired behavior. For the case of bug fixes, at the time of the PR, the
  non-regression tests should fail for the code base in the master branch
  and pass for the PR code.

* At least one paragraph of narrative documentation with links to
  references in the literature (with PDF links when possible) and
  the example. For more details on writing and building the
  documentation, see the :ref:`contribute_documentation` section.

* The documentation should also include expected time and space complexity
  of the algorithm and scalability, e.g. "this algorithm can scale to a
  large number of samples > 100000, but does not scale in dimensionality:
  n_features is expected to be lower than 100".

You can also check for common programming errors with the following tools:

* Code with a good unittest coverage (at least 80%, better 100%), check
  with::

    $ pip install pytest pytest-cov
    $ pytest --cov sklearn path/to/tests_for_package

  see also :ref:`testing_coverage`

* No flake8 warnings, check with::

    $ pip install flake8
    $ flake8 path/to/module.py

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the GitHub issue).

Also check out the :ref:`performance-howto` guide for more details on profiling
and Cython optimizations.

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

* Travis is used for testing on Linux platforms
* Appveyor is used for testing on Windows platforms
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
-------------

We are glad to accept any sort of documentation: function docstrings,
reStructuredText documents (like this one), tutorials, etc. reStructuredText
documents live in the source code repository under the ``doc/`` directory.

You can edit the documentation using any text editor, and then generate the
HTML output by typing ``make html`` from the ``doc/`` directory. Alternatively,
``make`` can be used to quickly generate the documentation without the example
gallery. The resulting HTML files will be placed in ``_build/html/stable`` and are viewable
in a web browser. See the ``README``file in the ``doc/`` directory for more information.


Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Building the documentation requires installing some additional packages::

    pip install sphinx sphinx-gallery numpydoc matplotlib Pillow pandas scikit-image joblib

To build the documentation, you need to be in the ``doc`` folder::

    cd doc

It also requires having the version of scikit-learn installed that corresponds
to the documentation, e.g.::

    pip install --editable ..

To generate the full web site, including the example gallery::

    make html

Generating the example gallery will run all our examples which takes a
while. To save some time, you can use:

- ``make html-noplot``: this will generate the documentation without the
  example gallery. This is useful when changing a docstring for example.
- ``EXAMPLES_PATTERN=your_regex_goes_here make html``: only the examples
  matching ``your_regex_goes_here`` will be run. This is particularly
  useful if you are modifying a few examples.

That should create all the documentation in the ``_build/html/stable``
directory.  Set the environment variable `NO_MATHJAX=1` if you intend to view
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

* For unwritten formatting rules, try to follow existing good works:

    * For "References" in docstrings, see the Silhouette Coefficient
      (:func:`sklearn.metrics.silhouette_score`).

* When editing reStructuredText (``.rst``) files, try to keep line length under
  80 characters when possible (exceptions include links and tables).

Generated documentation on CircleCI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
------------------------------------

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

.. note:: **Workflow to improve test coverage**

   To test code coverage, you need to install the `coverage
   <https://pypi.org/project/coverage/>`_ package in addition to pytest.

   1. Run 'make test-coverage'. The output lists for each file the line
      numbers that are not tested.

   2. Find a low hanging fruit, looking at which lines are not tested,
      write or adapt a test specifically for these lines.

   3. Loop.



Developers web site
-------------------

More information can be found on the `developer's wiki
<https://github.com/scikit-learn/scikit-learn/wiki>`_.


Issue Tracker Tags
------------------
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

Input validation
----------------

.. currentmodule:: sklearn.utils

The module :mod:`sklearn.utils` contains various functions for doing input
validation and conversion. Sometimes, ``np.asarray`` suffices for validation;
do *not* use ``np.asanyarray`` or ``np.atleast_2d``, since those let NumPy's
``np.matrix`` through, which has a different API
(e.g., ``*`` means dot product on ``np.matrix``,
but Hadamard product on ``np.ndarray``).

In other cases, be sure to call :func:`check_array` on any array-like argument
passed to a scikit-learn API function. The exact parameters to use depends
mainly on whether and which ``scipy.sparse`` matrices must be accepted.

For more information, refer to the :ref:`developers-utils` page.

Random Numbers
--------------

If your code depends on a random number generator, do not use
``numpy.random.random()`` or similar routines.  To ensure
repeatability in error checking, the routine should accept a keyword
``random_state`` and use this to construct a
``numpy.random.RandomState`` object.
See :func:`sklearn.utils.check_random_state` in :ref:`developers-utils`.

Here's a simple example of code using some of the above guidelines::

    from sklearn.utils import check_array, check_random_state

    def choose_random_sample(X, random_state=0):
        """
        Choose a random point from X

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            array representing the data
        random_state : RandomState or an int seed (0 by default)
            A random number generator instance to define the state of the
            random permutations generator.

        Returns
        -------
        x : numpy array, shape (n_features,)
            A random point selected from X
        """
        X = check_array(X)
        random_state = check_random_state(random_state)
        i = random_state.randint(X.shape[0])
        return X[i]

If you use randomness in an estimator instead of a freestanding function,
some additional guidelines apply.

First off, the estimator should take a ``random_state`` argument to its
``__init__`` with a default value of ``None``.
It should store that argument's value, **unmodified**,
in an attribute ``random_state``.
``fit`` can call ``check_random_state`` on that attribute
to get an actual random number generator.
If, for some reason, randomness is needed after ``fit``,
the RNG should be stored in an attribute ``random_state_``.
The following example should make this clear::

    class GaussianNoise(BaseEstimator, TransformerMixin):
        """This estimator ignores its input and returns random Gaussian noise.

        It also does not adhere to all scikit-learn conventions,
        but showcases how to handle randomness.
        """

        def __init__(self, n_components=100, random_state=None):
            self.random_state = random_state

        # the arguments are ignored anyway, so we make them optional
        def fit(self, X=None, y=None):
            self.random_state_ = check_random_state(self.random_state)

        def transform(self, X):
            n_samples = X.shape[0]
            return self.random_state_.randn(n_samples, n_components)

The reason for this setup is reproducibility:
when an estimator is ``fit`` twice to the same data,
it should produce an identical model both times,
hence the validation in ``fit``, not ``__init__``.

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
use the decorator ``deprecated`` on a property.
E.g., renaming an attribute ``labels_`` to ``classes_`` can be done as::

    @property
    @deprecated("Attribute labels_ was deprecated in version 0.13 and "
                "will be removed in 0.15. Use 'classes_' instead")
    def labels_(self):
        return self.classes_

If a parameter has to be deprecated, use ``DeprecationWarning`` appropriately.
In the following example, k is deprecated and renamed to n_clusters::

    import warnings

    def example_function(n_clusters=8, k='not_used'):
        if k != 'not_used':
            warnings.warn("'k' was renamed to n_clusters in version 0.13 and "
                          "will be removed in 0.15.", DeprecationWarning)
            n_clusters = k

When the change is in a class, we validate and raise warning in ``fit``::

  import warnings

  class ExampleEstimator(BaseEstimator):
      def __init__(self, n_clusters=8, k='not_used'):
          self.n_clusters = n_clusters
          self.k = k

      def fit(self, X, y):
          if self.k != 'not_used':
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

Python versions supported
-------------------------

Since scikit-learn 0.21, only Python 3.5 and newer is supported.

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

.. _api_overview:

APIs of scikit-learn objects
============================

To have a uniform API, we try to have a common basic API for all the
objects. In addition, to avoid the proliferation of framework code, we
try to adopt simple conventions and limit to a minimum the number of
methods an object must implement.

Elements of the scikit-learn API are described more definitively in the
:ref:`glossary`.

Different objects
-----------------

The main objects in scikit-learn are (one class can implement
multiple interfaces):

:Estimator:

    The base object, implements a ``fit`` method to learn from data, either::

      estimator = estimator.fit(data, targets)

    or::

      estimator = estimator.fit(data)

:Predictor:

    For supervised learning, or some unsupervised problems, implements::

      prediction = predictor.predict(data)

    Classification algorithms usually also offer a way to quantify certainty
    of a prediction, either using ``decision_function`` or ``predict_proba``::

      probability = predictor.predict_proba(data)

:Transformer:

    For filtering or modifying the data, in a supervised or unsupervised
    way, implements::

      new_data = transformer.transform(data)

    When fitting and transforming can be performed much more efficiently
    together than separately, implements::

      new_data = transformer.fit_transform(data)

:Model:

    A model that can give a `goodness of fit <https://en.wikipedia.org/wiki/Goodness_of_fit>`_
    measure or a likelihood of unseen data, implements (higher is better)::

      score = model.score(data)

Estimators
----------

The API has one predominant object: the estimator. A estimator is an
object that fits a model based on some training data and is capable of
inferring some properties on new data. It can be, for instance, a
classifier or a regressor. All estimators implement the fit method::

    estimator.fit(X, y)

All built-in estimators also have a ``set_params`` method, which sets
data-independent parameters (overriding previous parameter values passed
to ``__init__``).

All estimators in the main scikit-learn codebase should inherit from
``sklearn.base.BaseEstimator``.

Instantiation
^^^^^^^^^^^^^

This concerns the creation of an object. The object's ``__init__`` method
might accept constants as arguments that determine the estimator's behavior
(like the C constant in SVMs). It should not, however, take the actual training
data as an argument, as this is left to the ``fit()`` method::

    clf2 = SVC(C=2.3)
    clf3 = SVC([[1, 2], [2, 3]], [-1, 1]) # WRONG!


The arguments accepted by ``__init__`` should all be keyword arguments
with a default value. In other words, a user should be able to instantiate
an estimator without passing any arguments to it. The arguments should all
correspond to hyperparameters describing the model or the optimisation
problem the estimator tries to solve. These initial arguments (or parameters)
are always remembered by the estimator.
Also note that they should not be documented under the "Attributes" section,
but rather under the "Parameters" section for that estimator.

In addition, **every keyword argument accepted by** ``__init__`` **should
correspond to an attribute on the instance**. Scikit-learn relies on this to
find the relevant attributes to set on an estimator when doing model selection.

To summarize, an ``__init__`` should look like::

    def __init__(self, param1=1, param2=2):
        self.param1 = param1
        self.param2 = param2

There should be no logic, not even input validation,
and the parameters should not be changed.
The corresponding logic should be put where the parameters are used,
typically in ``fit``.
The following is wrong::

    def __init__(self, param1=1, param2=2, param3=3):
        # WRONG: parameters should not be modified
        if param1 > 1:
            param2 += 1
        self.param1 = param1
        # WRONG: the object's attributes should have exactly the name of
        # the argument in the constructor
        self.param3 = param2

The reason for postponing the validation is that the same validation
would have to be performed in ``set_params``,
which is used in algorithms like ``GridSearchCV``.

Fitting
^^^^^^^

The next thing you will probably want to do is to estimate some
parameters in the model. This is implemented in the ``fit()`` method.

The ``fit()`` method takes the training data as arguments, which can be one
array in the case of unsupervised learning, or two arrays in the case
of supervised learning.

Note that the model is fitted using X and y, but the object holds no
reference to X and y. There are, however, some exceptions to this, as in
the case of precomputed kernels where this data must be stored for use by
the predict method.

============= ======================================================
Parameters
============= ======================================================
X             array-like, shape (n_samples, n_features)

y             array, shape (n_samples,)

kwargs        optional data-dependent parameters.
============= ======================================================

``X.shape[0]`` should be the same as ``y.shape[0]``. If this requisite
is not met, an exception of type ``ValueError`` should be raised.

``y`` might be ignored in the case of unsupervised learning. However, to
make it possible to use the estimator as part of a pipeline that can
mix both supervised and unsupervised transformers, even unsupervised
estimators need to accept a ``y=None`` keyword argument in
the second position that is just ignored by the estimator.
For the same reason, ``fit_predict``, ``fit_transform``, ``score``
and ``partial_fit`` methods need to accept a ``y`` argument in
the second place if they are implemented.

The method should return the object (``self``). This pattern is useful
to be able to implement quick one liners in an IPython session such as::

  y_predicted = SVC(C=100).fit(X_train, y_train).predict(X_test)

Depending on the nature of the algorithm, ``fit`` can sometimes also
accept additional keywords arguments. However, any parameter that can
have a value assigned prior to having access to the data should be an
``__init__`` keyword argument. **fit parameters should be restricted
to directly data dependent variables**. For instance a Gram matrix or
an affinity matrix which are precomputed from the data matrix ``X`` are
data dependent. A tolerance stopping criterion ``tol`` is not directly
data dependent (although the optimal value according to some scoring
function probably is).

When ``fit`` is called, any previous call to ``fit`` should be ignored. In
general, calling ``estimator.fit(X1)`` and then ``estimator.fit(X2)`` should
be the same as only calling ``estimator.fit(X2)``. However, this may not be
true in practice when ``fit`` depends on some random process, see
:term:`random_state`. Another exception to this rule is when the
hyper-parameter ``warm_start`` is set to ``True`` for estimators that
support it. ``warm_start=True`` means that the previous state of the
trainable parameters of the estimator are reused instead of using the
default initialization strategy.

Estimated Attributes
^^^^^^^^^^^^^^^^^^^^

Attributes that have been estimated from the data must always have a name
ending with trailing underscore, for example the coefficients of
some regression estimator would be stored in a ``coef_`` attribute after
``fit`` has been called.

The estimated attributes are expected to be overridden when you call ``fit``
a second time.

Optional Arguments
^^^^^^^^^^^^^^^^^^

In iterative algorithms, the number of iterations should be specified by
an integer called ``n_iter``.

Pairwise Attributes
^^^^^^^^^^^^^^^^^^^

An estimator that accept ``X`` of shape ``(n_samples, n_samples)`` and defines
a :term:`_pairwise` property equal to ``True`` allows for cross-validation of
the dataset, e.g. when ``X`` is a precomputed kernel matrix. Specifically,
the :term:`_pairwise` property is used by ``utils.metaestimators._safe_split``
to slice rows and columns.

Rolling your own estimator
==========================
If you want to implement a new estimator that is scikit-learn-compatible,
whether it is just for you or for contributing it to scikit-learn, there are
several internals of scikit-learn that you should be aware of in addition to
the scikit-learn API outlined above. You can check whether your estimator
adheres to the scikit-learn interface and standards by running
:func:`utils.estimator_checks.check_estimator` on the class::

  >>> from sklearn.utils.estimator_checks import check_estimator
  >>> from sklearn.svm import LinearSVC
  >>> check_estimator(LinearSVC)  # passes

The main motivation to make a class compatible to the scikit-learn estimator
interface might be that you want to use it together with model evaluation and
selection tools such as :class:`model_selection.GridSearchCV` and
:class:`pipeline.Pipeline`.

Before detailing the required interface below, we describe two ways to achieve
the correct interface more easily.

.. topic:: Project template:

    We provide a `project template <https://github.com/scikit-learn-contrib/project-template/>`_
    which helps in the creation of Python packages containing scikit-learn compatible estimators.
    It provides:

    * an initial git repository with Python package directory structure
    * a template of a scikit-learn estimator
    * an initial test suite including use of ``check_estimator``
    * directory structures and scripts to compile documentation and example
      galleries
    * scripts to manage continuous integration (testing on Linux and Windows)
    * instructions from getting started to publishing on `PyPi <https://pypi.org/>`_

.. topic:: ``BaseEstimator`` and mixins:

    We tend to use "duck typing", so building an estimator which follows
    the API suffices for compatibility, without needing to inherit from or
    even import any scikit-learn classes.

    However, if a dependency on scikit-learn is acceptable in your code,
    you can prevent a lot of boilerplate code
    by deriving a class from ``BaseEstimator``
    and optionally the mixin classes in ``sklearn.base``.
    For example, below is a custom classifier, with more examples included
    in the scikit-learn-contrib
    `project template <https://github.com/scikit-learn-contrib/project-template/blob/master/skltemplate/_template.py>`__.

      >>> import numpy as np
      >>> from sklearn.base import BaseEstimator, ClassifierMixin
      >>> from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
      >>> from sklearn.utils.multiclass import unique_labels
      >>> from sklearn.metrics import euclidean_distances
      >>> class TemplateClassifier(BaseEstimator, ClassifierMixin):
      ...
      ...     def __init__(self, demo_param='demo'):
      ...         self.demo_param = demo_param
      ...
      ...     def fit(self, X, y):
      ...
      ...         # Check that X and y have correct shape
      ...         X, y = check_X_y(X, y)
      ...         # Store the classes seen during fit
      ...         self.classes_ = unique_labels(y)
      ...
      ...         self.X_ = X
      ...         self.y_ = y
      ...         # Return the classifier
      ...         return self
      ...
      ...     def predict(self, X):
      ...
      ...         # Check is fit had been called
      ...         check_is_fitted(self, ['X_', 'y_'])
      ...
      ...         # Input validation
      ...         X = check_array(X)
      ...
      ...         closest = np.argmin(euclidean_distances(X, self.X_), axis=1)
      ...         return self.y_[closest]


get_params and set_params
-------------------------
All scikit-learn estimators have ``get_params`` and ``set_params`` functions.
The ``get_params`` function takes no arguments and returns a dict of the
``__init__`` parameters of the estimator, together with their values.
It must take one keyword argument, ``deep``,
which receives a boolean value that determines
whether the method should return the parameters of sub-estimators
(for most estimators, this can be ignored).
The default value for ``deep`` should be true.

The ``set_params`` on the other hand takes as input a dict of the form
``'parameter': value`` and sets the parameter of the estimator using this dict.
Return value must be estimator itself.

While the ``get_params`` mechanism is not essential (see :ref:`cloning` below),
the ``set_params`` function is necessary as it is used to set parameters during
grid searches.

The easiest way to implement these functions, and to get a sensible
``__repr__`` method, is to inherit from ``sklearn.base.BaseEstimator``. If you
do not want to make your code dependent on scikit-learn, the easiest way to
implement the interface is::

    def get_params(self, deep=True):
        # suppose this estimator has parameters "alpha" and "recursive"
        return {"alpha": self.alpha, "recursive": self.recursive}

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self


Parameters and init
-------------------
As :class:`model_selection.GridSearchCV` uses ``set_params``
to apply parameter setting to estimators,
it is essential that calling ``set_params`` has the same effect
as setting parameters using the ``__init__`` method.
The easiest and recommended way to accomplish this is to
**not do any parameter validation in** ``__init__``.
All logic behind estimator parameters,
like translating string arguments into functions, should be done in ``fit``.

Also it is expected that parameters with trailing ``_`` are **not to be set
inside the** ``__init__`` **method**. All and only the public attributes set by
fit have a trailing ``_``. As a result the existence of parameters with
trailing ``_`` is used to check if the estimator has been fitted.

.. _cloning:

Cloning
-------
For use with the :mod:`model_selection` module,
an estimator must support the ``base.clone`` function to replicate an estimator.
This can be done by providing a ``get_params`` method.
If ``get_params`` is present, then ``clone(estimator)`` will be an instance of
``type(estimator)`` on which ``set_params`` has been called with clones of
the result of ``estimator.get_params()``.

Objects that do not provide this method will be deep-copied
(using the Python standard function ``copy.deepcopy``)
if ``safe=False`` is passed to ``clone``.

Pipeline compatibility
----------------------
For an estimator to be usable together with ``pipeline.Pipeline`` in any but the
last step, it needs to provide a ``fit`` or ``fit_transform`` function.
To be able to evaluate the pipeline on any data but the training set,
it also needs to provide a ``transform`` function.
There are no special requirements for the last step in a pipeline, except that
it has a ``fit`` function. All ``fit`` and ``fit_transform`` functions must
take arguments ``X, y``, even if y is not used. Similarly, for ``score`` to be
usable, the last step of the pipeline needs to have a ``score`` function that
accepts an optional ``y``.

Estimator types
---------------
Some common functionality depends on the kind of estimator passed.
For example, cross-validation in :class:`model_selection.GridSearchCV` and
:func:`model_selection.cross_val_score` defaults to being stratified when used
on a classifier, but not otherwise. Similarly, scorers for average precision
that take a continuous prediction need to call ``decision_function`` for classifiers,
but ``predict`` for regressors. This distinction between classifiers and regressors
is implemented using the ``_estimator_type`` attribute, which takes a string value.
It should be ``"classifier"`` for classifiers and ``"regressor"`` for
regressors and ``"clusterer"`` for clustering methods, to work as expected.
Inheriting from ``ClassifierMixin``, ``RegressorMixin`` or ``ClusterMixin``
will set the attribute automatically.  When a meta-estimator needs to distinguish
among estimator types, instead of checking ``_estimator_type`` directly, helpers
like :func:`base.is_classifier` should be used.

Working notes
-------------

For unresolved issues, TODOs, and remarks on ongoing work, developers are
advised to maintain notes on the `GitHub wiki
<https://github.com/scikit-learn/scikit-learn/wiki>`__.

Specific models
---------------

Classifiers should accept ``y`` (target) arguments to ``fit`` that are
sequences (lists, arrays) of either strings or integers.  They should not
assume that the class labels are a contiguous range of integers; instead, they
should store a list of classes in a ``classes_`` attribute or property.  The
order of class labels in this attribute should match the order in which
``predict_proba``, ``predict_log_proba`` and ``decision_function`` return their
values.  The easiest way to achieve this is to put::

    self.classes_, y = np.unique(y, return_inverse=True)

in ``fit``.  This returns a new ``y`` that contains class indexes, rather than
labels, in the range [0, ``n_classes``).

A classifier's ``predict`` method should return
arrays containing class labels from ``classes_``.
In a classifier that implements ``decision_function``,
this can be achieved with::

    def predict(self, X):
        D = self.decision_function(X)
        return self.classes_[np.argmax(D, axis=1)]

In linear models, coefficients are stored in an array called ``coef_``, and the
independent term is stored in ``intercept_``.  ``sklearn.linear_model.base``
contains a few base classes and mixins that implement common linear model
patterns.

The :mod:`sklearn.utils.multiclass` module contains useful functions
for working with multiclass and multilabel problems.

Estimator Tags
--------------
.. warning::

    The estimator tags are experimental and the API is subject to change.

Scikit-learn introduced estimator tags in version 0.21.  These are annotations
of estimators that allow programmatic inspection of their capabilities, such as
sparse matrix support, supported output types and supported methods.  The
estimator tags are a dictionary returned by the method ``_get_tags()``.  These
tags are used by the common tests and the :func:`sklearn.utils.estomator_checks.check_estimator` function to
decide what tests to run and what input data is appropriate. Tags can depends on
estimator parameters or even system architecture and can in general only be
determined at runtime.

The default value of all tags except for ``X_types`` is ``False``.

The current set of estimator tags are:

non_deterministic
    whether the estimator is not deterministic given a fixed ``random_state``

requires_positive_data - unused for now
    whether the estimator requires positive X.

no_validation
    whether the estimator skips input-validation. This is only meant for stateless and dummy transformers!

multioutput - unused for now
    whether a regressor supports multi-target outputs or a classifier supports multi-class multi-output.

multilabel
    whether the estimator supports multilabel output

stateless
    whether the estimator needs access to data for fitting. Even though
    an estimator is stateless, it might still need a call to ``fit`` for initialization.

allow_nan
    whether the estimator supports data with missing values encoded as np.NaN

poor_score
    whether the estimator fails to provide a "reasonable" test-set score, which
    currently for regression is an R2 of 0.5 on a subset of the boston housing
    dataset, and for classification an accuracy of 0.83 on
    ``make_blobs(n_samples=300, random_state=0)``. These datasets and values
    are based on current estimators in sklearn and might be replaced by
    something more systematic.

multioutput_only
    whether estimator supports only multi-output classification or regression.

_skip_test
    whether to skip common tests entirely. Don't use this unless you have a *very good* reason.

X_types
    Supported input types for X as list of strings. Tests are currently only run if '2darray' is contained
    in the list, signifying that the estimator takes continuous 2d numpy arrays as input. The default
    value is ['2darray']. Other possible types are ``'string'``, ``'sparse'``,
    ``'categorical'``, ``dict``, ``'1dlabels'`` and ``'2dlabels'``.
    The goals is that in the future the supported input type will determine the
    data used during testsing, in particular for ``'string'``, ``'sparse'`` and
    ``'categorical'`` data.  For now, the test for sparse data do not make use
    of the ``'sparse'`` tag.


In addition to the tags, estimators are also need to declare any non-optional
parameters to ``__init__`` in the ``_required_parameters`` class attribute,
which is a list or tuple.  If ``_required_parameters`` is only
``["estimator"]`` or ``["base_estimator"]``, then the estimator will be
instantiated with an instance of ``LinearDiscriminantAnalysis`` (or
``RidgeRegression`` if the estimator is a regressor) in the tests. The choice
of these two models is somewhat idiosyncratic but both should provide robust
closed-form solutions.

.. _reading-code:

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

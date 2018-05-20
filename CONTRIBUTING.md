
Contributing to scikit-learn
============================

**Note: This document is a 'getting started' summary for contributing code,
documentation, testing, and filing issues.** Visit the [**Contributing
page**](http://scikit-learn.org/dev/developers/contributing.html)
for the full contributor's guide. Please read it carefully to help make
the code review process go as smoothly as possible and maximize the
likelihood of your contribution being merged.

How to contribute
-----------------

The preferred workflow for contributing to scikit-learn is to fork the
[main repository](https://github.com/scikit-learn/scikit-learn) on
GitHub, clone, and develop on a branch. Steps:

1. Fork the [project repository](https://github.com/scikit-learn/scikit-learn)
   by clicking on the 'Fork' button near the top right of the page. This creates
   a copy of the code under your GitHub user account. For more details on
   how to fork a repository see [this guide](https://help.github.com/articles/fork-a-repo/).

2. Clone your fork of the scikit-learn repo from your GitHub account to your local disk:

   ```bash
   $ git clone git@github.com:YourLogin/scikit-learn.git
   $ cd scikit-learn
   ```

3. Create a ``feature`` branch to hold your development changes:

   ```bash
   $ git checkout -b my-feature
   ```

   Always use a ``feature`` branch. It's good practice to never work on the ``master`` branch!

4. Develop the feature on your feature branch. Add changed files using ``git add`` and then ``git commit`` files:

   ```bash
   $ git add modified_files
   $ git commit
   ```

   to record your changes in Git, then push the changes to your GitHub account with:

   ```bash
   $ git push -u origin my-feature
   ```

5. Follow [these instructions](https://help.github.com/articles/creating-a-pull-request-from-a-fork)
to create a pull request from your fork. This will send an email to the committers.

(If any of the above seems like magic to you, please look up the
[Git documentation](https://git-scm.com/documentation) on the web, or ask a friend or another contributor for help.)

Pull Request Checklist
----------------------

We recommended that your contribution complies with the
following rules before you submit a pull request:

-  Follow the
   [coding-guidelines](http://scikit-learn.org/dev/developers/contributing.html#coding-guidelines).

-  Use, when applicable, the validation tools and scripts in the
   `sklearn.utils` submodule.  A list of utility routines available
   for developers can be found in the
   [Utilities for Developers](http://scikit-learn.org/dev/developers/utilities.html#developers-utils)
   page.

-  Give your pull request a helpful title that summarises what your
   contribution does. In some cases `Fix <ISSUE TITLE>` is enough.
   `Fix #<ISSUE NUMBER>` is not enough.

-  Often pull requests resolve one or more other issues (or pull requests).
   If merging your pull request means that some other issues/PRs should
   be closed, you should
   [use keywords to create link to them](https://github.com/blog/1506-closing-issues-via-pull-requests/)
   (e.g., `Fixes #1234`; multiple issues/PRs are allowed as long as each one
   is preceded by a keyword). Upon merging, those issues/PRs will
   automatically be closed by GitHub. If your pull request is simply related
   to some other issues/PRs, create a link to them without using the keywords
   (e.g., `See also #1234`).

-  All public methods should have informative docstrings with sample
   usage presented as doctests when appropriate.

-  Please prefix the title of your pull request with `[MRG]` (Ready for
   Merge), if the contribution is complete and ready for a detailed review.
   Two core developers will review your code and change the prefix of the pull
   request to `[MRG + 1]` and `[MRG + 2]` on approval, making it eligible
   for merging. An incomplete contribution -- where you expect to do more work before
   receiving a full review -- should be prefixed `[WIP]` (to indicate a work
   in progress) and changed to `[MRG]` when it matures. WIPs may be useful
   to: indicate you are working on something to avoid duplicated work,
   request broad review of functionality or API, or seek collaborators.
   WIPs often benefit from the inclusion of a
   [task list](https://github.com/blog/1375-task-lists-in-gfm-issues-pulls-comments)
   in the PR description.

-  All other tests pass when everything is rebuilt from scratch. On
   Unix-like systems, check with (from the toplevel source folder):

      ```bash
      $ make
      ```

-  When adding additional functionality, provide at least one
   example script in the ``examples/`` folder. Have a look at other
   examples for reference. Examples should demonstrate why the new
   functionality is useful in practice and, if possible, compare it
   to other methods available in scikit-learn.

-  Documentation and high-coverage tests are necessary for enhancements to be
   accepted. Bug-fixes or new features should be provided with 
   [non-regression tests](https://en.wikipedia.org/wiki/Non-regression_testing).
   These tests verify the correct behavior of the fix or feature. In this
   manner, further modifications on the code base are granted to be consistent
   with the desired behavior.
   For the Bug-fixes case, at the time of the PR, this tests should fail for
   the code base in master and pass for the PR code.


-  At least one paragraph of narrative documentation with links to
   references in the literature (with PDF links when possible) and
   the example.

-  The documentation should also include expected time and space
   complexity of the algorithm and scalability, e.g. "this algorithm
   can scale to a large number of samples > 100000, but does not
   scale in dimensionality: n_features is expected to be lower than
   100".

You can also check for common programming errors with the following
tools:

-  Code with good unittest **coverage** (at least 80%), check with:

  ```bash
  $ pip install pytest pytest-cov
  $ pytest --cov sklearn path/to/tests_for_package
  ```

-  No flake8 warnings, check with:

  ```bash
  $ pip install flake8
  $ flake8 path/to/module.py
  ```

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the GitHub issue).

Filing bugs
-----------
We use GitHub issues to track all bugs and feature requests; feel free to
open an issue if you have found a bug or wish to see a feature implemented.

It is recommended to check that your issue complies with the
following rules before submitting:

-  Verify that your issue is not being currently addressed by other
   [issues](https://github.com/scikit-learn/scikit-learn/issues?q=)
   or [pull requests](https://github.com/scikit-learn/scikit-learn/pulls?q=).

-  If you are submitting an algorithm or feature request, please verify that
   the algorithm fulfills our
   [new algorithm requirements](http://scikit-learn.org/dev/faq.html#what-are-the-inclusion-criteria-for-new-algorithms).

-  Please ensure all code snippets and error messages are formatted in
   appropriate code blocks.
   See [Creating and highlighting code blocks](https://help.github.com/articles/creating-and-highlighting-code-blocks).

-  Please include your operating system type and version number, as well
   as your Python, scikit-learn, numpy, and scipy versions. This information
   can be found by running the following code snippet:

  ```python
  import platform; print(platform.platform())
  import sys; print("Python", sys.version)
  import numpy; print("NumPy", numpy.__version__)
  import scipy; print("SciPy", scipy.__version__)
  import sklearn; print("Scikit-Learn", sklearn.__version__)
  ```

-  Please be specific about what estimators and/or functions are involved
   and the shape of the data, as appropriate; please include a
   [reproducible](http://stackoverflow.com/help/mcve) code snippet
   or link to a [gist](https://gist.github.com). If an exception is raised,
   please provide the traceback.

New contributor tips
--------------------

A great way to start contributing to scikit-learn is to pick an item from the
list of
[good first issues](https://github.com/scikit-learn/scikit-learn/labels/good%20first%20issue). If
you have already contributed to scikit-learn look at
[Easy issues](https://github.com/scikit-learn/scikit-learn/labels/Easy)
instead. Resolving these issues allow you to start contributing to the project
without much prior knowledge. Your assistance in this area will be greatly
appreciated by the more experienced developers as it helps free up their time to
concentrate on other issues.

Documentation
-------------

We are glad to accept any sort of documentation: function docstrings,
reStructuredText documents (like this one), tutorials, etc.
reStructuredText documents live in the source code repository under the
doc/ directory.

You can edit the documentation using any text editor and then generate
the HTML output by typing ``make html`` from the doc/ directory.
Alternatively, ``make`` can be used to quickly generate the
documentation without the example gallery. The resulting HTML files will
be placed in ``_build/html/stable`` and are viewable in a web browser. See the
``README`` file in the ``doc/`` directory for more information.

For building the documentation, you will need
[sphinx](http://sphinx.pocoo.org/),
[matplotlib](http://matplotlib.org/), and
[pillow](http://pillow.readthedocs.io/en/latest/).

When you are writing documentation, it is important to keep a good
compromise between mathematical and algorithmic details, and give
intuition to the reader on what the algorithm does. It is best to always
start with a small paragraph with a hand-waving explanation of what the
method does to the data and a figure (coming from an example)
illustrating it.

Further Information
-------------------

Visit the [Contributing Code](http://scikit-learn.org/stable/developers/contributing.html#coding-guidelines)
section of the website for more information including conforming to the
API spec and profiling contributed code.


Contributing to scikit-learn
============================

**Note: This document is a quick overview of open source contribution.** Visit the official [**Contributing
page**](http://scikit-learn.org/dev/developers/contributing.html)
for more comprehensive instructions. Read carefully to increase your chances for a merged pull request.

How to contribute
-----------------

To contribute to scikit-learn, fork the
[main GitHub repository](https://github.com/scikit-learn/scikit-learn)
, clone, and develop on a new branch. Steps:

1. Fork the [project repository](https://github.com/scikit-learn/scikit-learn)
   by clicking on the 'Fork' button near the top right of the repo page. This creates
   a copy of scikit-learn on your GitHub account. For more details on forking, see [this guide](https://help.github.com/articles/fork-a-repo/).

2. Clone your fork of scikit-learn to your local disk:

   ```bash
   $ git clone git@github.com:YourLogin/scikit-learn.git
   $ cd scikit-learn
   ```

3. Create a ``feature`` branch within your clone:

   ```bash
   $ git checkout -b my-feature
   ```

   Always use a ``feature`` branch. Don't work on the ``master`` branch!

4. Develop the feature on the branch. Add changed files using ``git add``, and then ``git commit`` files:

   ```bash
   $ git add modified_files
   $ git commit
   ```

   to record your changes in Git. Then push the changes up to your GitHub account with:

   ```bash
   $ git push -u origin my-feature
   ```

5. Follow [these instructions](https://help.github.com/articles/creating-a-pull-request-from-a-fork)
to create a pull request from your fork. The project maintainers will receive an email about your pull request.

(If any of the above is confusing, please read the
[Git documentation](https://git-scm.com/documentation), or ask a friend or another contributor for help.)

Pull Request Checklist
----------------------

We recommended that your contribution comply with these
important rules before you submit a pull request:

-  Follow the
   [coding-guidelines](http://scikit-learn.org/dev/developers/contributing.html#coding-guidelines).

-  Use, when applicable, the validation tools and scripts in the
   `sklearn.utils` submodule.  A list of available utility routines can be found [here](http://scikit-learn.org/dev/developers/utilities.html#developers-utils)
   page.

-  Give your pull request a helpful title that summarises your
   contribution. In some cases `Fix <ISSUE TITLE>` is enough.
   `Fix #<ISSUE NUMBER>` is not enough.

-  Often pull requests resolve other issues (or pull requests).
   If merging your pull request means that some other issues/PRs should
   be closed, you should
   [use keywords to link to them](https://github.com/blog/1506-closing-issues-via-pull-requests/)
   (e.g., `Fixes #1234`; multiple issues/PRs are allowed as long as each one
   is preceded by a keyword). Upon merging, GitHub will automatically close those issues/PRs. If your pull request is simply related
   to some other issues/PRs, create a link to them without using the keywords
   (e.g., `See also #1234`).

-  All public methods should have informative docstrings with sample
   usage presented as doctests when appropriate.

-  Prefix the title of your pull request with `[MRG]` (Ready for
   Merge), if the contribution is complete and ready for a detailed review.
   Two core developers will review your code and change the prefix of the pull
   request to `[MRG + 1]` and `[MRG + 2]` on approval, making it eligible
   for merging. An incomplete contribution -- where you expect to do more work before
   receiving a full review -- should be prefixed `[WIP]` (to indicate a work
   in progress) and changed to `[MRG]` when it matures. WIPs may be useful
   to: (1) indicate you are working on something to avoid duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators.
   WIPs often benefit from the inclusion of a
   [task list](https://github.com/blog/1375-task-lists-in-gfm-issues-pulls-comments)
   in the PR description.

-  Confirm that all other tests pass when everything is rebuilt from scratch. On
   Unix-like systems, check with (from the toplevel source folder):

      ```bash
      $ make
      ```

-  When adding additional functionality, provide at least one
   example script in the ``examples/`` folder. Have a look at other
   examples for reference. Examples should demonstrate why the new
   functionality is useful in practice and, if possible, compare it
   to other methods available in scikit-learn.

-  Confirm that tests pass in the code of your pull request but fail in the codebase. Bug-fixes or new features should be provided with 
   [non-regression tests](https://en.wikipedia.org/wiki/Non-regression_testing).
   These tests verify the correct behavior of the fix or feature. This way, further modifications on the code base will be consistent
   with the desired behavior.

-  Add at least one paragraph of narrative documentation with links to
   references in the literature (with PDF links when possible) and
   the example.

-  Include expected time and space complexity of the algorithm and scalability in the documentation. E.g. "this algorithm
   can scale to a large number of samples > 100000, but does not
   scale in dimensionality: n_features is expected to be lower than
   100".

You can also check for common programming errors with tools:

-  Code with good unittest **coverage** (at least 80%), check with:

  ```bash
  $ pip install pytest pytest-cov
  $ pytest --cov sklearn path/to/tests_for_package
  ```

-  Check that there are no flake8 warnings:

  ```bash
  $ pip install flake8
  $ flake8 path/to/module.py
  ```

Bonus points for contributions that include a performance analysis with
a benchmark script and profiling output (please report on the mailing
list or on the GitHub issue).

Filing bugs
-----------
We use GitHub issues to track all bugs and feature requests; feel free to open an issue if you have found a bug or wish to see a feature implemented.

Check that your issue complies with the following rules:

-  Your issue should not being currently addressed by other
   [issues](https://github.com/scikit-learn/scikit-learn/issues?q=)
   or [pull requests](https://github.com/scikit-learn/scikit-learn/pulls?q=).

-  If you are submitting an algorithm or feature request, verify that the algorithm fulfills our
   [new algorithm requirements](http://scikit-learn.org/dev/faq.html#what-are-the-inclusion-criteria-for-new-algorithms).

-  Ensure all code snippets and error messages are formatted in
   appropriate code blocks.
   See [Creating and highlighting code blocks](https://help.github.com/articles/creating-and-highlighting-code-blocks).

-  Include your operating system type and version number, as well
   as your Python, scikit-learn, numpy, and scipy versions. This information
   can be found by running the following code snippet:

  ```python
  import platform; print(platform.platform())
  import sys; print("Python", sys.version)
  import numpy; print("NumPy", numpy.__version__)
  import scipy; print("SciPy", scipy.__version__)
  import sklearn; print("Scikit-Learn", sklearn.__version__)
  ```

-  Be specific about what estimators and/or functions are involved
   and the shape of the data, as appropriate; please include a
   [reproducible](http://stackoverflow.com/help/mcve) code snippet
   or link to a [gist](https://gist.github.com). If an exception is raised,
   please provide the traceback.

New contributor tips
--------------------

A great way to start contributing to scikit-learn is to pick an item from the list of
[good first issues](https://github.com/scikit-learn/scikit-learn/labels/good%20first%20issue). If
you have already contributed to scikit-learn look at
[Easy issues](https://github.com/scikit-learn/scikit-learn/labels/Easy)
instead. Resolving these issues allow you to start contributing to the project
without much prior knowledge. Your assistance in this area will be greatly
appreciated by the more experienced developers, as it helps free up their time to concentrate on other issues.

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

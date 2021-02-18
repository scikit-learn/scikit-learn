=======
Support
=======

There are several ways to get in touch with the developers.


.. _mailing_lists:

Mailing List
============

- The main mailing list is `scikit-learn
  <https://mail.python.org/mailman/listinfo/scikit-learn>`_.

- There is also a commit list `scikit-learn-commits
  <https://lists.sourceforge.net/lists/listinfo/scikit-learn-commits>`_,
  where updates to the main repository and test failures get notified.


.. _user_questions:

User questions
==============

- Some scikit-learn developers support users on StackOverflow using
  the `[scikit-learn] <https://stackoverflow.com/questions/tagged/scikit-learn>`_
  tag.

- For general theoretical or methodological Machine Learning questions
  `stack exchange <https://stats.stackexchange.com/>`_ is probably a more
  suitable venue.

In both cases please use a descriptive question in the title field (e.g.
no "Please help with scikit-learn!" as this is not a question) and put
details on what you tried to achieve, what were the expected results and
what you observed instead in the details field.

Code and data snippets are welcome. Minimalistic (up to ~20 lines long)
reproduction script very helpful.

Please describe the nature of your data and how you preprocessed it:
what is the number of samples, what is the number and type of features
(i.d. categorical or numerical) and for supervised learning tasks,
what target are your trying to predict: binary, multiclass (1 out of
``n_classes``) or multilabel (``k`` out of ``n_classes``) classification
or continuous variable regression.

User questions should **not be asked on the bug tracker**, as it crowds
the list of issues and makes the development of the project harder.

.. _bug_tracker:

Bug tracker
===========

If you think you've encountered a bug, please report it to the issue tracker:

https://github.com/scikit-learn/scikit-learn/issues

Don't forget to include:

  - steps (or better script) to reproduce,

  - expected outcome,

  - observed outcome or Python (or gdb) tracebacks

To help developers fix your bug faster, please link to a https://gist.github.com
holding a standalone minimalistic python script that reproduces your bug and
optionally a minimalistic subsample of your dataset (for instance, exported
as CSV files using ``numpy.savetxt``).

Note: Gists are Git cloneable repositories and thus you can use Git to
push datafiles to them.


.. _gitter:

Gitter
===

Some developers like to hang out on scikit-learn Gitter room:
https://gitter.im/scikit-learn/scikit-learn.


.. _documentation_resources:

Documentation resources
=======================

This documentation is relative to |release|. Documentation for
other versions can be found `here
<http://scikit-learn.org/dev/versions.html>`__.

Printable pdf documentation for old versions can be found `here
<https://sourceforge.net/projects/scikit-learn/files/documentation/>`_.

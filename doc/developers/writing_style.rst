.. _writing_style:

===================
Writing style guide
===================

This page describes how we write the scikit-learn documentation: the user guide,
docstrings, and examples. Most contributors, including most maintainers, are
not native English speakers, and neither are most readers. The goal of these
conventions is text that is clear, consistent, and easy to read for everyone.

You do not need to memorize this page. Many of the rules are checked by
:ref:`Vale <writing_style_vale>`, which points out problems in the lines you
change and links back to the relevant section here. Reviewers may also refer
to this page.

These conventions apply to new and edited text. There is no need to rewrite
existing text that you are not otherwise changing.

.. _writing_style_reader:

Addressing the reader
=====================

- Address the reader as *you*: "You can pass a callable to ``scoring``."
- Use *we* for recommendations made by the project ("We recommend scaling the
  data before fitting an SVM.") and for walking through an example together
  ("We now fit the model on the training set.").
- Avoid the impersonal *one*. Write "you can use" instead of "one can use", and
  "you should" instead of "one should".
- Do not switch between *you*, *we*, and *the user* within a paragraph.

.. _writing_style_contractions:

Contractions
============

Do not use contractions. Write "do not", "is not", "it is", and "cannot"
instead of "don't", "isn't", "it's", and "can't".

In examples, rephrase "let's" rather than writing "let us": "We now look at the
coefficients" instead of "Let's look at the coefficients".

.. _writing_style_abbreviations:

Abbreviations
=============

- Prefer "for example" to "e.g." and "that is" to "i.e.". The full words are
  easier to read, and the two abbreviations are often confused. A direct swap
  does not always read well: often the aside can be dropped, or the sentence
  rephrased.

  .. list-table::
     :header-rows: 1

     * - Instead of
       - Write
     * - (e.g., on large datasets)
       - (for example, on large datasets)
     * - The number of weak learners (i.e. regression trees)
       - The number of weak learners (regression trees)
     * - Default is `None`, i.e. the Ward algorithm is unstructured.
       - By default (`None`), the Ward algorithm is unstructured.

- If you do use "e.g." or "i.e.", write them with both periods.
- Write "also known as" instead of "aka", and "with respect to" instead of
  "w.r.t.".

.. _writing_style_terminology:

Terminology
===========

Use the spelling of a term that the :ref:`glossary` uses. In particular:

.. list-table::
   :header-rows: 1

   * - Write
     - Not
   * - dataset
     - data set
   * - data point
     - datapoint
   * - multiclass, multilabel, multioutput
     - multi-class, multi-label, multi-output
   * - hyperparameter
     - hyper-parameter
   * - cross-validation, cross-validated
     - cross validation, cross validated
   * - preprocessing, precomputed, subsample
     - pre-processing, pre-computed, sub-sample
   * - overfit, underfit
     - over-fit, under-fit
   * - meta-estimator
     - meta estimator, metaestimator
   * - one-vs-rest (OvR)
     - one-vs-all, one versus all (OvA)
   * - non-zero, non-linear, non-negative
     - non zero, non linear, non negative
   * - NumPy, SciPy, GitHub, scikit-learn
     - Numpy, Scipy, Github, Scikit-Learn

Both "mini-batch" and "minibatch" are in use. Use one of them consistently
within a page.

.. _writing_style_spelling:

Spelling
========

Use American spelling: "behavior", "neighbor", "optimize", "normalization",
"modeling", "labeled", and "learned" instead of "behaviour", "neighbour",
"optimise", "normalisation", "modelling", "labelled", and "learnt". Keep the
original spelling in titles of cited papers.

.. _writing_style_idioms:

Idiomatic English
=================

Some constructions are common in the documentation but are not idiomatic
English. They are usually direct translations from other languages.

.. list-table::
   :header-rows: 1

   * - Instead of
     - Write
   * - ``normalize`` allows to report ratios
     - ``normalize`` allows you to report ratios, or ``normalize`` reports ratios
   * - we recommend to use
     - we recommend using
   * - As other classifiers, ...
     - Like other classifiers, ...
   * - depending of
     - depending on
   * - the same than
     - the same as
   * - with regards to
     - with regard to
   * - an amount of samples, less samples
     - a number of samples, fewer samples
   * - keywords arguments, parameters names
     - keyword arguments, parameter names
   * - If None then ...
     - If None, then ...

.. _writing_style_concise:

Concise writing
===============

- Leave out phrases that announce a statement instead of making it: "please
  note that", "it should be noted that", "it is worth noting that". If a point
  is a caveat that readers must not miss, put it in a ``.. note::`` block.
- Prefer the short form: "to" instead of "in order to", "because" instead of
  "due to the fact that", "whether to" instead of "whether or not to", "you
  can" instead of "it is possible to".
- Keep sentences short. A sentence with more than about 40 words is usually
  easier to read as two.
- Do not write "as mentioned above" or "as discussed before". Readers often
  arrive at a section directly, and text inside a dropdown may be hidden. Link
  to the section with ``:ref:`` instead.
- Avoid words that suggest something is easy or evident, such as "obviously",
  "clearly", and "of course". What is obvious to the author is often not
  obvious to the reader.

.. _writing_style_vale:

Checking your writing with Vale
===============================

`Vale <https://vale.sh>`__ is a prose linter. Its configuration is in
``.vale.ini`` at the root of the repository, and the rules are in
``build_tools/vale/styles/sklearn``. Vale checks ``.rst`` files, docstrings,
and the text blocks of the examples.

If you have :ref:`set up pre-commit <pre_commit>`, Vale runs on every commit
and reports problems in the lines you changed. Each alert has a level:

- **error**: must be fixed. The commit is blocked.
- **warning**: should be fixed.
- **suggestion**: worth considering. Ignore it if the change would not improve
  the text.

To check the lines your branch changed compared to ``main`` in a given file,
run:

.. code-block:: console

  $ pre-commit run vale --files doc/modules/svm.rst

To see every alert in a file, including in lines you did not change, run the
script behind the hook directly. It needs ``pip install vale docutils``:

.. code-block:: console

  $ python build_tools/vale/run_vale.py --all-lines doc/modules/svm.rst

Vale only checks patterns. It cannot tell whether a sentence is clear or a
paragraph well structured, and it sometimes flags text that is fine. When an
alert is wrong, ignore it and mention it in your pull request so that we can
improve the rule.

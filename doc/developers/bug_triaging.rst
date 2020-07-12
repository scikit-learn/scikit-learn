.. _bug_triaging:

Bug triaging and issue curation
================================

The `issue tracker <https://github.com/scikit-learn/scikit-learn/issues>`_
is important to the communication in the project: it helps
developers identify major projects to work on, as well as to discuss
priorities. For this reason, it is important to curate it, adding labels
to issues and closing issues that are not necessary.

Working on issues to improve them
--------------------------------------

Improving issues increases their chances of being successfully resolved.
Guidelines on submitting good issues can be found :ref:`here
<filing_bugs>`. A third party can give useful feedback or even add
comments on the issue, while core-developpers or members of the triage
team can edit the issue description and title.

The following actions are typically useful:

  - documenting issues that are missing elements to reproduce the problem
    such as code samples

  - correcting incorrect use of code formatting

  - making sure that the title and description are explicit about the
    problem to be solved

  - linking to related issues or discussions while briefly describing how
    they are related, for instance "See also #xyz for a similar attempt
    at this" or "See also #xyz where the same thing happened in
    SomeEstimator" provides context and helps the discussion.

.. topic:: Fruitful discussions

   Online discussions may be harder than it seems at first glance, in
   particular given that a person new to open-source may have a very
   different understanding of the process than a seasonned maintainer.

   Overall, it is useful to stay positive and assume good will. `The
   following article
   <http://gael-varoquaux.info/programming/technical-discussions-are-hard-a-few-tips.html>`_
   explores how to lead online discussions in the context of open source.

Triaging operations for members of the core and triage teams
-------------------------------------------------------------

In addition to the above, members of the core team and the triage team
can do the following important tasks:

- Update labels for issues and PRs

- Follow up on stalled PRs, to see if they must be relabeled as
  stalled and needing help (this is typically very important in the context
  of sprints, where the risk is to create many unfinished PRs)

- Triage issues:

  - **close usage questions** and politely point the reporter to use
    Stack Overflow instead.

  - **close duplicate issues**, but only after checking that they are
    indeed duplicate. Ideally, the original submitter moves the
    discussion to the older, duplicate issue

  - **close issues that cannot be replicated**, after leaving time (at
    least a week) to add extra information

:ref:`Saved replies <saved_replies>` are useful to gain time and yet be
welcoming and polite when triaging.


.. topic:: Closing issues: a tough call

    When uncertain on whether an issue should be closed or not, it is
    best to strive for consensus with the original poster, and possibly
    to seek relevant expertise. However, when the issue is a usage
    question, or when it has been considered as unclear for many years it
    should be closed.

A typical workflow for triaging issues
----------------------------------------

The following workflow [*]_ is a good way to approach issue triaging:

1. Thank the reporter for opening an issue

   The issue tracker is many people’s first interaction with the
   scikit-learn project itself, beyond just using the library. As such,
   we want it to be a welcoming, pleasant experience.

2. Is this a usage question? If so close it with a polite message
   (:ref:`here is an example <saved_replies>`).

3. Is the necessary information provided?

   If crucial information (like the version of scikit-learn used), is
   missing feel free to ask for that and label the issue with "Needs
   info".

4. Is this a duplicate issue?

   We have many open issues. If a new issue seems to be a duplicate,
   point to the original issue. If it is a clear duplicate, or consensus
   is that it is redundant, close it. Make sure to still thank the
   reporter, and encourage them to chime in on the original issue, and
   perhaps try to fix it.

   If the new issue provides relevant information, such as a better or
   slightly different example, add it to the original issue as a comment
   or an edit to the original post.


5. Make sure that the title accurately reflects the issue. Edit it
   yourself if it's not clear.

6. Is the issue minimal and reproducible?

   For bug reports, we ask that the reporter provide a minimal
   reproducible example. See
   https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports
   for a good explanation. If the example is not reproducible, or if
   it's clearly not minimal, feel free to ask the reporter if they can
   provide and example or simplify the provided one. Do acknowledge that
   writing minimal reproducible examples is hard work. If the reporter
   is struggling, you can try to write one yourself.

   If a reproducible example is provided, but you see a simplification,
   edit the original post with your simpler reproducible example.

4. If a reproducible example can't be provided, add the “Bug: triage”
   label.

5. Add the relevant labels, such as "Documentation" when the issue is
   about documentation, "Bug" if it is clearly a bug, "Enhancement" if it
   is an enhancement request, ...

   If the issue is clearly defined and the fix seems relatively
   straightforward, label the issue as “Good first issue”.

   An additional useful step can be to tag the corresponding module e.g.
   `sklearn.linear_models` when relevant.

.. [*] Adapted from the pandas project https://dev.pandas.io/docs/development/maintaining.html

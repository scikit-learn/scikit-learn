.. _minimal-reproducer:

===========================================
Craft a minimal reproducer for scikit-learn
===========================================


Whether if you are submitting a bug report, designing a pytest for a feature you
are contributing to, or simply posting a question in a forum, being able to
craft minimal, reproducible examples (or minimal, workable examples) is the key
to communicating with the community.

There are very good guidelines on the internet such as `this StackOverflow
document <https://stackoverflow.com/help/mcve>`_ or `this blogpost by Matthew
Rocklin <https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_
on how to craft Minimal Complete Verifiable Examples. Our goal is not to be
repetitive with those references but rather to provide a step-by-step guide on
how to narrow down a bug until you have reached the shortest possible code to
reproduce it.

Understand the output of your code
----------------------------------
Try to answer this questions to ensure you are facing a bug:
    - Where, what, why

Identify the type of problem you are solving
--------------------------------------------
Is it classification, regression, clustering, etc.

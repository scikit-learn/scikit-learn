.. _faq:

===========================
Frequently Asked Questions
===========================

Here we try to give some answers to questions that regularly pop up on the mailing list.

What is the project name (a lot of people get it wrong)?
--------------------------------------------------------
scikit-learn, or sklearn for short, but not scikit or SciKit nor sci-kit learn. Also not scikits.learn or scikits-learn, which where previously used.

How do you pronounce the project name?
------------------------------------------
sy-kit learn. sci stands for science! 

Why scikit?
------------
There are multiple scikits, which are scientific toolboxes build around scipy.
You can find a list here:
Apart from scikit-learn, other popular ones are scikits.statsmodel and scikit-image.

How can I contribute to scikit-learn?
-----------------------------------------
See :ref:`contributing`.

Can I add this new algorithm that I (or someone else) just published?
-------------------------------------------------------------------------
No. As a rule we only add well-established algorithms. A rule of thumb is at least
3 years since publications, 1000+ cites and wide use and usefullness.

Can I add this classical algorithm from the 80s?
---------------------------------------------------
Depends. If there is a common usecase within the scope of scikit-learn, such
as classification, regression or clustering, where it outperforms methods
that are already implemented in scikit-learn, we will consider it.

Why did you remove HMMs from scikit-learn?
--------------------------------------------
See :ref:`adding_graphical_models`.

.. _adding_graphical_models:

Will you add graphical models or sequence prediction to scikit-learn?
------------------------------------------------------------------------
Not in the forseeable future. 
scikit-learn tries to provide a unified API for the basic tasks in machine
learning, with pipelines and meta-algorithms like grid search to tie
everything together. The required concepts, APIs, algorithms and
expertise required for stuctured learning are different from what
scikit-learn has to offer. If we started doing arbitrary structured
learning, we'd need to redesign the whole package and the project
would likely collapse under its own weight.

There are two project with API similar to scikit-learn that
do structured prediction:

* `pystruct <http://pystruct.github.io/>`_ handles general structured
  learning (focuses on SSVMs on arbitrary graph structures with
  approximate inference; defines the notion of sample as an instance of
  the graph structure)

* `seqlearn <http://larsmans.github.io/seqlearn/>`_ handles sequences only (focuses on
  exact inference; has HMMs, but mostly for the sake of completeness;
  treats a feature vector as a sample and uses an offset encoding for
  the dependencies between feature vectors)

===========
Naive Bayes
===========

.. currentmodule:: scikits.learn.naive_bayes


**Naive Bayes** algorithms are a set of supervised learning methods
based on applying Bayes' theorem with the "naive" assumption of independence
between every pair of features. Given a class variable :math:`c` and a
dependent set of feature variables :math:`f_1` through :math:`f_n`, Bayes'
theorem states the following relationship:

.. math::

   p(c \mid f_1,\dots,f_n) \propto p(c) p(\mid f_1,\dots,f_n \mid c)

Using the naive assumption this relationship is simplified:

.. math::

   p(c \mid f_1,\dots,f_n) \propto p(c) \prod_{i=1}^{n} p(f_i \mid c)

   \Downarrow

   \hat{c} = \arg\max_c p(c) \prod_{i=1}^{n} p(f_i \mid c),

where we used the Maximum a Posteriori estimator.

The different naive Bayes classifiers differ by the assumption on the
distribution of :math:`p(f_i \mid c)`:

In spite of their naive design and apparently over-simplified assumptions,
naive Bayes classifiers have worked quite well in many real-world situations,
famously document classification and spam filtering. They requires a small
amount of training data to estimate the necessary parameters.

The decoupling of the class conditional feature distributions means that each
distribution can be independently estimated as a one dimensional distribution.
This in turn helps to alleviate problems stemming from the curse of
dimensionality.


Gaussian Naive Bayes
--------------------

:class:`GaussianNB` implements the Gaussian Naive Bayes algorithm for
classification. The likelihood of the features is assumed to be gaussian:

.. math::

   p(f_i \mid c) &= \frac{1}{\sqrt{2\pi\sigma^2_c}} \exp^{-\frac{ (f_i - \mu_c)^2}{2\pi\sigma^2_c}}

The parameters of the distribution, :math:`\sigma_c` and :math:`\mu_c` are
estimated using maximum likelihood.

.. topic:: Examples:

 * :ref:`example_naive_bayes.py`,

Multinomial Naive Bayes
-----------------------

:class:`MultinomialNB` implements the Multinomial Naive Bayes algorithm for classification.
Multinomial Naive Bayes models the distribution of words in a document as a
multinomial. The distribution is parametrized by the vector
:math:`\overline{\theta_c} = (\theta_{c1},\ldots,\theta_{cn})` where :math:`c`
is the class of document, :math:`n` is the size of the vocabulary and :math:`\theta_{ci}`
is the probability of word :math:`i` appearing in a document of class :math:`c`.
The likelihood of document :math:`d` is,

.. math::

   p(d \mid \overline{\theta_c}) &= \frac{ (\sum_i f_i)! }{\prod_i f_i !} \prod_i(\theta_{ci})^{f_i}

where :math:`f_{i}` is the frequency count of word :math:`i`. It can be shown
that the maximum posterior probability is,

.. math::

   \hat{c} = \arg\max_c [ \log p(\overline{\theta_c}) + \sum_i f_i \log \theta_{ci} ]

The vector of parameters :math:`\overline{\theta_c}` is estimated by a smoothed
version of maximum likelihood,

.. math::

    \hat{\theta}_{ci} = \frac{ N_{ci} + \alpha_i }{N_c + \alpha }

where :math:`N_{ci}` is the number of times word :math:`i` appears in a document
of class :math:`c` and :math:`N_{c}` is the total count of words in a document
of class :math:`c`. The smoothness priors :math:`\alpha_i` and their sum 
:math:`\alpha` account for words not seen in the learning samples.

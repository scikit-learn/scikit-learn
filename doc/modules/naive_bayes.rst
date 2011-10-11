===========
Naive Bayes
===========

.. currentmodule:: sklearn.naive_bayes


Naive Bayes methods are a set of supervised learning algorithms
based on applying Bayes' theorem with the "naive" assumption of independence
between every pair of features. Given a class variable :math:`y` and a
dependent feature vector :math:`x_1` through :math:`x_n`,
Bayes' theorem states the following relationship:

.. math::

   P(y \mid x_1, \dots, x_n) = \frac{P(y) P(x_1, \dots x_n \mid y)}
                                    {P(x_1, \dots, x_n)}

Using the naive independence assumption that

.. math::

   P(x_i | y, x_1, \dots, x_{i-1}, x_{i+1}, \dots, x_n) = P(x_i | y),

for all :math:`i`, this relationship is simplified to

.. math::

   P(y \mid x_1, \dots, x_n) = \frac{P(y) \prod_{i=1}^{n} P(x_i \mid y)}
                                    {P(x_1, \dots, x_n)}

Since :math:`P(x_1, \dots, x_n)` is constant given the input,
we can use the following classification rule:

.. math::

   P(y \mid x_1, \dots, x_n) \propto P(y) \prod_{i=1}^{n} P(x_i \mid y)

   \Downarrow

   \hat{y} = \arg\max_y P(y) \prod_{i=1}^{n} P(f_i \mid y),

and we can use Maximum A Posteriori (MAP) estimation to estimate
:math:`P(y)` and :math:`P(x_i \mid y)`;
the former is then the relative frequency of class :math:`y`
in the training set.

The different Naive Bayes classifiers differ mainly by the assumptions they make
regarding the distribution of :math:`P(x_i \mid y)`.

In spite of their apparently over-simplified assumptions, Naive Bayes
classifiers have worked quite well in many real-world situations, famously
document classification and spam filtering. They requires a small amount
of training data to estimate the necessary parameters. (For theoretical
reasons why Naive Bayes works well, and on which types of data it does, see
the references below.)

Naive Bayes learners and classifiers can be extremely fast compared to more
sophisticated methods.
The decoupling of the class conditional feature distributions means that each
distribution can be independently estimated as a one dimensional distribution.
This in turn helps to alleviate problems stemming from the curse of
dimensionality.

.. topic:: References:

 * H. Zhang (2004). `The optimality of Naive Bayes.
   <http://www.cs.unb.ca/profs/hzhang/publications/FLAIRS04ZhangH.pdf>`_
   Proc. FLAIRS.


Gaussian Naive Bayes
--------------------

:class:`GaussianNB` implements the Gaussian Naive Bayes algorithm for
classification. The likelihood of the features is assumed to be Gaussian:

.. math::

   P(x_i \mid y) &= \frac{1}{\sqrt{2\pi\sigma^2_y}} \exp\left(-\frac{ (x_i - \mu_y)^2}{2\pi\sigma^2_y}\right)

The parameters :math:`\sigma_y` and :math:`\mu_y`
are estimated using maximum likelihood.

.. topic:: Examples:

 * :ref:`example_naive_bayes.py`


.. _multinomial_naive_bayes:

Multinomial Naive Bayes
-----------------------

:class:`MultinomialNB` implements the Naive Bayes algorithm for multinomially
distributed data, and is one of the two classic Naive Bayes variants used in
text classification (where the data are typically represented as word vector
counts, although tf-idf vectors are also known to work well in practice).
The distribution is parametrized by vectors
:math:`\theta_y = (\theta_{y1},\ldots,\theta_{yn})`
for each class :math:`y`, where :math:`n` is the number of features
(in text classification, the size of the vocabulary)
and :math:`\theta_{yi}` is the probability :math:`P(y \mid x_i)`
of feature :math:`i` appearing in a sample belonging to class :math:`y`.

The parameters :math:`\theta_y` is estimated by a smoothed
version of maximum likelihood, i.e. relative frequency counting:

.. math::

    \hat{\theta}_{yi} = \frac{ N_{yi} + \alpha}{N_y + \alpha n}

where :math:`N_{yi} = \sum_{x \in T} x_i` is
the number of times feature :math:`i` appears in a sample of class :math:`y`
in the training set :math:`T`,
and :math:`N_{y} = \sum_{i=1}^{|T|} N_{yi}` is the total count of
all features for class :math:`y`.

The smoothing priors :math:`\alpha \ge 0` accounts for
features not present in the learning samples and prevents zero probabilities
in further computations.
Setting :math:`\alpha = 1` is called Laplace smoothing,
while :math:`\alpha < 1` is called Lidstone smoothing.


.. _bernoulli_naive_bayes:

Bernoulli Naive Bayes
---------------------

:class:`BernoulliNB` implements the Naive Bayes training and classification
algorithms for data that is distributed according to multivariate Bernoulli
distributions. It requires samples to be represented as binary-valued/boolean
feature vectors; if handed any other kind of data, it binarizes it (depending
on the ``binarize`` parameter).

In the case of text classification, word occurrence vectors (rather than word
count vectors) may be used to train and use this classifier. ``BernoulliNB``
might perform better on some datasets, especially those with shorter documents,
because it explicitly penalizes the non-occurrence of words/features in a
dataset where ``MultinomialNB`` would only notice a zero count, but for text
classification ``MultinomialNB`` will generally be better. It is advisable to
evaluate both models, if time permits.

.. topic:: References:

 * C.D. Manning, P. Raghavan and H. Schütze (2008). Introduction to
   Information Retrieval. Cambridge University Press, pp. 234-265.

 * A. McCallum and K. Nigam (1998).
   `A comparison of event models for Naive Bayes text classification.
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.46.1529>`_
   Proc. AAAI/ICML-98 Workshop on Learning for Text Categorization, pp. 41-48.

 * V. Metsis, I. Androutsopoulos and G. Paliouras (2006).
   `Spam filtering with Naive Bayes -- Which Naive Bayes?
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.61.5542>`_
   3rd Conf. on Email and Anti-Spam (CEAS).


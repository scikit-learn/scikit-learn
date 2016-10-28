.. _naive_bayes:

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

   \hat{y} = \arg\max_y P(y) \prod_{i=1}^{n} P(x_i \mid y),

and we can use Maximum A Posteriori (MAP) estimation to estimate
:math:`P(y)` and :math:`P(x_i \mid y)`;
the former is then the relative frequency of class :math:`y`
in the training set.

The different naive Bayes classifiers differ mainly by the assumptions they
make regarding the distribution of :math:`P(x_i \mid y)`.

In spite of their apparently over-simplified assumptions, naive Bayes
classifiers have worked quite well in many real-world situations, famously
document classification and spam filtering. They require a small amount
of training data to estimate the necessary parameters. (For theoretical
reasons why naive Bayes works well, and on which types of data it does, see
the references below.)

Naive Bayes learners and classifiers can be extremely fast compared to more
sophisticated methods.
The decoupling of the class conditional feature distributions means that each
distribution can be independently estimated as a one dimensional distribution.
This in turn helps to alleviate problems stemming from the curse of
dimensionality.

On the flip side, although naive Bayes is known as a decent classifier,
it is known to be a bad estimator, so the probability outputs from
``predict_proba`` are not to be taken too seriously.

.. topic:: References:

 * H. Zhang (2004). `The optimality of Naive Bayes.
   <http://www.cs.unb.ca/~hzhang/publications/FLAIRS04ZhangH.pdf>`_
   Proc. FLAIRS.

.. _gaussian_naive_bayes:

Gaussian Naive Bayes
--------------------

:class:`GaussianNB` implements the Gaussian Naive Bayes algorithm for
classification. The likelihood of the features is assumed to be Gaussian:

.. math::

   P(x_i \mid y) &= \frac{1}{\sqrt{2\pi\sigma^2_y}} \exp\left(-\frac{(x_i - \mu_y)^2}{2\sigma^2_y}\right)

The parameters :math:`\sigma_y` and :math:`\mu_y`
are estimated using maximum likelihood.

    >>> from sklearn import datasets
    >>> iris = datasets.load_iris()
    >>> from sklearn.naive_bayes import GaussianNB
    >>> gnb = GaussianNB()
    >>> y_pred = gnb.fit(iris.data, iris.target).predict(iris.data)
    >>> print("Number of mislabeled points out of a total %d points : %d"
    ...       % (iris.data.shape[0],(iris.target != y_pred).sum()))
    Number of mislabeled points out of a total 150 points : 6

.. _multinomial_naive_bayes:

Multinomial Naive Bayes
-----------------------

:class:`MultinomialNB` implements the naive Bayes algorithm for multinomially
distributed data, and is one of the two classic naive Bayes variants used in
text classification (where the data are typically represented as word vector
counts, although tf-idf vectors are also known to work well in practice).
The distribution is parametrized by vectors
:math:`\theta_y = (\theta_{y1},\ldots,\theta_{yn})`
for each class :math:`y`, where :math:`n` is the number of features
(in text classification, the size of the vocabulary)
and :math:`\theta_{yi}` is the probability :math:`P(x_i \mid y)`
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

:class:`BernoulliNB` implements the naive Bayes training and classification
algorithms for data that is distributed according to multivariate Bernoulli
distributions; i.e., there may be multiple features but each one is assumed
to be a binary-valued (Bernoulli, boolean) variable.
Therefore, this class requires samples to be represented as binary-valued
feature vectors; if handed any other kind of data, a ``BernoulliNB`` instance
may binarize its input (depending on the ``binarize`` parameter).

The decision rule for Bernoulli naive Bayes is based on

.. math::

    P(x_i \mid y) = P(i \mid y) x_i + (1 - P(i \mid y)) (1 - x_i)

which differs from multinomial NB's rule
in that it explicitly penalizes the non-occurrence of a feature :math:`i`
that is an indicator for class :math:`y`,
where the multinomial variant would simply ignore a non-occurring feature.

In the case of text classification, word occurrence vectors (rather than word
count vectors) may be used to train and use this classifier. ``BernoulliNB``
might perform better on some datasets, especially those with shorter documents.
It is advisable to evaluate both models, if time permits.

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


Out-of-core naive Bayes model fitting
-------------------------------------

Naive Bayes models can be used to tackle large scale classification problems
for which the full training set might not fit in memory. To handle this case,
:class:`MultinomialNB`, :class:`BernoulliNB`, and :class:`GaussianNB`
expose a ``partial_fit`` method that can be used
incrementally as done with other classifiers as demonstrated in
:ref:`sphx_glr_auto_examples_applications_plot_out_of_core_classification.py`. All naive Bayes
classifiers support sample weighting.

Contrary to the ``fit`` method, the first call to ``partial_fit`` needs to be
passed the list of all the expected class labels.

For an overview of available strategies in scikit-learn, see also the
:ref:`out-of-core learning <scaling_strategies>` documentation.

.. note::

   The ``partial_fit`` method call of naive Bayes models introduces some
   computational overhead. It is recommended to use data chunk sizes that are as
   large as possible, that is as the available RAM allows.

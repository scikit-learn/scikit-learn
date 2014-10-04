.. _naive_bayes:

===========
朴素贝叶斯
===========

.. currentmodule:: sklearn.naive_bayes

朴素贝叶斯方法是一个基于贝叶斯定力的监督学习算法，其名字是由于采用了一个简单地假设：特征之间不存在相关性。对于一个类别变量 :math:`y` 和一系列特征 :math:`x_1` 至 :math:`x_n`  ，贝叶斯定理可以导出：

.. math::

   P(y \mid x_1, \dots, x_n) = \frac{P(y) P(x_1, \dots x_n \mid y)}
                                    {P(x_1, \dots, x_n)}

利用相互独立的假设：

.. math::

   P(x_i | y, x_1, \dots, x_{i-1}, x_{i+1}, \dots, x_n) = P(x_i | y),

对于所有的 :math:`i` 我们可以进一步简化

.. math::

   P(y \mid x_1, \dots, x_n) = \frac{P(y) \prod_{i=1}^{n} P(x_i \mid y)}
                                    {P(x_1, \dots, x_n)}

因为 :math:`P(x_1, \dots, x_n)` 是输入的常量，因此我们可以采用如下的分类规则：

.. math::

   P(y \mid x_1, \dots, x_n) \propto P(y) \prod_{i=1}^{n} P(x_i \mid y)

   \Downarrow

   \hat{y} = \arg\max_y P(y) \prod_{i=1}^{n} P(x_i \mid y),

其中，我们可以采用最大后验概率（Maximum A Posteriori MAP）来估计 :math:`P(y)` 和 :math:`P(x_i \mid y)` 。前者是 :math:`y` 类别在训练样本中的频率。

不同的朴素贝叶斯分类器的区别在于它们对于概率分布 :math:`P(x_i \mid y)` 的假设。

尽管其假设是过度简化的，但是朴素贝叶斯分类器在实际应用的表现很好，例如文档分类，垃圾邮件检测。他们只需要一小部分训练样本就可以统计出必要的系数（关于这点的理论讨论见参考）。

朴素贝叶斯分类器相比复杂的模型非常迅速。将分类的条件特征分离意味着假设每个分布都是独立的一维分布。这意味着可以减轻高维度带来的“诅咒”。

另一方面，尽管朴素贝叶斯的分类很有效，但是并不是很好的预测器。因此其输出地预测概率 ``predict_proba`` 并不可靠。


.. topic:: 参考:

 * H. Zhang (2004). `The optimality of Naive Bayes.
   <http://www.cs.unb.ca/profs/hzhang/publications/FLAIRS04ZhangH.pdf>`_
   Proc. FLAIRS.


高斯朴素贝叶斯
--------------------

:class:`GaussianNB` 采用高斯朴素贝叶斯算法进行分类。其假设特征的分布式高斯型的。

.. math::

   P(x_i \mid y) &= \frac{1}{\sqrt{2\pi\sigma^2_y}} \exp\left(-\frac{(x_i - \mu_y)^2}{2\sigma^2_y}\right)

其中参数 :math:`\sigma_y` 赫尔 :math:`\mu_y` 是通过最大似然方法进行统计的。

    >>> from sklearn import datasets
    >>> iris = datasets.load_iris()
    >>> from sklearn.naive_bayes import GaussianNB
    >>> gnb = GaussianNB()
    >>> y_pred = gnb.fit(iris.data, iris.target).predict(iris.data)
    >>> print("Number of mislabeled points out of a total %d points : %d"
    ...       % (iris.data.shape[0],(iris.target != y_pred).sum()))
    Number of mislabeled points out of a total 150 points : 6

.. _multinomial_naive_bayes:

多项式贝叶斯
-----------------------

:class:`MultinomialNB` 采用朴素贝叶斯算法针对多形式分布的概率。这是两个经典的处理文本分类的朴素贝叶斯算法。分布是有向量 :math:`\theta_y = (\theta_{y1},\ldots,\theta_{yn})` 来代表每一个类别 :math:`y` ， 其中 :math:`n` 是特征的数目，而 :math:`\theta_{yi}` 是 特征 :math:`i` 在样本被分类为 :math:`y` 中的概率 :math:`P(x_i \mid y)` 。

系数 :math:`\theta_y` 是通过平滑的最大似然概率方法，譬如相对频率计数：

.. math::

    \hat{\theta}_{yi} = \frac{ N_{yi} + \alpha}{N_y + \alpha n}

其中 :math:`N_{yi} = \sum_{x \in T} x_i` 是特征 :math:`i` 在样本 :math:`T` 中被分为特征 :math:`y` 的次数，而 :`N_{y} = \sum_{i=1}^{|T|} N_{yi}` 是样本被分为 :math:`y` 类的总数。

平滑的先验参数 :math:`\alpha \ge 0` 用来代表没有在样本中体现的特征，从而避免计算过程中出现0的情况。设置 :math:`\alpha = 1` 被称为拉普拉斯平滑，而 :math:`\alpha < 1` Lidstone平滑。


.. _bernoulli_naive_bayes:

伯努利朴素贝叶斯
---------------------

:class:`BernoulliNB` 通过朴素贝叶斯算法来解决数据服从多变量伯努利分布的情况，如某一个变量只有两个类别，且满足伯努利分布。因此这个算法要求样本有二值特征的特征向量，否则 ``BernoulliNB`` 将按照参数 ``binarize`` 进行分类。

伯努利朴素贝叶斯的分类规则依据如下；

.. math::

    P(x_i \mid y) = P(i \mid y) x_i + (1 - P(i \mid y)) (1 - x_i)

这与多项式朴素贝叶斯方法不同之处在于其对于没有出现的特征 :math:`i` 会加以惩罚，而多项式朴素贝叶斯不会。

对于文本分类，词频向量可以用来训练这个分类器。 ``BernoulliNB`` 有时会表现更好，尤其是文本较短的情况。如果时间允许，建议尝试两个方法。

.. topic:: 参考:

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


Out-of-core 朴素贝叶斯模型拟合
-------------------------------------

对于大数据的文类问题，朴素贝叶斯方法可能会超出内存的容量。对于这个问题 :class:`MultinomialNB` ， :class:`BernoulliNB` 和 :class:`GaussianNB` 可以通过采用 ``partial_fit`` 来避免（参见 :ref:`example_applications_plot_out_of_core_classification.py` ）。除了 :class:`GaussianNB` 其他两个分类器均支持取样权重。

区别于 ``fit`` ，当使用 ``partial_fit`` 需要传递所有的类别信息。

对于所有scikit-learn采用的策略，请参见 :ref:`out-of-core learning <scaling_strategies>` 。

note::

  ``partial_fit`` 会导致冗余。因此建议使用内存所允许的最大样本大小进行拟合。

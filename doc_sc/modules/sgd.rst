.. _sgd:

===========================
随机梯度下降法
===========================

.. currentmodule:: sklearn.linear_model

**随机梯度下降法 Stochastic Gradient Descent (SGD)** 是一个简单且有效的线性分类方法。其采用凸的损失函数，如线性 `Support Vector Machines <http://en.wikipedia.org/wiki/Support_vector_machine>`_ 和 `Logistic Regression <http://en.wikipedia.org/wiki/Logistic_regression>`_ 。尽管 SGD 已经在机器学习界存在了很长时间，目前其还是在大尺度学习领域的重要研究对象。

SGD在大尺度学习和稀疏问题上，尤其是文本分类和语言学习上有着成功的应用。由于数据的稀疏，常常会有10^5的样本，而其特征有10^5。

随机梯度下降法的优势在于:

    + 快速

    + 实现简单

劣势在于：

    + SGD 需要一些超参数，如复杂度参数，和重复步数。

    + SGD 对于特征的大小敏感。

分类
==============

.. warning::

  确定每次都对样本进行重新的排序，或者使用 ``shuffle=True`` 。

:class:`SGDClassifier` 类采用一个直接的随机梯度下降方法，并支持不同的损失函数和惩罚的分类。

.. figure:: ../auto_examples/linear_model/images/plot_sgd_separating_hyperplane_001.png
   :target: ../auto_examples/linear_model/plot_sgd_separating_hyperplane.html
   :align: center
   :scale: 75

如同其他的分类方法，SGD需要两个输入变量，大小为 [n_samples, n_features] 的训练样本X和大小为[n_samples]的目标分类::

    >>> from sklearn.linear_model import SGDClassifier
    >>> X = [[0., 0.], [1., 1.]]
    >>> y = [0, 1]
    >>> clf = SGDClassifier(loss="hinge", penalty="l2")
    >>> clf.fit(X, y)
    SGDClassifier(alpha=0.0001, class_weight=None, epsilon=0.1, eta0=0.0,
           fit_intercept=True, l1_ratio=0.15, learning_rate='optimal',
           loss='hinge', n_iter=5, n_jobs=1, penalty='l2', power_t=0.5,
           random_state=None, shuffle=False, verbose=0, warm_start=False)


拟合结束后，模型可以被用作预测::

    >>> clf.predict([[2., 2.]])
    array([1])

以上 SGD拟合一个线性模型，其成员 ``coef_`` 包含模型的系数::

    >>> clf.coef_
    array([[ 9.91080278,  9.91080278]])

成员 ``intercept_`` 为模型的截距（或者偏差）::

    >>> clf.intercept_                                    # doctest: +ELLIPSIS
    array([-9.990...])

是否使用截距可以通过 ``fit_intercept`` 来调节。

通过 :meth:`SGDClassifier.decision_function` 可以得到距离超平面的距离::

    >>> clf.decision_function([[2., 2.]])
    array([ 29.65318117])

具体的损失函数可以通过 ``loss`` 来调节。 :class:`SGDClassifier` 已包含如下的损失函数::

  * ``loss="hinge"``: (soft-margin) linear Support Vector Machine,
  * ``loss="modified_huber"``: smoothed hinge loss,
  * ``loss="log"``: logistic regression,
  * and all regression losses below.

前两个损失函数只会在样本违反了边界限制的事后才会更新模型的系数。这样可以使得拟合过程非常迅速，且模型比较疏散，即使采用L2的惩罚方式。

选择 ``loss="log"`` 或者 ``loss="modified_huber"`` 则可以使用 ``predict_proba`` 方法来估计每个样本 :math:`x` 预测值的概率 :math:`P(y|x)`::

    >>> clf = SGDClassifier(loss="log").fit(X, y)
    >>> clf.predict_proba([[1., 1.]])
    array([[ 0.0000005,  0.9999995]])

具体的惩罚方式可以通过 ``penalty`` 参数调节，SGD 包含如下方式::

  * ``penalty="l2"``: L2 norm penalty on ``coef_``.
  * ``penalty="l1"``: L1 norm penalty on ``coef_``.
  * ``penalty="elasticnet"``: Convex combination of L2 and L1;
    ``(1 - l1_ratio) * L2 + l1_ratio * L1``.

默认设置为 ``penalty="l2"`` 。 L1方法会导致稀疏的模型，使得大部分系数为0。弹性网络可以解决一定程度上L1带来的低效。参数 ``l1_ratio`` 用来调节L1与L2惩罚的相对权重。

:class:`SGDClassifier` 支持多类别分类，其策略是“一对所有”（OVA）。对 :math:`K` 分类中的每一个类别，构建一个分类器来区分其与其他 :math:`K-1` 个类别。在测试阶段，我们可以计算每一个预测的置信评分（如距离超平面的距离），并从中选择得分最高的分类。下图展示了本方法在iris数据中的应用。其中虚线表示三个分类器。背景颜色代表由分类器决定的分类空间。

.. figure:: ../auto_examples/linear_model/images/plot_sgd_iris_001.png
   :target: ../auto_examples/linear_model/plot_sgd_iris.html
   :align: center
   :scale: 75

对于多类别情况 ``coef_`` 是一个大小为 ``shape=[n_classes, n_features]`` 二维矩阵， ``intercept_`` 是一个大小为 ``shape=[n_classes]`` 的向量。 ``coef_`` 的第i行表示OVA的第i个分类器的权重向量。类别是按照升序排列（参见属性 ``classes_`` ）。注意，原则上， ``loss="log"`` 和 ``loss="modified_huber"`` 更适用于OVA分类，因为他们会计算概率模型。

:class:`SGDClassifier` 通过 ``class_weight`` 和 ``sample_weight`` 可以赋予类别或者样本更多的权重。参考下面的例子和函数 :meth:`SGDClassifier.fit` 的说明。

.. topic:: 示例:

 - :ref:`example_linear_model_plot_sgd_separating_hyperplane.py`,
 - :ref:`example_linear_model_plot_sgd_iris.py`
 - :ref:`example_linear_model_plot_sgd_weighted_samples.py`
 - :ref:`example_svm_plot_separating_hyperplane_unbalanced.py` (See the `Note`)

回归
==========

:class:`SGDRegressor` 采用一个直接的随机梯度下降方法来支持不同的损失函数和惩罚的分类来进行线性回归。 :class:`SGDRegressor` 适用于大数据情形 (> 10.000)，而其他情况，我们建议使用 :class:`Ridge` ， :class:`Lasso` 或者 :class:`ElasticNet` 。

具体的损失函数可以通过参数 ``loss`` 进行设定。 :class:`SGDRegressor` 包含以下几种情况：

  * ``loss="squared_loss"``: Ordinary least squares,
  * ``loss="huber"``: Huber loss for robust regression,
  * ``loss="epsilon_insensitive"``: linear Support Vector Regression.

Huber and epsilon-insensitive 损失函数可以被用作稳健回归。其中不敏感区域的大小由 ``epsilon`` 确定。其大小与样本特征的大小相关。


随机梯度下降方法方法在稀疏样本中应用
===========================================

.. note:: 这个针对稀疏样本的方法会以上的方法产生略微不同的结果。这是由于降低了对截距的学习速率。

对于稀疏数据，我们利用 `scipy.sparse <http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.html>`_. 构建稀疏矩阵。为了最大化效率，建议使用CSR矩阵格式 `scipy.sparse.csr_matrix <http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html>`_.

.. topic:: 示例:

 - :ref:`example_text_document_classification_20newsgroups.py`

计算复杂度
==========

SGD的主要优势在于计算效率，其计算复杂度正比于样本数目。对于大小为 (n, p) 的训练样本，其复杂度为 :math:`O(k n \bar p)` ，其中 :math:`k` 是重复次数， :math:`\bar p` 是单个样本非零特征的数目。

最近的理论研究表明，计算时间并不是严的随着样本的增大而增加。

使用技巧
=====================

  * 随机梯度下降方法对于数据特征的大小敏感，因此建议对在拟合之前对数据进行重整化。如将输入变量 X 整理到 [0,1] 或 [-1,+1] 区间，或者使其中值为0，方差为1。 注意同样的重整化需要应用到测试样本中。此不可以通过类 :class:`StandardScaler` 完成::

      from sklearn.preprocessing import StandardScaler
      scaler = StandardScaler()
      scaler.fit(X_train)  # Don't cheat - fit only on training data
      X_train = scaler.transform(X_train)
      X_test = scaler.transform(X_test)  # apply same transformation to test data

    如果样本属性有一些本征特质（如词频，特征），那么重整化是不必要的。

  * 找到适合的复杂度参数 :math:`\alpha` 。可以通过 :class:`GridSearchCV` 来自动寻找，通常取值范围为 ``10.0**-np.arange(1,7)`` 。

  * 经验上讲，我们发现 SGD 在观测10^6个样本后收敛。因此一个对于步数的合理估计为 ``n_iter = np.ceil(10**6 / n)`` ，其中 ``n`` 是样本大小。

  * 如果你用SGD来作PCA分析，那么通常将特征乘以一个系数 `c` 来使得其方差为1。


.. topic:: 参考:

 * `"Efficient BackProp" <yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf>`_
   Y. LeCun, L. Bottou, G. Orr, K. Müller - In Neural Networks: Tricks
   of the Trade 1998.

.. _sgd_mathematical_formulation:

数学基础
========================

给定一系列训练样本 :math:`(x_1, y_1), \ldots, (x_n, y_n)` 其中 :math:`x_i \in \mathbf{R}^n` 和 :math:`y_i \in \{-1,1\}` ， 我们的目标是寻找一个线性评分函数 :math:`f(x) = w^T x + b` 其中系数为 :math:`w \in \mathbf{R}^m` ，截距为 :math:`b \in \mathbf{R}` 。为了做预测，我们简单地判断 :math:`f(x)` 的符号。  一个常用的标准是最小化训练误差：

.. math::

    E(w,b) = \frac{1}{n}\sum_{i=1}^{n} L(y_i, f(x_i)) + \alpha R(w)

其中 :math:`L` 是损失函数，用以测量预测的准确程度。
:math:`R` 是惩罚函数来控制模型的复杂度。 :math:`\alpha > 0` 是一个大于
0的超系数。

:math:`L` 的不同选择会产生不同的分类方法：

   - Hinge: (soft-margin) Support Vector Machines.
   - Log:   Logistic Regression.
   - Least-Squares: Ridge Regression.
   - Epsilon-Insensitive: (soft-margin) Support Vector Regression.

以上损失函数所对应的误差上限展示如下::

.. figure:: ../auto_examples/linear_model/images/plot_sgd_loss_functions_001.png
   :align: center
   :scale: 75

常用的复杂度参数:math:`R` 选择如下：

   - L2 norm: :math:`R(w) := \frac{1}{2} \sum_{i=1}^{n} w_i^2`,
   - L1 norm: :math:`R(w) := \sum_{i=1}^{n} |w_i|`, 导致稀疏模型
   - 弹性网络: :math:`R(w) := \frac{\rho}{2} \sum_{i=1}^{n} w_i^2 +
     (1-\rho) \sum_{i=1}^{n} |w_i|`, 通过 :math:`\rho` 来调节 L1与L2的
     权重

下图展示了在系数空间中不同复杂度参数选择的影响 :math:`R(w) = 1`.

.. figure:: ../auto_examples/linear_model/images/plot_sgd_penalties_001.png
   :align: center
   :scale: 75

SGD
---

随机下降方法是一个无限制条件的优化方法。对比一般的梯度下降方法，SGD每次只通过一个样本数据来近似计算真实的梯度 :math:`E(w,b)` 。 :class:`SGDClassifier` 采用一阶SGD学习方法。此方法重复应用到训练样本上，对于每一个样本，其按照如下方法更新模型系数::

.. math::

    w \leftarrow w - \eta (\alpha \frac{\partial R(w)}{\partial w}
    + \frac{\partial L(w^T x_i + b, y_i)}{\partial w})

其中 :math:`\eta` 是学习效率，控制在系数空间的步长。 截距 :math:`b` 用同样的办法更新，但是不需要控制步长。

学习速率 :math:`\eta` 可以是个常数，也可以逐步降低。例如，对于分类，默认的学习效率按照如下方法更新 (``learning_rate='optimal'``)::

.. math::

    \eta^{(t)} = \frac {1}{\alpha  (t_0 + t)}

其中 :math:`t` 是时间步数，（一共会有 `n_samples * n_iter` 时间步）， :math:`t_0` 由乐观估计（参见 Léon Bottou）。预期的初始更新和预期的权重大小相一致（假设训练样本的方差为1）其具体的定义可以参见 :class:`BaseSGD` 中的  ``_init_t`` 。


对于回归问题，默认的学习效率如下 (``learning_rate='invscaling'``)

.. math::

    \eta^{(t)} = \frac{eta_0}{t^{power\_t}}

其中 :math:`eta_0` 和 :math:`power\_t` 是超系数，由 ``eta0`` 和 ``power_t`` 设定。

设定 ``learning_rate='constant'`` 和 ``eta0`` 来选择恒定学习效率。

模型的系数可以通过 ``coef_`` 和 ``intercept_`` 来获得::

     - ``coef_`` 对应权重 :math:`w`

     -  ``intercept_`` 对应截距 :math:`b`

.. topic:: 参考:

 * `"Solving large scale linear prediction problems using stochastic
   gradient descent algorithms"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.58.7377>`_
   T. Zhang - In Proceedings of ICML '04.

 * `"Regularization and variable selection via the elastic net"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.124.4696>`_
   H. Zou, T. Hastie - Journal of the Royal Statistical Society Series B,
   67 (2), 301-320.


实现细节
======================

SGD方法的实现是基于 Léon Bottou `Stochastic Gradient SVM
<http://leon.bottou.org/projects/sgd>`_  。 类似于 SvmSGD 权重向量是一个标量与向量的乘积，这样可以有效地在L2调控下更新权重。 对于稀疏特征向量，由于更频繁的更新，截距会采用一个更小的学习效率下（一般情况的1%）。训练样本会按照顺序进行分析，而学习效率会按照 Shalev-Shwartz et al. 2007 逐步降低。对于多类别分类的情况，我们采用“一对多”的方案。我们采用 Tsuruoka et al. 2009 的限定梯度方法来控制L1权重（以及弹性网络）。程序由Cython完成。

.. topic:: 参考:

 * `"Stochastic Gradient Descent" <http://leon.bottou.org/projects/sgd>`_ L. Bottou - Website, 2010.

 * `"The Tradeoffs of Large Scale Machine Learning" <http://leon.bottou.org/slides/largescale/lstut.pdf>`_ L. Bottou - Website, 2011.

 * `"Pegasos: Primal estimated sub-gradient solver for svm"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.74.8513>`_
   S. Shalev-Shwartz, Y. Singer, N. Srebro - In Proceedings of ICML '07.

 * `"Stochastic gradient descent training for l1-regularized log-linear models with cumulative penalty"
   <http://www.aclweb.org/anthology/P/P09/P09-1054.pdf>`_
   Y. Tsuruoka, J. Tsujii, S. Ananiadou -  In Proceedings of the AFNLP/ACL '09.

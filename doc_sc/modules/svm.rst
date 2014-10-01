.. _svm:

=======================
支持向量机
=======================

.. currentmodule:: sklearn.svm

**支持向量机 Support vector machines (SVMs)** 是一系列监督式学习的方法，被用 :ref:`分类 <svm_classification>` ， :ref:`回归 <svm_regression>` 和 :ref:`检测异常 <svm_outlier_detection>` 。

支持向量机的优势在于：

    - 对高维数据有效。

    - 当数据维度多于样本数目时，本方法依然有效。

    - 使用一部分训练数据进行特征识别（叫做support vector），因此在内存上也是有效的

    - 灵活性：不同的 :ref:`svm_kernels` 可以被定义为决策函数。程序已经包含了常用的函数核，用户也可以自定义函数核。

支持向量机的劣势在于：

    - 当数据维度多于样本数目时，计算效率较低。

    - SVM 并不能直接给出概率估计，因此需要通过5层交叉检验。(见下 :ref:`Scores and probabilities <scores_probabilities>`）。

scikit-learn中的支持向量机可以计算密集矩阵（ ``numpy.ndarray`` 及 ``numpy.asarray`` ） 或者稀疏矩阵 (``scipy.sparse``) 样本向量的情形。当用SVM预测稀疏矩阵时，需要先对数据进行拟合。为了最佳的执行效率，选择 C-ordered ``numpy.ndarray`` 密集矩阵和 ``scipy.sparse.csr_matrix`` 稀疏矩阵，并且 ``dtype=float64`` 。

.. _svm_classification:

分类
==============

:class:`SVC`, :class:`NuSVC` 和 :class:`LinearSVC` 是进行多类别分类的工具。


.. figure:: ../auto_examples/svm/images/plot_iris_001.png
   :target: ../auto_examples/svm/plot_iris.html
   :align: center


:class:`SVC` 和 :class:`NuSVC` 是类似的方法，但输入的参数略有不同，并基于不同的数学公式（参见 :ref:`svm_mathematical_formulation` ）。另一方面， :class:`LinearSVC` 是通过线性核来进行支持向量分类。注意:class:`LinearSVC` 并不接受参数 ``kernel`` ，因为其已经假设了线性，另外其也不包括一些 :class:`SVC` 和 :class:`NuSVC` 的成员，如 ``support_`` 。

如另外两种方法 :class:`SVC`, :class:`NuSVC` ， :class:`LinearSVC` 接受两个输入向量， ``[n_samples, n_features]`` 大小的向量 X包含训练数据，大小 ``[n_samples]`` 分类Y包含字符串或者整数::

    >>> from sklearn import svm
    >>> X = [[0, 0], [1, 1]]
    >>> y = [0, 1]
    >>> clf = svm.SVC()
    >>> clf.fit(X, y)  # doctest: +NORMALIZE_WHITESPACE
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
    gamma=0.0, kernel='rbf', max_iter=-1, probability=False, random_state=None,
    shrinking=True, tol=0.001, verbose=False)

拟合结束之后，模型可以用来作预测::

    >>> clf.predict([[2., 2.]])
    array([1])

SVM的决策函数依赖于训练数据的子集，称之为支持向量。这些支持向量的性质可以在成员 ``support_vectors_``, ``support_`` 和
``n_support`` 中找到::

    >>> # get support vectors
    >>> clf.support_vectors_
    array([[ 0.,  0.],
           [ 1.,  1.]])
    >>> # get indices of support vectors
    >>> clf.support_ # doctest: +ELLIPSIS
    array([0, 1]...)
    >>> # get number of support vectors for each class
    >>> clf.n_support_ # doctest: +ELLIPSIS
    array([1, 1]...)

.. _svm_multi_class:

多类别分类
--------------------------

:class:`SVC` 和 :class:`NuSVC` 采用“一对一”的模式 (Knerr et al., 1990) 进行多类别分类。 如果 ``n_class`` 是分类的数目，那么 ``n_class * (n_class - 1) / 2`` 分类器将被构造出来，每一个将训练数据分为两个分类::

    >>> X = [[0], [1], [2], [3]]
    >>> Y = [0, 1, 2, 3]
    >>> clf = svm.SVC()
    >>> clf.fit(X, Y) # doctest: +NORMALIZE_WHITESPACE
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
    gamma=0.0, kernel='rbf', max_iter=-1, probability=False, random_state=None,
    shrinking=True, tol=0.001, verbose=False)
    >>> dec = clf.decision_function([[1]])
    >>> dec.shape[1] # 4 classes: 4*3/2 = 6
    6

另一方面 :class:`LinearSVC` 采用“一对多”的分类模式。因此需要训练 ``n_class`` 个分类器。当仅有两个分类时，只会有一个分类器::

    >>> lin_clf = svm.LinearSVC()
    >>> lin_clf.fit(X, Y) # doctest: +NORMALIZE_WHITESPACE
    LinearSVC(C=1.0, class_weight=None, dual=True, fit_intercept=True,
    intercept_scaling=1, loss='l2', max_iter=1000, multi_class='ovr',
    penalty='l2', random_state=None, tol=0.0001, verbose=0)
    >>> dec = lin_clf.decision_function([[1]])
    >>> dec.shape[1]
    4

决策函数的具体形式请参见 :ref:`svm_mathematical_formulation`

注意 :class:`LinearSVC` 还采用了另一种由 Crammer and Singer提出的分类方法，可以通过参数 ``multi_class='crammer_singer'`` 设定。这个方法是稳定的，而一般的一对多方法并不能保证。在实际应用当中，一对多分类方法更常被采用，因为其结果基本一致，且运行时间更短。

对于一对多的 :class:`LinearSVC` 方法，其属性 ``coef_`` 和 ``intercept_`` 的维度为 ``[n_class, n_features]`` 和 ``[n_class]`` 。每一行的系数代表了 ``n_class`` 个一对多分类器。

对于一对一的 " :class:`SVC` 方法，其属性的结构有点复杂。对于线性核，属性 ``coef_`` 和 ``intercept_``与 :class:`LinearSVC` 的结构基本一致，除了 ``coef_`` 变为了 ``[n_class * (n_class - 1) / 2, n_features]`` 以对应相应的分类器。 0 到 n 是 "0 vs 1", "0 vs 2" , ... "0 vs n"，进而是 "1 vs 2", "1 vs 3", "1 vs n", . . . "n-1 vs n".

``dual_coef_`` 的形状是 ``[n_class-1, n_SV]`` 。其每一列对应分类中用到的支持向量。每一行对应相应分类器的系数。

下面的例子可以更清楚的解释这个问题：

考虑一个三类别分类的例子，其中类别0有三个支持向量 :math:`v^{0}_0, v^{1}_0, v^{2}_0` 而分类1和2各有两个支持向量 :math:`v^{0}_1, v^{1}_1` 和 :math:`v^{0}_2, v^{1}_2` 。对于每一个支持向量 :math:`v^{j}_i` 有两个系数，我们将用以区分类别i和k的系数记作 :math:`\alpha^{j}_{i,k}` 。那么 ``dual_coef_`` 的存储如下：

+------------------------+------------------------+------------------+
|:math:`\alpha^{0}_{0,1}`|:math:`\alpha^{0}_{0,2}`|Coefficients      |
+------------------------+------------------------+for SVs of class 0|
|:math:`\alpha^{1}_{0,1}`|:math:`\alpha^{1}_{0,2}`|                  |
+------------------------+------------------------+                  |
|:math:`\alpha^{2}_{0,1}`|:math:`\alpha^{2}_{0,2}`|                  |
+------------------------+------------------------+------------------+
|:math:`\alpha^{0}_{1,0}`|:math:`\alpha^{0}_{1,2}`|Coefficients      |
+------------------------+------------------------+for SVs of class 1|
|:math:`\alpha^{1}_{1,0}`|:math:`\alpha^{1}_{1,2}`|                  |
+------------------------+------------------------+------------------+
|:math:`\alpha^{0}_{2,0}`|:math:`\alpha^{0}_{2,1}`|Coefficients      |
+------------------------+------------------------+for SVs of class 2|
|:math:`\alpha^{1}_{2,0}`|:math:`\alpha^{1}_{2,1}`|                  |
+------------------------+------------------------+------------------+


.. _scores_probabilities:

评分和概率
------------------------

:class:`SVC` 方法中的决策函数 ``decision_function`` 赋予每一个样本的每一个分类一个评分（两个分类时，只有一个评分）。当函数的属性``probability`` 被设置为 ``True`` 时，会进一步估计分类的概率（通过成员函数 ``predict_proba`` 和 ``predict_log_proba`` 。对于两个分类的情况，概率是通过Platt scaling进行校准：通过逻辑回归SVM的评分，再进一步交叉检验。对于多种分类，请参考 Wu et al. (2004).

对于大数据，Platt scaling是耗时的。此外，概率估计可能与评分并不一致。因此最高评分的样本与最高概率的样本并不可能是一致的。（例如，一个样本被预测 ``predict`` 为某一个分类，其概率``predict_proba`` 可能小于50%。）Platt 方法也存在理论问题。如果需要一个置信评分，但是不需要概率，那么最好设置 ``probability=False`` 并选择 ``decision_function`` 而不是 ``predict_proba`` 。

.. topic:: 参考:

 * Wu, Lin and Weng,
   `"Probability estimates for multi-class classification by pairwise coupling"
   <http://www.csie.ntu.edu.tw/~cjlin/papers/svmprob/svmprob.pdf>`_.
   JMLR 5:975-1005, 2004.


非平衡问题
--------------------

对于需要强调某些分类或者样本，可以使用 ``class_weight`` 和 ``sample_weight`` 进行调节。

:class:`SVC` (注意没有 :class:`NuSVC`) 在 ``fit`` 方法中含有一个参数 ``class_weight`` 。它是一个字典类型，
``{class_label : value}``, 其中value是一个大于0的数值来强调系数 ``C`` 分类 ``class_label`` ，使其变为 ``C * value`` 。

.. figure:: ../auto_examples/svm/images/plot_separating_hyperplane_unbalanced_001.png
   :target: ../auto_examples/svm/plot_separating_hyperplane_unbalanced.html
   :align: center
   :scale: 75


:class:`SVC`, :class:`NuSVC`, :class:`SVR`, :class:`NuSVR` 和:class:`OneClassSVM` 也可以包含针对个别样本的权重。通过函数``fit`` 中的参数 ``sample_weight`` 。与 ``class_weight`` 类似，其使得系数 ``C`` 变为 ``C * sample_weight[i]`` 。

.. figure:: ../auto_examples/svm/images/plot_weighted_samples_001.png
   :target: ../auto_examples/svm/plot_weighted_samples.html
   :align: center
   :scale: 75


.. topic:: 示例:

 * :ref:`example_svm_plot_iris.py`,
 * :ref:`example_svm_plot_separating_hyperplane.py`,
 * :ref:`example_svm_plot_separating_hyperplane_unbalanced.py`
 * :ref:`example_svm_plot_svm_anova.py`,
 * :ref:`example_svm_plot_svm_nonlinear.py`
 * :ref:`example_svm_plot_weighted_samples.py`,


.. _svm_regression:

回归 
==========

支持向量分类的方法可以扩展到回归问题。此方法被成为支持向量回归。

支持向量分类只需要使用一部分的训练数据，因为其不在乎超出边际的样本。与此类似，支持向量回归同样也只依赖于一部分的训练数据，因为成本函数忽略那些接近模型预测值的样本。有两个支持向量回归的模型: :class:`SVR` 和 :class:`NuSVR`.

同分类一样，拟合函数需要两个输入变量X和y。不同之处在于y可以是浮点类型而不必是整数：

    >>> from sklearn import svm
    >>> X = [[0, 0], [2, 2]]
    >>> y = [0.5, 2.5]
    >>> clf = svm.SVR()
    >>> clf.fit(X, y) # doctest: +NORMALIZE_WHITESPACE
    SVR(C=1.0, cache_size=200, coef0=0.0, degree=3,
    epsilon=0.1, gamma=0.0, kernel='rbf', max_iter=-1, probability=False,
    random_state=None, shrinking=True, tol=0.001, verbose=False)
    >>> clf.predict([[1, 1]])
    array([ 1.5])


.. topic:: 示例:

 * :ref:`example_svm_plot_svm_regression.py`

.. _svm_outlier_detection:

密度分析，新奇检测
=======================================

单一分类的SVM模型可以用来作新奇检测。给定一个数据样本，它将探测到样本的边界，并判定新的数据是否属于原先的分类。 :class:`OneClassSVM` 用来实现这个目的。这是一个非监督式学习的类型。因此它的拟合方法只需要一个输入样本X，而不需要分类的标记y。

参见 :ref:`outlier_detection` 。

.. figure:: ../auto_examples/svm/images/plot_oneclass_001.png
   :target: ../auto_examples/svm/plot_oneclass.html
   :align: center
   :scale: 75


.. topic:: 示例:

 * :ref:`example_svm_plot_oneclass.py`
 * :ref:`example_applications_plot_species_distribution_modeling.py`


复杂度
==========

SVM是一个有效的工具，但是其计算和存储需要随着训练向量的增加而迅速增加。SMV的核心是个四阶计算问题，其求解用到了库 `libsvm`_ 其计算复杂度介于 :math:`O(n_{features} \times n_{samples}^2)` 与 :math:`O(n_{features} \times n_{samples}^3)` 之间（依赖于数据在缓存中的存储）。 如果数据很稀疏，那么 :math:`n_{features}` 应被替代为平均样本含有的观测特征数目。

另外注意，对于线性情况，:class:`LinearSVC` 采用库 `liblinear`_ ， 这要比基于 `libsvm`_ 的 :class:`SVC` 更为快速。


使用技巧
=====================


  * **避免数据拷贝**: 对于 :class:`SVC`, :class:`SVR`, :class:`NuSVC` 和 :class:`NuSVR` ，如果数据不是 C-ordered 连续，双精度，那么数据将被拷贝再转入C程序。因此运算前可以检查numpy 数组的 ``flags`` 属性。

     对于 :class:`LinearSVC` （和 :class:`LogisticRegression <sklearn.linear_model.LogisticRegression>` ）任何输入数组都会被转换到 `liblinear`_ 的数据格式。如果你希望拟合大数据而不拷贝数据，那么最好使用 :class:`SGDClassifier <sklearn.linear_model.SGDClassifier>` 。该方法的目标函数与 :class:`LinearSVC` 基本一致。

  * **核心缓存大小**: 对于 :class:`SVC`, :class:`SVR`, :class:`nuSVC` 和 :class:`NuSVR` ，其核心缓存大小会显著影响计算效率。如果你的内存比较大，那么建议设置 ``cache_size`` 到更大的值。其默认值为 200(MB)，建议提高到 500(MB) 或者 1000(MB) 。

  * **设置 C**: ``C`` 默认是一个 ``1`` 。通常是一个好的选择，如果你的数据噪声很大，那么适当降低取值可以让结果更稳定。

  * SVM不是一个比例不变的算法， 所以 **强烈建议将数据调整到适当比例** 。譬如将样本的特征取值整理到[0,1]或者[-1,1]，或者中值为零，方差为1。注意，同样的变换也要应用到测试向量上去，才可以得到有意义的结果。重整化的部分可以参考 :ref:`preprocessing` 。

  * :class:`NuSVC`/:class:`OneClassSVM`/:class:`NuSVR` 中的系数 ``nu`` 估算误差与支持向量的比例。

  * 在 :class:`SVC` ，如果数据是不平衡的（如很多正结果，很少负结果），那么设置 ``class_weight='auto'`` 或者尝试不同的 ``C`` 。

  * :class:`LinearSVC` 采用一个随机数来选择特征来拟合模型，所以每次的结果并不会非常一致。 此时可以降低容差参数 ``tol`` 。

  * 设定一阶复杂度参数``LinearSVC(loss='l2', penalty='l1', dual=False)`` 可以产生一个稀疏的结果。如，只有一部分的特征会贡献到决策函数中去。增加 ``C`` 会导致更多的特征被选择出来。避免空模型的最小的 ``C`` 可以通过函数 :func:`l1_min_c` 计算出来。

.. _svm_kernels:

函数核
================

函数核可以有如下选择：

  * 线性: :math:`\langle x, x'\rangle`.

  * 多项式: :math:`(\gamma \langle x, x'\rangle + r)^d`.
    :math:`d` 代表 ``degree``, :math:`r` 表示 ``coef0``.

  * rbf: :math:`\exp(-\gamma |x-x'|^2)`. :math:`\gamma` 代表 ``gamma`` ，且大于0。

  * sigmoid (:math:`\tanh(\gamma \langle x,x'\rangle + r)`)，其中 :math:`r`代表 ``coef0``.

不同的函数核需要在初始化时加以明确：

    >>> linear_svc = svm.SVC(kernel='linear')
    >>> linear_svc.kernel
    'linear'
    >>> rbf_svc = svm.SVC(kernel='rbf')
    >>> rbf_svc.kernel
    'rbf'


定制函数核
--------------

你可以定义选择自己的函数核，其可以是Python函数或者 Gram矩阵。

拥有自定义函数核的分类器的用法基本一致，除了：

    * 值 ``support_vectors_`` 为空，只有其索引被存在 ``support_``

    * 训练数据的指针被存在 ``fit()`` 函数中，而不是拷贝。因此任何对数据的更改将导致 ``predict()`` 函数的异常。


使用Python函数作为程序核
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

你可以将自定义的函数传递到 ``kernel`` 中。

自定义的函数核必须接受两个矩阵作为输入，并返回第三个矩阵

一下例子定义了一个函数核::

    >>> import numpy as np
    >>> from sklearn import svm
    >>> def my_kernel(x, y):
    ...     return np.dot(x, y.T)
    ...
    >>> clf = svm.SVC(kernel=my_kernel)

.. topic:: 示例:

 * :ref:`example_svm_plot_custom_kernel.py`.

使用 Gram 矩阵
~~~~~~~~~~~~~~~~~~~~~

设置 ``kernel='precomputed'`` 并将 Gram矩阵取代数据X作为 ``fit`` 函数的输入。此时对于测试样本，同样需要输入 Gram 矩阵。

    >>> import numpy as np
    >>> from sklearn import svm
    >>> X = np.array([[0, 0], [1, 1]])
    >>> y = [0, 1]
    >>> clf = svm.SVC(kernel='precomputed')
    >>> # linear kernel computation
    >>> gram = np.dot(X, X.T)
    >>> clf.fit(gram, y) # doctest: +NORMALIZE_WHITESPACE
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
    gamma=0.0, kernel='precomputed', max_iter=-1, probability=False,
    random_state=None, shrinking=True, tol=0.001, verbose=False)
    >>> # predict on training examples
    >>> clf.predict(gram)
    array([0, 1])

RBF 核系数
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

当通过 *径向基函数 Radial Basis Function* (RBF) 来训练SVM时，需要考虑两个系数 ``C`` 和 ``gamma`` 。 系数 ``C`` 用来平衡错误的分类和模型的简化。一个较低的 ``C`` 将使得判定面比较平滑，而较大的 ``C`` 尝试将所有的训练样本预测正确。  ``gamma`` 定义每一个单独样本对你和的贡献。较大的 ``gamma`` 对邻近的样本会有更大的影响。

适当的选择 ``C`` 和 ``gamma`` 是影响SVM表现的关键。一般建议使用 :class:`GridSearchCV` 来寻找合适的 ``C`` 和 ``gamma`` 。

.. topic:: 示例:

 * :ref:`example_svm_plot_rbf_parameters.py`

.. _svm_mathematical_formulation:

数学基础
========================

支持向量机在高维中构造一个超平面或者一系列超平面来进行分类，回归或者其他任务。直觉上说，一个好的分割到最近的样本有最大的间距。一般而言，越大的间距会有着更准确的分类结果。

.. figure:: ../auto_examples/svm/images/plot_separating_hyperplane_001.png
   :align: center
   :scale: 75

SVC
---

给定训练向量 :math:`x_i \in R^p`, i=1,..., n, 在两个分类，和一个向量 :math:`y \in R^n` 有 :math:`y_i \in \{1, -1\}`, SVC 试图解决如下问题：

.. math::

    \min_ {w, b, \zeta} \frac{1}{2} w^T w + C \sum_{i=1, n} \zeta_i


    \textrm {subject to } & y_i (w^T \phi (x_i) + b) \geq 1 - \zeta_i,\\
    & \zeta_i \geq 0, i=1, ..., n

或者

.. math::

   \min_{\alpha} \frac{1}{2} \alpha^T Q \alpha - e^T \alpha


   \textrm {subject to } & y^T \alpha = 0\\
   & 0 \leq \alpha_i \leq C, i=1, ..., l

其中 :math:`e` 是单位向量， :math:`C > 0` 是上界， :math:`Q` 一个 :math:`n` 乘 :math:`n` 半正定矩阵。 :math:`Q_{ij} \equiv K(x_i, x_j)` 和 :math:`\phi (x_i)^T \phi (x)` 是函数核。此处训练向量被函数 :math:`\phi` 影身到一个更高维的空间。


决策函数为：

.. math:: \operatorname{sgn}(\sum_{i=1}^n y_i \alpha_i K(x_i, x) + \rho)

.. note::

    由于SVM采用库 `libsvm`_  和 `liblinear`_ 沿用 ``C`` 作为控制参数，其他大部分分析选择用 ``alpha`` 。 两者关系如下 :math:`C = \frac{n\_samples}{alpha}` 。

.. TODO multiclass case ?/

系数 :math:`y_i \alpha_i` 存在 ``dual_coef_` 中，支持向量存在 ``support_vectors_`` ，独立想 :math:`-\rho` 存在 ``intercept_`` 。

.. topic:: 参考:

 * `"Automatic Capacity Tuning of Very Large VC-dimension Classifiers"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.17.7215>`_
   I Guyon, B Boser, V Vapnik - Advances in neural information
   processing 1993,


 * `"Support-vector networks"
   <http://www.springerlink.com/content/k238jx04hm87j80g/>`_
   C. Cortes, V. Vapnik, Machine Leaming, 20, 273-297 (1995)



NuSVC
-----

这里我们介绍一个新的参数 :math:`\nu` 来控制支持向量的数目和误差。这个参数 :math:`\nu \in (0, 1]` 调控误差的上界和支持向量的下界。

可以证明， :math:`\nu` 是一个对 :math:`C` 的重新表达，在数学上是等价的。


实现细节
======================

在内部实现上，我们采用 `libsvm`_ 和 `liblinear`_ 来处理所有的计算。这些库是通过C和Cython封装的。

.. _`libsvm`: http://www.csie.ntu.edu.tw/~cjlin/libsvm/
.. _`liblinear`: http://www.csie.ntu.edu.tw/~cjlin/liblinear/

.. topic:: 参考:

  关于库的具体细节，请参考

    - `LIBSVM: a library for Support Vector Machines
      <http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf>`_

    - `LIBLINEAR -- A Library for Large Linear Classification
      <http://www.csie.ntu.edu.tw/~cjlin/liblinear/>`_



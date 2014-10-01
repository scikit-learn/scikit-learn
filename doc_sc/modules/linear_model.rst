.. _linear_model:

=========================
广义线性模型 （GLM）
=========================

.. currentmodule:: sklearn.linear_model

下面将介绍一系列回归方法。回归方法主要针对目标值是一系列输入变量的线性
组合。数学上说，如果 :math:`\hat{y}` 是目标值

.. math::    \hat{y}(w, x) = w_0 + w_1 x_1 + ... + w_p x_p

在程序中，我们用 ``coef_`` 代表 :math:`w = (w_1,
..., w_p)` 用 ``intercept_`` 代表 :math:`w_0` 。

通过线性模型对数据进行分类的讨论请见 :ref:`Logistic_regression`.


.. _ordinary_least_squares:

最小二乘法 （OLS）
=======================

:class:`LinearRegression` 拟合线性模型的系数
:math:`w = (w_1, ..., w_p)` 来最小化观测数据与预测值的总方差。数学上
      讲，这是在解决如下问题：

.. math:: \underset{w}{min\,} {|| X w - y||_2}^2

.. figure:: ../auto_examples/linear_model/images/plot_ols_001.png
   :target: ../auto_examples/linear_model/plot_ols.html
   :align: center
   :scale: 50%

:class:`LinearRegression` 通过 `fit` 函数拟合输入的 X, y
并将系数 :math:`w` 存在变量 `coef\_` 中:

    >>> from sklearn import linear_model
    >>> clf = linear_model.LinearRegression()
    >>> clf.fit ([[0, 0], [1, 1], [2, 2]], [0, 1, 2])
    LinearRegression(copy_X=True, fit_intercept=True, normalize=False)
    >>> clf.coef_
    array([ 0.5,  0.5])

然而，最小二乘法的系数取决于模型参数的独立性。 当参数相关，那么
矩阵 :math:`X` 各列存在线性相关。因此矩阵近乎退化矩阵, 而观测值中的随
机误差将对最小二乘法的结果有较大影响。这个不谨慎的实验设计将导致这种 *
多重共线性* 。


.. topic:: 示例:

   * :ref:`example_linear_model_plot_ols.py`


最小二乘法计算复杂度
---------------------------------

本方法的计算复杂度取决于对矩阵X的奇异值分解。 对于矩阵 (n, p)，计算成本为 :math:`O(n p^2)`, 其中 :math:`n \geq p`.

.. _ridge_regression:

岭回归
================

:class:`Ridge` regression 通过降低复杂度，来解决一些
:ref:`ordinary_least_squares` 的问题。岭回归的系数将最小化以下方差和：


.. math::

   \underset{w}{min\,} {{|| X w - y||_2}^2 + \alpha {||w||_2}^2}

其中 :math:`\alpha \geq 0` 是一个复杂度参数，控制系数的复杂度。
:math:`\alpha` 越大，系数收敛越明显，也能有效地降低多重共线性的影响。

.. figure:: ../auto_examples/linear_model/images/plot_ridge_path_001.png
   :target: ../auto_examples/linear_model/plot_ridge_path.html
   :align: center
   :scale: 50%

如同其他线性模型 :class:`Ridge` 通过 `fit` 函数来读入数列 X, y 并将线
性模型系数 :math:`w` 存在类成员 `coef\_` 中::

    >>> from sklearn import linear_model
    >>> clf = linear_model.Ridge (alpha = .5)
    >>> clf.fit ([[0, 0], [0, 0], [1, 1]], [0, .1, 1]) # doctest: +NORMALIZE_WHITESPACE
    Ridge(alpha=0.5, copy_X=True, fit_intercept=True, max_iter=None,
          normalize=False, solver='auto', tol=0.001)
    >>> clf.coef_
    array([ 0.34545455,  0.34545455])
    >>> clf.intercept_ #doctest: +ELLIPSIS
    0.13636...


.. topic:: Examples:

   * :ref:`example_linear_model_plot_ridge_path.py`
   * :ref:`example_text_document_classification_20newsgroups.py`


岭回归计算复杂度
----------------

本方法的计算复杂度和 :ref:`ordinary_least_squares` 一致。

.. FIXME:
.. Not completely true: OLS is solved by an SVD, while Ridge is solved by
.. the method of normal equations (Cholesky), there is a big flop difference
.. between these


复杂度参数设置：广义交叉验证
------------------------------------------------------------------

:class:`RidgeCV` 采用交叉验证来调节岭回归的复杂度参数。此处的方法与
GridSearchCV 基本一致，除了采用默认的广义交叉验证（GCV）， 一个
高效的留一验证（LOOCV）::

    >>> from sklearn import linear_model
    >>> clf = linear_model.RidgeCV(alphas=[0.1, 1.0, 10.0])
    >>> clf.fit([[0, 0], [0, 0], [1, 1]], [0, .1, 1])       # doctest: +SKIP
    RidgeCV(alphas=[0.1, 1.0, 10.0], cv=None, fit_intercept=True, scoring=None,
        normalize=False)
    >>> clf.alpha_                                      # doctest: +SKIP
    0.1

.. topic:: 参考

    * "Notes on Regularized Least Squares", Rifkin & Lippert (`technical report
      <http://cbcl.mit.edu/projects/cbcl/publications/ps/MIT-CSAIL-TR-2007-025.pdf>`_,
      `course slides
      <http://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf>`_).


.. _lasso:

Lasso
=====

:class:`Lasso` 是岭回归的变种，其更有效地降低系数的数目。因此Lasso在
compressed sensing有着广泛的应用。在一定条件下，其计算结果能重
建非零系数（参见
:ref:`example_applications_plot_tomography_l1_reconstruction.py`）。

数学上讲，它包含一个线性模型的先验概率 :math:`\ell_1` 来进行调整。而目
标的最小化函数为：

.. math::  \underset{w}{min\,} { \frac{1}{2n_{samples}} ||X w - y||_2 ^ 2 + \alpha ||w||_1}

Lasso 的系数讲最小化方差和加上复杂度参数 :math:`\alpha ||w||_1` 。
:math:`\alpha` 是常数， :math:`||w||_1` 是 :math:`\ell_1`-norm 的系数
      向量

:class:`Lasso` 类采用的是梯度下降方法求解参数。另一种计算方法参见
       :ref:`least_angle_regression`::

    >>> clf = linear_model.Lasso(alpha = 0.1)
    >>> clf.fit([[0, 0], [1, 1]], [0, 1])
    Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
       normalize=False, positive=False, precompute='auto', tol=0.0001,
       warm_start=False)
    >>> clf.predict([[1, 1]])
    array([ 0.8])

此外，对于底层任务，函数 :func:`lasso_path` 在求解路径上的系数也很有用
处。

.. topic:: 示例:

  * :ref:`example_linear_model_plot_lasso_and_elasticnet.py`
  * :ref:`example_applications_plot_tomography_l1_reconstruction.py`


.. note:: **Lasso 用于特征选择**

      由于Lasso回归降低了模型的复杂度，因此可以用作特征选择。参见
      :ref:`l1_feature_selection`.

.. note:: **随机稀疏模型**

     对于特征选择和 sparse recovery，可以参见 :ref:`randomized_l1`.


复杂度参数设置
--------------------------------

`alpha` 用来控制模型的稀疏程度。

交叉检验
^^^^^^^^^^^^^^^^^^^^^^^

scikit-learn 通过交叉检验 :class:`LassoCV` 和 :class:`LassoLarsCV` 来
调节 `alpha` 系数的选择。 其中 :class:`LassoLarsCV` 是基于
:ref:`least_angle_regression` 算法，其解释如下：

对于存在共线性的高维数据，推荐使用 :class:`LassoCV` 。
而 :class:`LassoLarsCV` 在探索适用的 `alpha` 系数上更有优势，并且当样
本的数目小于观测数据的数目时，它比 :class:`LassoCV` 更为快速。

.. |lasso_cv_1| image:: ../auto_examples/linear_model/images/plot_lasso_model_selection_002.png
    :target: ../auto_examples/linear_model/plot_lasso_model_selection.html
    :scale: 48%

.. |lasso_cv_2| image:: ../auto_examples/linear_model/images/plot_lasso_model_selection_003.png
    :target: ../auto_examples/linear_model/plot_lasso_model_selection.html
    :scale: 48%

.. centered:: |lasso_cv_1| |lasso_cv_2|


基于信息量标准的模型选择
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

另一方面， :class:`LassoLarsIC` 通过采用 `赤池信息量准则
<http://zh.wikipedia.org/wiki/%E8%B5%A4%E6%B1%A0%E4%BF%A1%E6%81%AF%E9%87%8F%E5%87%86%E5%88%99>` 
(Akaike information criterion AIC) 和 贝叶斯信息量准则 (Bayes Information criterion BIC)。
这种寻找最佳 `alpha` 的方法相对于 k-fold 交叉检验更为快速，因为只需要
计算一次而不是 k+1 次 regularization path。然而，此方法需要适当的估计
问题的自由度。当问题的模型假设错误时，自由度的估计错误将会导致错误的结
果。

.. figure:: ../auto_examples/linear_model/images/plot_lasso_model_selection_001.png
    :target: ../auto_examples/linear_model/plot_lasso_model_selection.html
    :align: center
    :scale: 50%


.. topic:: 示例:

  * :ref:`example_linear_model_plot_lasso_model_selection.py`


弹性网络
===========
:class:`ElasticNet` 是线性模型配合一阶和二阶的复杂度来进行调节。
这个组合可以像 :class:`Lasso` 一样对稀疏模型进行拟合，同时保持 :class:`Ridge` 方法对于富再度的控制。 我们通过 `l1_ratio` 系数
来调节一阶和二阶系数的权重。

弹性网络对于不同特征间存在相关性的模型有更好的拟合结果。Lasso更倾向于
随机选择一个特征，而神经网络将选择所有的特征。

这个平衡Lasso和岭回归的优势在于神经网络在旋转变换下更为稳定。

此方法在于最小化目标函数如下：

.. math::

    \underset{w}{min\,} { \frac{1}{2n_{samples}} ||X w - y||_2 ^ 2 + \alpha \rho ||w||_1 +
    \frac{\alpha(1-\rho)}{2} ||w||_2 ^ 2}


.. figure:: ../auto_examples/linear_model/images/plot_lasso_coordinate_descent_path_001.png
   :target: ../auto_examples/linear_model/plot_lasso_coordinate_descent_path.html
   :align: center
   :scale: 50%

:class:`ElasticNetCV` 类可以通过交叉检验来调节系数 ``alpha``
       (:math:`\alpha`) 和 ``l1_ratio`` (:math:`\rho`)。

.. topic:: 示例:

  * :ref:`example_linear_model_plot_lasso_and_elasticnet.py`
  * :ref:`example_linear_model_plot_lasso_coordinate_descent_path.py`


.. _multi_task_lasso:

多任务Lasso
================

:class:`MultiTaskLasso` 类是同时计算多重回归问题的稀疏稀疏的线性模型：
       `y` 是一个二维数列，(n_samples, n_tasks)。 其中限制条件是所有的
       回归问题的特征是一样的，因此称之为多任务。

下图比较了通过简单的Lasso和 MultiTaskLasso 得到的非零系数 W。
对比可见MultiTaskLasso的系数更为统一。

.. |multi_task_lasso_1| image:: ../auto_examples/linear_model/images/plot_multi_task_lasso_support_001.png
    :target: ../auto_examples/linear_model/plot_multi_task_lasso_support.html
    :scale: 48%

.. |multi_task_lasso_2| image:: ../auto_examples/linear_model/images/plot_multi_task_lasso_support_002.png
    :target: ../auto_examples/linear_model/plot_multi_task_lasso_support.html
    :scale: 48%

.. centered:: |multi_task_lasso_1| |multi_task_lasso_2|

.. centered:: 时间序列模型拟合，任何有效的特征一直有效。

.. topic:: 示例:

  * :ref:`example_linear_model_plot_multi_task_lasso_support.py`

数学上讲，这是一个线性模型，通过
:math:`\ell_1` :math:`\ell_2` 来调节复杂度。
而最小化的目标函数是：

.. math::  \underset{w}{min\,} { \frac{1}{2n_{samples}} ||X W - Y||_2 ^ 2 + \alpha ||W||_{21}}

其中

.. math:: ||W||_21 = \sum_i \sqrt{\sum_j w_{ij}^2}


:class:`MultiTaskLasso` 梯度下降的方法来求解系数。

.. _least_angle_regression:

最小角回归 Least-angle regression
======================

最小角回归 (LARS) 是针对高维数据的回归算法，由Bradley Efron, Trevor Hastie, Iain
Johnstone and Robert Tibshirani 提出，参见 `wiki <http://en.wikipedia.org/wiki/Least-angle_regression>` 。

LARS的优势在于

  - 当 p >> n 时，计算有效 （当观测特征数目远大于样本数目）

  - 计算效率和 forward selection 一致，计算复杂度和最小二乘法一样。

  - 其计算结果包含求解路径，可以在交叉验证中重复使用

  - 当两个变量与输出变量相关性一致时，他们的系数变化也是类似的。此算法
    的行为与直觉相一致，而且更为稳定。

  - 可以很容易的修改成为Lasso。

LARS的劣势在于

  - 由于LARS是基于对残值的拟合。所以其对于噪声非常敏感。此点请参见 Weisberg
    在 Efron et al. (2004) Annals of Statistics 文章中的讨论。

LARS是通过 :class:`Lars` 类实现，或者其底层函数 :func:`lars_path`.


LARS Lasso
==========

:class:`LassoLars` 通过 LARS 算法的Lasso模型，而不依赖梯度下降方法。
this yields the exact solution, which is piecewise linear as a
function of the norm of its coefficients. （翻译不确定）

.. figure:: ../auto_examples/linear_model/images/plot_lasso_lars_001.png
   :target: ../auto_examples/linear_model/plot_lasso_lars.html
   :align: center
   :scale: 50%

::

   >>> from sklearn import linear_model
   >>> clf = linear_model.LassoLars(alpha=.1)
   >>> clf.fit([[0, 0], [1, 1]], [0, 1])  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
   LassoLars(alpha=0.1, copy_X=True, eps=..., fit_intercept=True,
        fit_path=True, max_iter=500, normalize=True, precompute='auto',
        verbose=False)
   >>> clf.coef_    # doctest: +ELLIPSIS
   array([ 0.717157...,  0.        ])

.. topic:: 示例:

 * :ref:`example_linear_model_plot_lasso_lars.py`

LARS算法提供系数的路径，并且几乎不依赖于复杂度参数。函数 :func:`lars_path`
用来获得系数路径。

数学形式
------------------------

本算法与 forward stepwise regression基本一致，除了每一步并不包含所有变
量。系数每一步向与梯度呈60度的方向前进。

LARS 的解并不是一个向量，而是一条曲线标记系数向量的一阶绝对值。
系数路径将存储在 ``coef_path_`` ，其大小为 (n_features,
max_features+1)。而且第一列为0。

.. topic:: 参考:

 * 算法的原始文献在 `Least Angle Regression
   <http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>`_
   by Hastie et al.


.. _omp:

正交匹配追寻 (OMP)
=================================
:class:`OrthogonalMatchingPursuit`类和函数 :func:`orthogonal_mp` 通过
OMP算法来求解线性模型，当未知系数的数目一直。

作为一个 forward feature selection 方法，如
:ref:`least_angle_regression` ，OMP将最优解近似为一系列非零元素：

.. math:: \text{arg\,min\,} ||y - X\gamma||_2^2 \text{ subject to } \
    ||\gamma||_0 \leq n_{nonzero\_coefs}

此外，OMP将针对特定的误差，而不是特定的系数数目。其数学表达如下：

.. math:: \text{arg\,min\,} ||\gamma||_0 \text{ subject to } ||y-X\gamma||_2^2 \
    \leq \text{tol}

OMP是基于贪婪算法，每一步寻找与当前残差最为相关的元素。这与简单地匹配
追寻相一致，但每一步更为有效。每一步的残差都是通过计算正交投影得到。


.. topic:: 示例:

 * :ref:`example_linear_model_plot_omp.py`

.. topic:: 参考:

 * http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

 * `Matching pursuits with time-frequency dictionaries
   <http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf>`_,
   S. G. Mallat, Z. Zhang,

贝叶斯回归
===================

贝叶斯回归技术被用来在计算系数时包含复杂度参数：复杂度参数不是提前预置，
而是随数据进行调整。

这是通过在模型中引入超参数（hyper parameter）来实现。参见 `uninformative priors
<http://en.wikipedia.org/wiki/Non-informative_prior#Uninformative_priors>`__

:math:`\ell_{2}` 复杂度参数在 `Ridge Regression`_ 等价于基于高斯先验分
布来寻找最大化后验概率的系数 :math:`w` 和精度 :math:`\lambda^{-1}` 。
区别于手动设置 `\lambda` ，此处我们可以将其作为来自数据的随机变量来处理。

作为完整的概率模型，输出变量 :math:`y` 将假设为围绕 :math:`X w` 的高斯
分布：

.. math::  p(y|X,w,\alpha) = \mathcal{N}(y|X w,\alpha)

:math:`\alpha` 如上所述将作为基于数据的随机变量。

贝叶斯回归的优势在于：

    - 它是基于现有数据的。

    - 它可以在系数拟合中包含复杂度参数。

劣势在于：

    - 推导模型比较耗时。


.. topic:: 参考

 * A good introduction to Bayesian methods is given in C. Bishop: Pattern
   Recognition and Machine learning

 * Original Algorithm is detailed in the  book `Bayesian learning for neural
   networks` by Radford M. Neal

.. _bayesian_ridge_regression:

贝叶斯岭回归
-------------------------

:class:`BayesianRidge` 类拟合回归问题的统计模型如下。
系数 :math:`w` 的先验分布是采用球对称高斯分布：

.. math:: p(w|\lambda) =
    \mathcal{N}(w|0,\lambda^{-1}\bold{I_{p}})

:math:`\alpha` 和 :math:`\lambda` 的先验分布是采用 `gamma
distributions <http://en.wikipedia.org/wiki/Gamma_distribution>`__， 
高斯分布的共轭先验分布。

这样的得到的模型被称为 *贝叶斯岭回归 Bayesian Ridge Regression* ，并与经典的
:class:`Ridge` 相类似。系数 :math:`w` ， :math:`\alpha` 和 
:math:`\lambda` 是在拟合过程中统一计算。仅剩的是超参数，系数
      :math:`\alpha` 和 :math:`\lambda`  的gamma先验分布
的参数。这些通常作为非信息 *non-informative* 。 最后系数是通过最大化 *边际似然值* 。

默认值 :math:`\alpha_1 = \alpha_2 =  \lambda_1 = \lambda_2 = 1.e^{-6}`.


.. figure:: ../auto_examples/linear_model/images/plot_bayesian_ridge_001.png
   :target: ../auto_examples/linear_model/plot_bayesian_ridge.html
   :align: center
   :scale: 50%


贝叶斯岭回归应用于回归拟合::

    >>> from sklearn import linear_model
    >>> X = [[0., 0.], [1., 1.], [2., 2.], [3., 3.]]
    >>> Y = [0., 1., 2., 3.]
    >>> clf = linear_model.BayesianRidge()
    >>> clf.fit(X, Y)
    BayesianRidge(alpha_1=1e-06, alpha_2=1e-06, compute_score=False, copy_X=True,
           fit_intercept=True, lambda_1=1e-06, lambda_2=1e-06, n_iter=300,
           normalize=False, tol=0.001, verbose=False)

拟合之后，模型可以用来预测新的观测结果::

    >>> clf.predict ([[1, 0.]])
    array([ 0.50000013])


模型的权重系数 :math:`w` 可以通过如下方法得到::

    >>> clf.coef_
    array([ 0.49999993,  0.49999993])

基于贝叶斯框架，这里得到的权重系数与 :ref:`ordinary_least_squares` 的
结果并不完全相同。但是贝叶斯岭回归对于ill-posed的问题更为稳定。

.. topic:: 示例:

 * :ref:`example_linear_model_plot_bayesian_ridge.py`

.. topic:: 参考:

  * More details can be found in the article `Bayesian Interpolation
    <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.27.9072&rep=rep1&type=pdf>`_
    by MacKay, David J. C.



自动相关决策  Automatic Relevance Determination - ARD
---------------------------------------

:class:`ARDRegression` 与 `贝叶斯岭回归`_ 类似，但产生更加稀疏的系数
:math:`w` [1]_  [2]_ 。 :class:`ARDRegression` 中 :math:`w` 采用另一种先验分布，不再假设高斯分
布是球对称的，而是一个椭球形的高斯分布。

因此对于每个系数 :math:`w_{i}` ，其服从高斯分布，中心为零，方差为
:math:`\lambda_{i}` ：

.. math:: p(w|\lambda) = \mathcal{N}(w|0,A^{-1})

其中 :math:`diag \; (A) = \lambda = \{\lambda_{1},...,\lambda_{p}\}`.

对比 `贝叶斯岭回归`_ ，每个系数 :math:`w_{i}`
用于其自己的方差 :math:`\lambda_i` 。所有 :math:`\lambda_i` 的先验分布
是由伽马分布的超系数 :math:`\lambda_1` 和 :math:`\lambda_2` 决定。

.. figure:: ../auto_examples/linear_model/images/plot_ard_001.png
   :target: ../auto_examples/linear_model/plot_ard.html
   :align: center
   :scale: 50%


.. topic:: 示例:

  * :ref:`example_linear_model_plot_ard.py`

.. topic:: 参考:

    .. [1] Christopher M. Bishop: Pattern Recognition and Machine Learning, Chapter 7.2.1

    .. [2] David Wipf and Srikantan Nagarajan: `A new view of automatic relevance determination. <http://books.nips.cc/papers/files/nips20/NIPS2007_0976.pdf>`_

.. _Logistic_regression:

逻辑回归 Logistic regression
===================

逻辑回归实际上是一个线性分类模型，而不是回归方法。
因此它试图降低成本函数而不是减少残差（如线性回归）。
逻辑回归也被称为 logit regression，最大熵分类 maximum-entropy classification (MaxEnt)
或者 log-linear classifier。

:class:`LogisticRegression` 类可以采用含一阶或者二阶复杂度参数的（L1，
L2 ）的逻辑回归。 一阶将产生稀疏的拟合系数。因此，
:func:`sklearn.svm.l1_min_c` 函数将计算最小的C值来避免所有特征的系数都
变为0。

.. topic:: 示例:

  * :ref:`example_linear_model_plot_logistic_l1_l2_sparsity.py`

  * :ref:`example_linear_model_plot_logistic_path.py`

.. note:: **特征选择与稀疏逻辑回归**

   因为逻辑回归配合L1可以产生稀疏的系数，所以可以用此方法选择特征。参
   见： :ref:`l1_feature_selection`.

随机梯度下降  Stochastic Gradient Descent - SGD
=================================

随机梯度下降方法是一个简单但是有效的拟合线性模型的方法。当样本或者特征
的数目很大时，这点尤为突出。
``partial_fit`` method allows only/out-of-core learning.

:class:`SGDClassifier` 与 :class:`SGDRegressor` 类提供了相应的函数对拟
合线性模型进行分类和回归。注意其可以采用不同的成本函数和惩罚函数。例如
当 ``loss="log"`` 时 :class:`SGDClassifier` 拟合逻辑回归模型，而当
``loss="hinge"`` 将拟合SVM。

.. topic:: 参考

 * :ref:`sgd`

感知器 Perceptron
==========

:class:`Perceptron` 类是另一个适用于大数据的简单地算法。在默认情况下：

    - 不需要设定学习速率

    - 不需要复杂度参数

    - 仅当出现错误时更新模型

最后一个特征说明此方法比SGD方法（loss=hinge）更为迅速而且结果更稀疏。

.. _passive_aggressive:

被动进取算法 Passive Aggressive Algorithms
=============================

被动进取算法是一系列大规模学习方法。与感知器类似，它们并不需要设定学习
速率。但是其需要一个复杂度系数 ``C`` 。

对于分类 :class:`PassiveAggressiveClassifier` 被设定为
``loss='hinge'`` (PA-I) 或者 ``loss='squared_hinge'`` (PA-II) 。对于回归
:class:`PassiveAggressiveRegressor` 可被设定为
``loss='epsilon_insensitive'`` (PA-I) 或者
``loss='squared_epsilon_insensitive'`` (PA-II).

.. topic:: 参考:


 * `"Online Passive-Aggressive Algorithms"
   <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>`_
   K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR 7 (2006)

对异常值的稳健性： RANSAC
==============================

RANSAC (RANdom SAmple Consensus)互动算法针对数据中的正常值来稳健地估计系数。

RANSAC 一非确定性方法来在一定程度上获得合理的结果。这取决于循环步骤的多少
（ `max_trials` 系数）。其通常用于线性或者非线性拟合，在摄影测绘的计
算机视觉中有广泛应用。

该算法将完整的输入数据分为还有误差的正常值和异常值（如错误的测量，或者
不合理的假设）。其结果将完全依赖于正常值。

.. figure:: ../auto_examples/linear_model/images/plot_ransac_001.png
   :target: ../auto_examples/linear_model/plot_ransac.html
   :align: center
   :scale: 50%

每一个循环进行如下：

1. 从整个样本中选择 `min_samples` 个随机样本，并检验数据完整性（参见
   `is_data_valid` ）。
2. 对该随机取样进行拟合（ `base_estimator.fit` ）并检查模型的正确性
   （参见 `is_model_valid` ）。
3. 将全部数据按照上一步的模型分为正常值和异常值
    (`base_estimator.predict(X) - y`) - 所有样本的误差小于 `residual_threshold`
    可以视作正常值。
4. 如果该模型的正常值的数目最多，则将该模型视作最佳模型。对于正常值一
    样多的模型，选择拟合度高的模型。

以上步骤循环若干次 (`max_trials`) 或者达到某个特定的停止条件
(`stop_n_inliers` 和 `stop_score`)。最后的模型为以上的最佳模型。

`is_data_valid` 和 `is_model_valid` 函数用来检测随机取样的结果是否是退
化的。 如果拟合的模型不需要识别退化情况，那么 `is_data_valid` 应该被用
来检查数据完整性，来提高计算效率。


.. topic:: 示例:

  * :ref:`example_linear_model_plot_ransac.py`

.. topic:: 参考:

 * http://en.wikipedia.org/wiki/RANSAC
 * `"Random Sample Consensus: A Paradigm for Model Fitting with Applications to
   Image Analysis and Automated Cartography"
   <http://www.cs.columbia.edu/~belhumeur/courses/compPhoto/ransac.pdf>`_
   Martin A. Fischler and Robert C. Bolles - SRI International (1981)
 * `"Performance Evaluation of RANSAC Family"
   <http://www.bmva.org/bmvc/2009/Papers/Paper355/Paper355.pdf>`_
   Sunglok Choi, Taemin Kim and Wonpil Yu - BMVC (2009)


.. _polynomial_regression:

多项式拟合 Polynomial regression
===================================================================

.. currentmodule:: sklearn.preprocessing

一个机器学习的常见模式是用线性模型来拟合非线性模型。这个思路能保持线性
模型的计算效率，而且增进其使用范围。

譬如，线性模型可以用来拟合 **非线性特征** 。 一个通常的二维线性回归模
型如下：

.. math::    \hat{y}(w, x) = w_0 + w_1 x_1 + w_2 x_2

如果我们想你和一个抛物面而不是平面，那么二阶多项式的表示如下：

.. math::    \hat{y}(w, x) = w_0 + w_1 x_1 + w_2 x_2 + w_3 x_1 x_2 + w_4 x_1^2 + w_5 x_2^2

而这有时还可以保持一个线性模型的形式。方法是采用新的变量：

.. math::  z = [x_1, x_2, x_1 x_2, x_1^2, x_2^2]

通过重新定义变量，我们的模型可以变换为：

.. math::    \hat{y}(w, x) = w_0 + w_1 z_1 + w_2 z_2 + w_3 z_3 + w_4 z_4 + w_5 z_5

我们可以得到 *多项式拟合* 与此前线性模型是同一种问题，（模型对于系数
:math:`w` 是线性的），也可以用同样的手段解决。因此通过变换到高维空间，
我们可以利用线性模型解决更多的问题。

以下这个例子，利用多项式特征来拟合一位数据：

.. figure:: ../auto_examples/linear_model/images/plot_polynomial_interpolation_001.png
   :target: ../auto_examples/linear_model/plot_polynomial_interpolation.html
   :align: center
   :scale: 50%

上例中用到了 :class:`PolynomialFeatures` 类的预处理。该类将输入的数据矩
阵转化为满足给定多项式阶数的新矩阵。其使用如下：

    >>> from sklearn.preprocessing import PolynomialFeatures
    >>> import numpy as np
    >>> X = np.arange(6).reshape(3, 2)
    >>> X
    array([[0, 1],
           [2, 3],
           [4, 5]])
    >>> poly = PolynomialFeatures(degree=2)
    >>> poly.fit_transform(X)
    array([[ 1,  0,  1,  0,  0,  1],
           [ 1,  2,  3,  4,  6,  9],
           [ 1,  4,  5, 16, 20, 25]])

输入变量 ``X`` ，:math:`[x_1, x_2]` ，被转换为
:math:`[1, x_1, x_2, x_1^2, x_1 x_2, x_2^2]` ，进而被用到任何线性模型
中。

这类变换可以通过 :ref:`Pipeline <pipeline>` 工具来进一步与拟合算法整合
为一个独立的模型，其使用方法如下：

    >>> from sklearn.preprocessing import PolynomialFeatures
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.pipeline import Pipeline
    >>> model = Pipeline([('poly', PolynomialFeatures(degree=3)),
    ...                   ('linear', LinearRegression(fit_intercept=False))])
    >>> # fit to an order-3 polynomial data
    >>> x = np.arange(5)
    >>> y = 3 - 2 * x + x ** 2 - x ** 3
    >>> model = model.fit(x[:, np.newaxis], y)
    >>> model.named_steps['linear'].coef_
    array([ 3., -2.,  1., -1.])

这个线性模型拟合多项式特征，并得出了准确的结果。

有时并不需要包含所有的高阶特征，而只需要那些交叉项 （ *interaction
features*  at :math:`d` distinct features）。这可以通过调节 :class:`PolynomialFeatures` 中的 ``interaction_only=True`` 。

例如当处理逻辑特征时， :math:`x_i^n = x_i` 因此无用。但
:math:`x_i x_j` 表达两个逻辑算子的共轭。此时我们可以通过线性分类来解决 异或 （XOR）
问题：

    >>> from sklearn.linear_model import Perceptron
    >>> from sklearn.preprocessing import PolynomialFeatures
    >>> X = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    >>> y = X[:, 0] ^ X[:, 1]
    >>> X = PolynomialFeatures(interaction_only=True).fit_transform(X)
    >>> X
    array([[1, 0, 0, 0],
           [1, 0, 1, 0],
           [1, 1, 0, 0],
           [1, 1, 1, 1]])
    >>> clf = Perceptron(fit_intercept=False, n_iter=10).fit(X, y)
    >>> clf.score(X, y)
    1.0

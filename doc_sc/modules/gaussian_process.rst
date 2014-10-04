

.. _gaussian_process:

==================
高斯过程
==================

.. currentmodule:: sklearn.gaussian_process

**机器学习中的高斯过程 (GPML)** 是一个广义的监督学习方法用以解决 *回归* 问题。其方法已经扩展到了 *概率分类* ，但是在目前的算法实现中，我们只采用 *回归* 的后处理。

机器学习中的高斯过程的优势有：

    - 预测是对样本的内插。 (至少对于正常相关性的模型）。

    - 预测的是概率 （高斯分布），所以可以去计算经验上的置信区间和超越概率评估，进而去在有兴趣的区间做进一步的拟合（如在线拟合，主动拟合）。

    - 灵活性: :ref:`现行回归模型 <linear_model>` 和 :ref:`相关性模型 <correlation_models>` 可以由用户选择。也可以使用稳定的自定义模型。

机器学习中的高斯过程的劣势有：

    - 这个不是稀疏的算法。其需要使用全部的样本，特征信息来进行预测。

    - 在高维情况中，计算效率较低。当维数高于二十左右时，其计算效率非常差。

    - 分类是藉由后期处理。首先需要解决回归问题来预测浮点型的数值 :math:`y` ，再根据其进行分类。

由于预测值得高斯属性，本方法有着广泛的应用：如全局优化，概率分类。

示例
========

回归示例简介
----------------------------------

假设我们需要寻找 :math:`g(x) = x \sin(x)` 的替代函数。那么我们需要在实验数据上拟合该函数。我们首先定义一个高斯过程模型，其回归模型和关联系数可以通过进一步的参数定义。进而我们可以将该模型对数据拟合，其中取决于参数的数目，我们可以需要通过最大似然法估计，或者采用给定参数。

.. figure:: ../auto_examples/gaussian_process/images/plot_gp_regression_001.png
   :target: ../auto_examples/gaussian_process/plot_gp_regression.html
   :align: center

::

    >>> import numpy as np
    >>> from sklearn import gaussian_process
    >>> def f(x):
    ...	    return x * np.sin(x)
    >>> X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
    >>> y = f(X).ravel()
    >>> x = np.atleast_2d(np.linspace(0, 10, 1000)).T
    >>> gp = gaussian_process.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=1e-1)
    >>> gp.fit(X, y)  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    GaussianProcess(beta0=None, corr=<function squared_exponential at 0x...>,
            normalize=True, nugget=array(2.22...-15),
            optimizer='fmin_cobyla', random_start=1, random_state=...
            regr=<function constant at 0x...>, storage_mode='full',
            theta0=array([[ 0.01]]), thetaL=array([[ 0.0001]]),
            thetaU=array([[ 0.1]]), verbose=False)
    >>> y_pred, sigma2_pred = gp.predict(x, eval_MSE=True)


拟合含噪声的数据
------------------

当拟合数据含有噪声的时候，高斯过程模型可以具体给定每个数据点的方差。 :class:`GaussianProcess` 通过参数 ``nugget`` 来添加样本数据的协方差矩阵的对角线。整体而言，这是一个 Tikhonov 规范化的过程。当协方差函数是平方指数型时，其代表的是相对方差，如下：

.. math::
   \mathrm{nugget}_i = \left[\frac{\sigma_i}{y_i}\right]^2

通过适当设置 ``nugget`` 和 ``corr`` 高斯过程可以用来可靠地从含噪声的数据中重构原始函数。

.. figure:: ../auto_examples/gaussian_process/images/plot_gp_regression_002.png
   :target: ../auto_examples/gaussian_process/plot_gp_regression.html
   :align: center

.. topic:: 其他示例

  * :ref:`example_gaussian_process_plot_gp_probabilistic_classification_after_regression.py`



数学基础
========================


基本假设
----------------------

假设我们对一个模拟实验进行建模，譬如数学函数：

.. math::

        g: & \mathbb{R}^{n_{\rm features}} \rightarrow \mathbb{R} \\
           & X \mapsto y = g(X)

GPML的起始假设是这个函数是一个高斯过程 :math:`G` 的条件取样路径，而该过程满足下面的条件

.. math::

        G(X) = f(X)^T \beta + Z(X)

其中 :math:`f(X)^T \beta` 是一个线性回归模型，而 :math:`Z(X)` 是一个期望为0的高斯过程，满足平稳的协方差：

.. math::

        C(X, X') = \sigma^2 R(|X - X'|)

:math:`\sigma^2` 是其方差， :math:`R` 是相关性函数且仅依赖于样本间的绝对距离（依赖于具体特征）。

基于上述公式我们可以看到GPML知识一个基本的最小二乘法的扩展：

.. math::

        g(X) \approx f(X)^T \beta

除了我们额外假设了样本的相关函数。的确，最小二乘法假设样本仅存在自相关，即：当且仅当 :math:`X = X'` 时 :math:`R(|X - X'|)=1` ，其余时候为0，也称谓狄拉克相关函数，有时候其在 kriging 文献中也被成为 *nugget* 相关模型。


最佳线性无偏差预测 (BLUP)
------------------------------------------

下面我们来推导基于观测的样本路径的 *最佳线性无偏差预期 best linear unbiased prediction* :math:`g` :

.. math::

    \hat{G}(X) = G(X | y_1 = g(X_1), ...,
                                y_{n_{\rm samples}} = g(X_{n_{\rm samples}}))

它的推导基于以下给定的属性：

- 其是线性的，（取样的线性组合）：

.. math::

    \hat{G}(X) \equiv a(X)^T y

- 其是无偏差的：

.. math::

    \mathbb{E}[G(X) - \hat{G}(X)] = 0

- 其是最佳的（从平均方差角度考虑）：

.. math::

    \hat{G}(X)^* = \arg \min\limits_{\hat{G}(X)} \;
                                            \mathbb{E}[(G(X) - \hat{G}(X))^2]

因此最优化的权重向量 :math:`a(X)` 是以下受限优化问题的解：

.. math::

    a(X)^* = \arg \min\limits_{a(X)} & \; \mathbb{E}[(G(X) - a(X)^T y)^2] \\
                       {\rm s. t.} & \; \mathbb{E}[G(X) - a(X)^T y] = 0

将此问题用拉格朗日方法重写，并满足一阶的优化条件，我们可以得到一个解析解。详细的推导请见参考。

最终，我们可以证明BLUP是一个高斯随机变量，其平均为：

.. math::

    \mu_{\hat{Y}}(X) = f(X)^T\,\hat{\beta} + r(X)^T\,\gamma

其方差为：

.. math::

    \sigma_{\hat{Y}}^2(X) = \sigma_{Y}^2\,
    ( 1
    - r(X)^T\,R^{-1}\,r(X)
    + u(X)^T\,(F^T\,R^{-1}\,F)^{-1}\,u(X)
    )

其中我们引入：

* 相关矩阵，每一个项都通过自相关函数定义，以及内置的参数 :math:`\theta`:

.. math::

    R_{i\,j} = R(|X_i - X_j|, \theta), \; i,\,j = 1, ..., m

* 交叉相关向量在取样以及在DOE中的预测（什么是DOE？）:

.. math::

    r_i = R(|X - X_i|, \theta), \; i = 1, ..., m

* 回归矩阵（如当 :math:`f` 是多项式的基时的 Vandermonde 矩阵）

.. math::

    F_{i\,j} = f_i(X_j), \; i = 1, ..., p, \, j = 1, ..., m

* 广义最小二乘法回归权重

.. math::

    \hat{\beta} =(F^T\,R^{-1}\,F)^{-1}\,F^T\,R^{-1}\,Y

* 向量

.. math::

    \gamma & = R^{-1}(Y - F\,\hat{\beta}) \\
    u(X) & = F^T\,R^{-1}\,r(X) - f(X)

值得注意的是高斯过程预测的概率解释是完全解析的，且基于基本的现行代数。进一步说，平均预期是可以是直接的线性组合（点积），方差则需要两个矩阵的求逆，而相关矩阵则需要进一步的Cholesky分解。


经验最佳线性无偏差预测 (EBLUP)
----------------------------------------------------

截止目前，自相关和回归模型都是假设给定的。在实际应用中，这两者是不可能提前知道的。因此我们需要依靠经验选择。参见 :ref:`correlation_models` 。

考虑到这些选择，我们需要估算在BLUP中剩余未知的参数。因此我们需要使用提供的观测样本并配合一些推测技巧。目前的程序是基于Matlab的DACE工具箱的原理，采用 *最大似然估计* 的方法。具体的方法请参见DACE的手册。这个方法是求解自相关系数的全局优化解。具体是通过 ``scipyoptimize`` 中的 ``fmin_cobyla`` 优化算法。当存在各向异性时，我们需要依靠 Welch's componentwise
optimization algorithm （见参考）。

对于机器学习中的高斯过程的理论介绍请看参考：

.. topic:: 参考：

    * `DACE, A Matlab Kriging Toolbox
      <http://www2.imm.dtu.dk/~hbn/dace/>`_ S Lophaven, HB Nielsen, J
      Sondergaard 2002


    * `Screening, predicting, and computer experiments
      <http://www.jstor.org/pss/1269548>`_ WJ Welch, RJ Buck, J Sacks,
      HP Wynn, TJ Mitchell, and MD Morris Technometrics 34(1) 15--25,
      1992


    * `Gaussian Processes for Machine Learning
      <http://www.gaussianprocess.org/gpml/chapters/RW.pdf>`_ CE
      Rasmussen, CKI Williams MIT Press, 2006 (Ed. T Diettrich)


    * `The design and analysis of computer experiments
      <http://www.stat.osu.edu/~comp_exp/book.html>`_ TJ Santner, BJ
      Williams, W Notz Springer, 2003



.. _correlation_models:

相关模型
==================

一般的相关模型与著名的支持向量机相一致。它们的基本假设是等价的。它们必须满足
Mercer条件，并且保持稳定。注意，相关模型的假设需要与原始实验的性质相一致。譬如：

* 如果原始实验是光滑的（无限可微），那么我们需要采用 *平方指数相关模型* 。
* 如果不是这样，那么我们应该采用 *指数相关模型* 。
* 如果有一个相关模型，其输入变量是可微的阶数，那么我们应该采用Matern相关模型，但是其尚未在本模块中实现 (TODO)。

关于选择相关模型的进一步讨论，请参见 Rasmussen & Williams 的书。

.. _regression_models:


回归模型
=================

通常回归模型涉及到0阶（常数），一阶和二阶多项式。但是我们可以以Python函数的形式定义，其输入变量是X，返回一个计算结果的向量。唯一的限制是函数的数目不能超过观测样本的数目，否则，回归问题将没有确定解。

实现细节
======================

目前的实现是基于对 DACE Matlab 工具箱的翻译。

.. topic:: 参考:

    * `DACE, A Matlab Kriging Toolbox
      <http://www2.imm.dtu.dk/~hbn/dace/>`_ S Lophaven, HB Nielsen, J
      Sondergaard 2002,

    * W.J. Welch, R.J. Buck, J. Sacks, H.P. Wynn, T.J. Mitchell, and M.D.
      Morris (1992). Screening, predicting, and computer experiments.
      Technometrics, 34(1) 15--25.

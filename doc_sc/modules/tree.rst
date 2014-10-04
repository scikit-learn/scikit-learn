.. _tree:

==============
决策树
==============

.. currentmodule:: sklearn.tree

**决策树 (DT)** 是一类无参数监督学习的 :ref:`分类
<tree_classification>` 和 :ref:`回归 <tree_regression>` 方法。其目标是通过从数据特征中学到的简单决策规则来构建模型进行预测。

在下例中，从数据中学习的决策树可以通过一系列“如果-那么-否则”的决策规则来近似正弦函数曲线。越深的树可以构建更复杂而又准确的模型。 

.. figure:: ../auto_examples/tree/images/plot_tree_regression_001.png
   :target: ../auto_examples/tree/plot_tree_regression.html
   :scale: 75
   :align: center

决策树的优势有：

    - 易于理解和图示。

    - 只需要少量的数据准备。其他技术往往要求数据的归一化，建立隐变量及去除空白值。注意，本模块不支持缺失数据。

    - 使用决策树的计算成本（如预测数据）与训练数据量呈指数关系。

    - 可以处理数值或者类别数据。其他技术通常只适用于一种情况。更多信息参见 :ref:`算法 <tree_algorithms>` 。

    - 可以处理多输出问题

    - 采用白盒方式（区别于黑盒）。对于模型中的一个给定情况，其解释可以很容易的从那些逻辑判断中得到，与此相对，一些黑盒方式（如神经网络）则很难阐释内部判断标准。

    - 可以通过统计测试来验证模型。因此，我们可以判断模型的可靠性。

    - 当一些假设不成立的事后，其预测能力依然出色。

决策树的劣势：

    - 决策树学习可能会产生一个过分复杂的模型，且并不能很好的刻画数据。这被称为过度拟合。一些机制，如修剪（目前不支持），设置每个节点的最小取样数，设置树的最大深度可以有效地避免这个问题。

    - 决策树可能不稳定。数据上的很小的变动可以导致完全不同的树结构。这个问题可以通过一个树集来减轻。

    - 寻求最佳树被证明是一个NP问题。所以使用的决策树学习算法常常是基于启发式算法，例如用贪婪算法来寻找局域最优。这类算法无法保证一个全局最优的决策树。这个问题可以通过训练多个决策树，每一个随机选择取样和特征来避免。

    - 决策树难于表达某些概念，譬如异或，对等，多路复合。

    - 当某些分类占主导时，决策树会产生偏差。因此在拟合之前，建议平衡数据。

.. _tree_classification:

分类
==============

:class:`DecisionTreeClassifier` 类可以对数据进行多类别分类。

如其他分类方法一样 :class:`DecisionTreeClassifier` 需要两个输入变量， 一个大小为 ``[n_samples, n_features]`` 的数列包含了训练取样的特征，另一个大小为 ``[n_samples]`` 的整数数列Y包含取样的分类::

    >>> from sklearn import tree
    >>> X = [[0, 0], [1, 1]]
    >>> Y = [0, 1]
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(X, Y)

拟合结束之后，该模型可以用来预测新的数据::

    >>> clf.predict([[2., 2.]])
    array([1])

:class:`DecisionTreeClassifier` 可以处理两类别的分类（ [-1, 1]）和多类别的分类（
[0, ..., K-1]）。

以 Iris 数据为例，我们可以构建以下的树模型::

    >>> from sklearn.datasets import load_iris
    >>> from sklearn import tree
    >>> iris = load_iris()
    >>> clf = tree.DecisionTreeClassifier()
    >>> clf = clf.fit(iris.data, iris.target)

当完成训练之后，我们可以利用函数  :func:`export_graphviz` 按照 `Graphviz <http://www.graphviz.org/>`_ 格式将树结构输出。下面是Iris数据的结果输出::

    >>> from sklearn.externals.six import StringIO
    >>> with open("iris.dot", 'w') as f:
    ...     f = tree.export_graphviz(clf, out_file=f)

之后我们可以通过 ``dot`` 工具来输出PDF文件（或者其他文件类型）： ``dot -Tpdf iris.dot -o iris.pdf`` 。

::

    >>> import os
    >>> os.unlink('iris.dot')

或者我们可以利用Python的模块 ``pydot`` 进行直接输出::

    >>> from sklearn.externals.six import StringIO  # doctest: +SKIP
    >>> import pydot # doctest: +SKIP
    >>> dot_data = StringIO() # doctest: +SKIP
    >>> tree.export_graphviz(clf, out_file=dot_data) # doctest: +SKIP
    >>> graph = pydot.graph_from_dot_data(dot_data.getvalue()) # doctest: +SKIP
    >>> graph.write_pdf("iris.pdf") # doctest: +SKIP

.. only:: html

    .. figure:: ../images/iris.svg
       :align: center

.. only:: latex

    .. figure:: ../images/iris.pdf
       :align: center

拟合结束之后，我们可以通过模型进行预测::

    >>> clf.predict(iris.data[0, :])
    array([0])

.. figure:: ../auto_examples/tree/images/plot_iris_001.png
   :target: ../auto_examples/tree/plot_iris.html
   :align: center
   :scale: 75

.. topic:: 示例:

 * :ref:`example_tree_plot_iris.py`


.. _tree_regression:

回归
==========

.. figure:: ../auto_examples/tree/images/plot_tree_regression_001.png
   :target: ../auto_examples/tree/plot_tree_regression.html
   :scale: 75
   :align: center

决策树可以通过 :class:`DecisionTreeRegressor` 类来解决回归问题。

和分类的设定一样，拟合函数需要两个输入数列X和Y。唯一的不同是Y在这里是浮点类型而不是整数类型。

    >>> from sklearn import tree
    >>> X = [[0, 0], [2, 2]]
    >>> y = [0.5, 2.5]
    >>> clf = tree.DecisionTreeRegressor()
    >>> clf = clf.fit(X, y)
    >>> clf.predict([[1, 1]])
    array([ 0.5])


.. topic:: 示例:

 * :ref:`example_tree_plot_tree_regression.py`


.. _tree_multioutput:

多输出问题
=====================

多输出问题是一个有个多个预测输出的监督学习问题，即输出Y是一个二维数列 ``[n_samples, n_outputs]`` 。

当输出结果间没有相关性的时候，一个非常简单的解决办法是构造n个独立的模型，每一个针对一个预测输出。然而由于输出预测是基于同样的输入变量，因此它们存在着相关性。因此一个更好的办法是构造一个模型来解决所有的输出预测。这有两方面好处，一是整体计算的时间降低了，因为只需要拟合一个模型，二是一般而言，预测的准确性会有所提高。

对于决策树，我们需要作如下的更改。

  - 在每一个终端节点（树叶）上存储n个输出分类结果，而不是一个
  - 决策条件要综合所有的在该节点上分类结果。

本模块在 :class:`DecisionTreeClassifier` 和 :class:`DecisionTreeRegressor` 均可以用来解决多输出问题。如果用决策树去拟合一个大小为 ``[n_samples, n_outputs]`` 的分类Y，那么输出的预测结果是：

  * 函数 ``predict`` 输出 ``n_outputs`` 个预测分类；

  * 函数 ``predict_proba`` 输出 ``n_output``  个数列来描述分类的概率。
    

在示例 :ref:`example_tree_plot_tree_regression_multioutput.py` 中我们展示如何使用多输出决策树回归。其中X是一个实数变量，而输出Y是X的正弦和余弦值。

.. figure:: ../auto_examples/tree/images/plot_tree_regression_multioutput_001.png
   :target: ../auto_examples/tree/plot_tree_regression_multioutput.html
   :scale: 75
   :align: center

在 :ref:`example_plot_multioutput_face_completion.py` 中我们展示如何使用多输出决策树分类，其中X是画面上半部分的像素，而输出是下半部分的输出。

.. figure:: ../auto_examples/images/plot_multioutput_face_completion_001.png
   :target: ../auto_examples/plot_multioutput_face_completion.html
   :scale: 75
   :align: center

.. topic:: 示例:

 * :ref:`example_tree_plot_tree_regression_multioutput.py`
 * :ref:`example_plot_multioutput_face_completion.py`

.. topic:: 参考:

 * M. Dumont et al,  `Fast multi-class image annotation with random subwindows
   and multiple output randomized trees
   <http://www.montefiore.ulg.ac.be/services/stochastic/pubs/2009/DMWG09/dumont-visapp09-shortpaper.pdf>`_, International Conference on
   Computer Vision Theory and Applications 2009

.. _tree_complexity:

计算复杂度
==========

一般而言，一个平衡二叉树的构造时间是 :math:`O(n_{samples}n_{features}\log(n_{samples}))` 而检索时间是 :math:`O(\log(n_{samples}))` 。尽管本算法试图构造一个平衡树，但往往树并不是平衡的。假设子树还是大致平衡的，那么搜索所有节点找到最佳特征的时间为 :math:`O(n_{features})` 。而在每个节点需要 :math:`O(n_{features}n_{samples}\log(n_{samples}))` 的时间。所以检索整个树的时间是 :math:`O(n_{features}n_{samples}^{2}\log(n_{samples}))` 。

Scikit-learn提供一个更有效的决策树构造算法。一个直白的算法需要重新计算每个决策点的特征分布（分类）或者平均值（回归）。通过提前将相关的样本进行排序，并记录一个分类计数，我们可以降低每个节点的计算复杂度到 :math:`O(n_{features}\log(n_{samples}))` ，这样，总的计算复杂度为 :math:`O(n_{features}n_{samples}\log(n_{samples}))` 。


实用技巧
=====================

  * 对于有很多特征的样本，决策树会过度拟合。一个合适的样本和特征的比例是重要的。一个有很少取样的高维决策树很容易过度拟合。

  * 在拟合决策树之前，对样本的重要特征进行选择会提高决策树的准确性。（参考 :ref:`PCA <PCA>` ， :ref:`ICA <ICA>` 或者 :ref:`feature_selection` ）。

  * 通过函数 ``export`` 来将决策树图像化。选择 ``max_depth=3`` 作为起始的树深度，并根据树的拟合程度逐渐增加深度。

  * 注意每增加一层树的深度，就需要增加一倍样本。通过调节 ``max_depth`` 来控制树的大小，避免过度拟合。

  * 使用 ``min_samples_split`` 或者 ``min_samples_leaf`` 来控制在每个叶节点的取样数目。一个过小的值常常意味着存在过度拟合。而一个过大的值会让决策树拟合不佳。尝试 ``min_samples_leaf=5`` 作为起始的估计。这两者的主要区别是， ``min_samples_leaf`` 保证每个叶节点的取样数目，而 ``min_samples_split`` 却会产生更小的叶节点。但是 ``min_samples_split`` 在文献中更为常见。

  * 平衡你的数据，否则会导致决策树偏向选择占多数的类别。平衡类别可以通过选择同样数目的各种类别样本，或者优化样本的权重（ ``sample_weight`` ）使得每个分类有同样的权重。另外注意，利用权重的修建参数，如 ``min_weight_fraction_leaf`` 也同样可以降低用 ``min_samples_leaf`` 情况下对多数类别的偏差。

  * 当样本是有权重的，使用 ``min_weight_fraction_leaf`` 可以保证每个叶节点都有一定的权重分配。

  * 所有的决策树内部都使用 ``np.float32`` 数列。如果训练数据不是这个格式，那么数据将拷贝。



.. _tree_algorithms:

树算法： ID3, C4.5, C5.0 和 CART
==========================================

一下来解释不同的决策树的算法，以及scikit-learn所采用的算法

ID3_ (Iterative Dichotomiser 3)在 1986 由 Ross Quinlan 提出的。这个算法采用贪婪算法来构造一个多路树，每一个节点会寻找获得最大信息收益的分类。树首先长到最大的状态，再通过修建过程来提高决策树对未知数据的概括能力。

C4.5是ID3的继承者，其通过采用数值处理为离散数据而不再要求特征是类别的数据。 C4.5将ID3的输出结果转换为一系列“如果-那么-否则”的条件规则。这些规则的准确性决定了那条规则将被首先应用。修建是通过测试去除规则的前提条件后能否提高准确性来判断的。、

C5.0是Quinlan最新的授权版本。其相比C4.5而言，使用更少的内存和构建更少的条件达到更高的准确性。

CART_ (Classification and Regression Trees) 与 C4.5 非常类似。其不同之处在于计算连续数值，但不输出决策规则。 CART通过特征和阈值来构建二叉树以提高节点的信息增量。

scikit-learn采用了CART的优化算法。

.. _ID3: http://en.wikipedia.org/wiki/ID3_algorithm
.. _CART: http://en.wikipedia.org/wiki/Predictive_analytics#Classification_and_regression_trees


.. _tree_mathematical_formulation:

数学基础
========================

给定训练向量 :math:`x_i \in R^n`, i=1,..., l 和一个类别向量
:math:`y \in R^l` ，决策树将会不断地将空间分割，知道同一类别的样本在一起。

我们用 :math:`Q` 代表在节点 :math:`m` 上的数据。每一个候选决策 :math:`\theta = (j, t_m)`  依赖于特征 :math:`j` 和阈值 :math:`t_m` 将样本分为 :math:`Q_{left}(\theta)` 和 :math:`Q_{right}(\theta)` 两个子集。

.. math::

    Q_{left}(\theta) = {(x, y) | x_j <= t_m}

    Q_{right}(\theta) = Q \setminus Q_{left}(\theta)

:math:`m` 的错误率有杂度函数 :math:`H()` 进行计算。其需要具体的问题（回归或分类）进行选择：

.. math::

   G(Q, \theta) = \frac{n_{left}}{N_m} H(Q_{left}(\theta))
   + \frac{n_{right}}{N_m} H(Q_{right}(\theta))

进而选择参数来降低错误率

.. math::

    \theta^* = \operatorname{argmin}_\theta  G(Q, \theta)

对于子集 :math:`Q_{left}(\theta^*)` 和 :math:`Q_{right}(\theta^*)` 知道达到最大的深度 :math:`N_m < \min_{samples}` 或者 :math:`N_m = 1` 。

分类标准
-----------------------

如果目标分类的类别由 0,1,...,K-1 表示，么对于属于节点 :math:`m` 有 :math:`N_m` 样本的集合 :math:`R_m` ，那么选择

.. math::

    p_{mk} = 1/ N_m \sum_{x_i \in R_m} I(y_i = k)

来表示k类别在节点 :math:`m` 中的比重。

一般的错误率是计算Gini：

.. math::

    H(X_m) = \sum_k p_{mk} (1 - p_{mk})

交叉熵：

.. math::

    H(X_m) = \sum_k p_{mk} \log(p_{mk})

或者错误分类：

.. math::

    H(X_m) = 1 - \max(p_{mk})

回归标准
-------------------

如果预测目标是连续数值，那么对于属于节点 :math:`m` 有 :math:`N_m` 样本的集合 :math:`R_m` 一般选择平均方差

.. math::

    c_m = \frac{1}{N_m} \sum_{i \in N_m} y_i

    H(X_m) = \frac{1}{N_m} \sum_{i \in N_m} (y_i - c_m)^2


.. topic:: 示例:

    * http://en.wikipedia.org/wiki/Decision_tree_learning

    * http://en.wikipedia.org/wiki/Predictive_analytics

    * L. Breiman, J. Friedman, R. Olshen, and C. Stone. Classification and
      Regression Trees. Wadsworth, Belmont, CA, 1984.

    * J.R. Quinlan. C4. 5: programs for machine learning. Morgan Kaufmann, 1993.

    * T. Hastie, R. Tibshirani and J. Friedman.
      Elements of Statistical Learning, Springer, 2009.

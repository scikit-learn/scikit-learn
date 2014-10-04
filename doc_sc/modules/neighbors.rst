.. _neighbors:

=================
最近邻算法
=================

.. sectionauthor:: Jake Vanderplas <vanderplas@astro.washington.edu>

.. currentmodule:: sklearn.neighbors

:mod:`sklearn.neighbors` 提供基于最近邻算法的非监督或监督学习方法。非监督的最近邻算法是很多学习方法的基础，例如，谱聚类（spectral clustering）和流刑学习（manifold learning）。监督的近邻学习分为两类： `classification`_ 针对数据有分离的类别和 `regression`_ 针对数据连续取值情况。

最邻近算法的原理是对于新的数据去寻找最邻近的一些训练取样，并以它们的类别作为新数据的类别。临近的训练取样数目依赖于用户的选择（k-近邻算法）或者根据局域的取样密度（半径依赖的近邻算法）。距离一般定义为标准的欧几里得距离。基于近邻算法的方法都被视为“非概括性”的机器学习方法，因为他们只是通过“记忆”所有的训练样本来作预测（往往会转换到一个更快编码的结构，如 :ref:`Ball Tree <ball_tree>` 或者 :ref:`KD Tree <kd_tree>` ）。

除却其简洁，最近邻算法还非常成功的在很多分类和回归问题中得以应用，如：数字识别，卫星图像识别。作为一个无参数方法，其常常在决策面不规则的情况下有着成功应用。

:mod:`sklearn.neighbors` 中的类可以处理 Numpy 数列或者`scipy.sparse` 稀疏矩阵。对于紧密矩阵，大部分的距离度规是支持的。对于稀疏矩阵，任意 Minkowski 度规都是支持的。

有很多机器学习的算法都是基于紧邻算法的，一个例子是 :ref:`kernel density estimation <kernel_density>` ，将在 :ref:`density estimation <density_estimation>` 中介绍。

.. _unsupervised_neighbors:

非监督最近邻算法
==============================

:class:`NearestNeighbors` 类采用非监督最紧邻算法。其作为一个统一的接口包含三个不同的算法： :class:`BallTree`, :class:`KDTree` 和一个暴力计算的方法基于 :mod:`sklearn.metrics.pairwise` 。 紧邻的搜索算法是通过参数 ``'algorithm'`` 进行选择（ ``['auto', 'ball_tree', 'kd_tree', 'brute']`` ）。作为默认值 ``'auto'`` 程序将自动选择最佳的方案。在 `Nearest Neighbor Algorithms`_ 将介绍不同方案的优劣。

    .. warning::

        对于最近邻算法，如果两个近邻 :math:`k+1` 与 :math:`k` 有着同样的距离，但是不同的分类，那么最终的分类将取决于取样的顺序。

搜索最近邻
-----------------------------
对于寻找两个数据样本间的最近邻，可以选择使用 :mod:`sklearn.neighbors` ::

    >>> from sklearn.neighbors import NearestNeighbors
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
    >>> distances, indices = nbrs.kneighbors(X)
    >>> indices                                           # doctest: +ELLIPSIS
    array([[0, 1],
           [1, 0],
           [2, 1],
           [3, 4],
           [4, 3],
           [5, 4]]...)
    >>> distances
    array([[ 0.        ,  1.        ],
           [ 0.        ,  1.        ],
           [ 0.        ,  1.41421356],
           [ 0.        ,  1.        ],
           [ 0.        ,  1.        ],
           [ 0.        ,  1.41421356]])

由于搜寻是基于训练取样，所以最近的取样就是其本身，而且距离为0。

而且我们同样可以高效地产生一个稀疏矩阵来表示最近邻的联系：

    >>> nbrs.kneighbors_graph(X).toarray()
    array([[ 1.,  1.,  0.,  0.,  0.,  0.],
           [ 1.,  1.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  1.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  1.,  1.]])

由于我们的数据在排序上相邻的取样在空间上也相邻，所以会产生这样局域对焦的矩阵形式。这种稀疏矩阵有一些列的应用，参见 :class:`sklearn.manifold.Isomap` ， :class:`sklearn.manifold.LocallyLinearEmbedding` 和 :class:`sklearn.cluster.SpectralClustering`.

KDTree 和 BallTree 类
---------------------------
另外，我们可以使用 :class:`KDTree` or :class:`BallTree` 类来直接寻找最近邻。这些是被 :class:`NearestNeighbors` 封装了的算法。球树和KD树有着同样的接口，这里我们通过KD树的例子来解释：

    >>> from sklearn.neighbors import KDTree
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> kdt = KDTree(X, leaf_size=30, metric='euclidean')
    >>> kdt.query(X, k=2, return_distance=False)          # doctest: +ELLIPSIS
    array([[0, 1],
           [1, 0],
           [2, 1],
           [3, 4],
           [4, 3],
           [5, 4]]...)

更多参数请参见 :class:`KDTree` and :class:`BallTree` 的文档，如搜索策略，变换的距离度规等。对于距离度规的选择，可以参考 :class:`DistanceMetric` 。

.. _classification:

最近邻分类
================================

最近邻分类算法是一种“基于例子的学习” 或者 “非概括性学习”。它并不尝试去构造一个统一的模型，而是储存训练样本的信息来进行预测。分类是通过计算近邻的简单多数的类别：测试样本的类别取决于最近邻的代表性类别。


scikit-learn 采用了两种不同的算法：
:class:`KNeighborsClassifier` 是基于 :math:`k` 个最近邻，其中 :math:`k` 是一个用户定义的整数。 :class:`RadiusNeighborsClassifier` 是基于一定半径之内的近邻，其中半径 :math:`r` 使用户定义的数值。

:class:`KNeighborsClassifier` 是两者中更为常用的方法。一个适当的 :math:`k` 很依赖于样本数据。一个大的 :math:`k` 可以压低噪声，但是会导致分类的边界不再清晰。

当数据不是均匀取样的，那么则应选择 :class:`RadiusNeighborsClassifier` 。用户定义一个固定的半径 :math:`r` ，当在取样低的区域，就会选取较少的近邻来完成分类。当数据维数较高是，此方法效率较低。这是由于所谓的“维数的诅咒”。

基本的紧邻方法采用统一的权重：对于待分类的数据，通过计算近邻的简单多数决定其分类。在有些情况，通过增加近邻的权重会有更好的效果。这可以通过调节参数 ``weights`` 来完成。默认的参数 ``weights = 'uniform'`` 选择统一的权重，而 ``weights = 'distance'`` 情况下权重反比于距离。用户还可以自定义权重。


.. |classification_1| image:: ../auto_examples/neighbors/images/plot_classification_001.png
   :target: ../auto_examples/neighbors/plot_classification.html
   :scale: 50

.. |classification_2| image:: ../auto_examples/neighbors/images/plot_classification_002.png
   :target: ../auto_examples/neighbors/plot_classification.html
   :scale: 50

.. centered:: |classification_1| |classification_2|

.. topic:: 示例:

  * :ref:`example_neighbors_plot_classification.py`

.. _regression:

最近邻回归
============================

最近邻回归被用在取样数据的分类时连续值的情况。预测值将取决于近邻的平均。

scikit-learn 采用两种不同的近邻算法：
:class:`KNeighborsRegressor` 是基于 :math:`k` 个近邻，其中 :math:`k` 是用户定义的整数。而 :class:`RadiusNeighborsRegressor` 是基于固定半径 :math:`r` 。

最基本的最近邻回归的权重算法与最近邻分类算法的权重一致。

.. figure:: ../auto_examples/neighbors/images/plot_regression_001.png
   :target: ../auto_examples/neighbors/plot_regression.html
   :align: center
   :scale: 75

下例展示了多输出的最紧邻回归算法： :ref:`example_plot_multioutput_face_completion.py` 。在这个例子中，输入X是脸部照片的上半部分的像素，而输出Y是下半部分的脸部像素。

.. figure:: ../auto_examples/images/plot_multioutput_face_completion_001.png
   :target: ../auto_examples/plot_multioutput_face_completion.html
   :scale: 75
   :align: center


.. topic:: 示例:

  * :ref:`example_neighbors_plot_regression.py`

  * :ref:`example_plot_multioutput_face_completion.py`


最近邻算法
===========================

.. _brute_force:

暴力计算
-----------

最近邻的快速求解一直是机器学习的热门研究领域。最直白的近邻搜索是通过直接计算所有数据对间的距离完成的。对于一个有 :math:`N` 个 :math:`D` 维样本的问题，其计算复杂度为 :math:`O[D N^2]` 。对于小样本问题，暴力计算还有具有优势的。然而当样布扥数目 :math:`N` 变大，很快此方法不再适用。
:mod:`sklearn.neighbors` 中通过 ``algorithm = 'brute'`` 来采用暴力计算。具体的函数为 :mod:`sklearn.metrics.pairwise`.

.. _kd_tree:

K-D 树
--------

为了提高计算效率，一系列基于树结构的算法被开发出来。大致上，这些算法通过树结构来降低计算间距的信息的任务。其基本想法是如果取样 :math:`A` 距离取样 :math:`B` 非常远，而 :math:`B` 与取样 :math:`C` 非常近，那么我们知道 :math:`A` 与 :math:`C` 的距离非常远，因此 *不再需要计算它们间的距离* 。这样我们就将计算效率从暴力计算的 :math:`N` 就提高了 :math:`O[D N \log(N)]` 。

*KD 树* （k维树的缩写）是一个利用这个信息的早期算法。它将二维四叉树和三维八叉树的想法进一步扩展到高维。KD树是一个二叉树，不断将系数空间沿数据维度分割，并将取样植入嵌套的空间之中。KD树的构造非常迅速，由于分割是基于数据轴的方向，不需要计算 :math:`D` 维的距离。一旦构造完成，一个数据的近邻可以通过 :math:`O[\log(N)]` 个距离计算得到。虽然KD树方法对于低维数据非常迅速（ :math:`D < 20` ），但对于高维并不有效。这也是“维度的诅咒”所导致的。在scikit-learn中KD树最近邻方法可以通过参数 ``algorithm = 'kd_tree'`` 选择，并通过 :class:`KDTree` 类进行计算。


.. topic:: 参考:

   * `"Multidimensional binary search trees used for associative searching"
     <http://dl.acm.org/citation.cfm?doid=361002.361007>`_,
     Bentley, J.L., Communications of the ACM (1975)


.. _ball_tree:

球树
---------

为了应对KD树在高维时的低效，*球树* 结构被发展出来。对比KD树沿笛卡尔坐标轴进行分割，球树在一系列嵌套的超球面上进行数据分割。这样使得树结构的构造成本增加，但是在高维情况下更为有效。

一个球树是将数据分入节点，每个节点的中心为 :math:`C` ，半径为 :math:`r` 。一个节点的潜在近邻可以通过三角不等式判定：

.. math::   |x+y| \leq |x| + |y|

通过这样的设定，测试数据点与中心点的距离可以用来估计与所有在该节点内的取样的距离。由于球树节点的球对称结构，其可以在高位情况下比KD数更有效率，虽然具体的表现很依赖于训练样本的结构。在scikit-learn中球树最近邻方法可以通过参数 ``algorithm = 'ball_tree'`` 来选择，并通过 :class:`sklearn.neighbors.BallTree` 计算。或者直接通过 :class:`BallTree` 类计算。

.. topic:: 参考:

   * `"Five balltree construction algorithms"
     <http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.91.8209>`_,
     Omohundro, S.M., International Computer Science Institute
     Technical Report (1989)

最近邻方法选择
-------------------------------------
最佳的算法是取决于具体的数据，以下提供几个考虑因素：

* 样本数目 :math:`N` （ ``n_samples`` ）和维数  :math:`D` （ ``n_features`` ）。

  * *暴力计算* 搜索时间正比于 :math:`O[D N]`
  * *球树* 搜索时间正比于 :math:`O[D \log(N)]`
  * *KD 树* 搜索时间正比于 :math:`D` ，但相对比较难准确估算。对于低维情况 :math:`D<20`  耗时为 :math:`O[D\log(N)]` 而搜索很有效。对于高维情况，耗时增加到 :math:`O[DN]` 因此构造树的花销会使得暴力计算更为迅速。

  对于比较小的样本（ :math:`N<30` ） :math:`\log(N)` 与 :math:`N` 差别不大，因此暴力计算可以更为有效。 :class:`KDTree` 和 :class:`BallTree` 对此提供了一个 *leaf size* 参数：当样本数目较小时，搜索将切换到暴力计算。这样使得这两种算法在小样本时有着较高的效率。

* 数据结构： 数据的 *本征维度* 和/或 *稀疏性* 。 本征维度是指数据所在的流形的维度 :math:`d \le D` ，其可以是线性或者非线性的嵌入到系数空间中。稀疏性是指数据填满系数空间的程度（注意区别稀疏矩阵的概念，样本数据的矩阵可以有非零项，但其对应的 **结构** 可以仍然是稀疏的）。

  * *暴力计算* 搜寻时间不因结构的不同而改变。
  * *球树* 和 *KD 树* 搜寻时间会显著的依赖于数据结构。大致上讲，稀疏的数据，低的本征维度会让搜索更快。由于KD树的内部表示是沿着系数的数轴，对于任意的数据结构，它不会有像球树那样显著的改善。、

  由于机器学习的数据往往有真很好的结构，所以很适合基于树结构的搜索。

* 寻找近邻的数目 :math:`k` 。

  * *暴力计算* 搜寻时间与 :math:`k` 基本无关。
  * *球树* 和 *KD 数* 的搜索时间会随着 :math:`k` 的增加而变慢。这是由两个因素导致的：1. 一个更大的 :math:`k` 会要求搜索更大的系数空间；2. 当 :math:`k > 1` 时，需要内部检索穿过的树。

  当 :math:`k` 与 :math:`N` 的大小相当时，树结构简化节点的优势被削弱。在这个情况下，暴力计算会更加有效。

* 需要搜寻的次数。球树和KD树都需要构造过程。当搜寻次数变多时，这部分的成本将变得不再重要。如果仅有少量数据需要搜索，那么构造过程将占用很大的时间。这时候，暴力计算将是个更好的选择。

目前，当 :math:`k < N/2` 且 ``'effective_metric_'`` 存在在 ``'kd_tree'`` 的 ``'VALID_METRICS'`` 列表中的时候， ``algorithm = 'auto'`` 会选择 ``'kd_tree'`` 。当 :math:`k < N/2` 且 ``'effective_metric_'`` 不存在在 ``'kd_tree'`` 的 ``'VALID_METRICS'`` 列表中的时候， ``algorithm = 'auto'`` 会选择 ``'ball_tree'`` 。而当 :math:`k >= N/2` 的时候，会选择 ``'brute'`` 。此处的选择是基于待搜索的数据数目与训练样本的数目相若，而且 ``leaf_size`` 是在默认值 ``30`` 。

``leaf_size`` 的作用
-----------------------

如前所述，当样本较小的时候，暴力计算的方法比树结构的搜索更为有效。这点在树结构中也有应用：在每个节点中的数据样本是通过暴力计算来进行搜索的。而这个节点内的数据数目由 ``leaf_size`` 来确定。其选择有着如下考量：

**构造时间**
  一个更大的 ``leaf_size`` 会减少需要构造的节点数目，因此构造时间更短。

**搜索时间**
  一个过大或者过小的 ``leaf_size`` 都会导致更长的搜索时间。当 ``leaf_size`` 接近于1的时候，遍历节点的时间会显著增加搜索时间。而当 ``leaf_size`` 接近样本数目的时候，搜索将变为暴力计算。因此默认的 ``leaf_size = 30`` 在这两者间的一个平衡。

**内存**
  当 ``leaf_size`` 增加，树结构占用的内存将减少。这点对于球树尤为突出。其需要存储 :math:`D` 维每个节点。这相当于需要存储 ``1 / leaf_size`` 乘以样本大小。

暴力计算中并不需要 ``leaf_size``  。


最近邻中心分类法
===========================

:class:`NearestCentroid` 分类法是一个简单的算法，其每个分类依赖于在中心的成员类别。实际上这相当于 :class:`sklearn.KMeans` 算法中更新分类的步骤。其也不需要任何参数，这也是其作为分类方法的一个基本特质。然而它却在非凸类别中表现不佳，或者当不同类别的方差差异较大的时候。对与样本类别方差不一致的情形，参见线性判别分析 (:class:`sklearn.lda.LDA`) 和二次判别分析 (:class:`sklearn.qda.QDA`) 。默认的 :class:`NearestCentroid` 使用非常简单：

    >>> from sklearn.neighbors.nearest_centroid import NearestCentroid
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = NearestCentroid()
    >>> clf.fit(X, y)
    NearestCentroid(metric='euclidean', shrink_threshold=None)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]


最近邻缩小中心
-------------------------

:class:`NearestCentroid` 分类中有一个参数 ``shrink_threshold`` 来应用缩小中心分类法。在应用中每一个中心的特征将被除以该分类在其中的方差，并进一步除以 ``shrink_threshold`` 。注意，当一个特征为0时，在上述变换下仍然为0。这个方法是在降低特征的贡献，但有时会很有效果，例如降低噪声。

在下面的例子中，通过一个很小的收缩值，我们将准确度从 0.81 提高到了 0.82。

.. |nearest_centroid_1| image:: ../auto_examples/neighbors/images/plot_nearest_centroid_001.png
   :target: ../auto_examples/neighbors/plot_classification.html
   :scale: 50

.. |nearest_centroid_2| image:: ../auto_examples/neighbors/images/plot_nearest_centroid_002.png
   :target: ../auto_examples/neighbors/plot_classification.html
   :scale: 50

.. centered:: |nearest_centroid_1| |nearest_centroid_2|

.. topic:: 示例:

  * :ref:`example_neighbors_plot_nearest_centroid.py`

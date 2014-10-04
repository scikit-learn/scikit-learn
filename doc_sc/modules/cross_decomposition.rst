.. _cross_decomposition:

===================
交叉分解
===================

.. currentmodule:: sklearn.cross_decomposition

交叉分解模块包含两类算法：部分最小二乘法（partial least squares PLS） 和正则相关分析（canonical correlation analysis CCA）。

这两个算法都是用来寻找两个多变量样本数据间的线性关系:  ``fit`` 函数的输入变量 ``X`` 和 ``Y`` 都是二维向量。

.. figure:: ../auto_examples/cross_decomposition/images/plot_compare_cross_decomposition_001.png
   :target: ../auto_examples/cross_decomposition/plot_compare_cross_decomposition.html
   :scale: 75%
   :align: center


交叉分解算法尝试寻找两个矩阵(X and Y)间的基本关系。他们是通过隐变量来刻画其间的协变关系。他们试图在X空间中寻找一个高维的方向来最大化的解释Y在沿该方向上的方差。PLS回归适用于样本的特征多于样本数目的时候，而标准的回归方法将会失效。

在本模块中包含 :class:`PLSRegression` ， :class:`PLSCanonical` ， :class:`CCA` 和 :class:`PLSSVD` 。


.. topic:: 参考：

   * JA Wegelin
     `A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block case <https://www.stat.washington.edu/research/reports/2000/tr371.pdf>`_

.. topic:: 文献：

    * :ref:`example_cross_decomposition_plot_compare_cross_decomposition.py`

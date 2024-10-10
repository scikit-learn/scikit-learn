Version 1.6 (2024-10-10)
========================





Dropping support for building with setuptools
---------------------------------------------

From scikit-learn 1.6 onwards, support for building with setuptools has been
removed. Meson is the only supported way to build scikit-learn, see
:ref:`Building from source <install_bleeding_edge>` for more details.

:pr:`29400` by :user:`Loïc Estève <lesteve>` #29634

sklearn.linear_model
--------------------

- |Fix| :class:`linear_model.LogisticRegressionCV` corrects sample weight handling
  for the calculation of test scores.
  :pr:`29419` by :user:`Shruti Nath <snath-xoc>`. #29419

- |Fix| :class:`linear_model.LassoCV` and :class:`linear_model.ElasticNetCV` now
  take sample weights into accounts to define the search grid for the internally tuned
  `alpha` hyper-parameter. :pr:`29442` by :user:`John Hopfensperger <s-banach> and
  :user:`Shruti Nath <snath-xoc>`. #29442

- |API| Deprecates `copy_X` in :class:`linear_model.TheilSenRegressor` as the parameter
  has no effect. `copy_X` will be removed in 1.8.
  :pr:`29105` by :user:`Adam Li <adam2392>`. #29105

sklearn.tree
------------

- |Feature| :class:`tree.ExtraTreeClassifier` and :class:`tree.ExtraTreeRegressor` now
  support missing-values in the data matrix ``X``. Missing-values are handled by
  randomly moving all of the samples to the left, or right child node as the tree is
  traversed.
  :pr:`27966` by :user:`Adam Li <adam2392>`. #27966

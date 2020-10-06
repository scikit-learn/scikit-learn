.. Places parent toc into the sidebar

:parenttoc: True

.. include:: includes/big_toc_css.rst

.. _visualizations:

==============
Visualizations
==============

Scikit-learn defines a simple API for creating visualizations for machine
learning. The key feature of this API is to allow for quick plotting and
visual adjustments without recalculation. In the following example, we plot a
ROC curve for a fitted support vector machine:

.. code-block:: python

    from sklearn.model_selection import train_test_split
    from sklearn.svm import SVC
    from sklearn.metrics import plot_roc_curve
    from sklearn.datasets import load_wine

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    svc = SVC(random_state=42)
    svc.fit(X_train, y_train)

    svc_disp = plot_roc_curve(svc, X_test, y_test)

.. figure:: auto_examples/miscellaneous/images/sphx_glr_plot_roc_curve_visualization_api_001.png
    :target: auto_examples/miscellaneous/plot_roc_curve_visualization_api.html
    :align: center
    :scale: 75%

The returned `svc_disp` object allows us to continue using the already computed
ROC curve for SVC in future plots. In this case, the `svc_disp` is a
:class:`~sklearn.metrics.RocCurveDisplay` that stores the computed values as
attributes called `roc_auc`, `fpr`, and `tpr`. Next, we train a random forest
classifier and plot the previously computed roc curve again by using the `plot`
method of the `Display` object.

.. code-block:: python

    import matplotlib.pyplot as plt
    from sklearn.ensemble import RandomForestClassifier

    rfc = RandomForestClassifier(random_state=42)
    rfc.fit(X_train, y_train)

    ax = plt.gca()
    rfc_disp = plot_roc_curve(rfc, X_test, y_test, ax=ax, alpha=0.8)
    svc_disp.plot(ax=ax, alpha=0.8)

.. figure:: auto_examples/miscellaneous/images/sphx_glr_plot_roc_curve_visualization_api_002.png
    :target: auto_examples/miscellaneous/plot_roc_curve_visualization_api.html
    :align: center
    :scale: 75%

Notice that we pass `alpha=0.8` to the plot functions to adjust the alpha
values of the curves.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_miscellaneous_plot_roc_curve_visualization_api.py`
    * :ref:`sphx_glr_auto_examples_miscellaneous_plot_partial_dependence_visualization_api.py`
    * :ref:`sphx_glr_auto_examples_miscellaneous_plot_display_object_visualization.py`

Available Plotting Utilities
============================

Functions
---------

.. currentmodule:: sklearn

.. autosummary::

   inspection.plot_partial_dependence
   metrics.plot_confusion_matrix
   metrics.plot_det_curve
   metrics.plot_precision_recall_curve
   metrics.plot_roc_curve


Display Objects
---------------

.. currentmodule:: sklearn

.. autosummary::

   inspection.PartialDependenceDisplay
   metrics.ConfusionMatrixDisplay
   metrics.DetCurveDisplay
   metrics.PrecisionRecallDisplay
   metrics.RocCurveDisplay

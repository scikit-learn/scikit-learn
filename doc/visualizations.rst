.. _visualizations:

==============
Visualizations
==============

Scikit-learn defines a simple API for creating visualizations for machine
learning. The key feature of this API is to allow for quick plotting and
visual adjustments without recalculation. We provide `Display` classes that
expose two methods for creating plots: `from_estimator` and
`from_predictions`.

The `from_estimator` method generates a `Display` object from a fitted estimator and
input data (`X`, `y`).
The `from_predictions` method creates a `Display` object from true and predicted
values (`y_test`, `y_pred`), which
is useful when you only want to compute the predictions once. Using `from_predictions`
avoids us to recompute the predictions, but does not automatically resolve some
ambiguities.

The `Display` object stores the computed values (e. g. metric values) required for
plotting with Matplotlib. These computed values are the results of some derivatives
after we pass the raw predictions to `from_predictions`, or we get them from
an estimator via `from_estimator`.

Additionally, the plot method allows adding to an existing plot by passing the existing
plots :class:`matplotlib.axes.Axes` to the `ax` parameter.

In the following example, we plot a ROC curve for a fitted support
vector machine using `from_estimator`:

.. plot::
   :context: close-figs
   :align: center

    from sklearn.model_selection import train_test_split
    from sklearn.svm import SVC
    from sklearn.metrics import RocCurveDisplay
    from sklearn.datasets import load_wine

    X, y = load_wine(return_X_y=True)
    y = y == 2  # make binary
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    svc = SVC(random_state=42)
    svc.fit(X_train, y_train)

    svc_disp = RocCurveDisplay.from_estimator(svc, X_test, y_test)

If you already have the prediction values, you could instead use
`from_predictions` to do the same thing:


.. plot::
   :context: close-figs
   :align: center

   from sklearn.model_selection import train_test_split
   from sklearn.svm import SVC
   from sklearn.metrics import RocCurveDisplay
   from sklearn.datasets import load_wine

   X, y = load_wine(return_X_y=True)
   y = y == 2  # make binary
   X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
   svc = SVC(random_state=42).fit(X_train, y_train)
   y_pred = svc.decision_function(X_test)

   svc_disp = RocCurveDisplay.from_predictions(y_test, y_pred)


The returned `svc_disp` object allows us to continue using the already computed
ROC curve for SVC in future plots. In this case, the `svc_disp` is a
:class:`~sklearn.metrics.RocCurveDisplay` that stores the computed values as
attributes called `roc_auc`, `fpr`, and `tpr`.

Next, we train a random forest classifier and plot
the previously computed ROC curve again by using the `plot` method of the
`Display` object.

.. plot::
   :context: close-figs
   :align: center

    import matplotlib.pyplot as plt
    from sklearn.ensemble import RandomForestClassifier

    rfc = RandomForestClassifier(n_estimators=10, random_state=42)
    rfc.fit(X_train, y_train)

    ax = plt.gca()
    rfc_disp = RocCurveDisplay.from_estimator(rfc, X_test, y_test, ax=ax, alpha=0.8)
    svc_disp.plot(ax=ax, alpha=0.8)

Notice that we pass `alpha=0.8` to the plot functions to adjust the alpha
values of the curves.


.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_miscellaneous_plot_roc_curve_visualization_api.py`
* :ref:`sphx_glr_auto_examples_miscellaneous_plot_partial_dependence_visualization_api.py`
* :ref:`sphx_glr_auto_examples_miscellaneous_plot_display_object_visualization.py`
* :ref:`sphx_glr_auto_examples_calibration_plot_compare_calibration.py`

Available Plotting Utilities
============================

Display Objects
---------------

.. currentmodule:: sklearn

.. autosummary::

   calibration.CalibrationDisplay
   inspection.PartialDependenceDisplay
   inspection.DecisionBoundaryDisplay
   metrics.ConfusionMatrixDisplay
   metrics.DetCurveDisplay
   metrics.PrecisionRecallDisplay
   metrics.PredictionErrorDisplay
   metrics.RocCurveDisplay
   model_selection.LearningCurveDisplay
   model_selection.ValidationCurveDisplay

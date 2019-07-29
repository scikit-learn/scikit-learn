.. include:: includes/big_toc_css.rst

.. _visualizers:

Visualizers
-----------

Scikit-learn defines a simple API for creating visualizations for machine
learning. The key features of this API is to allow for quick plotting and
visual adjustments without recalculation. In the following example, we plot a
ROC curve for a fitted support vector machine:

.. code-block:: python

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    svc = SVC(random_state=42)
    svc.fit(X_train, y_train)

    viz_svc = plot_roc_curve(svc, X_test, y_test)

The returned `viz_svc` object allows us to continue using the already computed
ROC curve for SVC in future plots. In this case, the `viz` is a
`RocCurveVisualizer` that stores the computed values as attributes called
`auc_roc`, `fpr`, and `tpr`. Next, we train a random forest classifier and plot
the roc curve on the same axes.

.. code-block:: python

    rfc = DecisionTreeClassifier(random_state=42)
    rfc.fit(X_train, y_train)

    fig, ax = plt.subplots()
    viz_rfc = plot_roc_curve(rfc, X_test, y_test, ax=ax)
    viz_svc.plot(ax=ax)

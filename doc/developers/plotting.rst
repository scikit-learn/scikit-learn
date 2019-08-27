.. _plotting_api:

================================
Developing with the Plotting API
================================

Scikit-learn defines a simple API for creating visualizations for machine
learning. The key features of this API is to run calculations once and to have
the flexibility to adjust the visualizations after the fact.

Plotting API Overview
---------------------

This logic is encapsulated into a display object where the computed data is
stored and the plotting is done in a `plot` method. The display object's
`__init__` method contains only the data needed to create the visualization.
The `plot` method takes in parameters that only have to do with visualization,
such as a matplotlib axes. The `plot` method will store the matplotlib artists
as attributes allowing for style adjustments through the display object. A
`plot_*` helper function accepts parameters to do the computation and the
parameters used for plotting. After the helper function creates the display
object with the computed values, it calls the display's plot method. Note that
the `plot` method defines attributes related to matplotlib, such as the line
artist. This allows for customizations after calling the `plot` method.

For example, the `RocCurveDisplay` defines the following methods and
attributes:

.. code-block:: python

   class RocCurveDisplay:
       def __init__(self, fpr, tpr, roc_auc, estimator_name):
           ...
           self.fpr = fpr
           self.tpr = tpr
           self.roc_auc = roc_auc
           self.estimator_name = estimator_name

       def plot(self, ax=None, name=None, **kwargs):
           ...
           self.line_ = ...
           self.ax_ = ax
           self.figure_ = ax.figure_

   def plot_roc_curve(estimator, X, y, pos_label=None, sample_weight=None,
                      drop_intermediate=True, response_method="auto",
                      name=None, ax=None, **kwargs):
       # do computation
       viz = RocCurveDisplay(fpr, tpr, roc_auc, 
                                estimator.__class__.__name__)
       return viz.plot(ax=ax, name=name, **kwargs)

Plotting with Multiple Axes
---------------------------

The visualization API handles plotting multiple axes in two ways. The Display
object's `plot` method and the `plot_*` helper function can take a list of axes
or a single axes. If a list of axes is provided, `plot` will check if the the
number of axes is consistent with the number of axes it needs and then draw on
those axes. 

When a single axes is passed in, that axes defines a space for the multiple
axes to be placed. In this case, matplotlib's
`gridspec.GridSpecFromSubplotSpec` can be used to split up the space::

   import matplotlib.pyplot as plt
   from matplotlib.gridspec import GridSpecFromSubplotSpec

   fig, ax = plt.subplots()
   gs = GridSpecFromSubplotSpec(2, 2, subplot_spec=ax.get_subplotspec())

   ax_top_left = fig.add_subplot(gs[0, 0])
   ax_top_right = fig.add_subplot(gs[0, 1])
   ax_bottom = fig.add_subplot(gs[1, :])

By default, the `ax` keyworld in `plot` is `None`. In this case, the single
axes is created and the gridspec api is used to create the regions to plot in.

For example, :func:`sklearn.inspection.plot_partial_dependence` plots multiple
lines and contours using this API. The axes that is passed in or created that
defines the space is saved as a `bounding_ax_` attribute. The individual axes
created are stored in a `axes_` ndarray, corresponding to the axes position on
the grid. Positions that are not used are set to `None`. Furthermore, the
matplotlib Artists are stored in `lines_` and `contours_` where the key is the
position on the grid. When a list of axes is passsed in, the `axes_`, `lines_`,
and `contours_` keys is single int corresponding to the position on the passed
in list of axes. 

Read more in :ref:`sphx_glr_auto_examples_plot_roc_curve_visualization_api.py`
and the :ref:`User Guide <visualizations>`.

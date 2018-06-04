from __future__ import division
from itertools import product
from distutils.version import LooseVersion
import warnings
from textwrap import dedent

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

from . import utils
from .palettes import color_palette, blend_palette
from .external.six import string_types
from .distributions import distplot, kdeplot,  _freedman_diaconis_bins


__all__ = ["FacetGrid", "PairGrid", "JointGrid", "pairplot", "jointplot"]


class Grid(object):
    """Base class for grids of subplots."""
    _margin_titles = False
    _legend_out = True

    def set(self, **kwargs):
        """Set attributes on each subplot Axes."""
        for ax in self.axes.flat:
            ax.set(**kwargs)

        return self

    def savefig(self, *args, **kwargs):
        """Save the figure."""
        kwargs = kwargs.copy()
        kwargs.setdefault("bbox_inches", "tight")
        self.fig.savefig(*args, **kwargs)

    def add_legend(self, legend_data=None, title=None, label_order=None,
                   **kwargs):
        """Draw a legend, maybe placing it outside axes and resizing the figure.

        Parameters
        ----------
        legend_data : dict, optional
            Dictionary mapping label names to matplotlib artist handles. The
            default reads from ``self._legend_data``.
        title : string, optional
            Title for the legend. The default reads from ``self._hue_var``.
        label_order : list of labels, optional
            The order that the legend entries should appear in. The default
            reads from ``self.hue_names`` or sorts the keys in ``legend_data``.
        kwargs : key, value pairings
            Other keyword arguments are passed to the underlying legend methods
            on the Figure or Axes object.

        Returns
        -------
        self : Grid instance
            Returns self for easy chaining.

        """
        # Find the data for the legend
        legend_data = self._legend_data if legend_data is None else legend_data
        if label_order is None:
            if self.hue_names is None:
                label_order = np.sort(list(legend_data.keys()))
            else:
                label_order = list(map(utils.to_utf8, self.hue_names))

        blank_handle = mpl.patches.Patch(alpha=0, linewidth=0)
        handles = [legend_data.get(l, blank_handle) for l in label_order]
        title = self._hue_var if title is None else title
        try:
            title_size = mpl.rcParams["axes.labelsize"] * .85
        except TypeError:  # labelsize is something like "large"
            title_size = mpl.rcParams["axes.labelsize"]

        # Set default legend kwargs
        kwargs.setdefault("scatterpoints", 1)

        if self._legend_out:
            # Draw a full-figure legend outside the grid
            figlegend = self.fig.legend(handles, label_order, "center right",
                                        **kwargs)
            self._legend = figlegend
            figlegend.set_title(title)

            # Set the title size a roundabout way to maintain
            # compatability with matplotlib 1.1
            prop = mpl.font_manager.FontProperties(size=title_size)
            figlegend._legend_title_box._text.set_font_properties(prop)

            # Draw the plot to set the bounding boxes correctly
            self.fig.draw(self.fig.canvas.get_renderer())

            # Calculate and set the new width of the figure so the legend fits
            legend_width = figlegend.get_window_extent().width / self.fig.dpi
            figure_width = self.fig.get_figwidth()
            self.fig.set_figwidth(figure_width + legend_width)

            # Draw the plot again to get the new transformations
            self.fig.draw(self.fig.canvas.get_renderer())

            # Now calculate how much space we need on the right side
            legend_width = figlegend.get_window_extent().width / self.fig.dpi
            space_needed = legend_width / (figure_width + legend_width)
            margin = .04 if self._margin_titles else .01
            self._space_needed = margin + space_needed
            right = 1 - self._space_needed

            # Place the subplot axes to give space for the legend
            self.fig.subplots_adjust(right=right)

        else:
            # Draw a legend in the first axis
            ax = self.axes.flat[0]
            leg = ax.legend(handles, label_order, loc="best", **kwargs)
            leg.set_title(title)

            # Set the title size a roundabout way to maintain
            # compatability with matplotlib 1.1
            prop = mpl.font_manager.FontProperties(size=title_size)
            leg._legend_title_box._text.set_font_properties(prop)

        return self

    def _clean_axis(self, ax):
        """Turn off axis labels and legend."""
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.legend_ = None
        return self

    def _update_legend_data(self, ax):
        """Extract the legend data from an axes object and save it."""
        handles, labels = ax.get_legend_handles_labels()
        data = {l: h for h, l in zip(handles, labels)}
        self._legend_data.update(data)

    def _get_palette(self, data, hue, hue_order, palette):
        """Get a list of colors for the hue variable."""
        if hue is None:
            palette = color_palette(n_colors=1)

        else:
            hue_names = utils.categorical_order(data[hue], hue_order)
            n_colors = len(hue_names)

            # By default use either the current color palette or HUSL
            if palette is None:
                current_palette = utils.get_color_cycle()
                if n_colors > len(current_palette):
                    colors = color_palette("husl", n_colors)
                else:
                    colors = color_palette(n_colors=n_colors)

            # Allow for palette to map from hue variable names
            elif isinstance(palette, dict):
                color_names = [palette[h] for h in hue_names]
                colors = color_palette(color_names, n_colors)

            # Otherwise act as if we just got a list of colors
            else:
                colors = color_palette(palette, n_colors)

            palette = color_palette(colors, n_colors)

        return palette


_facet_docs = dict(

    data=dedent("""\
    data : DataFrame
        Tidy ("long-form") dataframe where each column is a variable and each
        row is an observation.\
    """),
    col_wrap=dedent("""\
    col_wrap : int, optional
        "Wrap" the column variable at this width, so that the column facets
        span multiple rows. Incompatible with a ``row`` facet.\
    """),
    share_xy=dedent("""\
    share{x,y} : bool, optional
        If true, the facets will share y axes across columns and/or x axes
        across rows.\
    """),
    size=dedent("""\
    size : scalar, optional
        Height (in inches) of each facet. See also: ``aspect``.\
    """),
    aspect=dedent("""\
    aspect : scalar, optional
        Aspect ratio of each facet, so that ``aspect * size`` gives the width
        of each facet in inches.\
    """),
    palette=dedent("""\
    palette : palette name, list, or dict, optional
        Colors to use for the different levels of the ``hue`` variable. Should
        be something that can be interpreted by :func:`color_palette`, or a
        dictionary mapping hue levels to matplotlib colors.\
    """),
    legend_out=dedent("""\
    legend_out : bool, optional
        If ``True``, the figure size will be extended, and the legend will be
        drawn outside the plot on the center right.\
    """),
    margin_titles=dedent("""\
    margin_titles : bool, optional
        If ``True``, the titles for the row variable are drawn to the right of
        the last column. This option is experimental and may not work in all
        cases.\
    """),
    )


class FacetGrid(Grid):
    """Subplot grid for plotting conditional relationships."""

    def __init__(self, data, row=None, col=None, hue=None, col_wrap=None,
                 sharex=True, sharey=True, size=3, aspect=1, palette=None,
                 row_order=None, col_order=None, hue_order=None, hue_kws=None,
                 dropna=True, legend_out=True, despine=True,
                 margin_titles=False, xlim=None, ylim=None, subplot_kws=None,
                 gridspec_kws=None):

        MPL_GRIDSPEC_VERSION = LooseVersion('1.4')
        OLD_MPL = LooseVersion(mpl.__version__) < MPL_GRIDSPEC_VERSION

        # Determine the hue facet layer information
        hue_var = hue
        if hue is None:
            hue_names = None
        else:
            hue_names = utils.categorical_order(data[hue], hue_order)

        colors = self._get_palette(data, hue, hue_order, palette)

        # Set up the lists of names for the row and column facet variables
        if row is None:
            row_names = []
        else:
            row_names = utils.categorical_order(data[row], row_order)

        if col is None:
            col_names = []
        else:
            col_names = utils.categorical_order(data[col], col_order)

        # Additional dict of kwarg -> list of values for mapping the hue var
        hue_kws = hue_kws if hue_kws is not None else {}

        # Make a boolean mask that is True anywhere there is an NA
        # value in one of the faceting variables, but only if dropna is True
        none_na = np.zeros(len(data), np.bool)
        if dropna:
            row_na = none_na if row is None else data[row].isnull()
            col_na = none_na if col is None else data[col].isnull()
            hue_na = none_na if hue is None else data[hue].isnull()
            not_na = ~(row_na | col_na | hue_na)
        else:
            not_na = ~none_na

        # Compute the grid shape
        ncol = 1 if col is None else len(col_names)
        nrow = 1 if row is None else len(row_names)
        self._n_facets = ncol * nrow

        self._col_wrap = col_wrap
        if col_wrap is not None:
            if row is not None:
                err = "Cannot use `row` and `col_wrap` together."
                raise ValueError(err)
            ncol = col_wrap
            nrow = int(np.ceil(len(col_names) / col_wrap))
        self._ncol = ncol
        self._nrow = nrow

        # Calculate the base figure size
        # This can get stretched later by a legend
        figsize = (ncol * size * aspect, nrow * size)

        # Validate some inputs
        if col_wrap is not None:
            margin_titles = False

        # Build the subplot keyword dictionary
        subplot_kws = {} if subplot_kws is None else subplot_kws.copy()
        gridspec_kws = {} if gridspec_kws is None else gridspec_kws.copy()
        if xlim is not None:
            subplot_kws["xlim"] = xlim
        if ylim is not None:
            subplot_kws["ylim"] = ylim

        # Initialize the subplot grid
        if col_wrap is None:
            kwargs = dict(figsize=figsize, squeeze=False,
                          sharex=sharex, sharey=sharey,
                          subplot_kw=subplot_kws,
                          gridspec_kw=gridspec_kws)

            if OLD_MPL:
                _ = kwargs.pop('gridspec_kw', None)
                if gridspec_kws:
                    msg = "gridspec module only available in mpl >= {}"
                    warnings.warn(msg.format(MPL_GRIDSPEC_VERSION))

            fig, axes = plt.subplots(nrow, ncol, **kwargs)
            self.axes = axes

        else:
            # If wrapping the col variable we need to make the grid ourselves
            if gridspec_kws:
                warnings.warn("`gridspec_kws` ignored when using `col_wrap`")

            n_axes = len(col_names)
            fig = plt.figure(figsize=figsize)
            axes = np.empty(n_axes, object)
            axes[0] = fig.add_subplot(nrow, ncol, 1, **subplot_kws)
            if sharex:
                subplot_kws["sharex"] = axes[0]
            if sharey:
                subplot_kws["sharey"] = axes[0]
            for i in range(1, n_axes):
                axes[i] = fig.add_subplot(nrow, ncol, i + 1, **subplot_kws)
            self.axes = axes

            # Now we turn off labels on the inner axes
            if sharex:
                for ax in self._not_bottom_axes:
                    for label in ax.get_xticklabels():
                        label.set_visible(False)
                    ax.xaxis.offsetText.set_visible(False)
            if sharey:
                for ax in self._not_left_axes:
                    for label in ax.get_yticklabels():
                        label.set_visible(False)
                    ax.yaxis.offsetText.set_visible(False)

        # Set up the class attributes
        # ---------------------------

        # First the public API
        self.data = data
        self.fig = fig
        self.axes = axes

        self.row_names = row_names
        self.col_names = col_names
        self.hue_names = hue_names
        self.hue_kws = hue_kws

        # Next the private variables
        self._nrow = nrow
        self._row_var = row
        self._ncol = ncol
        self._col_var = col

        self._margin_titles = margin_titles
        self._col_wrap = col_wrap
        self._hue_var = hue_var
        self._colors = colors
        self._legend_out = legend_out
        self._legend = None
        self._legend_data = {}
        self._x_var = None
        self._y_var = None
        self._dropna = dropna
        self._not_na = not_na

        # Make the axes look good
        fig.tight_layout()
        if despine:
            self.despine()

    __init__.__doc__ = dedent("""\
        Initialize the matplotlib figure and FacetGrid object.

        The :class:`FacetGrid` is an object that links a Pandas DataFrame to
        a matplotlib figure with a particular structure.

        In particular, :class:`FacetGrid` is used to draw plots with multiple
        Axes where each Axes shows the same relationship conditioned on
        different levels of some variable. It's possible to condition on up to
        three variables by assigning variables to the rows and columns of the
        grid and using different colors for the plot elements.

        The general approach to plotting here is called "small multiples",
        where the same kind of plot is repeated multiple times, and the
        specific use of small multiples to display the same relationship
        conditioned on one ore more other variables is often called a "trellis
        plot".

        The basic workflow is to initialize the :class:`FacetGrid` object with
        the dataset and the variables that are used to structure the grid. Then
        one or more plotting functions can be applied to each subset by calling
        :meth:`FacetGrid.map` or :meth:`FacetGrid.map_dataframe`. Finally, the
        plot can be tweaked with other methods to do things like change the
        axis labels, use different ticks, or add a legend. See the detailed
        code examples below for more information.

        Parameters
        ----------
        {data}
        row, col, hue : strings
            Variables that define subsets of the data, which will be drawn on
            separate facets in the grid. See the ``*_order`` parameters to
            control the order of levels of this variable.
        {col_wrap}
        {share_xy}
        {size}
        {aspect}
        {palette}
        {{row,col,hue}}_order : lists, optional
            Order for the levels of the faceting variables. By default, this
            will be the order that the levels appear in ``data`` or, if the
            variables are pandas categoricals, the category order.
        hue_kws : dictionary of param -> list of values mapping
            Other keyword arguments to insert into the plotting call to let
            other plot attributes vary across levels of the hue variable (e.g.
            the markers in a scatterplot).
        {legend_out}
        despine : boolean, optional
            Remove the top and right spines from the plots.
        {margin_titles}
        {{x, y}}lim: tuples, optional
            Limits for each of the axes on each facet (only relevant when
            share{{x, y}} is True.
        subplot_kws : dict, optional
            Dictionary of keyword arguments passed to matplotlib subplot(s)
            methods.
        gridspec_kws : dict, optional
            Dictionary of keyword arguments passed to matplotlib's ``gridspec``
            module (via ``plt.subplots``). Requires matplotlib >= 1.4 and is
            ignored if ``col_wrap`` is not ``None``.

        See Also
        --------
        PairGrid : Subplot grid for plotting pairwise relationships.
        lmplot : Combine a regression plot and a :class:`FacetGrid`.
        factorplot : Combine a categorical plot and a :class:`FacetGrid`.

        Examples
        --------

        Initialize a 2x2 grid of facets using the tips dataset:

        .. plot::
            :context: close-figs

            >>> import seaborn as sns; sns.set(style="ticks", color_codes=True)
            >>> tips = sns.load_dataset("tips")
            >>> g = sns.FacetGrid(tips, col="time", row="smoker")

        Draw a univariate plot on each facet:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> g = sns.FacetGrid(tips, col="time",  row="smoker")
            >>> g = g.map(plt.hist, "total_bill")

        (Note that it's not necessary to re-catch the returned variable; it's
        the same object, but doing so in the examples makes dealing with the
        doctests somewhat less annoying).

        Pass additional keyword arguments to the mapped function:

        .. plot::
            :context: close-figs

            >>> import numpy as np
            >>> bins = np.arange(0, 65, 5)
            >>> g = sns.FacetGrid(tips, col="time",  row="smoker")
            >>> g = g.map(plt.hist, "total_bill", bins=bins, color="r")

        Plot a bivariate function on each facet:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="time",  row="smoker")
            >>> g = g.map(plt.scatter, "total_bill", "tip", edgecolor="w")

        Assign one of the variables to the color of the plot elements:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="time",  hue="smoker")
            >>> g = (g.map(plt.scatter, "total_bill", "tip", edgecolor="w")
            ...       .add_legend())

        Change the size and aspect ratio of each facet:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="day", size=4, aspect=.5)
            >>> g = g.map(plt.hist, "total_bill", bins=bins)

        Specify the order for plot elements:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="smoker", col_order=["Yes", "No"])
            >>> g = g.map(plt.hist, "total_bill", bins=bins, color="m")

        Use a different color palette:

        .. plot::
            :context: close-figs

            >>> kws = dict(s=50, linewidth=.5, edgecolor="w")
            >>> g = sns.FacetGrid(tips, col="sex", hue="time", palette="Set1",
            ...                   hue_order=["Dinner", "Lunch"])
            >>> g = (g.map(plt.scatter, "total_bill", "tip", **kws)
            ...      .add_legend())

        Use a dictionary mapping hue levels to colors:

        .. plot::
            :context: close-figs

            >>> pal = dict(Lunch="seagreen", Dinner="gray")
            >>> g = sns.FacetGrid(tips, col="sex", hue="time", palette=pal,
            ...                   hue_order=["Dinner", "Lunch"])
            >>> g = (g.map(plt.scatter, "total_bill", "tip", **kws)
            ...      .add_legend())

        Additionally use a different marker for the hue levels:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="sex", hue="time", palette=pal,
            ...                   hue_order=["Dinner", "Lunch"],
            ...                   hue_kws=dict(marker=["^", "v"]))
            >>> g = (g.map(plt.scatter, "total_bill", "tip", **kws)
            ...      .add_legend())

        "Wrap" a column variable with many levels into the rows:

        .. plot::
            :context: close-figs

            >>> attend = sns.load_dataset("attention")
            >>> g = sns.FacetGrid(attend, col="subject", col_wrap=5, size=1.5)
            >>> g = g.map(plt.plot, "solutions", "score", marker=".")

        Define a custom bivariate function to map onto the grid:

        .. plot::
            :context: close-figs

            >>> from scipy import stats
            >>> def qqplot(x, y, **kwargs):
            ...     _, xr = stats.probplot(x, fit=False)
            ...     _, yr = stats.probplot(y, fit=False)
            ...     plt.scatter(xr, yr, **kwargs)
            >>> g = sns.FacetGrid(tips, col="smoker", hue="sex")
            >>> g = (g.map(qqplot, "total_bill", "tip", **kws)
            ...       .add_legend())

        Define a custom function that uses a ``DataFrame`` object and accepts
        column names as positional variables:

        .. plot::
            :context: close-figs

            >>> import pandas as pd
            >>> df = pd.DataFrame(
            ...     data=np.random.randn(90, 4),
            ...     columns=pd.Series(list("ABCD"), name="walk"),
            ...     index=pd.date_range("2015-01-01", "2015-03-31",
            ...                         name="date"))
            >>> df = df.cumsum(axis=0).stack().reset_index(name="val")
            >>> def dateplot(x, y, **kwargs):
            ...     ax = plt.gca()
            ...     data = kwargs.pop("data")
            ...     data.plot(x=x, y=y, ax=ax, grid=False, **kwargs)
            >>> g = sns.FacetGrid(df, col="walk", col_wrap=2, size=3.5)
            >>> g = g.map_dataframe(dateplot, "date", "val")

        Use different axes labels after plotting:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="smoker", row="sex")
            >>> g = (g.map(plt.scatter, "total_bill", "tip", color="g", **kws)
            ...       .set_axis_labels("Total bill (US Dollars)", "Tip"))

        Set other attributes that are shared across the facetes:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="smoker", row="sex")
            >>> g = (g.map(plt.scatter, "total_bill", "tip", color="r", **kws)
            ...       .set(xlim=(0, 60), ylim=(0, 12),
            ...            xticks=[10, 30, 50], yticks=[2, 6, 10]))

        Use a different template for the facet titles:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="size", col_wrap=3)
            >>> g = (g.map(plt.hist, "tip", bins=np.arange(0, 13), color="c")
            ...       .set_titles("{{col_name}} diners"))

        Tighten the facets:

        .. plot::
            :context: close-figs

            >>> g = sns.FacetGrid(tips, col="smoker", row="sex",
            ...                   margin_titles=True)
            >>> g = (g.map(plt.scatter, "total_bill", "tip", color="m", **kws)
            ...       .set(xlim=(0, 60), ylim=(0, 12),
            ...            xticks=[10, 30, 50], yticks=[2, 6, 10])
            ...       .fig.subplots_adjust(wspace=.05, hspace=.05))

        """).format(**_facet_docs)

    def facet_data(self):
        """Generator for name indices and data subsets for each facet.

        Yields
        ------
        (i, j, k), data_ijk : tuple of ints, DataFrame
            The ints provide an index into the {row, col, hue}_names attribute,
            and the dataframe contains a subset of the full data corresponding
            to each facet. The generator yields subsets that correspond with
            the self.axes.flat iterator, or self.axes[i, j] when `col_wrap`
            is None.

        """
        data = self.data

        # Construct masks for the row variable
        if self._nrow == 1 or self._col_wrap is not None:
            row_masks = [np.repeat(True, len(self.data))]
        else:
            row_masks = [data[self._row_var] == n for n in self.row_names]

        # Construct masks for the column variable
        if self._ncol == 1:
            col_masks = [np.repeat(True, len(self.data))]
        else:
            col_masks = [data[self._col_var] == n for n in self.col_names]

        # Construct masks for the hue variable
        if len(self._colors) == 1:
            hue_masks = [np.repeat(True, len(self.data))]
        else:
            hue_masks = [data[self._hue_var] == n for n in self.hue_names]

        # Here is the main generator loop
        for (i, row), (j, col), (k, hue) in product(enumerate(row_masks),
                                                    enumerate(col_masks),
                                                    enumerate(hue_masks)):
            data_ijk = data[row & col & hue & self._not_na]
            yield (i, j, k), data_ijk

    def map(self, func, *args, **kwargs):
        """Apply a plotting function to each facet's subset of the data.

        Parameters
        ----------
        func : callable
            A plotting function that takes data and keyword arguments. It
            must plot to the currently active matplotlib Axes and take a
            `color` keyword argument. If faceting on the `hue` dimension,
            it must also take a `label` keyword argument.
        args : strings
            Column names in self.data that identify variables with data to
            plot. The data for each variable is passed to `func` in the
            order the variables are specified in the call.
        kwargs : keyword arguments
            All keyword arguments are passed to the plotting function.

        Returns
        -------
        self : object
            Returns self.

        """
        # If color was a keyword argument, grab it here
        kw_color = kwargs.pop("color", None)

        # Check for categorical plots without order information
        if func.__module__ == "seaborn.categorical":
            if "order" not in kwargs:
                warning = ("Using the {} function without specifying "
                           "`order` is likely to produce an incorrect "
                           "plot.".format(func.__name__))
                warnings.warn(warning)
            if len(args) == 3 and "hue_order" not in kwargs:
                warning = ("Using the {} function without specifying "
                           "`hue_order` is likely to produce an incorrect "
                           "plot.".format(func.__name__))
                warnings.warn(warning)

        # Iterate over the data subsets
        for (row_i, col_j, hue_k), data_ijk in self.facet_data():

            # If this subset is null, move on
            if not data_ijk.values.size:
                continue

            # Get the current axis
            ax = self.facet_axis(row_i, col_j)

            # Decide what color to plot with
            kwargs["color"] = self._facet_color(hue_k, kw_color)

            # Insert the other hue aesthetics if appropriate
            for kw, val_list in self.hue_kws.items():
                kwargs[kw] = val_list[hue_k]

            # Insert a label in the keyword arguments for the legend
            if self._hue_var is not None:
                kwargs["label"] = utils.to_utf8(self.hue_names[hue_k])

            # Get the actual data we are going to plot with
            plot_data = data_ijk[list(args)]
            if self._dropna:
                plot_data = plot_data.dropna()
            plot_args = [v for k, v in plot_data.iteritems()]

            # Some matplotlib functions don't handle pandas objects correctly
            if func.__module__ is not None:
                if func.__module__.startswith("matplotlib"):
                    plot_args = [v.values for v in plot_args]

            # Draw the plot
            self._facet_plot(func, ax, plot_args, kwargs)

        # Finalize the annotations and layout
        self._finalize_grid(args[:2])

        return self

    def map_dataframe(self, func, *args, **kwargs):
        """Like `map` but passes args as strings and inserts data in kwargs.

        This method is suitable for plotting with functions that accept a
        long-form DataFrame as a `data` keyword argument and access the
        data in that DataFrame using string variable names.

        Parameters
        ----------
        func : callable
            A plotting function that takes data and keyword arguments. Unlike
            the `map` method, a function used here must "understand" Pandas
            objects. It also must plot to the currently active matplotlib Axes
            and take a `color` keyword argument. If faceting on the `hue`
            dimension, it must also take a `label` keyword argument.
        args : strings
            Column names in self.data that identify variables with data to
            plot. The data for each variable is passed to `func` in the
            order the variables are specified in the call.
        kwargs : keyword arguments
            All keyword arguments are passed to the plotting function.

        Returns
        -------
        self : object
            Returns self.

        """

        # If color was a keyword argument, grab it here
        kw_color = kwargs.pop("color", None)

        # Iterate over the data subsets
        for (row_i, col_j, hue_k), data_ijk in self.facet_data():

            # If this subset is null, move on
            if not data_ijk.values.size:
                continue

            # Get the current axis
            ax = self.facet_axis(row_i, col_j)

            # Decide what color to plot with
            kwargs["color"] = self._facet_color(hue_k, kw_color)

            # Insert the other hue aesthetics if appropriate
            for kw, val_list in self.hue_kws.items():
                kwargs[kw] = val_list[hue_k]

            # Insert a label in the keyword arguments for the legend
            if self._hue_var is not None:
                kwargs["label"] = self.hue_names[hue_k]

            # Stick the facet dataframe into the kwargs
            if self._dropna:
                data_ijk = data_ijk.dropna()
            kwargs["data"] = data_ijk

            # Draw the plot
            self._facet_plot(func, ax, args, kwargs)

        # Finalize the annotations and layout
        self._finalize_grid(args[:2])

        return self

    def _facet_color(self, hue_index, kw_color):

        color = self._colors[hue_index]
        if kw_color is not None:
            return kw_color
        elif color is not None:
            return color

    def _facet_plot(self, func, ax, plot_args, plot_kwargs):

        # Draw the plot
        func(*plot_args, **plot_kwargs)

        # Sort out the supporting information
        self._update_legend_data(ax)
        self._clean_axis(ax)

    def _finalize_grid(self, axlabels):
        """Finalize the annotations and layout."""
        self.set_axis_labels(*axlabels)
        self.set_titles()
        self.fig.tight_layout()

    def facet_axis(self, row_i, col_j):
        """Make the axis identified by these indices active and return it."""

        # Calculate the actual indices of the axes to plot on
        if self._col_wrap is not None:
            ax = self.axes.flat[col_j]
        else:
            ax = self.axes[row_i, col_j]

        # Get a reference to the axes object we want, and make it active
        plt.sca(ax)
        return ax

    def despine(self, **kwargs):
        """Remove axis spines from the facets."""
        utils.despine(self.fig, **kwargs)
        return self

    def set_axis_labels(self, x_var=None, y_var=None):
        """Set axis labels on the left column and bottom row of the grid."""
        if x_var is not None:
            self._x_var = x_var
            self.set_xlabels(x_var)
        if y_var is not None:
            self._y_var = y_var
            self.set_ylabels(y_var)
        return self

    def set_xlabels(self, label=None, **kwargs):
        """Label the x axis on the bottom row of the grid."""
        if label is None:
            label = self._x_var
        for ax in self._bottom_axes:
            ax.set_xlabel(label, **kwargs)
        return self

    def set_ylabels(self, label=None, **kwargs):
        """Label the y axis on the left column of the grid."""
        if label is None:
            label = self._y_var
        for ax in self._left_axes:
            ax.set_ylabel(label, **kwargs)
        return self

    def set_xticklabels(self, labels=None, step=None, **kwargs):
        """Set x axis tick labels on the bottom row of the grid."""
        for ax in self._bottom_axes:
            if labels is None:
                labels = [l.get_text() for l in ax.get_xticklabels()]
                if step is not None:
                    xticks = ax.get_xticks()[::step]
                    labels = labels[::step]
                    ax.set_xticks(xticks)
            ax.set_xticklabels(labels, **kwargs)
        return self

    def set_yticklabels(self, labels=None, **kwargs):
        """Set y axis tick labels on the left column of the grid."""
        for ax in self._left_axes:
            if labels is None:
                labels = [l.get_text() for l in ax.get_yticklabels()]
            ax.set_yticklabels(labels, **kwargs)
        return self

    def set_titles(self, template=None, row_template=None,  col_template=None,
                   **kwargs):
        """Draw titles either above each facet or on the grid margins.

        Parameters
        ----------
        template : string
            Template for all titles with the formatting keys {col_var} and
            {col_name} (if using a `col` faceting variable) and/or {row_var}
            and {row_name} (if using a `row` faceting variable).
        row_template:
            Template for the row variable when titles are drawn on the grid
            margins. Must have {row_var} and {row_name} formatting keys.
        col_template:
            Template for the row variable when titles are drawn on the grid
            margins. Must have {col_var} and {col_name} formatting keys.

        Returns
        -------
        self: object
            Returns self.

        """
        args = dict(row_var=self._row_var, col_var=self._col_var)
        kwargs["size"] = kwargs.pop("size", mpl.rcParams["axes.labelsize"])

        # Establish default templates
        if row_template is None:
            row_template = "{row_var} = {row_name}"
        if col_template is None:
            col_template = "{col_var} = {col_name}"
        if template is None:
            if self._row_var is None:
                template = col_template
            elif self._col_var is None:
                template = row_template
            else:
                template = " | ".join([row_template, col_template])

        row_template = utils.to_utf8(row_template)
        col_template = utils.to_utf8(col_template)
        template = utils.to_utf8(template)

        if self._margin_titles:
            if self.row_names is not None:
                # Draw the row titles on the right edge of the grid
                for i, row_name in enumerate(self.row_names):
                    ax = self.axes[i, -1]
                    args.update(dict(row_name=row_name))
                    title = row_template.format(**args)
                    bgcolor = self.fig.get_facecolor()
                    ax.annotate(title, xy=(1.02, .5), xycoords="axes fraction",
                                rotation=270, ha="left", va="center",
                                backgroundcolor=bgcolor, **kwargs)

            if self.col_names is not None:
                # Draw the column titles  as normal titles
                for j, col_name in enumerate(self.col_names):
                    args.update(dict(col_name=col_name))
                    title = col_template.format(**args)
                    self.axes[0, j].set_title(title, **kwargs)

            return self

        # Otherwise title each facet with all the necessary information
        if (self._row_var is not None) and (self._col_var is not None):
            for i, row_name in enumerate(self.row_names):
                for j, col_name in enumerate(self.col_names):
                    args.update(dict(row_name=row_name, col_name=col_name))
                    title = template.format(**args)
                    self.axes[i, j].set_title(title, **kwargs)
        elif self.row_names is not None and len(self.row_names):
            for i, row_name in enumerate(self.row_names):
                args.update(dict(row_name=row_name))
                title = template.format(**args)
                self.axes[i, 0].set_title(title, **kwargs)
        elif self.col_names is not None and len(self.col_names):
            for i, col_name in enumerate(self.col_names):
                args.update(dict(col_name=col_name))
                title = template.format(**args)
                # Index the flat array so col_wrap works
                self.axes.flat[i].set_title(title, **kwargs)
        return self

    @property
    def ax(self):
        """Easy access to single axes."""
        if self.axes.shape == (1, 1):
            return self.axes[0, 0]
        else:
            raise AttributeError

    @property
    def _inner_axes(self):
        """Return a flat array of the inner axes."""
        if self._col_wrap is None:
            return self.axes[:-1, 1:].flat
        else:
            axes = []
            n_empty = self._nrow * self._ncol - self._n_facets
            for i, ax in enumerate(self.axes):
                append = (i % self._ncol and
                          i < (self._ncol * (self._nrow - 1)) and
                          i < (self._ncol * (self._nrow - 1) - n_empty))
                if append:
                    axes.append(ax)
            return np.array(axes, object).flat

    @property
    def _left_axes(self):
        """Return a flat array of the left column of axes."""
        if self._col_wrap is None:
            return self.axes[:, 0].flat
        else:
            axes = []
            for i, ax in enumerate(self.axes):
                if not i % self._ncol:
                    axes.append(ax)
            return np.array(axes, object).flat

    @property
    def _not_left_axes(self):
        """Return a flat array of axes that aren't on the left column."""
        if self._col_wrap is None:
            return self.axes[:, 1:].flat
        else:
            axes = []
            for i, ax in enumerate(self.axes):
                if i % self._ncol:
                    axes.append(ax)
            return np.array(axes, object).flat

    @property
    def _bottom_axes(self):
        """Return a flat array of the bottom row of axes."""
        if self._col_wrap is None:
            return self.axes[-1, :].flat
        else:
            axes = []
            n_empty = self._nrow * self._ncol - self._n_facets
            for i, ax in enumerate(self.axes):
                append = (i >= (self._ncol * (self._nrow - 1)) or
                          i >= (self._ncol * (self._nrow - 1) - n_empty))
                if append:
                    axes.append(ax)
            return np.array(axes, object).flat

    @property
    def _not_bottom_axes(self):
        """Return a flat array of axes that aren't on the bottom row."""
        if self._col_wrap is None:
            return self.axes[:-1, :].flat
        else:
            axes = []
            n_empty = self._nrow * self._ncol - self._n_facets
            for i, ax in enumerate(self.axes):
                append = (i < (self._ncol * (self._nrow - 1)) and
                          i < (self._ncol * (self._nrow - 1) - n_empty))
                if append:
                    axes.append(ax)
            return np.array(axes, object).flat


class PairGrid(Grid):
    """Subplot grid for plotting pairwise relationships in a dataset."""

    def __init__(self, data, hue=None, hue_order=None, palette=None,
                 hue_kws=None, vars=None, x_vars=None, y_vars=None,
                 diag_sharey=True, size=2.5, aspect=1,
                 despine=True, dropna=True):
        """Initialize the plot figure and PairGrid object.

        Parameters
        ----------
        data : DataFrame
            Tidy (long-form) dataframe where each column is a variable and
            each row is an observation.
        hue : string (variable name), optional
            Variable in ``data`` to map plot aspects to different colors.
        hue_order : list of strings
            Order for the levels of the hue variable in the palette
        palette : dict or seaborn color palette
            Set of colors for mapping the ``hue`` variable. If a dict, keys
            should be values  in the ``hue`` variable.
        hue_kws : dictionary of param -> list of values mapping
            Other keyword arguments to insert into the plotting call to let
            other plot attributes vary across levels of the hue variable (e.g.
            the markers in a scatterplot).
        vars : list of variable names, optional
            Variables within ``data`` to use, otherwise use every column with
            a numeric datatype.
        {x, y}_vars : lists of variable names, optional
            Variables within ``data`` to use separately for the rows and
            columns of the figure; i.e. to make a non-square plot.
        size : scalar, optional
            Height (in inches) of each facet.
        aspect : scalar, optional
            Aspect * size gives the width (in inches) of each facet.
        despine : boolean, optional
            Remove the top and right spines from the plots.
        dropna : boolean, optional
            Drop missing values from the data before plotting.

        See Also
        --------
        pairplot : Easily drawing common uses of :class:`PairGrid`.
        FacetGrid : Subplot grid for plotting conditional relationships.

        Examples
        --------

        Draw a scatterplot for each pairwise relationship:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns; sns.set()
            >>> iris = sns.load_dataset("iris")
            >>> g = sns.PairGrid(iris)
            >>> g = g.map(plt.scatter)

        Show a univariate distribution on the diagonal:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris)
            >>> g = g.map_diag(plt.hist)
            >>> g = g.map_offdiag(plt.scatter)

        (It's not actually necessary to catch the return value every time,
        as it is the same object, but it makes it easier to deal with the
        doctests).

        Color the points using a categorical variable:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, hue="species")
            >>> g = g.map_diag(plt.hist)
            >>> g = g.map_offdiag(plt.scatter)
            >>> g = g.add_legend()

        Use a different style to show multiple histograms:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, hue="species")
            >>> g = g.map_diag(plt.hist, histtype="step", linewidth=3)
            >>> g = g.map_offdiag(plt.scatter)
            >>> g = g.add_legend()

        Plot a subset of variables

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, vars=["sepal_length", "sepal_width"])
            >>> g = g.map(plt.scatter)

        Pass additional keyword arguments to the functions

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris)
            >>> g = g.map_diag(plt.hist, edgecolor="w")
            >>> g = g.map_offdiag(plt.scatter, edgecolor="w", s=40)

        Use different variables for the rows and columns:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris,
            ...                  x_vars=["sepal_length", "sepal_width"],
            ...                  y_vars=["petal_length", "petal_width"])
            >>> g = g.map(plt.scatter)

        Use different functions on the upper and lower triangles:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris)
            >>> g = g.map_upper(plt.scatter)
            >>> g = g.map_lower(sns.kdeplot, cmap="Blues_d")
            >>> g = g.map_diag(sns.kdeplot, lw=3, legend=False)

        Use different colors and markers for each categorical level:

        .. plot::
            :context: close-figs

            >>> g = sns.PairGrid(iris, hue="species", palette="Set2",
            ...                  hue_kws={"marker": ["o", "s", "D"]})
            >>> g = g.map(plt.scatter, linewidths=1, edgecolor="w", s=40)
            >>> g = g.add_legend()

        """

        # Sort out the variables that define the grid
        if vars is not None:
            x_vars = list(vars)
            y_vars = list(vars)
        elif (x_vars is not None) or (y_vars is not None):
            if (x_vars is None) or (y_vars is None):
                raise ValueError("Must specify `x_vars` and `y_vars`")
        else:
            numeric_cols = self._find_numeric_cols(data)
            x_vars = numeric_cols
            y_vars = numeric_cols

        if np.isscalar(x_vars):
            x_vars = [x_vars]
        if np.isscalar(y_vars):
            y_vars = [y_vars]

        self.x_vars = list(x_vars)
        self.y_vars = list(y_vars)
        self.square_grid = self.x_vars == self.y_vars

        # Create the figure and the array of subplots
        figsize = len(x_vars) * size * aspect, len(y_vars) * size

        fig, axes = plt.subplots(len(y_vars), len(x_vars),
                                 figsize=figsize,
                                 sharex="col", sharey="row",
                                 squeeze=False)

        self.fig = fig
        self.axes = axes
        self.data = data

        # Save what we are going to do with the diagonal
        self.diag_sharey = diag_sharey
        self.diag_axes = None

        # Label the axes
        self._add_axis_labels()

        # Sort out the hue variable
        self._hue_var = hue
        if hue is None:
            self.hue_names = ["_nolegend_"]
            self.hue_vals = pd.Series(["_nolegend_"] * len(data),
                                      index=data.index)
        else:
            hue_names = utils.categorical_order(data[hue], hue_order)
            if dropna:
                # Filter NA from the list of unique hue names
                hue_names = list(filter(pd.notnull, hue_names))
            self.hue_names = hue_names
            self.hue_vals = data[hue]

        # Additional dict of kwarg -> list of values for mapping the hue var
        self.hue_kws = hue_kws if hue_kws is not None else {}

        self.palette = self._get_palette(data, hue, hue_order, palette)
        self._legend_data = {}

        # Make the plot look nice
        if despine:
            utils.despine(fig=fig)
        fig.tight_layout()

    def map(self, func, **kwargs):
        """Plot with the same function in every subplot.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """
        kw_color = kwargs.pop("color", None)
        for i, y_var in enumerate(self.y_vars):
            for j, x_var in enumerate(self.x_vars):
                hue_grouped = self.data.groupby(self.hue_vals)
                for k, label_k in enumerate(self.hue_names):

                    # Attempt to get data for this level, allowing for empty
                    try:
                        data_k = hue_grouped.get_group(label_k)
                    except KeyError:
                        data_k = pd.DataFrame(columns=self.data.columns,
                                              dtype=np.float)

                    ax = self.axes[i, j]
                    plt.sca(ax)

                    # Insert the other hue aesthetics if appropriate
                    for kw, val_list in self.hue_kws.items():
                        kwargs[kw] = val_list[k]

                    color = self.palette[k] if kw_color is None else kw_color
                    func(data_k[x_var], data_k[y_var],
                         label=label_k, color=color, **kwargs)

                self._clean_axis(ax)
                self._update_legend_data(ax)

        if kw_color is not None:
            kwargs["color"] = kw_color
        self._add_axis_labels()

        return self

    def map_diag(self, func, **kwargs):
        """Plot with a univariate function on each diagonal subplot.

        Parameters
        ----------
        func : callable plotting function
            Must take an x array as a positional arguments and draw onto the
            "currently active" matplotlib Axes. There is a special case when
            using a ``hue`` variable and ``plt.hist``; the histogram will be
            plotted with stacked bars.

        """
        # Add special diagonal axes for the univariate plot
        if self.square_grid and self.diag_axes is None:
            diag_axes = []
            for i, (var, ax) in enumerate(zip(self.x_vars,
                                              np.diag(self.axes))):
                if i and self.diag_sharey:
                    diag_ax = ax._make_twin_axes(sharex=ax,
                                                 sharey=diag_axes[0],
                                                 frameon=False)
                else:
                    diag_ax = ax._make_twin_axes(sharex=ax, frameon=False)
                diag_ax.set_axis_off()
                diag_axes.append(diag_ax)
            self.diag_axes = np.array(diag_axes, np.object)

        # Plot on each of the diagonal axes
        fixed_color = kwargs.pop("color", None)
        for i, var in enumerate(self.x_vars):
            ax = self.diag_axes[i]
            hue_grouped = self.data[var].groupby(self.hue_vals)

            # Special-case plt.hist with stacked bars
            if func is plt.hist:
                plt.sca(ax)

                vals = []
                for label in self.hue_names:
                    # Attempt to get data for this level, allowing for empty
                    try:
                        vals.append(np.asarray(hue_grouped.get_group(label)))
                    except KeyError:
                        vals.append(np.array([]))

                color = self.palette if fixed_color is None else fixed_color

                if "histtype" in kwargs:
                    func(vals, color=color, **kwargs)
                else:
                    func(vals, color=color, histtype="barstacked", **kwargs)

            else:
                plt.sca(ax)

                for k, label_k in enumerate(self.hue_names):

                    # Attempt to get data for this level, allowing for empty
                    try:
                        data_k = hue_grouped.get_group(label_k)
                    except KeyError:
                        data_k = np.array([])

                    if fixed_color is None:
                        color = self.palette[k]
                    else:
                        color = fixed_color

                    func(data_k, label=label_k, color=color, **kwargs)

            self._clean_axis(ax)

        self._add_axis_labels()

        return self

    def map_lower(self, func, **kwargs):
        """Plot with a bivariate function on the lower diagonal subplots.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """
        kw_color = kwargs.pop("color", None)
        for i, j in zip(*np.tril_indices_from(self.axes, -1)):
            hue_grouped = self.data.groupby(self.hue_vals)
            for k, label_k in enumerate(self.hue_names):

                # Attempt to get data for this level, allowing for empty
                try:
                    data_k = hue_grouped.get_group(label_k)
                except KeyError:
                    data_k = pd.DataFrame(columns=self.data.columns,
                                          dtype=np.float)

                ax = self.axes[i, j]
                plt.sca(ax)

                x_var = self.x_vars[j]
                y_var = self.y_vars[i]

                # Insert the other hue aesthetics if appropriate
                for kw, val_list in self.hue_kws.items():
                    kwargs[kw] = val_list[k]

                color = self.palette[k] if kw_color is None else kw_color
                func(data_k[x_var], data_k[y_var], label=label_k,
                     color=color, **kwargs)

            self._clean_axis(ax)
            self._update_legend_data(ax)

        if kw_color is not None:
            kwargs["color"] = kw_color
        self._add_axis_labels()

        return self

    def map_upper(self, func, **kwargs):
        """Plot with a bivariate function on the upper diagonal subplots.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """
        kw_color = kwargs.pop("color", None)
        for i, j in zip(*np.triu_indices_from(self.axes, 1)):

            hue_grouped = self.data.groupby(self.hue_vals)

            for k, label_k in enumerate(self.hue_names):

                # Attempt to get data for this level, allowing for empty
                try:
                    data_k = hue_grouped.get_group(label_k)
                except KeyError:
                    data_k = pd.DataFrame(columns=self.data.columns,
                                          dtype=np.float)

                ax = self.axes[i, j]
                plt.sca(ax)

                x_var = self.x_vars[j]
                y_var = self.y_vars[i]

                # Insert the other hue aesthetics if appropriate
                for kw, val_list in self.hue_kws.items():
                    kwargs[kw] = val_list[k]

                color = self.palette[k] if kw_color is None else kw_color
                func(data_k[x_var], data_k[y_var], label=label_k,
                     color=color, **kwargs)

            self._clean_axis(ax)
            self._update_legend_data(ax)

        if kw_color is not None:
            kwargs["color"] = kw_color

        return self

    def map_offdiag(self, func, **kwargs):
        """Plot with a bivariate function on the off-diagonal subplots.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """

        self.map_lower(func, **kwargs)
        self.map_upper(func, **kwargs)
        return self

    def _add_axis_labels(self):
        """Add labels to the left and bottom Axes."""
        for ax, label in zip(self.axes[-1, :], self.x_vars):
            ax.set_xlabel(label)
        for ax, label in zip(self.axes[:, 0], self.y_vars):
            ax.set_ylabel(label)

    def _find_numeric_cols(self, data):
        """Find which variables in a DataFrame are numeric."""
        # This can't be the best way to do this, but  I do not
        # know what the best way might be, so this seems ok
        numeric_cols = []
        for col in data:
            try:
                data[col].astype(np.float)
                numeric_cols.append(col)
            except (ValueError, TypeError):
                pass
        return numeric_cols


class JointGrid(object):
    """Grid for drawing a bivariate plot with marginal univariate plots."""

    def __init__(self, x, y, data=None, size=6, ratio=5, space=.2,
                 dropna=True, xlim=None, ylim=None):
        """Set up the grid of subplots.

        Parameters
        ----------
        x, y : strings or vectors
            Data or names of variables in ``data``.
        data : DataFrame, optional
            DataFrame when ``x`` and ``y`` are variable names.
        size : numeric
            Size of each side of the figure in inches (it will be square).
        ratio : numeric
            Ratio of joint axes size to marginal axes height.
        space : numeric, optional
            Space between the joint and marginal axes
        dropna : bool, optional
            If True, remove observations that are missing from `x` and `y`.
        {x, y}lim : two-tuples, optional
            Axis limits to set before plotting.

        See Also
        --------
        jointplot : High-level interface for drawing bivariate plots with
                    several different default plot kinds.

        Examples
        --------

        Initialize the figure but don't draw any plots onto it:

        .. plot::
            :context: close-figs

            >>> import seaborn as sns; sns.set(style="ticks", color_codes=True)
            >>> tips = sns.load_dataset("tips")
            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips)

        Add plots using default parameters:

        .. plot::
            :context: close-figs

            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot(sns.regplot, sns.distplot)

        Draw the join and marginal plots separately, which allows finer-level
        control other parameters:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter, color=".5", edgecolor="white")
            >>> g = g.plot_marginals(sns.distplot, kde=False, color=".5")

        Draw the two marginal plots separately:

        .. plot::
            :context: close-figs

            >>> import numpy as np
            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter, color="m", edgecolor="white")
            >>> _ = g.ax_marg_x.hist(tips["total_bill"], color="b", alpha=.6,
            ...                      bins=np.arange(0, 60, 5))
            >>> _ = g.ax_marg_y.hist(tips["tip"], color="r", alpha=.6,
            ...                      orientation="horizontal",
            ...                      bins=np.arange(0, 12, 1))

        Add an annotation with a statistic summarizing the bivariate
        relationship:

        .. plot::
            :context: close-figs

            >>> from scipy import stats
            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter,
            ...                  color="g", s=40, edgecolor="white")
            >>> g = g.plot_marginals(sns.distplot, kde=False, color="g")
            >>> g = g.annotate(stats.pearsonr)

        Use a custom function and formatting for the annotation

        .. plot::
            :context: close-figs

            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter,
            ...                  color="g", s=40, edgecolor="white")
            >>> g = g.plot_marginals(sns.distplot, kde=False, color="g")
            >>> rsquare = lambda a, b: stats.pearsonr(a, b)[0] ** 2
            >>> g = g.annotate(rsquare, template="{stat}: {val:.2f}",
            ...                stat="$R^2$", loc="upper left", fontsize=12)

        Remove the space between the joint and marginal axes:

        .. plot::
            :context: close-figs

            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips, space=0)
            >>> g = g.plot_joint(sns.kdeplot, cmap="Blues_d")
            >>> g = g.plot_marginals(sns.kdeplot, shade=True)

        Draw a smaller plot with relatively larger marginal axes:

        .. plot::
            :context: close-figs

            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips,
            ...                   size=5, ratio=2)
            >>> g = g.plot_joint(sns.kdeplot, cmap="Reds_d")
            >>> g = g.plot_marginals(sns.kdeplot, color="r", shade=True)

        Set limits on the axes:

        .. plot::
            :context: close-figs

            >>> g = sns.JointGrid(x="total_bill", y="tip", data=tips,
            ...                   xlim=(0, 50), ylim=(0, 8))
            >>> g = g.plot_joint(sns.kdeplot, cmap="Purples_d")
            >>> g = g.plot_marginals(sns.kdeplot, color="m", shade=True)

        """
        # Set up the subplot grid
        f = plt.figure(figsize=(size, size))
        gs = plt.GridSpec(ratio + 1, ratio + 1)

        ax_joint = f.add_subplot(gs[1:, :-1])
        ax_marg_x = f.add_subplot(gs[0, :-1], sharex=ax_joint)
        ax_marg_y = f.add_subplot(gs[1:, -1], sharey=ax_joint)

        self.fig = f
        self.ax_joint = ax_joint
        self.ax_marg_x = ax_marg_x
        self.ax_marg_y = ax_marg_y

        # Turn off tick visibility for the measure axis on the marginal plots
        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)

        # Turn off the ticks on the density axis for the marginal plots
        plt.setp(ax_marg_x.yaxis.get_majorticklines(), visible=False)
        plt.setp(ax_marg_x.yaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_y.xaxis.get_majorticklines(), visible=False)
        plt.setp(ax_marg_y.xaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_x.get_yticklabels(), visible=False)
        plt.setp(ax_marg_y.get_xticklabels(), visible=False)
        ax_marg_x.yaxis.grid(False)
        ax_marg_y.xaxis.grid(False)

        # Possibly extract the variables from a DataFrame
        if data is not None:
            x = data.get(x, x)
            y = data.get(y, y)

        for var in [x, y]:
            if isinstance(var, string_types):
                err = "Could not interpret input '{}'".format(var)
                raise ValueError(err)

        # Possibly drop NA
        if dropna:
            not_na = pd.notnull(x) & pd.notnull(y)
            x = x[not_na]
            y = y[not_na]

        # Find the names of the variables
        if hasattr(x, "name"):
            xlabel = x.name
            ax_joint.set_xlabel(xlabel)
        if hasattr(y, "name"):
            ylabel = y.name
            ax_joint.set_ylabel(ylabel)

        # Convert the x and y data to arrays for plotting
        self.x = np.asarray(x)
        self.y = np.asarray(y)

        if xlim is not None:
            ax_joint.set_xlim(xlim)
        if ylim is not None:
            ax_joint.set_ylim(ylim)

        # Make the grid look nice
        utils.despine(f)
        utils.despine(ax=ax_marg_x, left=True)
        utils.despine(ax=ax_marg_y, bottom=True)
        f.tight_layout()
        f.subplots_adjust(hspace=space, wspace=space)

    def plot(self, joint_func, marginal_func, annot_func=None):
        """Shortcut to draw the full plot.

        Use `plot_joint` and `plot_marginals` directly for more control.

        Parameters
        ----------
        joint_func, marginal_func: callables
            Functions to draw the bivariate and univariate plots.

        Returns
        -------
        self : JointGrid instance
            Returns `self`.

        """
        self.plot_marginals(marginal_func)
        self.plot_joint(joint_func)
        if annot_func is not None:
            self.annotate(annot_func)
        return self

    def plot_joint(self, func, **kwargs):
        """Draw a bivariate plot of `x` and `y`.

        Parameters
        ----------
        func : plotting callable
            This must take two 1d arrays of data as the first two
            positional arguments, and it must plot on the "current" axes.
        kwargs : key, value mappings
            Keyword argument are passed to the plotting function.

        Returns
        -------
        self : JointGrid instance
            Returns `self`.

        """
        plt.sca(self.ax_joint)
        func(self.x, self.y, **kwargs)

        return self

    def plot_marginals(self, func, **kwargs):
        """Draw univariate plots for `x` and `y` separately.

        Parameters
        ----------
        func : plotting callable
            This must take a 1d array of data as the first positional
            argument, it must plot on the "current" axes, and it must
            accept a "vertical" keyword argument to orient the measure
            dimension of the plot vertically.
        kwargs : key, value mappings
            Keyword argument are passed to the plotting function.

        Returns
        -------
        self : JointGrid instance
            Returns `self`.

        """
        kwargs["vertical"] = False
        plt.sca(self.ax_marg_x)
        func(self.x, **kwargs)

        kwargs["vertical"] = True
        plt.sca(self.ax_marg_y)
        func(self.y, **kwargs)

        return self

    def annotate(self, func, template=None, stat=None, loc="best", **kwargs):
        """Annotate the plot with a statistic about the relationship.

        Parameters
        ----------
        func : callable
            Statistical function that maps the x, y vectors either to (val, p)
            or to val.
        template : string format template, optional
            The template must have the format keys "stat" and "val";
            if `func` returns a p value, it should also have the key "p".
        stat : string, optional
            Name to use for the statistic in the annotation, by default it
            uses the name of `func`.
        loc : string or int, optional
            Matplotlib legend location code; used to place the annotation.
        kwargs : key, value mappings
            Other keyword arguments are passed to `ax.legend`, which formats
            the annotation.

        Returns
        -------
        self : JointGrid instance.
            Returns `self`.

        """
        default_template = "{stat} = {val:.2g}; p = {p:.2g}"

        # Call the function and determine the form of the return value(s)
        out = func(self.x, self.y)
        try:
            val, p = out
        except TypeError:
            val, p = out, None
            default_template, _ = default_template.split(";")

        # Set the default template
        if template is None:
            template = default_template

        # Default to name of the function
        if stat is None:
            stat = func.__name__

        # Format the annotation
        if p is None:
            annotation = template.format(stat=stat, val=val)
        else:
            annotation = template.format(stat=stat, val=val, p=p)

        # Draw an invisible plot and use the legend to draw the annotation
        # This is a bit of a hack, but `loc=best` works nicely and is not
        # easily abstracted.
        phantom, = self.ax_joint.plot(self.x, self.y, linestyle="", alpha=0)
        self.ax_joint.legend([phantom], [annotation], loc=loc, **kwargs)
        phantom.remove()

        return self

    def set_axis_labels(self, xlabel="", ylabel="", **kwargs):
        """Set the axis labels on the bivariate axes.

        Parameters
        ----------
        xlabel, ylabel : strings
            Label names for the x and y variables.
        kwargs : key, value mappings
            Other keyword arguments are passed to the set_xlabel or
            set_ylabel.

        Returns
        -------
        self : JointGrid instance
            returns `self`

        """
        self.ax_joint.set_xlabel(xlabel, **kwargs)
        self.ax_joint.set_ylabel(ylabel, **kwargs)
        return self

    def savefig(self, *args, **kwargs):
        """Wrap figure.savefig defaulting to tight bounding box."""
        kwargs.setdefault("bbox_inches", "tight")
        self.fig.savefig(*args, **kwargs)


def pairplot(data, hue=None, hue_order=None, palette=None,
             vars=None, x_vars=None, y_vars=None,
             kind="scatter", diag_kind="hist", markers=None,
             size=2.5, aspect=1, dropna=True,
             plot_kws=None, diag_kws=None, grid_kws=None):
    """Plot pairwise relationships in a dataset.

    By default, this function will create a grid of Axes such that each
    variable in ``data`` will by shared in the y-axis across a single row and
    in the x-axis across a single column. The diagonal Axes are treated
    differently, drawing a plot to show the univariate distribution of the data
    for the variable in that column.

    It is also possible to show a subset of variables or plot different
    variables on the rows and columns.

    This is a high-level interface for :class:`PairGrid` that is intended to
    make it easy to draw a few common styles. You should use :class:`PairGrid`
    directly if you need more flexibility.

    Parameters
    ----------
    data : DataFrame
        Tidy (long-form) dataframe where each column is a variable and
        each row is an observation.
    hue : string (variable name), optional
        Variable in ``data`` to map plot aspects to different colors.
    hue_order : list of strings
        Order for the levels of the hue variable in the palette
    palette : dict or seaborn color palette
        Set of colors for mapping the ``hue`` variable. If a dict, keys
        should be values  in the ``hue`` variable.
    vars : list of variable names, optional
        Variables within ``data`` to use, otherwise use every column with
        a numeric datatype.
    {x, y}_vars : lists of variable names, optional
        Variables within ``data`` to use separately for the rows and
        columns of the figure; i.e. to make a non-square plot.
    kind : {'scatter', 'reg'}, optional
        Kind of plot for the non-identity relationships.
    diag_kind : {'hist', 'kde'}, optional
        Kind of plot for the diagonal subplots.
    markers : single matplotlib marker code or list, optional
        Either the marker to use for all datapoints or a list of markers with
        a length the same as the number of levels in the hue variable so that
        differently colored points will also have different scatterplot
        markers.
    size : scalar, optional
        Height (in inches) of each facet.
    aspect : scalar, optional
        Aspect * size gives the width (in inches) of each facet.
    dropna : boolean, optional
        Drop missing values from the data before plotting.
    {plot, diag, grid}_kws : dicts, optional
        Dictionaries of keyword arguments.

    Returns
    -------
    grid : PairGrid
        Returns the underlying ``PairGrid`` instance for further tweaking.

    See Also
    --------
    PairGrid : Subplot grid for more flexible plotting of pairwise
               relationships.

    Examples
    --------

    Draw scatterplots for joint relationships and histograms for univariate
    distributions:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set(style="ticks", color_codes=True)
        >>> iris = sns.load_dataset("iris")
        >>> g = sns.pairplot(iris)

    Show different levels of a categorical variable by the color of plot
    elements:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, hue="species")

    Use a different color palette:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, hue="species", palette="husl")

    Use different markers for each level of the hue variable:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, hue="species", markers=["o", "s", "D"])

    Plot a subset of variables:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, vars=["sepal_width", "sepal_length"])

    Draw larger plots:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, size=3,
        ...                  vars=["sepal_width", "sepal_length"])

    Plot different variables in the rows and columns:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris,
        ...                  x_vars=["sepal_width", "sepal_length"],
        ...                  y_vars=["petal_width", "petal_length"])

    Use kernel density estimates for univariate plots:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, diag_kind="kde")

    Fit linear regression models to the scatter plots:

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, kind="reg")

    Pass keyword arguments down to the underlying functions (it may be easier
    to use :class:`PairGrid` directly):

    .. plot::
        :context: close-figs

        >>> g = sns.pairplot(iris, diag_kind="kde", markers="+",
        ...                  plot_kws=dict(s=50, edgecolor="b", linewidth=1),
        ...                  diag_kws=dict(shade=True))

    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError(
            "'data' must be pandas DataFrame object, not: {typefound}".format(
                typefound=type(data)))

    if plot_kws is None:
        plot_kws = {}
    if diag_kws is None:
        diag_kws = {}
    if grid_kws is None:
        grid_kws = {}

    # Set up the PairGrid
    diag_sharey = diag_kind == "hist"
    grid = PairGrid(data, vars=vars, x_vars=x_vars, y_vars=y_vars, hue=hue,
                    hue_order=hue_order, palette=palette,
                    diag_sharey=diag_sharey,
                    size=size, aspect=aspect, dropna=dropna, **grid_kws)

    # Add the markers here as PairGrid has figured out how many levels of the
    # hue variable are needed and we don't want to duplicate that process
    if markers is not None:
        if grid.hue_names is None:
            n_markers = 1
        else:
            n_markers = len(grid.hue_names)
        if not isinstance(markers, list):
            markers = [markers] * n_markers
        if len(markers) != n_markers:
            raise ValueError(("markers must be a singleton or a list of "
                              "markers for each level of the hue variable"))
        grid.hue_kws = {"marker": markers}

    # Maybe plot on the diagonal
    if grid.square_grid:
        if diag_kind == "hist":
            grid.map_diag(plt.hist, **diag_kws)
        elif diag_kind == "kde":
            diag_kws["legend"] = False
            grid.map_diag(kdeplot, **diag_kws)

    # Maybe plot on the off-diagonals
    if grid.square_grid and diag_kind is not None:
        plotter = grid.map_offdiag
    else:
        plotter = grid.map

    if kind == "scatter":
        plot_kws.setdefault("edgecolor", "white")
        plotter(plt.scatter, **plot_kws)
    elif kind == "reg":
        from .regression import regplot  # Avoid circular import
        plotter(regplot, **plot_kws)

    # Add a legend
    if hue is not None:
        grid.add_legend()

    return grid


def jointplot(x, y, data=None, kind="scatter", stat_func=stats.pearsonr,
              color=None, size=6, ratio=5, space=.2,
              dropna=True, xlim=None, ylim=None,
              joint_kws=None, marginal_kws=None, annot_kws=None, **kwargs):
    """Draw a plot of two variables with bivariate and univariate graphs.

    This function provides a convenient interface to the :class:`JointGrid`
    class, with several canned plot kinds. This is intended to be a fairly
    lightweight wrapper; if you need more flexibility, you should use
    :class:`JointGrid` directly.

    Parameters
    ----------
    x, y : strings or vectors
        Data or names of variables in ``data``.
    data : DataFrame, optional
        DataFrame when ``x`` and ``y`` are variable names.
    kind : { "scatter" | "reg" | "resid" | "kde" | "hex" }, optional
        Kind of plot to draw.
    stat_func : callable or None, optional
        Function used to calculate a statistic about the relationship and
        annotate the plot. Should map `x` and `y` either to a single value
        or to a (value, p) tuple. Set to ``None`` if you don't want to
        annotate the plot.
    color : matplotlib color, optional
        Color used for the plot elements.
    size : numeric, optional
        Size of the figure (it will be square).
    ratio : numeric, optional
        Ratio of joint axes size to marginal axes height.
    space : numeric, optional
        Space between the joint and marginal axes
    dropna : bool, optional
        If True, remove observations that are missing from ``x`` and ``y``.
    {x, y}lim : two-tuples, optional
        Axis limits to set before plotting.
    {joint, marginal, annot}_kws : dicts, optional
        Additional keyword arguments for the plot components.
    kwargs : key, value pairings
        Additional keyword arguments are passed to the function used to
        draw the plot on the joint Axes, superseding items in the
        ``joint_kws`` dictionary.

    Returns
    -------
    grid : :class:`JointGrid`
        :class:`JointGrid` object with the plot on it.

    See Also
    --------
    JointGrid : The Grid class used for drawing this plot. Use it directly if
                you need more flexibility.

    Examples
    --------

    Draw a scatterplot with marginal histograms:

    .. plot::
        :context: close-figs

        >>> import numpy as np, pandas as pd; np.random.seed(0)
        >>> import seaborn as sns; sns.set(style="white", color_codes=True)
        >>> tips = sns.load_dataset("tips")
        >>> g = sns.jointplot(x="total_bill", y="tip", data=tips)

    Add regression and kernel density fits:

    .. plot::
        :context: close-figs

        >>> g = sns.jointplot("total_bill", "tip", data=tips, kind="reg")

    Replace the scatterplot with a joint histogram using hexagonal bins:

    .. plot::
        :context: close-figs

        >>> g = sns.jointplot("total_bill", "tip", data=tips, kind="hex")

    Replace the scatterplots and histograms with density estimates and align
    the marginal Axes tightly with the joint Axes:

    .. plot::
        :context: close-figs

        >>> iris = sns.load_dataset("iris")
        >>> g = sns.jointplot("sepal_width", "petal_length", data=iris,
        ...                   kind="kde", space=0, color="g")

    Use a different statistic for the annotation:

    .. plot::
        :context: close-figs

        >>> from scipy.stats import spearmanr
        >>> g = sns.jointplot("size", "total_bill", data=tips,
        ...                   stat_func=spearmanr, color="m")

    Draw a scatterplot, then add a joint density estimate:

    .. plot::
        :context: close-figs

        >>> g = (sns.jointplot("sepal_length", "sepal_width",
        ...                    data=iris, color="k")
        ...         .plot_joint(sns.kdeplot, zorder=0, n_levels=6))

    Pass vectors in directly without using Pandas, then name the axes:

    .. plot::
        :context: close-figs

        >>> x, y = np.random.randn(2, 300)
        >>> g = (sns.jointplot(x, y, kind="hex", stat_func=None)
        ...         .set_axis_labels("x", "y"))

    Draw a smaller figure with more space devoted to the marginal plots:

    .. plot::
        :context: close-figs

        >>> g = sns.jointplot("total_bill", "tip", data=tips,
        ...                   size=5, ratio=3, color="g")

    Pass keyword arguments down to the underlying plots:

    .. plot::
        :context: close-figs

        >>> g = sns.jointplot("petal_length", "sepal_length", data=iris,
        ...                   marginal_kws=dict(bins=15, rug=True),
        ...                   annot_kws=dict(stat="r"),
        ...                   s=40, edgecolor="w", linewidth=1)

    """
    # Set up empty default kwarg dicts
    if joint_kws is None:
        joint_kws = {}
    joint_kws.update(kwargs)
    if marginal_kws is None:
        marginal_kws = {}
    if annot_kws is None:
        annot_kws = {}

    # Make a colormap based off the plot color
    if color is None:
        color = color_palette()[0]
    color_rgb = mpl.colors.colorConverter.to_rgb(color)
    colors = [utils.set_hls_values(color_rgb, l=l)
              for l in np.linspace(1, 0, 12)]
    cmap = blend_palette(colors, as_cmap=True)

    # Initialize the JointGrid object
    grid = JointGrid(x, y, data, dropna=dropna,
                     size=size, ratio=ratio, space=space,
                     xlim=xlim, ylim=ylim)

    # Plot the data using the grid
    if kind == "scatter":

        joint_kws.setdefault("color", color)
        grid.plot_joint(plt.scatter, **joint_kws)

        marginal_kws.setdefault("kde", False)
        marginal_kws.setdefault("color", color)
        grid.plot_marginals(distplot, **marginal_kws)

    elif kind.startswith("hex"):

        x_bins = min(_freedman_diaconis_bins(grid.x), 50)
        y_bins = min(_freedman_diaconis_bins(grid.y), 50)
        gridsize = int(np.mean([x_bins, y_bins]))

        joint_kws.setdefault("gridsize", gridsize)
        joint_kws.setdefault("cmap", cmap)
        grid.plot_joint(plt.hexbin, **joint_kws)

        marginal_kws.setdefault("kde", False)
        marginal_kws.setdefault("color", color)
        grid.plot_marginals(distplot, **marginal_kws)

    elif kind.startswith("kde"):

        joint_kws.setdefault("shade", True)
        joint_kws.setdefault("cmap", cmap)
        grid.plot_joint(kdeplot, **joint_kws)

        marginal_kws.setdefault("shade", True)
        marginal_kws.setdefault("color", color)
        grid.plot_marginals(kdeplot, **marginal_kws)

    elif kind.startswith("reg"):

        from .regression import regplot

        marginal_kws.setdefault("color", color)
        grid.plot_marginals(distplot, **marginal_kws)

        joint_kws.setdefault("color", color)
        grid.plot_joint(regplot, **joint_kws)

    elif kind.startswith("resid"):

        from .regression import residplot

        joint_kws.setdefault("color", color)
        grid.plot_joint(residplot, **joint_kws)

        x, y = grid.ax_joint.collections[0].get_offsets().T
        marginal_kws.setdefault("color", color)
        marginal_kws.setdefault("kde", False)
        distplot(x, ax=grid.ax_marg_x, **marginal_kws)
        distplot(y, vertical=True, fit=stats.norm, ax=grid.ax_marg_y,
                 **marginal_kws)
        stat_func = None
    else:
        msg = "kind must be either 'scatter', 'reg', 'resid', 'kde', or 'hex'"
        raise ValueError(msg)

    if stat_func is not None:
        grid.annotate(stat_func, **annot_kws)

    return grid

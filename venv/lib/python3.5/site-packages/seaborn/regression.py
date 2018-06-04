"""Plotting functions for linear models (broadly construed)."""
from __future__ import division
import copy
import itertools
from textwrap import dedent
import numpy as np
import pandas as pd
from scipy.spatial import distance
import matplotlib as mpl
import matplotlib.pyplot as plt

import warnings

try:
    import statsmodels
    assert statsmodels
    _has_statsmodels = True
except ImportError:
    _has_statsmodels = False

from .external.six import string_types
from .external.six.moves import range

from . import utils
from . import algorithms as algo
from .palettes import color_palette
from .axisgrid import FacetGrid, _facet_docs


__all__ = ["lmplot", "regplot", "residplot"]


class _LinearPlotter(object):
    """Base class for plotting relational data in tidy format.

    To get anything useful done you'll have to inherit from this, but setup
    code that can be abstracted out should be put here.

    """
    def establish_variables(self, data, **kws):
        """Extract variables from data or use directly."""
        self.data = data

        # Validate the inputs
        any_strings = any([isinstance(v, string_types) for v in kws.values()])
        if any_strings and data is None:
            raise ValueError("Must pass `data` if using named variables.")

        # Set the variables
        for var, val in kws.items():
            if isinstance(val, string_types):
                setattr(self, var, data[val])
            else:
                setattr(self, var, val)

    def dropna(self, *vars):
        """Remove observations with missing data."""
        vals = [getattr(self, var) for var in vars]
        vals = [v for v in vals if v is not None]
        not_na = np.all(np.column_stack([pd.notnull(v) for v in vals]), axis=1)
        for var in vars:
            val = getattr(self, var)
            if val is not None:
                setattr(self, var, val[not_na])

    def plot(self, ax):
        raise NotImplementedError


class _RegressionPlotter(_LinearPlotter):
    """Plotter for numeric independent variables with regression model.

    This does the computations and drawing for the `regplot` function, and
    is thus also used indirectly by `lmplot`.
    """
    def __init__(self, x, y, data=None, x_estimator=None, x_bins=None,
                 x_ci="ci", scatter=True, fit_reg=True, ci=95, n_boot=1000,
                 units=None, order=1, logistic=False, lowess=False,
                 robust=False, logx=False, x_partial=None, y_partial=None,
                 truncate=False, dropna=True, x_jitter=None, y_jitter=None,
                 color=None, label=None):

        # Set member attributes
        self.x_estimator = x_estimator
        self.ci = ci
        self.x_ci = ci if x_ci == "ci" else x_ci
        self.n_boot = n_boot
        self.scatter = scatter
        self.fit_reg = fit_reg
        self.order = order
        self.logistic = logistic
        self.lowess = lowess
        self.robust = robust
        self.logx = logx
        self.truncate = truncate
        self.x_jitter = x_jitter
        self.y_jitter = y_jitter
        self.color = color
        self.label = label

        # Validate the regression options:
        if sum((order > 1, logistic, robust, lowess, logx)) > 1:
            raise ValueError("Mutually exclusive regression options.")

        # Extract the data vals from the arguments or passed dataframe
        self.establish_variables(data, x=x, y=y, units=units,
                                 x_partial=x_partial, y_partial=y_partial)

        # Drop null observations
        if dropna:
            self.dropna("x", "y", "units", "x_partial", "y_partial")

        # Regress nuisance variables out of the data
        if self.x_partial is not None:
            self.x = self.regress_out(self.x, self.x_partial)
        if self.y_partial is not None:
            self.y = self.regress_out(self.y, self.y_partial)

        # Possibly bin the predictor variable, which implies a point estimate
        if x_bins is not None:
            self.x_estimator = np.mean if x_estimator is None else x_estimator
            x_discrete, x_bins = self.bin_predictor(x_bins)
            self.x_discrete = x_discrete
        else:
            self.x_discrete = self.x

        # Save the range of the x variable for the grid later
        self.x_range = self.x.min(), self.x.max()

    @property
    def scatter_data(self):
        """Data where each observation is a point."""
        x_j = self.x_jitter
        if x_j is None:
            x = self.x
        else:
            x = self.x + np.random.uniform(-x_j, x_j, len(self.x))

        y_j = self.y_jitter
        if y_j is None:
            y = self.y
        else:
            y = self.y + np.random.uniform(-y_j, y_j, len(self.y))

        return x, y

    @property
    def estimate_data(self):
        """Data with a point estimate and CI for each discrete x value."""
        x, y = self.x_discrete, self.y
        vals = sorted(np.unique(x))
        points, cis = [], []

        for val in vals:

            # Get the point estimate of the y variable
            _y = y[x == val]
            est = self.x_estimator(_y)
            points.append(est)

            # Compute the confidence interval for this estimate
            if self.x_ci is None:
                cis.append(None)
            else:
                units = None
                if self.x_ci == "sd":
                    sd = np.std(_y)
                    _ci = est - sd, est + sd
                else:
                    if self.units is not None:
                        units = self.units[x == val]
                    boots = algo.bootstrap(_y, func=self.x_estimator,
                                           n_boot=self.n_boot, units=units)
                    _ci = utils.ci(boots, self.x_ci)
                cis.append(_ci)

        return vals, points, cis

    def fit_regression(self, ax=None, x_range=None, grid=None):
        """Fit the regression model."""
        # Create the grid for the regression
        if grid is None:
            if self.truncate:
                x_min, x_max = self.x_range
            else:
                if ax is None:
                    x_min, x_max = x_range
                else:
                    x_min, x_max = ax.get_xlim()
            grid = np.linspace(x_min, x_max, 100)
        ci = self.ci

        # Fit the regression
        if self.order > 1:
            yhat, yhat_boots = self.fit_poly(grid, self.order)
        elif self.logistic:
            from statsmodels.genmod.generalized_linear_model import GLM
            from statsmodels.genmod.families import Binomial
            yhat, yhat_boots = self.fit_statsmodels(grid, GLM,
                                                    family=Binomial())
        elif self.lowess:
            ci = None
            grid, yhat = self.fit_lowess()
        elif self.robust:
            from statsmodels.robust.robust_linear_model import RLM
            yhat, yhat_boots = self.fit_statsmodels(grid, RLM)
        elif self.logx:
            yhat, yhat_boots = self.fit_logx(grid)
        else:
            yhat, yhat_boots = self.fit_fast(grid)

        # Compute the confidence interval at each grid point
        if ci is None:
            err_bands = None
        else:
            err_bands = utils.ci(yhat_boots, ci, axis=0)

        return grid, yhat, err_bands

    def fit_fast(self, grid):
        """Low-level regression and prediction using linear algebra."""
        X, y = np.c_[np.ones(len(self.x)), self.x], self.y
        grid = np.c_[np.ones(len(grid)), grid]
        reg_func = lambda _x, _y: np.linalg.pinv(_x).dot(_y)
        yhat = grid.dot(reg_func(X, y))
        if self.ci is None:
            return yhat, None

        beta_boots = algo.bootstrap(X, y, func=reg_func,
                                    n_boot=self.n_boot, units=self.units).T
        yhat_boots = grid.dot(beta_boots).T
        return yhat, yhat_boots

    def fit_poly(self, grid, order):
        """Regression using numpy polyfit for higher-order trends."""
        x, y = self.x, self.y
        reg_func = lambda _x, _y: np.polyval(np.polyfit(_x, _y, order), grid)
        yhat = reg_func(x, y)
        if self.ci is None:
            return yhat, None

        yhat_boots = algo.bootstrap(x, y, func=reg_func,
                                    n_boot=self.n_boot, units=self.units)
        return yhat, yhat_boots

    def fit_statsmodels(self, grid, model, **kwargs):
        """More general regression function using statsmodels objects."""
        X, y = np.c_[np.ones(len(self.x)), self.x], self.y
        grid = np.c_[np.ones(len(grid)), grid]
        import statsmodels.genmod.generalized_linear_model as glm

        def reg_func(_x, _y):
            try:
                yhat = model(_y, _x, **kwargs).fit().predict(grid)
            except glm.PerfectSeparationError:
                yhat = np.empty(len(grid))
                yhat.fill(np.nan)
            return yhat

        yhat = reg_func(X, y)
        if self.ci is None:
            return yhat, None

        yhat_boots = algo.bootstrap(X, y, func=reg_func,
                                    n_boot=self.n_boot, units=self.units)
        return yhat, yhat_boots

    def fit_lowess(self):
        """Fit a locally-weighted regression, which returns its own grid."""
        from statsmodels.nonparametric.smoothers_lowess import lowess
        grid, yhat = lowess(self.y, self.x).T
        return grid, yhat

    def fit_logx(self, grid):
        """Fit the model in log-space."""
        X, y = np.c_[np.ones(len(self.x)), self.x], self.y
        grid = np.c_[np.ones(len(grid)), np.log(grid)]

        def reg_func(_x, _y):
            _x = np.c_[_x[:, 0], np.log(_x[:, 1])]
            return np.linalg.pinv(_x).dot(_y)

        yhat = grid.dot(reg_func(X, y))
        if self.ci is None:
            return yhat, None

        beta_boots = algo.bootstrap(X, y, func=reg_func,
                                    n_boot=self.n_boot, units=self.units).T
        yhat_boots = grid.dot(beta_boots).T
        return yhat, yhat_boots

    def bin_predictor(self, bins):
        """Discretize a predictor by assigning value to closest bin."""
        x = self.x
        if np.isscalar(bins):
            percentiles = np.linspace(0, 100, bins + 2)[1:-1]
            bins = np.c_[utils.percentiles(x, percentiles)]
        else:
            bins = np.c_[np.ravel(bins)]

        dist = distance.cdist(np.c_[x], bins)
        x_binned = bins[np.argmin(dist, axis=1)].ravel()

        return x_binned, bins.ravel()

    def regress_out(self, a, b):
        """Regress b from a keeping a's original mean."""
        a_mean = a.mean()
        a = a - a_mean
        b = b - b.mean()
        b = np.c_[b]
        a_prime = a - b.dot(np.linalg.pinv(b).dot(a))
        return (a_prime + a_mean).reshape(a.shape)

    def plot(self, ax, scatter_kws, line_kws):
        """Draw the full plot."""
        # Insert the plot label into the correct set of keyword arguments
        if self.scatter:
            scatter_kws["label"] = self.label
        else:
            line_kws["label"] = self.label

        # Use the current color cycle state as a default
        if self.color is None:
            lines, = plt.plot(self.x.mean(), self.y.mean())
            color = lines.get_color()
            lines.remove()
        else:
            color = self.color

        # Ensure that color is hex to avoid matplotlib weidness
        color = mpl.colors.rgb2hex(mpl.colors.colorConverter.to_rgb(color))

        # Let color in keyword arguments override overall plot color
        scatter_kws.setdefault("color", color)
        line_kws.setdefault("color", color)

        # Draw the constituent plots
        if self.scatter:
            self.scatterplot(ax, scatter_kws)
        if self.fit_reg:
            self.lineplot(ax, line_kws)

        # Label the axes
        if hasattr(self.x, "name"):
            ax.set_xlabel(self.x.name)
        if hasattr(self.y, "name"):
            ax.set_ylabel(self.y.name)

    def scatterplot(self, ax, kws):
        """Draw the data."""
        # Treat the line-based markers specially, explicitly setting larger
        # linewidth than is provided by the seaborn style defaults.
        # This would ideally be handled better in matplotlib (i.e., distinguish
        # between edgewidth for solid glyphs and linewidth for line glyphs
        # but this should do for now.
        line_markers = ["1", "2", "3", "4", "+", "x", "|", "_"]
        if self.x_estimator is None:
            if "marker" in kws and kws["marker"] in line_markers:
                lw = mpl.rcParams["lines.linewidth"]
            else:
                lw = mpl.rcParams["lines.markeredgewidth"]
            kws.setdefault("linewidths", lw)

            if not hasattr(kws['color'], 'shape') or kws['color'].shape[1] < 4:
                kws.setdefault("alpha", .8)

            x, y = self.scatter_data
            ax.scatter(x, y, **kws)
        else:
            # TODO abstraction
            ci_kws = {"color": kws["color"]}
            ci_kws["linewidth"] = mpl.rcParams["lines.linewidth"] * 1.75
            kws.setdefault("s", 50)

            xs, ys, cis = self.estimate_data
            if [ci for ci in cis if ci is not None]:
                for x, ci in zip(xs, cis):
                    ax.plot([x, x], ci, **ci_kws)
            ax.scatter(xs, ys, **kws)

    def lineplot(self, ax, kws):
        """Draw the model."""
        xlim = ax.get_xlim()

        # Fit the regression model
        grid, yhat, err_bands = self.fit_regression(ax)

        # Get set default aesthetics
        fill_color = kws["color"]
        lw = kws.pop("lw", mpl.rcParams["lines.linewidth"] * 1.5)
        kws.setdefault("linewidth", lw)

        # Draw the regression line and confidence interval
        ax.plot(grid, yhat, **kws)
        if err_bands is not None:
            ax.fill_between(grid, *err_bands, facecolor=fill_color, alpha=.15)
        ax.set_xlim(*xlim)


_regression_docs = dict(

    model_api=dedent("""\
    There are a number of mutually exclusive options for estimating the
    regression model: ``order``, ``logistic``, ``lowess``, ``robust``, and
    ``logx``. See the parameter docs for more information on these options.\
    """),

    regplot_vs_lmplot=dedent("""\
    Understanding the difference between :func:`regplot` and :func:`lmplot` can
    be a bit tricky. In fact, they are closely related, as :func:`lmplot` uses
    :func:`regplot` internally and takes most of its parameters. However,
    :func:`regplot` is an axes-level function, so it draws directly onto an
    axes (either the currently active axes or the one provided by the ``ax``
    parameter), while :func:`lmplot` is a figure-level function and creates its
    own figure, which is managed through a :class:`FacetGrid`. This has a few
    consequences, namely that :func:`regplot` can happily coexist in a figure
    with other kinds of plots and will follow the global matplotlib color
    cycle. In contrast, :func:`lmplot` needs to occupy an entire figure, and
    the size and color cycle are controlled through function parameters,
    ignoring the global defaults.\
    """),

    x_estimator=dedent("""\
    x_estimator : callable that maps vector -> scalar, optional
        Apply this function to each unique value of ``x`` and plot the
        resulting estimate. This is useful when ``x`` is a discrete variable.
        If ``x_ci`` is given, this estimate will be bootstrapped and a
        confidence interval will be drawn.\
    """),
    x_bins=dedent("""\
    x_bins : int or vector, optional
        Bin the ``x`` variable into discrete bins and then estimate the central
        tendency and a confidence interval. This binning only influences how
        the scatterplot is drawn; the regression is still fit to the original
        data.  This parameter is interpreted either as the number of
        evenly-sized (not necessary spaced) bins or the positions of the bin
        centers. When this parameter is used, it implies that the default of
        ``x_estimator`` is ``numpy.mean``.\
    """),
    x_ci=dedent("""\
    x_ci : "ci", "sd", int in [0, 100] or None, optional
        Size of the confidence interval used when plotting a central tendency
        for discrete values of ``x``. If ``"ci"``, defer to the value of the
        ``ci`` parameter. If ``"sd"``, skip bootstrappig and show the standard
        deviation of the observations in each bin.\
    """),
    scatter=dedent("""\
    scatter : bool, optional
        If ``True``, draw a scatterplot with the underlying observations (or
        the ``x_estimator`` values).\
    """),
    fit_reg=dedent("""\
    fit_reg : bool, optional
        If ``True``, estimate and plot a regression model relating the ``x``
        and ``y`` variables.\
    """),
    ci=dedent("""\
    ci : int in [0, 100] or None, optional
        Size of the confidence interval for the regression estimate. This will
        be drawn using translucent bands around the regression line. The
        confidence interval is estimated using a bootstrap; for large
        datasets, it may be advisable to avoid that computation by setting
        this parameter to None.\
    """),
    n_boot=dedent("""\
    n_boot : int, optional
        Number of bootstrap resamples used to estimate the ``ci``. The default
        value attempts to balance time and stability; you may want to increase
        this value for "final" versions of plots.\
    """),
    units=dedent("""\
    units : variable name in ``data``, optional
        If the ``x`` and ``y`` observations are nested within sampling units,
        those can be specified here. This will be taken into account when
        computing the confidence intervals by performing a multilevel bootstrap
        that resamples both units and observations (within unit). This does not
        otherwise influence how the regression is estimated or drawn.\
    """),
    order=dedent("""\
    order : int, optional
        If ``order`` is greater than 1, use ``numpy.polyfit`` to estimate a
        polynomial regression.\
    """),
    logistic=dedent("""\
    logistic : bool, optional
        If ``True``, assume that ``y`` is a binary variable and use
        ``statsmodels`` to estimate a logistic regression model. Note that this
        is substantially more computationally intensive than linear regression,
        so you may wish to decrease the number of bootstrap resamples
        (``n_boot``) or set ``ci`` to None.\
    """),
    lowess=dedent("""\
    lowess : bool, optional
        If ``True``, use ``statsmodels`` to estimate a nonparametric lowess
        model (locally weighted linear regression). Note that confidence
        intervals cannot currently be drawn for this kind of model.\
    """),
    robust=dedent("""\
    robust : bool, optional
        If ``True``, use ``statsmodels`` to estimate a robust regression. This
        will de-weight outliers. Note that this is substantially more
        computationally intensive than standard linear regression, so you may
        wish to decrease the number of bootstrap resamples (``n_boot``) or set
        ``ci`` to None.\
    """),
    logx=dedent("""\
    logx : bool, optional
        If ``True``, estimate a linear regression of the form y ~ log(x), but
        plot the scatterplot and regression model in the input space. Note that
        ``x`` must be positive for this to work.\
    """),
    xy_partial=dedent("""\
    {x,y}_partial : strings in ``data`` or matrices
        Confounding variables to regress out of the ``x`` or ``y`` variables
        before plotting.\
    """),
    truncate=dedent("""\
    truncate : bool, optional
        By default, the regression line is drawn to fill the x axis limits
        after the scatterplot is drawn. If ``truncate`` is ``True``, it will
        instead by bounded by the data limits.\
    """),
    xy_jitter=dedent("""\
    {x,y}_jitter : floats, optional
        Add uniform random noise of this size to either the ``x`` or ``y``
        variables. The noise is added to a copy of the data after fitting the
        regression, and only influences the look of the scatterplot. This can
        be helpful when plotting variables that take discrete values.\
    """),
    scatter_line_kws=dedent("""\
    {scatter,line}_kws : dictionaries
        Additional keyword arguments to pass to ``plt.scatter`` and
        ``plt.plot``.\
    """),
    )
_regression_docs.update(_facet_docs)


def lmplot(x, y, data, hue=None, col=None, row=None, palette=None,
           col_wrap=None, size=5, aspect=1, markers="o", sharex=True,
           sharey=True, hue_order=None, col_order=None, row_order=None,
           legend=True, legend_out=True, x_estimator=None, x_bins=None,
           x_ci="ci", scatter=True, fit_reg=True, ci=95, n_boot=1000,
           units=None, order=1, logistic=False, lowess=False, robust=False,
           logx=False, x_partial=None, y_partial=None, truncate=False,
           x_jitter=None, y_jitter=None, scatter_kws=None, line_kws=None):

    # Reduce the dataframe to only needed columns
    need_cols = [x, y, hue, col, row, units, x_partial, y_partial]
    cols = np.unique([a for a in need_cols if a is not None]).tolist()
    data = data[cols]

    # Initialize the grid
    facets = FacetGrid(data, row, col, hue, palette=palette,
                       row_order=row_order, col_order=col_order,
                       hue_order=hue_order, size=size, aspect=aspect,
                       col_wrap=col_wrap, sharex=sharex, sharey=sharey,
                       legend_out=legend_out)

    # Add the markers here as FacetGrid has figured out how many levels of the
    # hue variable are needed and we don't want to duplicate that process
    if facets.hue_names is None:
        n_markers = 1
    else:
        n_markers = len(facets.hue_names)
    if not isinstance(markers, list):
        markers = [markers] * n_markers
    if len(markers) != n_markers:
        raise ValueError(("markers must be a singeton or a list of markers "
                          "for each level of the hue variable"))
    facets.hue_kws = {"marker": markers}

    # Hack to set the x limits properly, which needs to happen here
    # because the extent of the regression estimate is determined
    # by the limits of the plot
    if sharex:
        for ax in facets.axes.flat:
            ax.scatter(data[x], np.ones(len(data)) * data[y].mean()).remove()

    # Draw the regression plot on each facet
    regplot_kws = dict(
        x_estimator=x_estimator, x_bins=x_bins, x_ci=x_ci,
        scatter=scatter, fit_reg=fit_reg, ci=ci, n_boot=n_boot, units=units,
        order=order, logistic=logistic, lowess=lowess, robust=robust,
        logx=logx, x_partial=x_partial, y_partial=y_partial, truncate=truncate,
        x_jitter=x_jitter, y_jitter=y_jitter,
        scatter_kws=scatter_kws, line_kws=line_kws,
        )
    facets.map_dataframe(regplot, x, y, **regplot_kws)

    # Add a legend
    if legend and (hue is not None) and (hue not in [col, row]):
        facets.add_legend()
    return facets


lmplot.__doc__ = dedent("""\
    Plot data and regression model fits across a FacetGrid.

    This function combines :func:`regplot` and :class:`FacetGrid`. It is
    intended as a convenient interface to fit regression models across
    conditional subsets of a dataset.

    When thinking about how to assign variables to different facets, a general
    rule is that it makes sense to use ``hue`` for the most important
    comparison, followed by ``col`` and ``row``. However, always think about
    your particular dataset and the goals of the visualization you are
    creating.

    {model_api}

    The parameters to this function span most of the options in
    :class:`FacetGrid`, although there may be occasional cases where you will
    want to use that class and :func:`regplot` directly.

    Parameters
    ----------
    x, y : strings, optional
        Input variables; these should be column names in ``data``.
    {data}
    hue, col, row : strings
        Variables that define subsets of the data, which will be drawn on
        separate facets in the grid. See the ``*_order`` parameters to control
        the order of levels of this variable.
    {palette}
    {col_wrap}
    {size}
    {aspect}
    markers : matplotlib marker code or list of marker codes, optional
        Markers for the scatterplot. If a list, each marker in the list will be
        used for each level of the ``hue`` variable.
    {share_xy}
    {{hue,col,row}}_order : lists, optional
        Order for the levels of the faceting variables. By default, this will
        be the order that the levels appear in ``data`` or, if the variables
        are pandas categoricals, the category order.
    legend : bool, optional
        If ``True`` and there is a ``hue`` variable, add a legend.
    {legend_out}
    {x_estimator}
    {x_bins}
    {x_ci}
    {scatter}
    {fit_reg}
    {ci}
    {n_boot}
    {units}
    {order}
    {logistic}
    {lowess}
    {robust}
    {logx}
    {xy_partial}
    {truncate}
    {xy_jitter}
    {scatter_line_kws}

    See Also
    --------
    regplot : Plot data and a conditional model fit.
    FacetGrid : Subplot grid for plotting conditional relationships.
    pairplot : Combine :func:`regplot` and :class:`PairGrid` (when used with
               ``kind="reg"``).

    Notes
    -----

    {regplot_vs_lmplot}

    Examples
    --------

    These examples focus on basic regression model plots to exhibit the
    various faceting options; see the :func:`regplot` docs for demonstrations
    of the other options for plotting the data and models. There are also
    other examples for how to manipulate plot using the returned object on
    the :class:`FacetGrid` docs.

    Plot a simple linear relationship between two variables:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set(color_codes=True)
        >>> tips = sns.load_dataset("tips")
        >>> g = sns.lmplot(x="total_bill", y="tip", data=tips)

    Condition on a third variable and plot the levels in different colors:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", hue="smoker", data=tips)

    Use different markers as well as colors so the plot will reproduce to
    black-and-white more easily:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", hue="smoker", data=tips,
        ...                markers=["o", "x"])

    Use a different color palette:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", hue="smoker", data=tips,
        ...                palette="Set1")

    Map ``hue`` levels to colors with a dictionary:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", hue="smoker", data=tips,
        ...                palette=dict(Yes="g", No="m"))

    Plot the levels of the third variable across different columns:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", col="smoker", data=tips)

    Change the size and aspect ratio of the facets:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="size", y="total_bill", hue="day", col="day",
        ...                data=tips, aspect=.4, x_jitter=.1)

    Wrap the levels of the column variable into multiple rows:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", col="day", hue="day",
        ...                data=tips, col_wrap=2, size=3)

    Condition on two variables to make a full grid:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", row="sex", col="time",
        ...                data=tips, size=3)

    Use methods on the returned :class:`FacetGrid` instance to further tweak
    the plot:

    .. plot::
        :context: close-figs

        >>> g = sns.lmplot(x="total_bill", y="tip", row="sex", col="time",
        ...                data=tips, size=3)
        >>> g = (g.set_axis_labels("Total bill (US Dollars)", "Tip")
        ...       .set(xlim=(0, 60), ylim=(0, 12),
        ...            xticks=[10, 30, 50], yticks=[2, 6, 10])
        ...       .fig.subplots_adjust(wspace=.02))



    """).format(**_regression_docs)


def regplot(x, y, data=None, x_estimator=None, x_bins=None, x_ci="ci",
            scatter=True, fit_reg=True, ci=95, n_boot=1000, units=None,
            order=1, logistic=False, lowess=False, robust=False,
            logx=False, x_partial=None, y_partial=None,
            truncate=False, dropna=True, x_jitter=None, y_jitter=None,
            label=None, color=None, marker="o",
            scatter_kws=None, line_kws=None, ax=None):

    plotter = _RegressionPlotter(x, y, data, x_estimator, x_bins, x_ci,
                                 scatter, fit_reg, ci, n_boot, units,
                                 order, logistic, lowess, robust, logx,
                                 x_partial, y_partial, truncate, dropna,
                                 x_jitter, y_jitter, color, label)

    if ax is None:
        ax = plt.gca()

    scatter_kws = {} if scatter_kws is None else copy.copy(scatter_kws)
    scatter_kws["marker"] = marker
    line_kws = {} if line_kws is None else copy.copy(line_kws)
    plotter.plot(ax, scatter_kws, line_kws)
    return ax

regplot.__doc__ = dedent("""\
    Plot data and a linear regression model fit.

    {model_api}

    Parameters
    ----------
    x, y: string, series, or vector array
        Input variables. If strings, these should correspond with column names
        in ``data``. When pandas objects are used, axes will be labeled with
        the series name.
    {data}
    {x_estimator}
    {x_bins}
    {x_ci}
    {scatter}
    {fit_reg}
    {ci}
    {n_boot}
    {units}
    {order}
    {logistic}
    {lowess}
    {robust}
    {logx}
    {xy_partial}
    {truncate}
    {xy_jitter}
    label : string
        Label to apply to ether the scatterplot or regression line (if
        ``scatter`` is ``False``) for use in a legend.
    color : matplotlib color
        Color to apply to all plot elements; will be superseded by colors
        passed in ``scatter_kws`` or ``line_kws``.
    marker : matplotlib marker code
        Marker to use for the scatterplot glyphs.
    {scatter_line_kws}
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Returns
    -------
    ax : matplotlib Axes
        The Axes object containing the plot.

    See Also
    --------
    lmplot : Combine :func:`regplot` and :class:`FacetGrid` to plot multiple
             linear relationships in a dataset.
    jointplot : Combine :func:`regplot` and :class:`JointGrid` (when used with
                ``kind="reg"``).
    pairplot : Combine :func:`regplot` and :class:`PairGrid` (when used with
               ``kind="reg"``).
    residplot : Plot the residuals of a linear regression model.

    Notes
    -----

    {regplot_vs_lmplot}


    It's also easy to combine combine :func:`regplot` and :class:`JointGrid` or
    :class:`PairGrid` through the :func:`jointplot` and :func:`pairplot`
    functions, although these do not directly accept all of :func:`regplot`'s
    parameters.

    Examples
    --------

    Plot the relationship between two variables in a DataFrame:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set(color_codes=True)
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.regplot(x="total_bill", y="tip", data=tips)

    Plot with two variables defined as numpy arrays; use a different color:

    .. plot::
        :context: close-figs

        >>> import numpy as np; np.random.seed(8)
        >>> mean, cov = [4, 6], [(1.5, .7), (.7, 1)]
        >>> x, y = np.random.multivariate_normal(mean, cov, 80).T
        >>> ax = sns.regplot(x=x, y=y, color="g")

    Plot with two variables defined as pandas Series; use a different marker:

    .. plot::
        :context: close-figs

        >>> import pandas as pd
        >>> x, y = pd.Series(x, name="x_var"), pd.Series(y, name="y_var")
        >>> ax = sns.regplot(x=x, y=y, marker="+")

    Use a 68% confidence interval, which corresponds with the standard error
    of the estimate:

    .. plot::
        :context: close-figs

        >>> ax = sns.regplot(x=x, y=y, ci=68)

    Plot with a discrete ``x`` variable and add some jitter:

    .. plot::
        :context: close-figs

        >>> ax = sns.regplot(x="size", y="total_bill", data=tips, x_jitter=.1)

    Plot with a discrete ``x`` variable showing means and confidence intervals
    for unique values:

    .. plot::
        :context: close-figs

        >>> ax = sns.regplot(x="size", y="total_bill", data=tips,
        ...                  x_estimator=np.mean)

    Plot with a continuous variable divided into discrete bins:

    .. plot::
        :context: close-figs

        >>> ax = sns.regplot(x=x, y=y, x_bins=4)

    Fit a higher-order polynomial regression and truncate the model prediction:

    .. plot::
        :context: close-figs

        >>> ans = sns.load_dataset("anscombe")
        >>> ax = sns.regplot(x="x", y="y", data=ans.loc[ans.dataset == "II"],
        ...                  scatter_kws={{"s": 80}},
        ...                  order=2, ci=None, truncate=True)

    Fit a robust regression and don't plot a confidence interval:

    .. plot::
        :context: close-figs

        >>> ax = sns.regplot(x="x", y="y", data=ans.loc[ans.dataset == "III"],
        ...                  scatter_kws={{"s": 80}},
        ...                  robust=True, ci=None)

    Fit a logistic regression; jitter the y variable and use fewer bootstrap
    iterations:

    .. plot::
        :context: close-figs

        >>> tips["big_tip"] = (tips.tip / tips.total_bill) > .175
        >>> ax = sns.regplot(x="total_bill", y="big_tip", data=tips,
        ...                  logistic=True, n_boot=500, y_jitter=.03)

    Fit the regression model using log(x) and truncate the model prediction:

    .. plot::
        :context: close-figs

        >>> ax = sns.regplot(x="size", y="total_bill", data=tips,
        ...                  x_estimator=np.mean, logx=True, truncate=True)

    """).format(**_regression_docs)


def residplot(x, y, data=None, lowess=False, x_partial=None, y_partial=None,
              order=1, robust=False, dropna=True, label=None, color=None,
              scatter_kws=None, line_kws=None, ax=None):
    """Plot the residuals of a linear regression.

    This function will regress y on x (possibly as a robust or polynomial
    regression) and then draw a scatterplot of the residuals. You can
    optionally fit a lowess smoother to the residual plot, which can
    help in determining if there is structure to the residuals.

    Parameters
    ----------
    x : vector or string
        Data or column name in `data` for the predictor variable.
    y : vector or string
        Data or column name in `data` for the response variable.
    data : DataFrame, optional
        DataFrame to use if `x` and `y` are column names.
    lowess : boolean, optional
        Fit a lowess smoother to the residual scatterplot.
    {x, y}_partial : matrix or string(s) , optional
        Matrix with same first dimension as `x`, or column name(s) in `data`.
        These variables are treated as confounding and are removed from
        the `x` or `y` variables before plotting.
    order : int, optional
        Order of the polynomial to fit when calculating the residuals.
    robust : boolean, optional
        Fit a robust linear regression when calculating the residuals.
    dropna : boolean, optional
        If True, ignore observations with missing data when fitting and
        plotting.
    label : string, optional
        Label that will be used in any plot legends.
    color : matplotlib color, optional
        Color to use for all elements of the plot.
    {scatter, line}_kws : dictionaries, optional
        Additional keyword arguments passed to scatter() and plot() for drawing
        the components of the plot.
    ax : matplotlib axis, optional
        Plot into this axis, otherwise grab the current axis or make a new
        one if not existing.

    Returns
    -------
    ax: matplotlib axes
        Axes with the regression plot.

    See Also
    --------
    regplot : Plot a simple linear regression model.
    jointplot (with kind="resid"): Draw a residplot with univariate
                                   marginal distrbutions.

    """
    plotter = _RegressionPlotter(x, y, data, ci=None,
                                 order=order, robust=robust,
                                 x_partial=x_partial, y_partial=y_partial,
                                 dropna=dropna, color=color, label=label)

    if ax is None:
        ax = plt.gca()

    # Calculate the residual from a linear regression
    _, yhat, _ = plotter.fit_regression(grid=plotter.x)
    plotter.y = plotter.y - yhat

    # Set the regression option on the plotter
    if lowess:
        plotter.lowess = True
    else:
        plotter.fit_reg = False

    # Plot a horizontal line at 0
    ax.axhline(0, ls=":", c=".2")

    # Draw the scatterplot
    scatter_kws = {} if scatter_kws is None else scatter_kws
    line_kws = {} if line_kws is None else line_kws
    plotter.plot(ax, scatter_kws, line_kws)
    return ax


def coefplot(formula, data, groupby=None, intercept=False, ci=95,
             palette="husl"):
    """Plot the coefficients from a linear model.

    Parameters
    ----------
    formula : string
        patsy formula for ols model
    data : dataframe
        data for the plot; formula terms must appear in columns
    groupby : grouping object, optional
        object to group data with to fit conditional models
    intercept : bool, optional
        if False, strips the intercept term before plotting
    ci : float, optional
        size of confidence intervals
    palette : seaborn color palette, optional
        palette for the horizonal plots

    """
    msg = (
        "The `coefplot` function has been deprecated and will be removed "
        "in a future version."
    )
    warnings.warn(msg, UserWarning)
    if not _has_statsmodels:
        raise ImportError("The `coefplot` function requires statsmodels")
    import statsmodels.formula.api as sf

    alpha = 1 - ci / 100
    if groupby is None:
        coefs = sf.ols(formula, data).fit().params
        cis = sf.ols(formula, data).fit().conf_int(alpha)
    else:
        grouped = data.groupby(groupby)
        coefs = grouped.apply(lambda d: sf.ols(formula, d).fit().params).T
        cis = grouped.apply(lambda d: sf.ols(formula, d).fit().conf_int(alpha))

    # Possibly ignore the intercept
    if not intercept:
        coefs = coefs.ix[1:]

    n_terms = len(coefs)

    # Plot seperately depending on groupby
    w, h = mpl.rcParams["figure.figsize"]
    hsize = lambda n: n * (h / 2)
    wsize = lambda n: n * (w / (4 * (n / 5)))
    if groupby is None:
        colors = itertools.cycle(color_palette(palette, n_terms))
        f, ax = plt.subplots(1, 1, figsize=(wsize(n_terms), hsize(1)))
        for i, term in enumerate(coefs.index):
            color = next(colors)
            low, high = cis.ix[term]
            ax.plot([i, i], [low, high], c=color,
                    solid_capstyle="round", lw=2.5)
            ax.plot(i, coefs.ix[term], "o", c=color, ms=8)
        ax.set_xlim(-.5, n_terms - .5)
        ax.axhline(0, ls="--", c="dimgray")
        ax.set_xticks(range(n_terms))
        ax.set_xticklabels(coefs.index)

    else:
        n_groups = len(coefs.columns)
        f, axes = plt.subplots(n_terms, 1, sharex=True,
                               figsize=(wsize(n_groups), hsize(n_terms)))
        if n_terms == 1:
            axes = [axes]
        colors = itertools.cycle(color_palette(palette, n_groups))
        for ax, term in zip(axes, coefs.index):
            for i, group in enumerate(coefs.columns):
                color = next(colors)
                low, high = cis.ix[(group, term)]
                ax.plot([i, i], [low, high], c=color,
                        solid_capstyle="round", lw=2.5)
                ax.plot(i, coefs.loc[term, group], "o", c=color, ms=8)
            ax.set_xlim(-.5, n_groups - .5)
            ax.axhline(0, ls="--", c="dimgray")
            ax.set_title(term)
        ax.set_xlabel(groupby)
        ax.set_xticks(range(n_groups))
        ax.set_xticklabels(coefs.columns)


def interactplot(x1, x2, y, data=None, filled=False, cmap="RdBu_r",
                 colorbar=True, levels=30, logistic=False,
                 contour_kws=None, scatter_kws=None, ax=None, **kwargs):
    """Visualize a continuous two-way interaction with a contour plot.

    Parameters
    ----------
    x1, x2, y, strings or array-like
        Either the two independent variables and the dependent variable,
        or keys to extract them from `data`
    data : DataFrame
        Pandas DataFrame with the data in the columns.
    filled : bool
        Whether to plot with filled or unfilled contours
    cmap : matplotlib colormap
        Colormap to represent yhat in the countour plot.
    colorbar : bool
        Whether to draw the colorbar for interpreting the color values.
    levels : int or sequence
        Number or position of contour plot levels.
    logistic : bool
        Fit a logistic regression model instead of linear regression.
    contour_kws : dictionary
        Keyword arguments for contour[f]().
    scatter_kws : dictionary
        Keyword arguments for plot().
    ax : matplotlib axis
        Axis to draw plot in.

    Returns
    -------
    ax : Matplotlib axis
        Axis with the contour plot.

    """
    msg = (
        "The `interactplot` function has been deprecated and will be removed "
        "in a future version."
    )
    warnings.warn(msg, UserWarning)
    if not _has_statsmodels:
        raise ImportError("The `interactplot` function requires statsmodels")
    from statsmodels.regression.linear_model import OLS
    from statsmodels.genmod.generalized_linear_model import GLM
    from statsmodels.genmod.families import Binomial

    # Handle the form of the data
    if data is not None:
        x1 = data[x1]
        x2 = data[x2]
        y = data[y]
    if hasattr(x1, "name"):
        xlabel = x1.name
    else:
        xlabel = None
    if hasattr(x2, "name"):
        ylabel = x2.name
    else:
        ylabel = None
    if hasattr(y, "name"):
        clabel = y.name
    else:
        clabel = None
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    y = np.asarray(y)

    # Initialize the scatter keyword dictionary
    if scatter_kws is None:
        scatter_kws = {}
    if not ("color" in scatter_kws or "c" in scatter_kws):
        scatter_kws["color"] = "#222222"
    if "alpha" not in scatter_kws:
        scatter_kws["alpha"] = 0.75

    # Intialize the contour keyword dictionary
    if contour_kws is None:
        contour_kws = {}

    # Initialize the axis
    if ax is None:
        ax = plt.gca()

    # Plot once to let matplotlib sort out the axis limits
    ax.plot(x1, x2, "o", **scatter_kws)

    # Find the plot limits
    x1min, x1max = ax.get_xlim()
    x2min, x2max = ax.get_ylim()

    # Make the grid for the contour plot
    x1_points = np.linspace(x1min, x1max, 100)
    x2_points = np.linspace(x2min, x2max, 100)
    xx1, xx2 = np.meshgrid(x1_points, x2_points)

    # Fit the model with an interaction
    X = np.c_[np.ones(x1.size), x1, x2, x1 * x2]
    if logistic:
        lm = GLM(y, X, family=Binomial()).fit()
    else:
        lm = OLS(y, X).fit()

    # Evaluate the model on the grid
    eval = np.vectorize(lambda x1_, x2_: lm.predict([1, x1_, x2_, x1_ * x2_]))
    yhat = eval(xx1, xx2)

    # Default color limits put the midpoint at mean(y)
    y_bar = y.mean()
    c_min = min(np.percentile(y, 2), yhat.min())
    c_max = max(np.percentile(y, 98), yhat.max())
    delta = max(c_max - y_bar, y_bar - c_min)
    c_min, cmax = y_bar - delta, y_bar + delta
    contour_kws.setdefault("vmin", c_min)
    contour_kws.setdefault("vmax", c_max)

    # Draw the contour plot
    func_name = "contourf" if filled else "contour"
    contour = getattr(ax, func_name)
    c = contour(xx1, xx2, yhat, levels, cmap=cmap, **contour_kws)

    # Draw the scatter again so it's visible
    ax.plot(x1, x2, "o", **scatter_kws)

    # Draw a colorbar, maybe
    if colorbar:
        bar = plt.colorbar(c)

    # Label the axes
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if clabel is not None and colorbar:
        clabel = "P(%s)" % clabel if logistic else clabel
        bar.set_label(clabel, labelpad=15, rotation=270)

    return ax

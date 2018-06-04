from __future__ import division
from textwrap import dedent
import colorsys
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib as mpl
from matplotlib.collections import PatchCollection
import matplotlib.patches as Patches
import matplotlib.pyplot as plt
import warnings

from .external.six import string_types
from .external.six.moves import range

from . import utils
from .utils import iqr, categorical_order, remove_na
from .algorithms import bootstrap
from .palettes import color_palette, husl_palette, light_palette, dark_palette
from .axisgrid import FacetGrid, _facet_docs


__all__ = ["boxplot", "violinplot", "stripplot", "swarmplot", "lvplot",
           "pointplot", "barplot", "countplot", "factorplot"]


class _CategoricalPlotter(object):

    width = .8
    default_palette = "light"

    def establish_variables(self, x=None, y=None, hue=None, data=None,
                            orient=None, order=None, hue_order=None,
                            units=None):
        """Convert input specification into a common representation."""
        # Option 1:
        # We are plotting a wide-form dataset
        # -----------------------------------
        if x is None and y is None:

            # Do a sanity check on the inputs
            if hue is not None:
                error = "Cannot use `hue` without `x` or `y`"
                raise ValueError(error)

            # No hue grouping with wide inputs
            plot_hues = None
            hue_title = None
            hue_names = None

            # No statistical units with wide inputs
            plot_units = None

            # We also won't get a axes labels here
            value_label = None
            group_label = None

            # Option 1a:
            # The input data is a Pandas DataFrame
            # ------------------------------------

            if isinstance(data, pd.DataFrame):

                # Order the data correctly
                if order is None:
                    order = []
                    # Reduce to just numeric columns
                    for col in data:
                        try:
                            data[col].astype(np.float)
                            order.append(col)
                        except ValueError:
                            pass
                plot_data = data[order]
                group_names = order
                group_label = data.columns.name

                # Convert to a list of arrays, the common representation
                iter_data = plot_data.iteritems()
                plot_data = [np.asarray(s, np.float) for k, s in iter_data]

            # Option 1b:
            # The input data is an array or list
            # ----------------------------------

            else:

                # We can't reorder the data
                if order is not None:
                    error = "Input data must be a pandas object to reorder"
                    raise ValueError(error)

                # The input data is an array
                if hasattr(data, "shape"):
                    if len(data.shape) == 1:
                        if np.isscalar(data[0]):
                            plot_data = [data]
                        else:
                            plot_data = list(data)
                    elif len(data.shape) == 2:
                        nr, nc = data.shape
                        if nr == 1 or nc == 1:
                            plot_data = [data.ravel()]
                        else:
                            plot_data = [data[:, i] for i in range(nc)]
                    else:
                        error = ("Input `data` can have no "
                                 "more than 2 dimensions")
                        raise ValueError(error)

                # Check if `data` is None to let us bail out here (for testing)
                elif data is None:
                    plot_data = [[]]

                # The input data is a flat list
                elif np.isscalar(data[0]):
                    plot_data = [data]

                # The input data is a nested list
                # This will catch some things that might fail later
                # but exhaustive checks are hard
                else:
                    plot_data = data

                # Convert to a list of arrays, the common representation
                plot_data = [np.asarray(d, np.float) for d in plot_data]

                # The group names will just be numeric indices
                group_names = list(range((len(plot_data))))

            # Figure out the plotting orientation
            orient = "h" if str(orient).startswith("h") else "v"

        # Option 2:
        # We are plotting a long-form dataset
        # -----------------------------------

        else:

            # See if we need to get variables from `data`
            if data is not None:
                x = data.get(x, x)
                y = data.get(y, y)
                hue = data.get(hue, hue)
                units = data.get(units, units)

            # Validate the inputs
            for input in [x, y, hue, units]:
                if isinstance(input, string_types):
                    err = "Could not interpret input '{}'".format(input)
                    raise ValueError(err)

            # Figure out the plotting orientation
            orient = self.infer_orient(x, y, orient)

            # Option 2a:
            # We are plotting a single set of data
            # ------------------------------------
            if x is None or y is None:

                # Determine where the data are
                vals = y if x is None else x

                # Put them into the common representation
                plot_data = [np.asarray(vals)]

                # Get a label for the value axis
                if hasattr(vals, "name"):
                    value_label = vals.name
                else:
                    value_label = None

                # This plot will not have group labels or hue nesting
                groups = None
                group_label = None
                group_names = []
                plot_hues = None
                hue_names = None
                hue_title = None
                plot_units = None

            # Option 2b:
            # We are grouping the data values by another variable
            # ---------------------------------------------------
            else:

                # Determine which role each variable will play
                if orient == "v":
                    vals, groups = y, x
                else:
                    vals, groups = x, y

                # Get the categorical axis label
                group_label = None
                if hasattr(groups, "name"):
                    group_label = groups.name

                # Get the order on the categorical axis
                group_names = categorical_order(groups, order)

                # Group the numeric data
                plot_data, value_label = self._group_longform(vals, groups,
                                                              group_names)

                # Now handle the hue levels for nested ordering
                if hue is None:
                    plot_hues = None
                    hue_title = None
                    hue_names = None
                else:

                    # Get the order of the hue levels
                    hue_names = categorical_order(hue, hue_order)

                    # Group the hue data
                    plot_hues, hue_title = self._group_longform(hue, groups,
                                                                group_names)

                # Now handle the units for nested observations
                if units is None:
                    plot_units = None
                else:
                    plot_units, _ = self._group_longform(units, groups,
                                                         group_names)

        # Assign object attributes
        # ------------------------
        self.orient = orient
        self.plot_data = plot_data
        self.group_label = group_label
        self.value_label = value_label
        self.group_names = group_names
        self.plot_hues = plot_hues
        self.hue_title = hue_title
        self.hue_names = hue_names
        self.plot_units = plot_units

    def _group_longform(self, vals, grouper, order):
        """Group a long-form variable by another with correct order."""
        # Ensure that the groupby will work
        if not isinstance(vals, pd.Series):
            vals = pd.Series(vals)

        # Group the val data
        grouped_vals = vals.groupby(grouper)
        out_data = []
        for g in order:
            try:
                g_vals = np.asarray(grouped_vals.get_group(g))
            except KeyError:
                g_vals = np.array([])
            out_data.append(g_vals)

        # Get the vals axis label
        label = vals.name

        return out_data, label

    def establish_colors(self, color, palette, saturation):
        """Get a list of colors for the main component of the plots."""
        if self.hue_names is None:
            n_colors = len(self.plot_data)
        else:
            n_colors = len(self.hue_names)

        # Determine the main colors
        if color is None and palette is None:
            # Determine whether the current palette will have enough values
            # If not, we'll default to the husl palette so each is distinct
            current_palette = utils.get_color_cycle()
            if n_colors <= len(current_palette):
                colors = color_palette(n_colors=n_colors)
            else:
                colors = husl_palette(n_colors, l=.7)

        elif palette is None:
            # When passing a specific color, the interpretation depends
            # on whether there is a hue variable or not.
            # If so, we will make a blend palette so that the different
            # levels have some amount of variation.
            if self.hue_names is None:
                colors = [color] * n_colors
            else:
                if self.default_palette == "light":
                    colors = light_palette(color, n_colors)
                elif self.default_palette == "dark":
                    colors = dark_palette(color, n_colors)
                else:
                    raise RuntimeError("No default palette specified")
        else:

            # Let `palette` be a dict mapping level to color
            if isinstance(palette, dict):
                if self.hue_names is None:
                    levels = self.group_names
                else:
                    levels = self.hue_names
                palette = [palette[l] for l in levels]

            colors = color_palette(palette, n_colors)

        # Desaturate a bit because these are patches
        if saturation < 1:
            colors = color_palette(colors, desat=saturation)

        # Conver the colors to a common representations
        rgb_colors = color_palette(colors)

        # Determine the gray color to use for the lines framing the plot
        light_vals = [colorsys.rgb_to_hls(*c)[1] for c in rgb_colors]
        l = min(light_vals) * .6
        gray = mpl.colors.rgb2hex((l, l, l))

        # Assign object attributes
        self.colors = rgb_colors
        self.gray = gray

    def infer_orient(self, x, y, orient=None):
        """Determine how the plot should be oriented based on the data."""
        orient = str(orient)

        def is_categorical(s):
            try:
                # Correct way, but doesnt exist in older Pandas
                try:
                    return pd.api.types.is_categorical_dtype(s)
                except AttributeError:
                    return pd.core.common.is_categorical_dtype(s)
            except AttributeError:
                # Also works, but feels hackier
                return str(s.dtype) == "categorical"

        def is_not_numeric(s):
            try:
                np.asarray(s, dtype=np.float)
            except ValueError:
                return True
            return False

        no_numeric = "Neither the `x` nor `y` variable appears to be numeric."

        if orient.startswith("v"):
            return "v"
        elif orient.startswith("h"):
            return "h"
        elif x is None:
            return "v"
        elif y is None:
            return "h"
        elif is_categorical(y):
            if is_categorical(x):
                raise ValueError(no_numeric)
            else:
                return "h"
        elif is_not_numeric(y):
            if is_not_numeric(x):
                raise ValueError(no_numeric)
            else:
                return "h"
        else:
            return "v"

    @property
    def hue_offsets(self):
        """A list of center positions for plots when hue nesting is used."""
        n_levels = len(self.hue_names)
        if self.dodge:
            each_width = self.width / n_levels
            offsets = np.linspace(0, self.width - each_width, n_levels)
            offsets -= offsets.mean()
        else:
            offsets = np.zeros(n_levels)

        return offsets

    @property
    def nested_width(self):
        """A float with the width of plot elements when hue nesting is used."""
        if self.dodge:
            width = self.width / len(self.hue_names) * .98
        else:
            width = self.width
        return width

    def annotate_axes(self, ax):
        """Add descriptive labels to an Axes object."""
        if self.orient == "v":
            xlabel, ylabel = self.group_label, self.value_label
        else:
            xlabel, ylabel = self.value_label, self.group_label

        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)

        if self.orient == "v":
            ax.set_xticks(np.arange(len(self.plot_data)))
            ax.set_xticklabels(self.group_names)
        else:
            ax.set_yticks(np.arange(len(self.plot_data)))
            ax.set_yticklabels(self.group_names)

        if self.orient == "v":
            ax.xaxis.grid(False)
            ax.set_xlim(-.5, len(self.plot_data) - .5)
        else:
            ax.yaxis.grid(False)
            ax.set_ylim(-.5, len(self.plot_data) - .5)

        if self.hue_names is not None:
            leg = ax.legend(loc="best")
            if self.hue_title is not None:
                leg.set_title(self.hue_title)

                # Set the title size a roundabout way to maintain
                # compatability with matplotlib 1.1
                try:
                    title_size = mpl.rcParams["axes.labelsize"] * .85
                except TypeError:  # labelsize is something like "large"
                    title_size = mpl.rcParams["axes.labelsize"]
                prop = mpl.font_manager.FontProperties(size=title_size)
                leg._legend_title_box._text.set_font_properties(prop)

    def add_legend_data(self, ax, color, label):
        """Add a dummy patch object so we can get legend data."""
        rect = plt.Rectangle([0, 0], 0, 0,
                             linewidth=self.linewidth / 2,
                             edgecolor=self.gray,
                             facecolor=color,
                             label=label)
        ax.add_patch(rect)


class _BoxPlotter(_CategoricalPlotter):

    def __init__(self, x, y, hue, data, order, hue_order,
                 orient, color, palette, saturation,
                 width, dodge, fliersize, linewidth):

        self.establish_variables(x, y, hue, data, orient, order, hue_order)
        self.establish_colors(color, palette, saturation)

        self.dodge = dodge
        self.width = width
        self.fliersize = fliersize

        if linewidth is None:
            linewidth = mpl.rcParams["lines.linewidth"]
        self.linewidth = linewidth

    def draw_boxplot(self, ax, kws):
        """Use matplotlib to draw a boxplot on an Axes."""
        vert = self.orient == "v"

        props = {}
        for obj in ["box", "whisker", "cap", "median", "flier"]:
            props[obj] = kws.pop(obj + "props", {})

        for i, group_data in enumerate(self.plot_data):

            if self.plot_hues is None:

                # Handle case where there is data at this level
                if group_data.size == 0:
                    continue

                # Draw a single box or a set of boxes
                # with a single level of grouping
                box_data = remove_na(group_data)

                # Handle case where there is no non-null data
                if box_data.size == 0:
                    continue

                artist_dict = ax.boxplot(box_data,
                                         vert=vert,
                                         patch_artist=True,
                                         positions=[i],
                                         widths=self.width,
                                         **kws)
                color = self.colors[i]
                self.restyle_boxplot(artist_dict, color, props)
            else:
                # Draw nested groups of boxes
                offsets = self.hue_offsets
                for j, hue_level in enumerate(self.hue_names):

                    # Add a legend for this hue level
                    if not i:
                        self.add_legend_data(ax, self.colors[j], hue_level)

                    # Handle case where there is data at this level
                    if group_data.size == 0:
                        continue

                    hue_mask = self.plot_hues[i] == hue_level
                    box_data = remove_na(group_data[hue_mask])

                    # Handle case where there is no non-null data
                    if box_data.size == 0:
                        continue

                    center = i + offsets[j]
                    artist_dict = ax.boxplot(box_data,
                                             vert=vert,
                                             patch_artist=True,
                                             positions=[center],
                                             widths=self.nested_width,
                                             **kws)
                    self.restyle_boxplot(artist_dict, self.colors[j], props)
                    # Add legend data, but just for one set of boxes

    def restyle_boxplot(self, artist_dict, color, props):
        """Take a drawn matplotlib boxplot and make it look nice."""
        for box in artist_dict["boxes"]:
            box.update(dict(facecolor=color,
                            zorder=.9,
                            edgecolor=self.gray,
                            linewidth=self.linewidth))
            box.update(props["box"])
        for whisk in artist_dict["whiskers"]:
            whisk.update(dict(color=self.gray,
                              linewidth=self.linewidth,
                              linestyle="-"))
            whisk.update(props["whisker"])
        for cap in artist_dict["caps"]:
            cap.update(dict(color=self.gray,
                            linewidth=self.linewidth))
            cap.update(props["cap"])
        for med in artist_dict["medians"]:
            med.update(dict(color=self.gray,
                            linewidth=self.linewidth))
            med.update(props["median"])
        for fly in artist_dict["fliers"]:
            fly.update(dict(markerfacecolor=self.gray,
                            marker="d",
                            markeredgecolor=self.gray,
                            markersize=self.fliersize))
            fly.update(props["flier"])

    def plot(self, ax, boxplot_kws):
        """Make the plot."""
        self.draw_boxplot(ax, boxplot_kws)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()


class _ViolinPlotter(_CategoricalPlotter):

    def __init__(self, x, y, hue, data, order, hue_order,
                 bw, cut, scale, scale_hue, gridsize,
                 width, inner, split, dodge, orient, linewidth,
                 color, palette, saturation):

        self.establish_variables(x, y, hue, data, orient, order, hue_order)
        self.establish_colors(color, palette, saturation)
        self.estimate_densities(bw, cut, scale, scale_hue, gridsize)

        self.gridsize = gridsize
        self.width = width
        self.dodge = dodge

        if inner is not None:
            if not any([inner.startswith("quart"),
                        inner.startswith("box"),
                        inner.startswith("stick"),
                        inner.startswith("point")]):
                err = "Inner style '{}' not recognized".format(inner)
                raise ValueError(err)
        self.inner = inner

        if split and self.hue_names is not None and len(self.hue_names) != 2:
            msg = "There must be exactly two hue levels to use `split`.'"
            raise ValueError(msg)
        self.split = split

        if linewidth is None:
            linewidth = mpl.rcParams["lines.linewidth"]
        self.linewidth = linewidth

    def estimate_densities(self, bw, cut, scale, scale_hue, gridsize):
        """Find the support and density for all of the data."""
        # Initialize data structures to keep track of plotting data
        if self.hue_names is None:
            support = []
            density = []
            counts = np.zeros(len(self.plot_data))
            max_density = np.zeros(len(self.plot_data))
        else:
            support = [[] for _ in self.plot_data]
            density = [[] for _ in self.plot_data]
            size = len(self.group_names), len(self.hue_names)
            counts = np.zeros(size)
            max_density = np.zeros(size)

        for i, group_data in enumerate(self.plot_data):

            # Option 1: we have a single level of grouping
            # --------------------------------------------

            if self.plot_hues is None:

                # Strip missing datapoints
                kde_data = remove_na(group_data)

                # Handle special case of no data at this level
                if kde_data.size == 0:
                    support.append(np.array([]))
                    density.append(np.array([1.]))
                    counts[i] = 0
                    max_density[i] = 0
                    continue

                # Handle special case of a single unique datapoint
                elif np.unique(kde_data).size == 1:
                    support.append(np.unique(kde_data))
                    density.append(np.array([1.]))
                    counts[i] = 1
                    max_density[i] = 0
                    continue

                # Fit the KDE and get the used bandwidth size
                kde, bw_used = self.fit_kde(kde_data, bw)

                # Determine the support grid and get the density over it
                support_i = self.kde_support(kde_data, bw_used, cut, gridsize)
                density_i = kde.evaluate(support_i)

                # Update the data structures with these results
                support.append(support_i)
                density.append(density_i)
                counts[i] = kde_data.size
                max_density[i] = density_i.max()

            # Option 2: we have nested grouping by a hue variable
            # ---------------------------------------------------

            else:
                for j, hue_level in enumerate(self.hue_names):

                    # Handle special case of no data at this category level
                    if not group_data.size:
                        support[i].append(np.array([]))
                        density[i].append(np.array([1.]))
                        counts[i, j] = 0
                        max_density[i, j] = 0
                        continue

                    # Select out the observations for this hue level
                    hue_mask = self.plot_hues[i] == hue_level

                    # Strip missing datapoints
                    kde_data = remove_na(group_data[hue_mask])

                    # Handle special case of no data at this level
                    if kde_data.size == 0:
                        support[i].append(np.array([]))
                        density[i].append(np.array([1.]))
                        counts[i, j] = 0
                        max_density[i, j] = 0
                        continue

                    # Handle special case of a single unique datapoint
                    elif np.unique(kde_data).size == 1:
                        support[i].append(np.unique(kde_data))
                        density[i].append(np.array([1.]))
                        counts[i, j] = 1
                        max_density[i, j] = 0
                        continue

                    # Fit the KDE and get the used bandwidth size
                    kde, bw_used = self.fit_kde(kde_data, bw)

                    # Determine the support grid and get the density over it
                    support_ij = self.kde_support(kde_data, bw_used,
                                                  cut, gridsize)
                    density_ij = kde.evaluate(support_ij)

                    # Update the data structures with these results
                    support[i].append(support_ij)
                    density[i].append(density_ij)
                    counts[i, j] = kde_data.size
                    max_density[i, j] = density_ij.max()

        # Scale the height of the density curve.
        # For a violinplot the density is non-quantitative.
        # The objective here is to scale the curves relative to 1 so that
        # they can be multiplied by the width parameter during plotting.

        if scale == "area":
            self.scale_area(density, max_density, scale_hue)

        elif scale == "width":
            self.scale_width(density)

        elif scale == "count":
            self.scale_count(density, counts, scale_hue)

        else:
            raise ValueError("scale method '{}' not recognized".format(scale))

        # Set object attributes that will be used while plotting
        self.support = support
        self.density = density

    def fit_kde(self, x, bw):
        """Estimate a KDE for a vector of data with flexible bandwidth."""
        # Allow for the use of old scipy where `bw` is fixed
        try:
            kde = stats.gaussian_kde(x, bw)
        except TypeError:
            kde = stats.gaussian_kde(x)
            if bw != "scott":  # scipy default
                msg = ("Ignoring bandwidth choice, "
                       "please upgrade scipy to use a different bandwidth.")
                warnings.warn(msg, UserWarning)

        # Extract the numeric bandwidth from the KDE object
        bw_used = kde.factor

        # At this point, bw will be a numeric scale factor.
        # To get the actual bandwidth of the kernel, we multiple by the
        # unbiased standard deviation of the data, which we will use
        # elsewhere to compute the range of the support.
        bw_used = bw_used * x.std(ddof=1)

        return kde, bw_used

    def kde_support(self, x, bw, cut, gridsize):
        """Define a grid of support for the violin."""
        support_min = x.min() - bw * cut
        support_max = x.max() + bw * cut
        return np.linspace(support_min, support_max, gridsize)

    def scale_area(self, density, max_density, scale_hue):
        """Scale the relative area under the KDE curve.

        This essentially preserves the "standard" KDE scaling, but the
        resulting maximum density will be 1 so that the curve can be
        properly multiplied by the violin width.

        """
        if self.hue_names is None:
            for d in density:
                if d.size > 1:
                    d /= max_density.max()
        else:
            for i, group in enumerate(density):
                for d in group:
                    if scale_hue:
                        max = max_density[i].max()
                    else:
                        max = max_density.max()
                    if d.size > 1:
                        d /= max

    def scale_width(self, density):
        """Scale each density curve to the same height."""
        if self.hue_names is None:
            for d in density:
                d /= d.max()
        else:
            for group in density:
                for d in group:
                    d /= d.max()

    def scale_count(self, density, counts, scale_hue):
        """Scale each density curve by the number of observations."""
        if self.hue_names is None:
            for count, d in zip(counts, density):
                d /= d.max()
                d *= count / counts.max()
        else:
            for i, group in enumerate(density):
                for j, d in enumerate(group):
                    count = counts[i, j]
                    if scale_hue:
                        scaler = count / counts[i].max()
                    else:
                        scaler = count / counts.max()
                    d /= d.max()
                    d *= scaler

    @property
    def dwidth(self):

        if self.hue_names is None or not self.dodge:
            return self.width / 2
        elif self.split:
            return self.width / 2
        else:
            return self.width / (2 * len(self.hue_names))

    def draw_violins(self, ax):
        """Draw the violins onto `ax`."""
        fill_func = ax.fill_betweenx if self.orient == "v" else ax.fill_between
        for i, group_data in enumerate(self.plot_data):

            kws = dict(edgecolor=self.gray, linewidth=self.linewidth)

            # Option 1: we have a single level of grouping
            # --------------------------------------------

            if self.plot_hues is None:

                support, density = self.support[i], self.density[i]

                # Handle special case of no observations in this bin
                if support.size == 0:
                    continue

                # Handle special case of a single observation
                elif support.size == 1:
                    val = np.asscalar(support)
                    d = np.asscalar(density)
                    self.draw_single_observation(ax, i, val, d)
                    continue

                # Draw the violin for this group
                grid = np.ones(self.gridsize) * i
                fill_func(support,
                          grid - density * self.dwidth,
                          grid + density * self.dwidth,
                          facecolor=self.colors[i],
                          **kws)

                # Draw the interior representation of the data
                if self.inner is None:
                    continue

                # Get a nan-free vector of datapoints
                violin_data = remove_na(group_data)

                # Draw box and whisker information
                if self.inner.startswith("box"):
                    self.draw_box_lines(ax, violin_data, support, density, i)

                # Draw quartile lines
                elif self.inner.startswith("quart"):
                    self.draw_quartiles(ax, violin_data, support, density, i)

                # Draw stick observations
                elif self.inner.startswith("stick"):
                    self.draw_stick_lines(ax, violin_data, support, density, i)

                # Draw point observations
                elif self.inner.startswith("point"):
                    self.draw_points(ax, violin_data, i)

            # Option 2: we have nested grouping by a hue variable
            # ---------------------------------------------------

            else:
                offsets = self.hue_offsets
                for j, hue_level in enumerate(self.hue_names):

                    support, density = self.support[i][j], self.density[i][j]
                    kws["facecolor"] = self.colors[j]

                    # Add legend data, but just for one set of violins
                    if not i:
                        self.add_legend_data(ax, self.colors[j], hue_level)

                    # Handle the special case where we have no observations
                    if support.size == 0:
                        continue

                    # Handle the special case where we have one observation
                    elif support.size == 1:
                        val = np.asscalar(support)
                        d = np.asscalar(density)
                        if self.split:
                            d = d / 2
                        at_group = i + offsets[j]
                        self.draw_single_observation(ax, at_group, val, d)
                        continue

                    # Option 2a: we are drawing a single split violin
                    # -----------------------------------------------

                    if self.split:

                        grid = np.ones(self.gridsize) * i
                        if j:
                            fill_func(support,
                                      grid,
                                      grid + density * self.dwidth,
                                      **kws)
                        else:
                            fill_func(support,
                                      grid - density * self.dwidth,
                                      grid,
                                      **kws)

                        # Draw the interior representation of the data
                        if self.inner is None:
                            continue

                        # Get a nan-free vector of datapoints
                        hue_mask = self.plot_hues[i] == hue_level
                        violin_data = remove_na(group_data[hue_mask])

                        # Draw quartile lines
                        if self.inner.startswith("quart"):
                            self.draw_quartiles(ax, violin_data,
                                                support, density, i,
                                                ["left", "right"][j])

                        # Draw stick observations
                        elif self.inner.startswith("stick"):
                            self.draw_stick_lines(ax, violin_data,
                                                  support, density, i,
                                                  ["left", "right"][j])

                        # The box and point interior plots are drawn for
                        # all data at the group level, so we just do that once
                        if not j:
                            continue

                        # Get the whole vector for this group level
                        violin_data = remove_na(group_data)

                        # Draw box and whisker information
                        if self.inner.startswith("box"):
                            self.draw_box_lines(ax, violin_data,
                                                support, density, i)

                        # Draw point observations
                        elif self.inner.startswith("point"):
                            self.draw_points(ax, violin_data, i)

                    # Option 2b: we are drawing full nested violins
                    # -----------------------------------------------

                    else:
                        grid = np.ones(self.gridsize) * (i + offsets[j])
                        fill_func(support,
                                  grid - density * self.dwidth,
                                  grid + density * self.dwidth,
                                  **kws)

                        # Draw the interior representation
                        if self.inner is None:
                            continue

                        # Get a nan-free vector of datapoints
                        hue_mask = self.plot_hues[i] == hue_level
                        violin_data = remove_na(group_data[hue_mask])

                        # Draw box and whisker information
                        if self.inner.startswith("box"):
                            self.draw_box_lines(ax, violin_data,
                                                support, density,
                                                i + offsets[j])

                        # Draw quartile lines
                        elif self.inner.startswith("quart"):
                            self.draw_quartiles(ax, violin_data,
                                                support, density,
                                                i + offsets[j])

                        # Draw stick observations
                        elif self.inner.startswith("stick"):
                            self.draw_stick_lines(ax, violin_data,
                                                  support, density,
                                                  i + offsets[j])

                        # Draw point observations
                        elif self.inner.startswith("point"):
                            self.draw_points(ax, violin_data, i + offsets[j])

    def draw_single_observation(self, ax, at_group, at_quant, density):
        """Draw a line to mark a single observation."""
        d_width = density * self.dwidth
        if self.orient == "v":
            ax.plot([at_group - d_width, at_group + d_width],
                    [at_quant, at_quant],
                    color=self.gray,
                    linewidth=self.linewidth)
        else:
            ax.plot([at_quant, at_quant],
                    [at_group - d_width, at_group + d_width],
                    color=self.gray,
                    linewidth=self.linewidth)

    def draw_box_lines(self, ax, data, support, density, center):
        """Draw boxplot information at center of the density."""
        # Compute the boxplot statistics
        q25, q50, q75 = np.percentile(data, [25, 50, 75])
        whisker_lim = 1.5 * iqr(data)
        h1 = np.min(data[data >= (q25 - whisker_lim)])
        h2 = np.max(data[data <= (q75 + whisker_lim)])

        # Draw a boxplot using lines and a point
        if self.orient == "v":
            ax.plot([center, center], [h1, h2],
                    linewidth=self.linewidth,
                    color=self.gray)
            ax.plot([center, center], [q25, q75],
                    linewidth=self.linewidth * 3,
                    color=self.gray)
            ax.scatter(center, q50,
                       zorder=3,
                       color="white",
                       edgecolor=self.gray,
                       s=np.square(self.linewidth * 2))
        else:
            ax.plot([h1, h2], [center, center],
                    linewidth=self.linewidth,
                    color=self.gray)
            ax.plot([q25, q75], [center, center],
                    linewidth=self.linewidth * 3,
                    color=self.gray)
            ax.scatter(q50, center,
                       zorder=3,
                       color="white",
                       edgecolor=self.gray,
                       s=np.square(self.linewidth * 2))

    def draw_quartiles(self, ax, data, support, density, center, split=False):
        """Draw the quartiles as lines at width of density."""
        q25, q50, q75 = np.percentile(data, [25, 50, 75])

        self.draw_to_density(ax, center, q25, support, density, split,
                             linewidth=self.linewidth,
                             dashes=[self.linewidth * 1.5] * 2)
        self.draw_to_density(ax, center, q50, support, density, split,
                             linewidth=self.linewidth,
                             dashes=[self.linewidth * 3] * 2)
        self.draw_to_density(ax, center, q75, support, density, split,
                             linewidth=self.linewidth,
                             dashes=[self.linewidth * 1.5] * 2)

    def draw_points(self, ax, data, center):
        """Draw individual observations as points at middle of the violin."""
        kws = dict(s=np.square(self.linewidth * 2),
                   color=self.gray,
                   edgecolor=self.gray)

        grid = np.ones(len(data)) * center

        if self.orient == "v":
            ax.scatter(grid, data, **kws)
        else:
            ax.scatter(data, grid, **kws)

    def draw_stick_lines(self, ax, data, support, density,
                         center, split=False):
        """Draw individual observations as sticks at width of density."""
        for val in data:
            self.draw_to_density(ax, center, val, support, density, split,
                                 linewidth=self.linewidth * .5)

    def draw_to_density(self, ax, center, val, support, density, split, **kws):
        """Draw a line orthogonal to the value axis at width of density."""
        idx = np.argmin(np.abs(support - val))
        width = self.dwidth * density[idx] * .99

        kws["color"] = self.gray

        if self.orient == "v":
            if split == "left":
                ax.plot([center - width, center], [val, val], **kws)
            elif split == "right":
                ax.plot([center, center + width], [val, val], **kws)
            else:
                ax.plot([center - width, center + width], [val, val], **kws)
        else:
            if split == "left":
                ax.plot([val, val], [center - width, center], **kws)
            elif split == "right":
                ax.plot([val, val], [center, center + width], **kws)
            else:
                ax.plot([val, val], [center - width, center + width], **kws)

    def plot(self, ax):
        """Make the violin plot."""
        self.draw_violins(ax)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()


class _CategoricalScatterPlotter(_CategoricalPlotter):

    default_palette = "dark"

    @property
    def point_colors(self):
        """Return a color for each scatter point based on group and hue."""
        colors = []
        for i, group_data in enumerate(self.plot_data):

            # Initialize the array for this group level
            group_colors = np.empty((group_data.size, 3))

            if self.plot_hues is None:

                # Use the same color for all points at this level
                group_color = self.colors[i]
                group_colors[:] = group_color

            else:

                # Color the points based on  the hue level
                for j, level in enumerate(self.hue_names):
                    hue_color = self.colors[j]
                    if group_data.size:
                        group_colors[self.plot_hues[i] == level] = hue_color

            colors.append(group_colors)

        return colors

    def add_legend_data(self, ax):
        """Add empty scatterplot artists with labels for the legend."""
        if self.hue_names is not None:
            for rgb, label in zip(self.colors, self.hue_names):
                ax.scatter([], [],
                           color=mpl.colors.rgb2hex(rgb),
                           label=label,
                           s=60)


class _StripPlotter(_CategoricalScatterPlotter):
    """1-d scatterplot with categorical organization."""
    def __init__(self, x, y, hue, data, order, hue_order,
                 jitter, dodge, orient, color, palette):
        """Initialize the plotter."""
        self.establish_variables(x, y, hue, data, orient, order, hue_order)
        self.establish_colors(color, palette, 1)

        # Set object attributes
        self.dodge = dodge
        self.width = .8

        if jitter == 1:  # Use a good default for `jitter = True`
            jlim = 0.1
        else:
            jlim = float(jitter)
        if self.hue_names is not None and dodge:
            jlim /= len(self.hue_names)
        self.jitterer = stats.uniform(-jlim, jlim * 2).rvs

    def draw_stripplot(self, ax, kws):
        """Draw the points onto `ax`."""
        # Set the default zorder to 2.1, so that the points
        # will be drawn on top of line elements (like in a boxplot)
        for i, group_data in enumerate(self.plot_data):
            if self.plot_hues is None or not self.dodge:

                if self.hue_names is None:
                    hue_mask = np.ones(group_data.size, np.bool)
                else:
                    hue_mask = np.array([h in self.hue_names
                                         for h in self.plot_hues[i]], np.bool)
                    # Broken on older numpys
                    # hue_mask = np.in1d(self.plot_hues[i], self.hue_names)

                strip_data = group_data[hue_mask]

                # Plot the points in centered positions
                cat_pos = np.ones(strip_data.size) * i
                cat_pos += self.jitterer(len(strip_data))
                kws.update(c=self.point_colors[i][hue_mask])
                if self.orient == "v":
                    ax.scatter(cat_pos, strip_data, **kws)
                else:
                    ax.scatter(strip_data, cat_pos, **kws)

            else:
                offsets = self.hue_offsets
                for j, hue_level in enumerate(self.hue_names):
                    hue_mask = self.plot_hues[i] == hue_level
                    strip_data = group_data[hue_mask]

                    # Plot the points in centered positions
                    center = i + offsets[j]
                    cat_pos = np.ones(strip_data.size) * center
                    cat_pos += self.jitterer(len(strip_data))
                    kws.update(c=self.point_colors[i][hue_mask])
                    if self.orient == "v":
                        ax.scatter(cat_pos, strip_data, **kws)
                    else:
                        ax.scatter(strip_data, cat_pos, **kws)

    def plot(self, ax, kws):
        """Make the plot."""
        self.draw_stripplot(ax, kws)
        self.add_legend_data(ax)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()


class _SwarmPlotter(_CategoricalScatterPlotter):

    def __init__(self, x, y, hue, data, order, hue_order,
                 dodge, orient, color, palette):
        """Initialize the plotter."""
        self.establish_variables(x, y, hue, data, orient, order, hue_order)
        self.establish_colors(color, palette, 1)

        # Set object attributes
        self.dodge = dodge
        self.width = .8

    def could_overlap(self, xy_i, swarm, d):
        """Return a list of all swarm points that could overlap with target.

        Assumes that swarm is a sorted list of all points below xy_i.
        """
        _, y_i = xy_i
        neighbors = []
        for xy_j in reversed(swarm):
            _, y_j = xy_j
            if (y_i - y_j) < d:
                neighbors.append(xy_j)
            else:
                break
        return np.array(list(reversed(neighbors)))

    def position_candidates(self, xy_i, neighbors, d):
        """Return a list of (x, y) coordinates that might be valid."""
        candidates = [xy_i]
        x_i, y_i = xy_i
        left_first = True
        for x_j, y_j in neighbors:
            dy = y_i - y_j
            dx = np.sqrt(d ** 2 - dy ** 2) * 1.05
            cl, cr = (x_j - dx, y_i), (x_j + dx, y_i)
            if left_first:
                new_candidates = [cl, cr]
            else:
                new_candidates = [cr, cl]
            candidates.extend(new_candidates)
            left_first = not left_first
        return np.array(candidates)

    def first_non_overlapping_candidate(self, candidates, neighbors, d):
        """Remove candidates from the list if they overlap with the swarm."""

        # IF we have no neighbours, all candidates are good.
        if len(neighbors) == 0:
            return candidates[0]

        neighbors_x = neighbors[:, 0]
        neighbors_y = neighbors[:, 1]

        d_square = d ** 2

        for xy_i in candidates:
            x_i, y_i = xy_i

            dx = neighbors_x - x_i
            dy = neighbors_y - y_i

            sq_distances = np.power(dx, 2.0) + np.power(dy, 2.0)

            # good candidate does not overlap any of neighbors
            # which means that squared distance between candidate
            # and any of the neighbours has to be at least
            # square of the diameter
            good_candidate = np.all(sq_distances >= d_square)

            if good_candidate:
                return xy_i

        # If `position_candidates` works well
        # this should never happen
        raise Exception('No non-overlapping candidates found. '
                        'This should not happen.')

    def beeswarm(self, orig_xy, d):
        """Adjust x position of points to avoid overlaps."""
        # In this method, ``x`` is always the categorical axis
        # Center of the swarm, in point coordinates
        midline = orig_xy[0, 0]

        # Start the swarm with the first point
        swarm = [orig_xy[0]]

        # Loop over the remaining points
        for xy_i in orig_xy[1:]:

            # Find the points in the swarm that could possibly
            # overlap with the point we are currently placing
            neighbors = self.could_overlap(xy_i, swarm, d)

            # Find positions that would be valid individually
            # with respect to each of the swarm neighbors
            candidates = self.position_candidates(xy_i, neighbors, d)

            # Sort candidates by their centrality
            offsets = np.abs(candidates[:, 0] - midline)
            candidates = candidates[np.argsort(offsets)]

            # Find the first candidate that doesn't overlap any neighbours
            new_xy_i = self.first_non_overlapping_candidate(candidates,
                                                            neighbors, d)

            # Place it into the swarm
            swarm.append(new_xy_i)

        return np.array(swarm)

    def add_gutters(self, points, center, width):
        """Stop points from extending beyond their territory."""
        half_width = width / 2
        low_gutter = center - half_width
        off_low = points < low_gutter
        if off_low.any():
            points[off_low] = low_gutter
        high_gutter = center + half_width
        off_high = points > high_gutter
        if off_high.any():
            points[off_high] = high_gutter
        return points

    def swarm_points(self, ax, points, center, width, s, **kws):
        """Find new positions on the categorical axis for each point."""
        # Convert from point size (area) to diameter
        default_lw = mpl.rcParams["patch.linewidth"]
        lw = kws.get("linewidth", kws.get("lw", default_lw))
        dpi = ax.figure.dpi
        d = (np.sqrt(s) + lw) * (dpi / 72)

        # Transform the data coordinates to point coordinates.
        # We'll figure out the swarm positions in the latter
        # and then convert back to data coordinates and replot
        orig_xy = ax.transData.transform(points.get_offsets())

        # Order the variables so that x is the caegorical axis
        if self.orient == "h":
            orig_xy = orig_xy[:, [1, 0]]

        # Do the beeswarm in point coordinates
        new_xy = self.beeswarm(orig_xy, d)

        # Transform the point coordinates back to data coordinates
        if self.orient == "h":
            new_xy = new_xy[:, [1, 0]]
        new_x, new_y = ax.transData.inverted().transform(new_xy).T

        # Add gutters
        if self.orient == "v":
            self.add_gutters(new_x, center, width)
        else:
            self.add_gutters(new_y, center, width)

        # Reposition the points so they do not overlap
        points.set_offsets(np.c_[new_x, new_y])

    def draw_swarmplot(self, ax, kws):
        """Plot the data."""
        s = kws.pop("s")

        centers = []
        swarms = []

        # Set the categorical axes limits here for the swarm math
        if self.orient == "v":
            ax.set_xlim(-.5, len(self.plot_data) - .5)
        else:
            ax.set_ylim(-.5, len(self.plot_data) - .5)

        # Plot each swarm
        for i, group_data in enumerate(self.plot_data):

            if self.plot_hues is None or not self.dodge:

                width = self.width

                if self.hue_names is None:
                    hue_mask = np.ones(group_data.size, np.bool)
                else:
                    hue_mask = np.array([h in self.hue_names
                                         for h in self.plot_hues[i]], np.bool)
                    # Broken on older numpys
                    # hue_mask = np.in1d(self.plot_hues[i], self.hue_names)

                swarm_data = group_data[hue_mask]

                # Sort the points for the beeswarm algorithm
                sorter = np.argsort(swarm_data)
                swarm_data = swarm_data[sorter]
                point_colors = self.point_colors[i][hue_mask][sorter]

                # Plot the points in centered positions
                cat_pos = np.ones(swarm_data.size) * i
                kws.update(c=point_colors)
                if self.orient == "v":
                    points = ax.scatter(cat_pos, swarm_data, s=s, **kws)
                else:
                    points = ax.scatter(swarm_data, cat_pos, s=s, **kws)

                centers.append(i)
                swarms.append(points)

            else:
                offsets = self.hue_offsets
                width = self.nested_width

                for j, hue_level in enumerate(self.hue_names):
                    hue_mask = self.plot_hues[i] == hue_level
                    swarm_data = group_data[hue_mask]

                    # Sort the points for the beeswarm algorithm
                    sorter = np.argsort(swarm_data)
                    swarm_data = swarm_data[sorter]
                    point_colors = self.point_colors[i][hue_mask][sorter]

                    # Plot the points in centered positions
                    center = i + offsets[j]
                    cat_pos = np.ones(swarm_data.size) * center
                    kws.update(c=point_colors)
                    if self.orient == "v":
                        points = ax.scatter(cat_pos, swarm_data, s=s, **kws)
                    else:
                        points = ax.scatter(swarm_data, cat_pos, s=s, **kws)

                    centers.append(center)
                    swarms.append(points)

        # Update the position of each point on the categorical axis
        # Do this after plotting so that the numerical axis limits are correct
        for center, swarm in zip(centers, swarms):
            if swarm.get_offsets().size:
                self.swarm_points(ax, swarm, center, width, s, **kws)

    def plot(self, ax, kws):
        """Make the full plot."""
        self.draw_swarmplot(ax, kws)
        self.add_legend_data(ax)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()


class _CategoricalStatPlotter(_CategoricalPlotter):

    @property
    def nested_width(self):
        """A float with the width of plot elements when hue nesting is used."""
        if self.dodge:
            width = self.width / len(self.hue_names)
        else:
            width = self.width
        return width

    def estimate_statistic(self, estimator, ci, n_boot):

        if self.hue_names is None:
            statistic = []
            confint = []
        else:
            statistic = [[] for _ in self.plot_data]
            confint = [[] for _ in self.plot_data]

        for i, group_data in enumerate(self.plot_data):

            # Option 1: we have a single layer of grouping
            # --------------------------------------------

            if self.plot_hues is None:

                if self.plot_units is None:
                    stat_data = remove_na(group_data)
                    unit_data = None
                else:
                    unit_data = self.plot_units[i]
                    have = pd.notnull(np.c_[group_data, unit_data]).all(axis=1)
                    stat_data = group_data[have]
                    unit_data = unit_data[have]

                # Estimate a statistic from the vector of data
                if not stat_data.size:
                    statistic.append(np.nan)
                else:
                    statistic.append(estimator(stat_data))

                # Get a confidence interval for this estimate
                if ci is not None:

                    if stat_data.size < 2:
                        confint.append([np.nan, np.nan])
                        continue

                    if ci == "sd":

                        estimate = estimator(stat_data)
                        sd = np.std(stat_data)
                        confint.append((estimate - sd, estimate + sd))

                    else:

                        boots = bootstrap(stat_data, func=estimator,
                                          n_boot=n_boot,
                                          units=unit_data)
                        confint.append(utils.ci(boots, ci))

            # Option 2: we are grouping by a hue layer
            # ----------------------------------------

            else:
                for j, hue_level in enumerate(self.hue_names):

                    if not self.plot_hues[i].size:
                        statistic[i].append(np.nan)
                        if ci is not None:
                            confint[i].append((np.nan, np.nan))
                        continue

                    hue_mask = self.plot_hues[i] == hue_level
                    if self.plot_units is None:
                        stat_data = remove_na(group_data[hue_mask])
                        unit_data = None
                    else:
                        group_units = self.plot_units[i]
                        have = pd.notnull(
                            np.c_[group_data, group_units]
                            ).all(axis=1)
                        stat_data = group_data[hue_mask & have]
                        unit_data = group_units[hue_mask & have]

                    # Estimate a statistic from the vector of data
                    if not stat_data.size:
                        statistic[i].append(np.nan)
                    else:
                        statistic[i].append(estimator(stat_data))

                    # Get a confidence interval for this estimate
                    if ci is not None:

                        if stat_data.size < 2:
                            confint[i].append([np.nan, np.nan])
                            continue

                        if ci == "sd":

                            estimate = estimator(stat_data)
                            sd = np.std(stat_data)
                            confint[i].append((estimate - sd, estimate + sd))

                        else:

                            boots = bootstrap(stat_data, func=estimator,
                                              n_boot=n_boot,
                                              units=unit_data)
                            confint[i].append(utils.ci(boots, ci))

        # Save the resulting values for plotting
        self.statistic = np.array(statistic)
        self.confint = np.array(confint)

    def draw_confints(self, ax, at_group, confint, colors,
                      errwidth=None, capsize=None, **kws):

        if errwidth is not None:
            kws.setdefault("lw", errwidth)
        else:
            kws.setdefault("lw", mpl.rcParams["lines.linewidth"] * 1.8)

        for at, (ci_low, ci_high), color in zip(at_group,
                                                confint,
                                                colors):
            if self.orient == "v":
                ax.plot([at, at], [ci_low, ci_high], color=color, **kws)
                if capsize is not None:
                    ax.plot([at - capsize / 2, at + capsize / 2],
                            [ci_low, ci_low], color=color, **kws)
                    ax.plot([at - capsize / 2, at + capsize / 2],
                            [ci_high, ci_high], color=color, **kws)
            else:
                ax.plot([ci_low, ci_high], [at, at], color=color, **kws)
                if capsize is not None:
                    ax.plot([ci_low, ci_low],
                            [at - capsize / 2, at + capsize / 2],
                            color=color, **kws)
                    ax.plot([ci_high, ci_high],
                            [at - capsize / 2, at + capsize / 2],
                            color=color, **kws)


class _BarPlotter(_CategoricalStatPlotter):
    """Show point estimates and confidence intervals with bars."""

    def __init__(self, x, y, hue, data, order, hue_order,
                 estimator, ci, n_boot, units,
                 orient, color, palette, saturation, errcolor,
                 errwidth, capsize, dodge):
        """Initialize the plotter."""
        self.establish_variables(x, y, hue, data, orient,
                                 order, hue_order, units)
        self.establish_colors(color, palette, saturation)
        self.estimate_statistic(estimator, ci, n_boot)

        self.dodge = dodge

        self.errcolor = errcolor
        self.errwidth = errwidth
        self.capsize = capsize

    def draw_bars(self, ax, kws):
        """Draw the bars onto `ax`."""
        # Get the right matplotlib function depending on the orientation
        barfunc = ax.bar if self.orient == "v" else ax.barh
        barpos = np.arange(len(self.statistic))

        if self.plot_hues is None:

            # Draw the bars
            barfunc(barpos, self.statistic, self.width,
                    color=self.colors, align="center", **kws)

            # Draw the confidence intervals
            errcolors = [self.errcolor] * len(barpos)
            self.draw_confints(ax,
                               barpos,
                               self.confint,
                               errcolors,
                               self.errwidth,
                               self.capsize)

        else:

            for j, hue_level in enumerate(self.hue_names):

                # Draw the bars
                offpos = barpos + self.hue_offsets[j]
                barfunc(offpos, self.statistic[:, j], self.nested_width,
                        color=self.colors[j], align="center",
                        label=hue_level, **kws)

                # Draw the confidence intervals
                if self.confint.size:
                    confint = self.confint[:, j]
                    errcolors = [self.errcolor] * len(offpos)
                    self.draw_confints(ax,
                                       offpos,
                                       confint,
                                       errcolors,
                                       self.errwidth,
                                       self.capsize)

    def plot(self, ax, bar_kws):
        """Make the plot."""
        self.draw_bars(ax, bar_kws)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()


class _PointPlotter(_CategoricalStatPlotter):

    default_palette = "dark"

    """Show point estimates and confidence intervals with (joined) points."""
    def __init__(self, x, y, hue, data, order, hue_order,
                 estimator, ci, n_boot, units,
                 markers, linestyles, dodge, join, scale,
                 orient, color, palette, errwidth=None, capsize=None):
        """Initialize the plotter."""
        self.establish_variables(x, y, hue, data, orient,
                                 order, hue_order, units)
        self.establish_colors(color, palette, 1)
        self.estimate_statistic(estimator, ci, n_boot)

        # Override the default palette for single-color plots
        if hue is None and color is None and palette is None:
            self.colors = [color_palette()[0]] * len(self.colors)

        # Don't join single-layer plots with different colors
        if hue is None and palette is not None:
            join = False

        # Use a good default for `dodge=True`
        if dodge is True and self.hue_names is not None:
            dodge = .025 * len(self.hue_names)

        # Make sure we have a marker for each hue level
        if isinstance(markers, string_types):
            markers = [markers] * len(self.colors)
        self.markers = markers

        # Make sure we have a line style for each hue level
        if isinstance(linestyles, string_types):
            linestyles = [linestyles] * len(self.colors)
        self.linestyles = linestyles

        # Set the other plot components
        self.dodge = dodge
        self.join = join
        self.scale = scale
        self.errwidth = errwidth
        self.capsize = capsize

    @property
    def hue_offsets(self):
        """Offsets relative to the center position for each hue level."""
        if self.dodge:
            offset = np.linspace(0, self.dodge, len(self.hue_names))
            offset -= offset.mean()
        else:
            offset = np.zeros(len(self.hue_names))
        return offset

    def draw_points(self, ax):
        """Draw the main data components of the plot."""
        # Get the center positions on the categorical axis
        pointpos = np.arange(len(self.statistic))

        # Get the size of the plot elements
        lw = mpl.rcParams["lines.linewidth"] * 1.8 * self.scale
        mew = lw * .75
        markersize = np.pi * np.square(lw) * 2

        if self.plot_hues is None:

            # Draw lines joining each estimate point
            if self.join:
                color = self.colors[0]
                ls = self.linestyles[0]
                if self.orient == "h":
                    ax.plot(self.statistic, pointpos,
                            color=color, ls=ls, lw=lw)
                else:
                    ax.plot(pointpos, self.statistic,
                            color=color, ls=ls, lw=lw)

            # Draw the confidence intervals
            self.draw_confints(ax, pointpos, self.confint, self.colors,
                               self.errwidth, self.capsize)

            # Draw the estimate points
            marker = self.markers[0]
            hex_colors = [mpl.colors.rgb2hex(c) for c in self.colors]
            if self.orient == "h":
                x, y = self.statistic, pointpos
            else:
                x, y = pointpos, self.statistic
            ax.scatter(x, y,
                       linewidth=mew, marker=marker, s=markersize,
                       c=hex_colors, edgecolor=hex_colors)

        else:

            offsets = self.hue_offsets
            for j, hue_level in enumerate(self.hue_names):

                # Determine the values to plot for this level
                statistic = self.statistic[:, j]

                # Determine the position on the categorical and z axes
                offpos = pointpos + offsets[j]
                z = j + 1

                # Draw lines joining each estimate point
                if self.join:
                    color = self.colors[j]
                    ls = self.linestyles[j]
                    if self.orient == "h":
                        ax.plot(statistic, offpos, color=color,
                                zorder=z, ls=ls, lw=lw)
                    else:
                        ax.plot(offpos, statistic, color=color,
                                zorder=z, ls=ls, lw=lw)

                # Draw the confidence intervals
                if self.confint.size:
                    confint = self.confint[:, j]
                    errcolors = [self.colors[j]] * len(offpos)
                    self.draw_confints(ax, offpos, confint, errcolors,
                                       self.errwidth, self.capsize,
                                       zorder=z)

                # Draw the estimate points
                n_points = len(remove_na(offpos))
                marker = self.markers[j]
                hex_color = mpl.colors.rgb2hex(self.colors[j])
                if n_points:
                    point_colors = [hex_color for _ in range(n_points)]
                else:
                    point_colors = hex_color
                if self.orient == "h":
                    x, y = statistic, offpos
                else:
                    x, y = offpos, statistic
                if not len(remove_na(statistic)):
                    x, y = [], []
                ax.scatter(x, y, label=hue_level,
                           c=point_colors, edgecolor=point_colors,
                           linewidth=mew, marker=marker, s=markersize,
                           zorder=z)

    def plot(self, ax):
        """Make the plot."""
        self.draw_points(ax)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()


class _LVPlotter(_CategoricalPlotter):

    def __init__(self, x, y, hue, data, order, hue_order,
                 orient, color, palette, saturation,
                 width, dodge, k_depth, linewidth, scale, outlier_prop):

        # TODO assigning variables for None is unceccesary
        if width is None:
            width = .8
        self.width = width

        self.dodge = dodge

        if saturation is None:
            saturation = .75
        self.saturation = saturation

        if k_depth is None:
            k_depth = 'proportion'
        self.k_depth = k_depth

        if linewidth is None:
            linewidth = mpl.rcParams["lines.linewidth"]
        self.linewidth = linewidth

        if scale is None:
            scale = 'exponential'
        self.scale = scale

        self.outlier_prop = outlier_prop

        self.establish_variables(x, y, hue, data, orient, order, hue_order)
        self.establish_colors(color, palette, saturation)

    def _lv_box_ends(self, vals, k_depth='proportion', outlier_prop=None):
        """Get the number of data points and calculate `depth` of
        letter-value plot."""
        vals = np.asarray(vals)
        vals = vals[np.isfinite(vals)]
        n = len(vals)
        # If p is not set, calculate it so that 8 points are outliers
        if not outlier_prop:
            # Conventional boxplots assume this proportion of the data are
            # outliers.
            p = 0.007
        else:
            if ((outlier_prop > 1.) or (outlier_prop < 0.)):
                raise ValueError('outlier_prop not in range [0, 1]!')
            p = outlier_prop
        # Select the depth, i.e. number of boxes to draw, based on the method
        k_dict = {'proportion': (np.log2(n)) - int(np.log2(n*p)) + 1,
                  'tukey': (np.log2(n)) - 3,
                  'trustworthy': (np.log2(n) -
                                  np.log2(2*stats.norm.ppf((1-p))**2)) + 1}
        k = k_dict[k_depth]
        try:
            k = int(k)
        except ValueError:
            k = 1
        # If the number happens to be less than 0, set k to 0
        if k < 1.:
            k = 1
        # Calculate the upper box ends
        upper = [100*(1 - 0.5**(i+2)) for i in range(k, -1, -1)]
        # Calculate the lower box ends
        lower = [100*(0.5**(i+2)) for i in range(k, -1, -1)]
        # Stitch the box ends together
        percentile_ends = [(i, j) for i, j in zip(lower, upper)]
        box_ends = [np.percentile(vals, q) for q in percentile_ends]
        return box_ends, k

    def _lv_outliers(self, vals, k):
        """Find the outliers based on the letter value depth."""
        perc_ends = (100*(0.5**(k+2)), 100*(1 - 0.5**(k+2)))
        edges = np.percentile(vals, perc_ends)
        lower_out = vals[np.where(vals < edges[0])[0]]
        upper_out = vals[np.where(vals > edges[1])[0]]
        return np.concatenate((lower_out, upper_out))

    def _width_functions(self, width_func):
        # Dictionary of functions for computing the width of the boxes
        width_functions = {'linear': lambda h, i, k: (i + 1.) / k,
                           'exponential': lambda h, i, k: 2**(-k+i-1),
                           'area': lambda h, i, k: (1 - 2**(-k+i-2)) / h}
        return width_functions[width_func]

    def _lvplot(self, box_data, positions,
                color=[255. / 256., 185. / 256., 0.],
                vert=True, widths=1, k_depth='proportion',
                ax=None, outlier_prop=None, scale='exponential',
                **kws):

        x = positions[0]
        box_data = np.asarray(box_data)

        # If we only have one data point, plot a line
        if len(box_data) == 1:
            kws.update({'color': self.gray, 'linestyle': '-'})
            ys = [box_data[0], box_data[0]]
            xs = [x - widths / 2, x + widths / 2]
            if vert:
                xx, yy = xs, ys
            else:
                xx, yy = ys, xs
            ax.plot(xx, yy, **kws)
        else:
            # Get the number of data points and calculate "depth" of
            # letter-value plot
            box_ends, k = self._lv_box_ends(box_data, k_depth=k_depth,
                                            outlier_prop=outlier_prop)

            # Anonymous functions for calculating the width and height
            # of the letter value boxes
            width = self._width_functions(scale)

            # Function to find height of boxes
            def height(b):
                return b[1] - b[0]

            # Functions to construct the letter value boxes
            def vert_perc_box(x, b, i, k, w):
                rect = Patches.Rectangle((x - widths*w / 2, b[0]),
                                         widths*w,
                                         height(b), fill=True)
                return rect

            def horz_perc_box(x, b, i, k, w):
                rect = Patches.Rectangle((b[0], x - widths*w / 2),
                                         height(b), widths*w,
                                         fill=True)
                return rect

            # Scale the width of the boxes so the biggest starts at 1
            w_area = np.array([width(height(b), i, k)
                               for i, b in enumerate(box_ends)])
            w_area = w_area / np.max(w_area)

            # Calculate the medians
            y = np.median(box_data)

            # Calculate the outliers and plot
            outliers = self._lv_outliers(box_data, k)

            if vert:
                boxes = [vert_perc_box(x, b[0], i, k, b[1])
                         for i, b in enumerate(zip(box_ends, w_area))]

                # Plot the medians
                ax.plot([x - widths / 2, x + widths / 2], [y, y],
                        c='.15', alpha=.45, **kws)

                ax.scatter(np.repeat(x, len(outliers)), outliers,
                           marker='d', c=mpl.colors.rgb2hex(color), **kws)
            else:
                boxes = [horz_perc_box(x, b[0], i, k, b[1])
                         for i, b in enumerate(zip(box_ends, w_area))]

                # Plot the medians
                ax.plot([y, y], [x - widths / 2, x + widths / 2],
                        c='.15', alpha=.45, **kws)

                ax.scatter(outliers, np.repeat(x, len(outliers)),
                           marker='d', c=color, **kws)

            # Construct a color map from the input color
            rgb = [[1, 1, 1], list(color)]
            cmap = mpl.colors.LinearSegmentedColormap.from_list('new_map', rgb)
            collection = PatchCollection(boxes, cmap=cmap)

            # Set the color gradation
            collection.set_array(np.array(np.linspace(0, 1, len(boxes))))

            # Plot the boxes
            ax.add_collection(collection)

    def draw_letter_value_plot(self, ax, kws):
        """Use matplotlib to draw a letter value plot on an Axes."""
        vert = self.orient == "v"

        for i, group_data in enumerate(self.plot_data):

            if self.plot_hues is None:

                # Handle case where there is data at this level
                if group_data.size == 0:
                    continue

                # Draw a single box or a set of boxes
                # with a single level of grouping
                box_data = remove_na(group_data)

                # Handle case where there is no non-null data
                if box_data.size == 0:
                    continue

                color = self.colors[i]

                artist_dict = self._lvplot(box_data,
                                           positions=[i],
                                           color=color,
                                           vert=vert,
                                           widths=self.width,
                                           k_depth=self.k_depth,
                                           ax=ax,
                                           scale=self.scale,
                                           outlier_prop=self.outlier_prop,
                                           **kws)

            else:
                # Draw nested groups of boxes
                offsets = self.hue_offsets
                for j, hue_level in enumerate(self.hue_names):

                    # Add a legend for this hue level
                    if not i:
                        self.add_legend_data(ax, self.colors[j], hue_level)

                    # Handle case where there is data at this level
                    if group_data.size == 0:
                        continue

                    hue_mask = self.plot_hues[i] == hue_level
                    box_data = remove_na(group_data[hue_mask])

                    # Handle case where there is no non-null data
                    if box_data.size == 0:
                        continue

                    color = self.colors[j]
                    center = i + offsets[j]
                    artist_dict = self._lvplot(box_data,
                                               positions=[center],
                                               color=color,
                                               vert=vert,
                                               widths=self.nested_width,
                                               k_depth=self.k_depth,
                                               ax=ax,
                                               scale=self.scale,
                                               outlier_prop=self.outlier_prop,
                                               **kws)

    def plot(self, ax, boxplot_kws):
        """Make the plot."""
        self.draw_letter_value_plot(ax, boxplot_kws)
        self.annotate_axes(ax)
        if self.orient == "h":
            ax.invert_yaxis()

_categorical_docs = dict(

    # Shared narrative docs
    main_api_narrative=dedent("""\
    Input data can be passed in a variety of formats, including:

    - Vectors of data represented as lists, numpy arrays, or pandas Series
      objects passed directly to the ``x``, ``y``, and/or ``hue`` parameters.
    - A "long-form" DataFrame, in which case the ``x``, ``y``, and ``hue``
      variables will determine how the data are plotted.
    - A "wide-form" DataFrame, such that each numeric column will be plotted.
    - Anything accepted by ``plt.boxplot`` (e.g. a 2d array or list of vectors)

    In most cases, it is possible to use numpy or Python objects, but pandas
    objects are preferable because the associated names will be used to
    annotate the axes. Additionally, you can use Categorical types for the
    grouping variables to control the order of plot elements.\
    """),

    # Shared function parameters
    input_params=dedent("""\
    x, y, hue : names of variables in ``data`` or vector data, optional
        Inputs for plotting long-form data. See examples for interpretation.\
        """),
    string_input_params=dedent("""\
    x, y, hue : names of variables in ``data``
        Inputs for plotting long-form data. See examples for interpretation.\
        """),
    categorical_data=dedent("""\
    data : DataFrame, array, or list of arrays, optional
        Dataset for plotting. If ``x`` and ``y`` are absent, this is
        interpreted as wide-form. Otherwise it is expected to be long-form.\
    """),
    long_form_data=dedent("""\
    data : DataFrame
        Long-form (tidy) dataset for plotting. Each column should correspond
        to a variable, and each row should correspond to an observation.\
    """),
    order_vars=dedent("""\
    order, hue_order : lists of strings, optional
        Order to plot the categorical levels in, otherwise the levels are
        inferred from the data objects.\
        """),
    stat_api_params=dedent("""\
    estimator : callable that maps vector -> scalar, optional
        Statistical function to estimate within each categorical bin.
    ci : float or "sd" or None, optional
        Size of confidence intervals to draw around estimated values.  If
        "sd", skip bootstrapping and draw the standard deviation of the
        observations. If ``None``, no bootstrapping will be performed, and
        error bars will not be drawn.
    n_boot : int, optional
        Number of bootstrap iterations to use when computing confidence
        intervals.
    units : name of variable in ``data`` or vector data, optional
        Identifier of sampling units, which will be used to perform a
        multilevel bootstrap and account for repeated measures design.\
    """),
    orient=dedent("""\
    orient : "v" | "h", optional
        Orientation of the plot (vertical or horizontal). This is usually
        inferred from the dtype of the input variables, but can be used to
        specify when the "categorical" variable is a numeric or when plotting
        wide-form data.\
    """),
    color=dedent("""\
    color : matplotlib color, optional
        Color for all of the elements, or seed for a gradient palette.
    """),
    palette=dedent("""\
    palette : palette name, list, or dict, optional
        Color palette that maps either the grouping variable or the hue
        variable. If the palette is a dictionary, keys should be names of
        levels and values should be matplotlib colors.\
    """),
    saturation=dedent("""\
    saturation : float, optional
        Proportion of the original saturation to draw colors at. Large patches
        often look better with slightly desaturated colors, but set this to
        ``1`` if you want the plot colors to perfectly match the input color
        spec.\
    """),
    capsize=dedent("""\
         capsize : float, optional
             Width of the "caps" on error bars.
         """),
    errwidth=dedent("""\
         errwidth : float, optional
             Thickness of error bar lines (and caps).\
         """),
    width=dedent("""\
    width : float, optional
        Width of a full element when not using hue nesting, or width of all the
        elements for one level of the major grouping variable.\
    """),
    dodge=dedent("""\
    dodge : bool, optional
        When hue nesting is used, whether elements should be shifted along the
        categorical axis.\
    """),
    linewidth=dedent("""\
    linewidth : float, optional
        Width of the gray lines that frame the plot elements.\
    """),
    ax_in=dedent("""\
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.\
    """),
    ax_out=dedent("""\
    ax : matplotlib Axes
        Returns the Axes object with the boxplot drawn onto it.\
    """),

    # Shared see also
    boxplot=dedent("""\
    boxplot : A traditional box-and-whisker plot with a similar API.\
    """),
    violinplot=dedent("""\
    violinplot : A combination of boxplot and kernel density estimation.\
    """),
    stripplot=dedent("""\
    stripplot : A scatterplot where one variable is categorical. Can be used
                in conjunction with other plots to show each observation.\
    """),
    swarmplot=dedent("""\
    swarmplot : A categorical scatterplot where the points do not overlap. Can
                be used with other plots to show each observation.\
    """),
    barplot=dedent("""\
    barplot : Show point estimates and confidence intervals using bars.\
    """),
    countplot=dedent("""\
    countplot : Show the counts of observations in each categorical bin.\
    """),
    pointplot=dedent("""\
    pointplot : Show point estimates and confidence intervals using scatterplot
                glyphs.\
    """),
    factorplot=dedent("""\
    factorplot : Combine categorical plots and a class:`FacetGrid`.\
    """),
    lvplot=dedent("""\
    lvplot : An extension of the boxplot for long-tailed and large data sets.
    """),

    )

_categorical_docs.update(_facet_docs)


def boxplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
            orient=None, color=None, palette=None, saturation=.75,
            width=.8, dodge=True, fliersize=5, linewidth=None,
            whis=1.5, notch=False, ax=None, **kwargs):

    plotter = _BoxPlotter(x, y, hue, data, order, hue_order,
                          orient, color, palette, saturation,
                          width, dodge, fliersize, linewidth)

    if ax is None:
        ax = plt.gca()
    kwargs.update(dict(whis=whis, notch=notch))

    plotter.plot(ax, kwargs)
    return ax

boxplot.__doc__ = dedent("""\
    Draw a box plot to show distributions with respect to categories.

    A box plot (or box-and-whisker plot) shows the distribution of quantitative
    data in a way that facilitates comparisons between variables or across
    levels of a categorical variable. The box shows the quartiles of the
    dataset while the whiskers extend to show the rest of the distribution,
    except for points that are determined to be "outliers" using a method
    that is a function of the inter-quartile range.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    {orient}
    {color}
    {palette}
    {saturation}
    {width}
    {dodge}
    fliersize : float, optional
        Size of the markers used to indicate outlier observations.
    {linewidth}
    whis : float, optional
        Proportion of the IQR past the low and high quartiles to extend the
        plot whiskers. Points outside this range will be identified as
        outliers.
    notch : boolean, optional
        Whether to "notch" the box to indicate a confidence interval for the
        median. There are several other parameters that can control how the
        notches are drawn; see the ``plt.boxplot`` help for more information
        on them.
    {ax_in}
    kwargs : key, value mappings
        Other keyword arguments are passed through to ``plt.boxplot`` at draw
        time.

    Returns
    -------
    {ax_out}

    See Also
    --------
    {violinplot}
    {stripplot}
    {swarmplot}

    Examples
    --------

    Draw a single horizontal boxplot:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("whitegrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.boxplot(x=tips["total_bill"])

    Draw a vertical boxplot grouped by a categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="day", y="total_bill", data=tips)

    Draw a boxplot with nested grouping by two categorical variables:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="day", y="total_bill", hue="smoker",
        ...                  data=tips, palette="Set3")

    Draw a boxplot with nested grouping when some bins are empty:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="day", y="total_bill", hue="time",
        ...                  data=tips, linewidth=2.5)

    Control box order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="time", y="tip", data=tips,
        ...                  order=["Dinner", "Lunch"])

    Draw a boxplot for each numeric variable in a DataFrame:

    .. plot::
        :context: close-figs

        >>> iris = sns.load_dataset("iris")
        >>> ax = sns.boxplot(data=iris, orient="h", palette="Set2")

    Use ``hue`` without changing box position or width:

    .. plot::
        :context: close-figs

        >>> tips["weekend"] = tips["day"].isin(["Sat", "Sun"])
        >>> ax = sns.boxplot(x="day", y="total_bill", hue="weekend",
        ...                  data=tips, dodge=False)

    Use :func:`swarmplot` to show the datapoints on top of the boxes:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="day", y="total_bill", data=tips)
        >>> ax = sns.swarmplot(x="day", y="total_bill", data=tips, color=".25")

    Use :func:`factorplot` to combine a :func:`boxplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="box",
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def violinplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
               bw="scott", cut=2, scale="area", scale_hue=True, gridsize=100,
               width=.8, inner="box", split=False, dodge=True, orient=None,
               linewidth=None, color=None, palette=None, saturation=.75,
               ax=None, **kwargs):

    plotter = _ViolinPlotter(x, y, hue, data, order, hue_order,
                             bw, cut, scale, scale_hue, gridsize,
                             width, inner, split, dodge, orient, linewidth,
                             color, palette, saturation)

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax)
    return ax

violinplot.__doc__ = dedent("""\
    Draw a combination of boxplot and kernel density estimate.

    A violin plot plays a similar role as a box and whisker plot. It shows the
    distribution of quantitative data across several levels of one (or more)
    categorical variables such that those distributions can be compared. Unlike
    a box plot, in which all of the plot components correspond to actual
    datapoints, the violin plot features a kernel density estimation of the
    underlying distribution.

    This can be an effective and attractive way to show multiple distributions
    of data at once, but keep in mind that the estimation procedure is
    influenced by the sample size, and violins for relatively small samples
    might look misleadingly smooth.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    bw : {{'scott', 'silverman', float}}, optional
        Either the name of a reference rule or the scale factor to use when
        computing the kernel bandwidth. The actual kernel size will be
        determined by multiplying the scale factor by the standard deviation of
        the data within each bin.
    cut : float, optional
        Distance, in units of bandwidth size, to extend the density past the
        extreme datapoints. Set to 0 to limit the violin range within the range
        of the observed data (i.e., to have the same effect as ``trim=True`` in
        ``ggplot``.
    scale : {{"area", "count", "width"}}, optional
        The method used to scale the width of each violin. If ``area``, each
        violin will have the same area. If ``count``, the width of the violins
        will be scaled by the number of observations in that bin. If ``width``,
        each violin will have the same width.
    scale_hue : bool, optional
        When nesting violins using a ``hue`` variable, this parameter
        determines whether the scaling is computed within each level of the
        major grouping variable (``scale_hue=True``) or across all the violins
        on the plot (``scale_hue=False``).
    gridsize : int, optional
        Number of points in the discrete grid used to compute the kernel
        density estimate.
    {width}
    inner : {{"box", "quartile", "point", "stick", None}}, optional
        Representation of the datapoints in the violin interior. If ``box``,
        draw a miniature boxplot. If ``quartiles``, draw the quartiles of the
        distribution.  If ``point`` or ``stick``, show each underlying
        datapoint. Using ``None`` will draw unadorned violins.
    split : bool, optional
        When using hue nesting with a variable that takes two levels, setting
        ``split`` to True will draw half of a violin for each level. This can
        make it easier to directly compare the distributions.
    {dodge}
    {orient}
    {linewidth}
    {color}
    {palette}
    {saturation}
    {ax_in}

    Returns
    -------
    {ax_out}

    See Also
    --------
    {boxplot}
    {stripplot}
    {swarmplot}

    Examples
    --------

    Draw a single horizontal violinplot:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("whitegrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.violinplot(x=tips["total_bill"])

    Draw a vertical violinplot grouped by a categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", data=tips)

    Draw a violinplot with nested grouping by two categorical variables:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="smoker",
        ...                     data=tips, palette="muted")

    Draw split violins to compare the across the hue variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="smoker",
        ...                     data=tips, palette="muted", split=True)

    Control violin order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="time", y="tip", data=tips,
        ...                     order=["Dinner", "Lunch"])

    Scale the violin width by the number of observations in each bin:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="sex",
        ...                     data=tips, palette="Set2", split=True,
        ...                     scale="count")

    Draw the quartiles as horizontal lines instead of a mini-box:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="sex",
        ...                     data=tips, palette="Set2", split=True,
        ...                     scale="count", inner="quartile")

    Show each observation with a stick inside the violin:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="sex",
        ...                     data=tips, palette="Set2", split=True,
        ...                     scale="count", inner="stick")

    Scale the density relative to the counts across all bins:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="sex",
        ...                     data=tips, palette="Set2", split=True,
        ...                     scale="count", inner="stick", scale_hue=False)

    Use a narrow bandwidth to reduce the amount of smoothing:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", hue="sex",
        ...                     data=tips, palette="Set2", split=True,
        ...                     scale="count", inner="stick",
        ...                     scale_hue=False, bw=.2)

    Draw horizontal violins:

    .. plot::
        :context: close-figs

        >>> planets = sns.load_dataset("planets")
        >>> ax = sns.violinplot(x="orbital_period", y="method",
        ...                     data=planets[planets.orbital_period < 1000],
        ...                     scale="width", palette="Set3")

    Don't let density extend past extreme values in the data:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="orbital_period", y="method",
        ...                     data=planets[planets.orbital_period < 1000],
        ...                     cut=0, scale="width", palette="Set3")

    Use ``hue`` without changing violin position or width:

    .. plot::
        :context: close-figs

        >>> tips["weekend"] = tips["day"].isin(["Sat", "Sun"])
        >>> ax = sns.violinplot(x="day", y="total_bill", hue="weekend",
        ...                     data=tips, dodge=False)

    Use :func:`factorplot` to combine a :func:`violinplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="violin", split=True,
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def stripplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
              jitter=False, dodge=False, orient=None, color=None, palette=None,
              size=5, edgecolor="gray", linewidth=0, ax=None, **kwargs):

    if "split" in kwargs:
        dodge = kwargs.pop("split")
        msg = "The `split` parameter has been renamed to `dodge`."
        warnings.warn(msg, UserWarning)

    plotter = _StripPlotter(x, y, hue, data, order, hue_order,
                            jitter, dodge, orient, color, palette)
    if ax is None:
        ax = plt.gca()

    kwargs.setdefault("zorder", 3)
    size = kwargs.get("s", size)
    if linewidth is None:
        linewidth = size / 10
    if edgecolor == "gray":
        edgecolor = plotter.gray
    kwargs.update(dict(s=size ** 2,
                       edgecolor=edgecolor,
                       linewidth=linewidth))

    plotter.plot(ax, kwargs)
    return ax


stripplot.__doc__ = dedent("""\
    Draw a scatterplot where one variable is categorical.

    A strip plot can be drawn on its own, but it is also a good complement
    to a box or violin plot in cases where you want to show all observations
    along with some representation of the underlying distribution.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    jitter : float, ``True``/``1`` is special-cased, optional
        Amount of jitter (only along the categorical axis) to apply. This
        can be useful when you have many points and they overlap, so that
        it is easier to see the distribution. You can specify the amount
        of jitter (half the width of the uniform random variable support),
        or just use ``True`` for a good default.
    split : bool, optional
        When using ``hue`` nesting, setting this to ``True`` will separate
        the strips for different hue levels along the categorical axis.
        Otherwise, the points for each level will be plotted on top of
        each other.
    {orient}
    {color}
    {palette}
    size : float, optional
        Diameter of the markers, in points. (Although ``plt.scatter`` is used
        to draw the points, the ``size`` argument here takes a "normal"
        markersize and not size^2 like ``plt.scatter``.
    edgecolor : matplotlib color, "gray" is special-cased, optional
        Color of the lines around each point. If you pass ``"gray"``, the
        brightness is determined by the color palette used for the body
        of the points.
    {linewidth}
    {ax_in}

    Returns
    -------
    {ax_out}

    See Also
    --------
    {swarmplot}
    {boxplot}
    {violinplot}

    Examples
    --------

    Draw a single horizontal strip plot:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("whitegrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.stripplot(x=tips["total_bill"])

    Group the strips by a categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="day", y="total_bill", data=tips)

    Add jitter to bring out the distribution of values:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="day", y="total_bill", data=tips, jitter=True)

    Use a smaller amount of jitter:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="day", y="total_bill", data=tips, jitter=0.05)

    Draw horizontal strips:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="total_bill", y="day", data=tips,
        ...                    jitter=True)

    Draw outlines around the points:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="total_bill", y="day", data=tips,
        ...                    jitter=True, linewidth=1)

    Nest the strips within a second categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="sex", y="total_bill", hue="day",
        ...                    data=tips, jitter=True)

    Draw each level of the ``hue`` variable at different locations on the
    major categorical axis:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="day", y="total_bill", hue="smoker",
        ...                    data=tips, jitter=True,
        ...                    palette="Set2", dodge=True)

    Control strip order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.stripplot(x="time", y="tip", data=tips,
        ...                    order=["Dinner", "Lunch"])

    Draw strips with large points and different aesthetics:

    .. plot::
        :context: close-figs

        >>> ax =  sns.stripplot("day", "total_bill", "smoker", data=tips,
        ...                    palette="Set2", size=20, marker="D",
        ...                    edgecolor="gray", alpha=.25)

    Draw strips of observations on top of a box plot:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="tip", y="day", data=tips, whis=np.inf)
        >>> ax = sns.stripplot(x="tip", y="day", data=tips,
        ...                    jitter=True, color=".3")

    Draw strips of observations on top of a violin plot:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", data=tips,
        ...                     inner=None, color=".8")
        >>> ax = sns.stripplot(x="day", y="total_bill", data=tips, jitter=True)

    Use :func:`factorplot` to combine a :func:`stripplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="strip",
        ...                    jitter=True,
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def swarmplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
              dodge=False, orient=None, color=None, palette=None,
              size=5, edgecolor="gray", linewidth=0, ax=None, **kwargs):

    if "split" in kwargs:
        dodge = kwargs.pop("split")
        msg = "The `split` parameter has been renamed to `dodge`."
        warnings.warn(msg, UserWarning)

    plotter = _SwarmPlotter(x, y, hue, data, order, hue_order,
                            dodge, orient, color, palette)
    if ax is None:
        ax = plt.gca()

    kwargs.setdefault("zorder", 3)
    size = kwargs.get("s", size)
    if linewidth is None:
        linewidth = size / 10
    if edgecolor == "gray":
        edgecolor = plotter.gray
    kwargs.update(dict(s=size ** 2,
                       edgecolor=edgecolor,
                       linewidth=linewidth))

    plotter.plot(ax, kwargs)
    return ax


swarmplot.__doc__ = dedent("""\
    Draw a categorical scatterplot with non-overlapping points.

    This function is similar to :func:`stripplot`, but the points are adjusted
    (only along the categorical axis) so that they don't overlap. This gives a
    better representation of the distribution of values, although it does not
    scale as well to large numbers of observations (both in terms of the
    ability to show all the points and in terms of the computation needed
    to arrange them).

    This style of plot is often called a "beeswarm".

    A swarm plot can be drawn on its own, but it is also a good complement
    to a box or violin plot in cases where you want to show all observations
    along with some representation of the underlying distribution.

    Note that arranging the points properly requires an accurate transformation
    between data and point coordinates. This means that non-default axis limits
    should be set *before* drawing the swarm plot.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    split : bool, optional
        When using ``hue`` nesting, setting this to ``True`` will separate
        the strips for different hue levels along the categorical axis.
        Otherwise, the points for each level will be plotted in one swarm.
    {orient}
    {color}
    {palette}
    size : float, optional
        Diameter of the markers, in points. (Although ``plt.scatter`` is used
        to draw the points, the ``size`` argument here takes a "normal"
        markersize and not size^2 like ``plt.scatter``.
    edgecolor : matplotlib color, "gray" is special-cased, optional
        Color of the lines around each point. If you pass ``"gray"``, the
        brightness is determined by the color palette used for the body
        of the points.
    {linewidth}
    {ax_in}

    Returns
    -------
    {ax_out}

    See Also
    --------
    {boxplot}
    {violinplot}
    {stripplot}
    {factorplot}

    Examples
    --------

    Draw a single horizontal swarm plot:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("whitegrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.swarmplot(x=tips["total_bill"])

    Group the swarms by a categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.swarmplot(x="day", y="total_bill", data=tips)

    Draw horizontal swarms:

    .. plot::
        :context: close-figs

        >>> ax = sns.swarmplot(x="total_bill", y="day", data=tips)

    Color the points using a second categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.swarmplot(x="day", y="total_bill", hue="sex", data=tips)

    Split each level of the ``hue`` variable along the categorical axis:

    .. plot::
        :context: close-figs

        >>> ax = sns.swarmplot(x="day", y="total_bill", hue="smoker",
        ...                    data=tips, palette="Set2", dodge=True)

    Control swarm order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.swarmplot(x="time", y="tip", data=tips,
        ...                    order=["Dinner", "Lunch"])

    Plot using larger points:

    .. plot::
        :context: close-figs

        >>> ax = sns.swarmplot(x="time", y="tip", data=tips, size=6)

    Draw swarms of observations on top of a box plot:

    .. plot::
        :context: close-figs

        >>> ax = sns.boxplot(x="tip", y="day", data=tips, whis=np.inf)
        >>> ax = sns.swarmplot(x="tip", y="day", data=tips, color=".2")

    Draw swarms of observations on top of a violin plot:

    .. plot::
        :context: close-figs

        >>> ax = sns.violinplot(x="day", y="total_bill", data=tips, inner=None)
        >>> ax = sns.swarmplot(x="day", y="total_bill", data=tips,
        ...                    color="white", edgecolor="gray")

    Use :func:`factorplot` to combine a :func:`swarmplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="swarm",
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def barplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
            estimator=np.mean, ci=95, n_boot=1000, units=None,
            orient=None, color=None, palette=None, saturation=.75,
            errcolor=".26", errwidth=None, capsize=None, dodge=True,
            ax=None, **kwargs):

    plotter = _BarPlotter(x, y, hue, data, order, hue_order,
                          estimator, ci, n_boot, units,
                          orient, color, palette, saturation,
                          errcolor, errwidth, capsize, dodge)

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax, kwargs)
    return ax


barplot.__doc__ = dedent("""\
    Show point estimates and confidence intervals as rectangular bars.

    A bar plot represents an estimate of central tendency for a numeric
    variable with the height of each rectangle and provides some indication of
    the uncertainty around that estimate using error bars. Bar plots include 0
    in the quantitative axis range, and they are a good choice when 0 is a
    meaningful value for the quantitative variable, and you want to make
    comparisons against it.

    For datasets where 0 is not a meaningful value, a point plot will allow you
    to focus on differences between levels of one or more categorical
    variables.

    It is also important to keep in mind that a bar plot shows only the mean
    (or other estimator) value, but in many cases it may be more informative to
    show the distribution of values at each level of the categorical variables.
    In that case, other approaches such as a box or violin plot may be more
    appropriate.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    {stat_api_params}
    {orient}
    {color}
    {palette}
    {saturation}
    errcolor : matplotlib color
        Color for the lines that represent the confidence interval.
    {errwidth}
    {capsize}
    {dodge}
    {ax_in}
    kwargs : key, value mappings
        Other keyword arguments are passed through to ``plt.bar`` at draw
        time.

    Returns
    -------
    {ax_out}

    See Also
    --------
    {countplot}
    {pointplot}
    {factorplot}

    Examples
    --------

    Draw a set of vertical bar plots grouped by a categorical variable:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("whitegrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.barplot(x="day", y="total_bill", data=tips)

    Draw a set of vertical bars with nested grouping by a two variables:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot(x="day", y="total_bill", hue="sex", data=tips)

    Draw a set of horizontal bars:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot(x="tip", y="day", data=tips)

    Control bar order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot(x="time", y="tip", data=tips,
        ...                  order=["Dinner", "Lunch"])

    Use median as the estimate of central tendency:

    .. plot::
        :context: close-figs

        >>> from numpy import median
        >>> ax = sns.barplot(x="day", y="tip", data=tips, estimator=median)

    Show the standard error of the mean with the error bars:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot(x="day", y="tip", data=tips, ci=68)

    Show standard deviation of observations instead of a confidence interval:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot(x="day", y="tip", data=tips, ci="sd")

    Add "caps" to the error bars:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot(x="day", y="tip", data=tips, capsize=.2)

    Use a different color palette for the bars:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot("size", y="total_bill", data=tips,
        ...                  palette="Blues_d")

    Use ``hue`` without changing bar position or width:

    .. plot::
        :context: close-figs

        >>> tips["weekend"] = tips["day"].isin(["Sat", "Sun"])
        >>> ax = sns.barplot(x="day", y="total_bill", hue="weekend",
        ...                  data=tips, dodge=False)

    Plot all bars in a single color:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot("size", y="total_bill", data=tips,
        ...                  color="salmon", saturation=.5)

    Use ``plt.bar`` keyword arguments to further change the aesthetic:

    .. plot::
        :context: close-figs

        >>> ax = sns.barplot("day", "total_bill", data=tips,
        ...                  linewidth=2.5, facecolor=(1, 1, 1, 0),
        ...                  errcolor=".2", edgecolor=".2")

    Use :func:`factorplot` to combine a :func:`barplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="bar",
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def pointplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
              estimator=np.mean, ci=95, n_boot=1000, units=None,
              markers="o", linestyles="-", dodge=False, join=True, scale=1,
              orient=None, color=None, palette=None, errwidth=None,
              capsize=None, ax=None, **kwargs):

    plotter = _PointPlotter(x, y, hue, data, order, hue_order,
                            estimator, ci, n_boot, units,
                            markers, linestyles, dodge, join, scale,
                            orient, color, palette, errwidth, capsize)

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax)
    return ax


pointplot.__doc__ = dedent("""\
    Show point estimates and confidence intervals using scatter plot glyphs.

    A point plot represents an estimate of central tendency for a numeric
    variable by the position of scatter plot points and provides some
    indication of the uncertainty around that estimate using error bars.

    Point plots can be more useful than bar plots for focusing comparisons
    between different levels of one or more categorical variables. They are
    particularly adept at showing interactions: how the relationship between
    levels of one categorical variable changes across levels of a second
    categorical variable. The lines that join each point from the same ``hue``
    level allow interactions to be judged by differences in slope, which is
    easier for the eyes than comparing the heights of several groups of points
    or bars.

    It is important to keep in mind that a point plot shows only the mean (or
    other estimator) value, but in many cases it may be more informative to
    show the distribution of values at each level of the categorical variables.
    In that case, other approaches such as a box or violin plot may be more
    appropriate.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    {stat_api_params}
    markers : string or list of strings, optional
        Markers to use for each of the ``hue`` levels.
    linestyles : string or list of strings, optional
        Line styles to use for each of the ``hue`` levels.
    dodge : bool or float, optional
        Amount to separate the points for each level of the ``hue`` variable
        along the categorical axis.
    join : bool, optional
        If ``True``, lines will be drawn between point estimates at the same
        ``hue`` level.
    scale : float, optional
        Scale factor for the plot elements.
    {orient}
    {color}
    {palette}
    {errwidth}
    {capsize}
    {ax_in}

    Returns
    -------
    {ax_out}

    See Also
    --------
    {barplot}
    {factorplot}

    Examples
    --------

    Draw a set of vertical point plots grouped by a categorical variable:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("darkgrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.pointplot(x="time", y="total_bill", data=tips)

    Draw a set of vertical points with nested grouping by a two variables:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="time", y="total_bill", hue="smoker",
        ...                    data=tips)

    Separate the points for different hue levels along the categorical axis:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="time", y="total_bill", hue="smoker",
        ...                    data=tips, dodge=True)

    Use a different marker and line style for the hue levels:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="time", y="total_bill", hue="smoker",
        ...                    data=tips,
        ...                    markers=["o", "x"],
        ...                    linestyles=["-", "--"])

    Draw a set of horizontal points:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="tip", y="day", data=tips)

    Don't draw a line connecting each point:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="tip", y="day", data=tips, join=False)

    Use a different color for a single-layer plot:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot("time", y="total_bill", data=tips,
        ...                    color="#bb3f3f")

    Use a different color palette for the points:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="time", y="total_bill", hue="smoker",
        ...                    data=tips, palette="Set2")

    Control point order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="time", y="tip", data=tips,
        ...                    order=["Dinner", "Lunch"])

    Use median as the estimate of central tendency:

    .. plot::
        :context: close-figs

        >>> from numpy import median
        >>> ax = sns.pointplot(x="day", y="tip", data=tips, estimator=median)

    Show the standard error of the mean with the error bars:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="day", y="tip", data=tips, ci=68)

    Show standard deviation of observations instead of a confidence interval:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="day", y="tip", data=tips, ci="sd")

    Add "caps" to the error bars:

    .. plot::
        :context: close-figs

        >>> ax = sns.pointplot(x="day", y="tip", data=tips, capsize=.2)

    Use :func:`factorplot` to combine a :func:`barplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="point",
        ...                    dodge=True,
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def countplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
              orient=None, color=None, palette=None, saturation=.75,
              dodge=True, ax=None, **kwargs):

    estimator = len
    ci = None
    n_boot = 0
    units = None
    errcolor = None
    errwidth = None
    capsize = None

    if x is None and y is not None:
        orient = "h"
        x = y
    elif y is None and x is not None:
        orient = "v"
        y = x
    elif x is not None and y is not None:
        raise TypeError("Cannot pass values for both `x` and `y`")
    else:
        raise TypeError("Must pass values for either `x` or `y`")

    plotter = _BarPlotter(x, y, hue, data, order, hue_order,
                          estimator, ci, n_boot, units,
                          orient, color, palette, saturation,
                          errcolor, errwidth, capsize, dodge)

    plotter.value_label = "count"

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax, kwargs)
    return ax


countplot.__doc__ = dedent("""\
    Show the counts of observations in each categorical bin using bars.

    A count plot can be thought of as a histogram across a categorical, instead
    of quantitative, variable. The basic API and options are identical to those
    for :func:`barplot`, so you can compare counts across nested variables.

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    {orient}
    {color}
    {palette}
    {saturation}
    {dodge}
    {ax_in}
    kwargs : key, value mappings
        Other keyword arguments are passed to ``plt.bar``.

    Returns
    -------
    {ax_out}

    See Also
    --------
    {barplot}
    {factorplot}

    Examples
    --------

    Show value counts for a single categorical variable:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set(style="darkgrid")
        >>> titanic = sns.load_dataset("titanic")
        >>> ax = sns.countplot(x="class", data=titanic)

    Show value counts for two categorical variables:

    .. plot::
        :context: close-figs

        >>> ax = sns.countplot(x="class", hue="who", data=titanic)

    Plot the bars horizontally:

    .. plot::
        :context: close-figs

        >>> ax = sns.countplot(y="class", hue="who", data=titanic)

    Use a different color palette:

    .. plot::
        :context: close-figs

        >>> ax = sns.countplot(x="who", data=titanic, palette="Set3")

    Use ``plt.bar`` keyword arguments for a different look:

    .. plot::
        :context: close-figs

        >>> ax = sns.countplot(x="who", data=titanic,
        ...                    facecolor=(0, 0, 0, 0),
        ...                    linewidth=5,
        ...                    edgecolor=sns.color_palette("dark", 3))

    Use :func:`factorplot` to combine a :func:`countplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="class", hue="who", col="survived",
        ...                    data=titanic, kind="count",
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)


def factorplot(x=None, y=None, hue=None, data=None, row=None, col=None,
               col_wrap=None, estimator=np.mean, ci=95, n_boot=1000,
               units=None, order=None, hue_order=None, row_order=None,
               col_order=None, kind="point", size=4, aspect=1,
               orient=None, color=None, palette=None,
               legend=True, legend_out=True, sharex=True, sharey=True,
               margin_titles=False, facet_kws=None, **kwargs):

    # Determine the plotting function
    try:
        plot_func = globals()[kind + "plot"]
    except KeyError:
        err = "Plot kind '{}' is not recognized".format(kind)
        raise ValueError(err)

    # Alias the input variables to determine categorical order and palette
    # correctly in the case of a count plot
    if kind == "count":
        if x is None and y is not None:
            x_, y_, orient = y, y, "h"
        elif y is None and x is not None:
            x_, y_, orient = x, x, "v"
        else:
            raise ValueError("Either `x` or `y` must be None for count plots")
    else:
        x_, y_ = x, y

    # Determine the order for the whole dataset, which will be used in all
    # facets to ensure representation of all data in the final plot
    p = _CategoricalPlotter()
    p.establish_variables(x_, y_, hue, data, orient, order, hue_order)
    order = p.group_names
    hue_order = p.hue_names

    # Determine the palette to use
    # (FacetGrid will pass a value for ``color`` to the plotting function
    # so we need to define ``palette`` to get default behavior for the
    # categorical functions
    p.establish_colors(color, palette, 1)
    if kind != "point" or hue is not None:
        palette = p.colors

    # Determine keyword arguments for the facets
    facet_kws = {} if facet_kws is None else facet_kws
    facet_kws.update(
        data=data, row=row, col=col,
        row_order=row_order, col_order=col_order,
        col_wrap=col_wrap, size=size, aspect=aspect,
        sharex=sharex, sharey=sharey,
        legend_out=legend_out, margin_titles=margin_titles,
        dropna=False,
        )

    # Determine keyword arguments for the plotting function
    plot_kws = dict(
        order=order, hue_order=hue_order,
        orient=orient, color=color, palette=palette,
        )
    plot_kws.update(kwargs)

    if kind in ["bar", "point"]:
        plot_kws.update(
            estimator=estimator, ci=ci, n_boot=n_boot, units=units,
            )

    # Initialize the facets
    g = FacetGrid(**facet_kws)

    # Draw the plot onto the facets
    g.map_dataframe(plot_func, x, y, hue, **plot_kws)

    # Special case axis labels for a count type plot
    if kind == "count":
        if x is None:
            g.set_axis_labels(x_var="count")
        if y is None:
            g.set_axis_labels(y_var="count")

    if legend and (hue is not None) and (hue not in [x, row, col]):
        hue_order = list(map(utils.to_utf8, hue_order))
        g.add_legend(title=hue, label_order=hue_order)

    return g


factorplot.__doc__ = dedent("""\
    Draw a categorical plot onto a FacetGrid.

    The default plot that is shown is a point plot, but other seaborn
    categorical plots can be chosen with the ``kind`` parameter, including
    box plots, violin plots, bar plots, or strip plots.

    It is important to choose how variables get mapped to the plot structure
    such that the most important comparisons are easiest to make. As a general
    rule, it is easier to compare positions that are closer together, so the
    ``hue`` variable should be used for the most important comparisons. For
    secondary comparisons, try to share the quantitative axis (so, use ``col``
    for vertical plots and ``row`` for horizontal plots). Note that, although
    it is possible to make rather complex plots using this function, in many
    cases you may be better served by created several smaller and more focused
    plots than by trying to stuff many comparisons into one figure.

    After plotting, the :class:`FacetGrid` with the plot is returned and can
    be used directly to tweak supporting plot details or add other layers.

    Note that, unlike when using the underlying plotting functions directly,
    data must be passed in a long-form DataFrame with variables specified by
    passing strings to ``x``, ``y``, ``hue``, and other parameters.

    As in the case with the underlying plot functions, if variables have a
    ``categorical`` data type, the correct orientation of the plot elements,
    the levels of the categorical variables, and their order will be inferred
    from the objects. Otherwise you may have to use the function parameters
    (``orient``, ``order``, ``hue_order``, etc.) to set up the plot correctly.

    Parameters
    ----------
    {string_input_params}
    {long_form_data}
    row, col : names of variables in ``data``, optional
        Categorical variables that will determine the faceting of the grid.
    {col_wrap}
    {stat_api_params}
    {order_vars}
    row_order, col_order : lists of strings, optional
        Order to organize the rows and/or columns of the grid in, otherwise the
        orders are inferred from the data objects.
    kind : {{``point``, ``bar``, ``count``, ``box``, ``violin``, ``strip``}}
        The kind of plot to draw.
    {size}
    {aspect}
    {orient}
    {color}
    {palette}
    legend : bool, optional
        If ``True`` and there is a ``hue`` variable, draw a legend on the plot.
    {legend_out}
    {share_xy}
    {margin_titles}
    facet_kws : dict, optional
        Dictionary of other keyword arguments to pass to :class:`FacetGrid`.
    kwargs : key, value pairings
        Other keyword arguments are passed through to the underlying plotting
        function.

    Returns
    -------
    g : :class:`FacetGrid`
        Returns the :class:`FacetGrid` object with the plot on it for further
        tweaking.

    Examples
    --------

    Draw a single facet to use the :class:`FacetGrid` legend placement:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set(style="ticks")
        >>> exercise = sns.load_dataset("exercise")
        >>> g = sns.factorplot(x="time", y="pulse", hue="kind", data=exercise)

    Use a different plot kind to visualize the same data:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="time", y="pulse", hue="kind",
        ...                    data=exercise, kind="violin")

    Facet along the columns to show a third categorical variable:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="time", y="pulse", hue="kind",
        ...                    col="diet", data=exercise)

    Use a different size and aspect ratio for the facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="time", y="pulse", hue="kind",
        ...                    col="diet", data=exercise,
        ...                    size=5, aspect=.8)

    Make many column facets and wrap them into the rows of the grid:

    .. plot::
        :context: close-figs

        >>> titanic = sns.load_dataset("titanic")
        >>> g = sns.factorplot("alive", col="deck", col_wrap=4,
        ...                    data=titanic[titanic.deck.notnull()],
        ...                    kind="count", size=2.5, aspect=.8)

    Plot horizontally and pass other keyword arguments to the plot function:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="age", y="embark_town",
        ...                    hue="sex", row="class",
        ...                    data=titanic[titanic.embark_town.notnull()],
        ...                    orient="h", size=2, aspect=3.5, palette="Set3",
        ...                    kind="violin", dodge=True, cut=0, bw=.2)

    Use methods on the returned :class:`FacetGrid` to tweak the presentation:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="who", y="survived", col="class",
        ...                    data=titanic, saturation=.5,
        ...                    kind="bar", ci=None, aspect=.6)
        >>> (g.set_axis_labels("", "Survival Rate")
        ...   .set_xticklabels(["Men", "Women", "Children"])
        ...   .set_titles("{{col_name}} {{col_var}}")
        ...   .set(ylim=(0, 1))
        ...   .despine(left=True))  #doctest: +ELLIPSIS
        <seaborn.axisgrid.FacetGrid object at 0x...>

    """).format(**_categorical_docs)


def lvplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
           orient=None, color=None, palette=None, saturation=.75,
           width=.8, dodge=True, k_depth='proportion', linewidth=None,
           scale='exponential', outlier_prop=None, ax=None, **kwargs):

    plotter = _LVPlotter(x, y, hue, data, order, hue_order,
                         orient, color, palette, saturation,
                         width, dodge, k_depth, linewidth, scale, outlier_prop)

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax, kwargs)
    return ax

lvplot.__doc__ = dedent("""\
    Draw a letter value plot to show distributions of large datasets.

    Letter value (LV) plots are non-parametric estimates of the distribution of
    a dataset, similar to boxplots. LV plots are also similar to violin plots
    but without the need to fit a kernel density estimate. Thus, LV plots are
    fast to generate, directly interpretable in terms of the distribution of
    data, and easy to understand. For a more extensive explanation of letter
    value plots and their properties, see Hadley Wickham's excellent paper on
    the topic:

    http://vita.had.co.nz/papers/letter-value-plot.html

    {main_api_narrative}

    Parameters
    ----------
    {input_params}
    {categorical_data}
    {order_vars}
    {orient}
    {color}
    {palette}
    {saturation}
    {width}
    {dodge}
    k_depth : "proportion" | "tukey" | "trustworthy", optional
        The number of boxes, and by extension number of percentiles, to draw.
        All methods are detailed in Wickham's paper. Each makes different
        assumptions about the number of outliers and leverages different
        statistical properties.
    {linewidth}
    scale : "linear" | "exonential" | "area"
        Method to use for the width of the letter value boxes. All give similar
        results visually. "linear" reduces the width by a constant linear
        factor, "exponential" uses the proportion of data not covered, "area"
        is proportional to the percentage of data covered.
    outlier_prop : float, optional
        Proportion of data believed to be outliers. Is used in conjuction with
        k_depth to determine the number of percentiles to draw. Defaults to
        0.007 as a proportion of outliers. Should be in range [0, 1].
    {ax_in}
    kwargs : key, value mappings
        Other keyword arguments are passed through to ``plt.plot`` and
        ``plt.scatter`` at draw time.

    Returns
    -------
    {ax_out}

    See Also
    --------
    {violinplot}
    {boxplot}

    Examples
    --------

    Draw a single horizontal letter value plot:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns
        >>> sns.set_style("whitegrid")
        >>> tips = sns.load_dataset("tips")
        >>> ax = sns.lvplot(x=tips["total_bill"])

    Draw a vertical letter value plot grouped by a categorical variable:

    .. plot::
        :context: close-figs

        >>> ax = sns.lvplot(x="day", y="total_bill", data=tips)

    Draw a letter value plot with nested grouping by two categorical variables:

    .. plot::
        :context: close-figs

        >>> ax = sns.lvplot(x="day", y="total_bill", hue="smoker",
        ...                 data=tips, palette="Set3")

    Draw a letter value plot with nested grouping when some bins are empty:

    .. plot::
        :context: close-figs

        >>> ax = sns.lvplot(x="day", y="total_bill", hue="time",
        ...                 data=tips, linewidth=2.5)

    Control box order by passing an explicit order:

    .. plot::
        :context: close-figs

        >>> ax = sns.lvplot(x="time", y="tip", data=tips,
        ...                 order=["Dinner", "Lunch"])

    Draw a letter value plot for each numeric variable in a DataFrame:

    .. plot::
        :context: close-figs

        >>> iris = sns.load_dataset("iris")
        >>> ax = sns.lvplot(data=iris, orient="h", palette="Set2")

    Use :func:`stripplot` to show the datapoints on top of the boxes:

    .. plot::
        :context: close-figs

        >>> ax = sns.lvplot(x="day", y="total_bill", data=tips)
        >>> ax = sns.stripplot(x="day", y="total_bill", data=tips,
        ...                    size=4, jitter=True, color="gray")

    Use :func:`factorplot` to combine a :func:`lvplot` and a
    :class:`FacetGrid`. This allows grouping within additional categorical
    variables. Using :func:`factorplot` is safer than using :class:`FacetGrid`
    directly, as it ensures synchronization of variable order across facets:

    .. plot::
        :context: close-figs

        >>> g = sns.factorplot(x="sex", y="total_bill",
        ...                    hue="smoker", col="time",
        ...                    data=tips, kind="lv",
        ...                    size=4, aspect=.7);

    """).format(**_categorical_docs)

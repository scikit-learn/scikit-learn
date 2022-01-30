import tempfile
import copy

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy

import numpy.testing as npt
try:
    import pandas.testing as pdt
except ImportError:
    import pandas.util.testing as pdt
import pytest

from .. import matrix as mat
from .. import color_palette
from .._testing import assert_colors_equal

try:
    import fastcluster

    assert fastcluster
    _no_fastcluster = False
except ImportError:
    _no_fastcluster = True


# Copied from master onto v0.11 here to fix break introduced by
# cherry pick commit 49fbd353

class TestHeatmap:
    rs = np.random.RandomState(sum(map(ord, "heatmap")))

    x_norm = rs.randn(4, 8)
    letters = pd.Series(["A", "B", "C", "D"], name="letters")
    df_norm = pd.DataFrame(x_norm, index=letters)

    x_unif = rs.rand(20, 13)
    df_unif = pd.DataFrame(x_unif)

    default_kws = dict(vmin=None, vmax=None, cmap=None, center=None,
                       robust=False, annot=False, fmt=".2f", annot_kws=None,
                       cbar=True, cbar_kws=None, mask=None)

    def test_ndarray_input(self):

        p = mat._HeatMapper(self.x_norm, **self.default_kws)
        npt.assert_array_equal(p.plot_data, self.x_norm)
        pdt.assert_frame_equal(p.data, pd.DataFrame(self.x_norm))

        npt.assert_array_equal(p.xticklabels, np.arange(8))
        npt.assert_array_equal(p.yticklabels, np.arange(4))

        assert p.xlabel == ""
        assert p.ylabel == ""

    def test_df_input(self):

        p = mat._HeatMapper(self.df_norm, **self.default_kws)
        npt.assert_array_equal(p.plot_data, self.x_norm)
        pdt.assert_frame_equal(p.data, self.df_norm)

        npt.assert_array_equal(p.xticklabels, np.arange(8))
        npt.assert_array_equal(p.yticklabels, self.letters.values)

        assert p.xlabel == ""
        assert p.ylabel == "letters"

    def test_df_multindex_input(self):

        df = self.df_norm.copy()
        index = pd.MultiIndex.from_tuples([("A", 1), ("B", 2),
                                           ("C", 3), ("D", 4)],
                                          names=["letter", "number"])
        index.name = "letter-number"
        df.index = index

        p = mat._HeatMapper(df, **self.default_kws)

        combined_tick_labels = ["A-1", "B-2", "C-3", "D-4"]
        npt.assert_array_equal(p.yticklabels, combined_tick_labels)
        assert p.ylabel == "letter-number"

        p = mat._HeatMapper(df.T, **self.default_kws)

        npt.assert_array_equal(p.xticklabels, combined_tick_labels)
        assert p.xlabel == "letter-number"

    @pytest.mark.parametrize("dtype", [float, np.int64, object])
    def test_mask_input(self, dtype):
        kws = self.default_kws.copy()

        mask = self.x_norm > 0
        kws['mask'] = mask
        data = self.x_norm.astype(dtype)
        p = mat._HeatMapper(data, **kws)
        plot_data = np.ma.masked_where(mask, data)

        npt.assert_array_equal(p.plot_data, plot_data)

    def test_mask_limits(self):
        """Make sure masked cells are not used to calculate extremes"""

        kws = self.default_kws.copy()

        mask = self.x_norm > 0
        kws['mask'] = mask
        p = mat._HeatMapper(self.x_norm, **kws)

        assert p.vmax == np.ma.array(self.x_norm, mask=mask).max()
        assert p.vmin == np.ma.array(self.x_norm, mask=mask).min()

        mask = self.x_norm < 0
        kws['mask'] = mask
        p = mat._HeatMapper(self.x_norm, **kws)

        assert p.vmin == np.ma.array(self.x_norm, mask=mask).min()
        assert p.vmax == np.ma.array(self.x_norm, mask=mask).max()

    def test_default_vlims(self):

        p = mat._HeatMapper(self.df_unif, **self.default_kws)
        assert p.vmin == self.x_unif.min()
        assert p.vmax == self.x_unif.max()

    def test_robust_vlims(self):

        kws = self.default_kws.copy()
        kws["robust"] = True
        p = mat._HeatMapper(self.df_unif, **kws)

        assert p.vmin == np.percentile(self.x_unif, 2)
        assert p.vmax == np.percentile(self.x_unif, 98)

    def test_custom_sequential_vlims(self):

        kws = self.default_kws.copy()
        kws["vmin"] = 0
        kws["vmax"] = 1
        p = mat._HeatMapper(self.df_unif, **kws)

        assert p.vmin == 0
        assert p.vmax == 1

    def test_custom_diverging_vlims(self):

        kws = self.default_kws.copy()
        kws["vmin"] = -4
        kws["vmax"] = 5
        kws["center"] = 0
        p = mat._HeatMapper(self.df_norm, **kws)

        assert p.vmin == -4
        assert p.vmax == 5

    def test_array_with_nans(self):

        x1 = self.rs.rand(10, 10)
        nulls = np.zeros(10) * np.nan
        x2 = np.c_[x1, nulls]

        m1 = mat._HeatMapper(x1, **self.default_kws)
        m2 = mat._HeatMapper(x2, **self.default_kws)

        assert m1.vmin == m2.vmin
        assert m1.vmax == m2.vmax

    def test_mask(self):

        df = pd.DataFrame(data={'a': [1, 1, 1],
                                'b': [2, np.nan, 2],
                                'c': [3, 3, np.nan]})

        kws = self.default_kws.copy()
        kws["mask"] = np.isnan(df.values)

        m = mat._HeatMapper(df, **kws)

        npt.assert_array_equal(np.isnan(m.plot_data.data),
                               m.plot_data.mask)

    def test_custom_cmap(self):

        kws = self.default_kws.copy()
        kws["cmap"] = "BuGn"
        p = mat._HeatMapper(self.df_unif, **kws)
        assert p.cmap == mpl.cm.BuGn

    def test_centered_vlims(self):

        kws = self.default_kws.copy()
        kws["center"] = .5

        p = mat._HeatMapper(self.df_unif, **kws)

        assert p.vmin == self.df_unif.values.min()
        assert p.vmax == self.df_unif.values.max()

    def test_default_colors(self):

        vals = np.linspace(.2, 1, 9)
        cmap = mpl.cm.binary
        ax = mat.heatmap([vals], cmap=cmap)
        fc = ax.collections[0].get_facecolors()
        cvals = np.linspace(0, 1, 9)
        npt.assert_array_almost_equal(fc, cmap(cvals), 2)

    def test_custom_vlim_colors(self):

        vals = np.linspace(.2, 1, 9)
        cmap = mpl.cm.binary
        ax = mat.heatmap([vals], vmin=0, cmap=cmap)
        fc = ax.collections[0].get_facecolors()
        npt.assert_array_almost_equal(fc, cmap(vals), 2)

    def test_custom_center_colors(self):

        vals = np.linspace(.2, 1, 9)
        cmap = mpl.cm.binary
        ax = mat.heatmap([vals], center=.5, cmap=cmap)
        fc = ax.collections[0].get_facecolors()
        npt.assert_array_almost_equal(fc, cmap(vals), 2)

    def test_cmap_with_properties(self):

        kws = self.default_kws.copy()
        cmap = copy.copy(mpl.cm.get_cmap("BrBG"))
        cmap.set_bad("red")
        kws["cmap"] = cmap
        hm = mat._HeatMapper(self.df_unif, **kws)
        npt.assert_array_equal(
            cmap(np.ma.masked_invalid([np.nan])),
            hm.cmap(np.ma.masked_invalid([np.nan])))

        kws["center"] = 0.5
        hm = mat._HeatMapper(self.df_unif, **kws)
        npt.assert_array_equal(
            cmap(np.ma.masked_invalid([np.nan])),
            hm.cmap(np.ma.masked_invalid([np.nan])))

        kws = self.default_kws.copy()
        cmap = copy.copy(mpl.cm.get_cmap("BrBG"))
        cmap.set_under("red")
        kws["cmap"] = cmap
        hm = mat._HeatMapper(self.df_unif, **kws)
        npt.assert_array_equal(cmap(-np.inf), hm.cmap(-np.inf))

        kws["center"] = .5
        hm = mat._HeatMapper(self.df_unif, **kws)
        npt.assert_array_equal(cmap(-np.inf), hm.cmap(-np.inf))

        kws = self.default_kws.copy()
        cmap = copy.copy(mpl.cm.get_cmap("BrBG"))
        cmap.set_over("red")
        kws["cmap"] = cmap
        hm = mat._HeatMapper(self.df_unif, **kws)
        npt.assert_array_equal(cmap(-np.inf), hm.cmap(-np.inf))

        kws["center"] = .5
        hm = mat._HeatMapper(self.df_unif, **kws)
        npt.assert_array_equal(cmap(np.inf), hm.cmap(np.inf))

    def test_tickabels_off(self):
        kws = self.default_kws.copy()
        kws['xticklabels'] = False
        kws['yticklabels'] = False
        p = mat._HeatMapper(self.df_norm, **kws)
        assert p.xticklabels == []
        assert p.yticklabels == []

    def test_custom_ticklabels(self):
        kws = self.default_kws.copy()
        xticklabels = list('iheartheatmaps'[:self.df_norm.shape[1]])
        yticklabels = list('heatmapsarecool'[:self.df_norm.shape[0]])
        kws['xticklabels'] = xticklabels
        kws['yticklabels'] = yticklabels
        p = mat._HeatMapper(self.df_norm, **kws)
        assert p.xticklabels == xticklabels
        assert p.yticklabels == yticklabels

    def test_custom_ticklabel_interval(self):

        kws = self.default_kws.copy()
        xstep, ystep = 2, 3
        kws['xticklabels'] = xstep
        kws['yticklabels'] = ystep
        p = mat._HeatMapper(self.df_norm, **kws)

        nx, ny = self.df_norm.T.shape
        npt.assert_array_equal(p.xticks, np.arange(0, nx, xstep) + .5)
        npt.assert_array_equal(p.yticks, np.arange(0, ny, ystep) + .5)
        npt.assert_array_equal(p.xticklabels,
                               self.df_norm.columns[0:nx:xstep])
        npt.assert_array_equal(p.yticklabels,
                               self.df_norm.index[0:ny:ystep])

    def test_heatmap_annotation(self):

        ax = mat.heatmap(self.df_norm, annot=True, fmt=".1f",
                         annot_kws={"fontsize": 14})
        for val, text in zip(self.x_norm.flat, ax.texts):
            assert text.get_text() == "{:.1f}".format(val)
            assert text.get_fontsize() == 14

    def test_heatmap_annotation_overwrite_kws(self):

        annot_kws = dict(color="0.3", va="bottom", ha="left")
        ax = mat.heatmap(self.df_norm, annot=True, fmt=".1f",
                         annot_kws=annot_kws)
        for text in ax.texts:
            assert text.get_color() == "0.3"
            assert text.get_ha() == "left"
            assert text.get_va() == "bottom"

    def test_heatmap_annotation_with_mask(self):

        df = pd.DataFrame(data={'a': [1, 1, 1],
                                'b': [2, np.nan, 2],
                                'c': [3, 3, np.nan]})
        mask = np.isnan(df.values)
        df_masked = np.ma.masked_where(mask, df)
        ax = mat.heatmap(df, annot=True, fmt='.1f', mask=mask)
        assert len(df_masked.compressed()) == len(ax.texts)
        for val, text in zip(df_masked.compressed(), ax.texts):
            assert "{:.1f}".format(val) == text.get_text()

    def test_heatmap_annotation_mesh_colors(self):

        ax = mat.heatmap(self.df_norm, annot=True)
        mesh = ax.collections[0]
        assert len(mesh.get_facecolors()) == self.df_norm.values.size

        plt.close("all")

    def test_heatmap_annotation_other_data(self):
        annot_data = self.df_norm + 10

        ax = mat.heatmap(self.df_norm, annot=annot_data, fmt=".1f",
                         annot_kws={"fontsize": 14})

        for val, text in zip(annot_data.values.flat, ax.texts):
            assert text.get_text() == "{:.1f}".format(val)
            assert text.get_fontsize() == 14

    def test_heatmap_annotation_with_limited_ticklabels(self):
        ax = mat.heatmap(self.df_norm, fmt=".2f", annot=True,
                         xticklabels=False, yticklabels=False)
        for val, text in zip(self.x_norm.flat, ax.texts):
            assert text.get_text() == "{:.2f}".format(val)

    def test_heatmap_cbar(self):

        f = plt.figure()
        mat.heatmap(self.df_norm)
        assert len(f.axes) == 2
        plt.close(f)

        f = plt.figure()
        mat.heatmap(self.df_norm, cbar=False)
        assert len(f.axes) == 1
        plt.close(f)

        f, (ax1, ax2) = plt.subplots(2)
        mat.heatmap(self.df_norm, ax=ax1, cbar_ax=ax2)
        assert len(f.axes) == 2
        plt.close(f)

    @pytest.mark.xfail(mpl.__version__ == "3.1.1",
                       reason="matplotlib 3.1.1 bug")
    def test_heatmap_axes(self):

        ax = mat.heatmap(self.df_norm)

        xtl = [int(l.get_text()) for l in ax.get_xticklabels()]
        assert xtl == list(self.df_norm.columns)
        ytl = [l.get_text() for l in ax.get_yticklabels()]
        assert ytl == list(self.df_norm.index)

        assert ax.get_xlabel() == ""
        assert ax.get_ylabel() == "letters"

        assert ax.get_xlim() == (0, 8)
        assert ax.get_ylim() == (4, 0)

    def test_heatmap_ticklabel_rotation(self):

        f, ax = plt.subplots(figsize=(2, 2))
        mat.heatmap(self.df_norm, xticklabels=1, yticklabels=1, ax=ax)

        for t in ax.get_xticklabels():
            assert t.get_rotation() == 0

        for t in ax.get_yticklabels():
            assert t.get_rotation() == 90

        plt.close(f)

        df = self.df_norm.copy()
        df.columns = [str(c) * 10 for c in df.columns]
        df.index = [i * 10 for i in df.index]

        f, ax = plt.subplots(figsize=(2, 2))
        mat.heatmap(df, xticklabels=1, yticklabels=1, ax=ax)

        for t in ax.get_xticklabels():
            assert t.get_rotation() == 90

        for t in ax.get_yticklabels():
            assert t.get_rotation() == 0

        plt.close(f)

    def test_heatmap_inner_lines(self):

        c = (0, 0, 1, 1)
        ax = mat.heatmap(self.df_norm, linewidths=2, linecolor=c)
        mesh = ax.collections[0]
        assert mesh.get_linewidths()[0] == 2
        assert tuple(mesh.get_edgecolor()[0]) == c

    def test_square_aspect(self):

        ax = mat.heatmap(self.df_norm, square=True)
        obs_aspect = ax.get_aspect()
        # mpl>3.3 returns 1 for setting "equal" aspect
        # so test for the two possible equal outcomes
        assert obs_aspect == "equal" or obs_aspect == 1

    def test_mask_validation(self):

        mask = mat._matrix_mask(self.df_norm, None)
        assert mask.shape == self.df_norm.shape
        assert mask.values.sum() == 0

        with pytest.raises(ValueError):
            bad_array_mask = self.rs.randn(3, 6) > 0
            mat._matrix_mask(self.df_norm, bad_array_mask)

        with pytest.raises(ValueError):
            bad_df_mask = pd.DataFrame(self.rs.randn(4, 8) > 0)
            mat._matrix_mask(self.df_norm, bad_df_mask)

    def test_missing_data_mask(self):

        data = pd.DataFrame(np.arange(4, dtype=float).reshape(2, 2))
        data.loc[0, 0] = np.nan
        mask = mat._matrix_mask(data, None)
        npt.assert_array_equal(mask, [[True, False], [False, False]])

        mask_in = np.array([[False, True], [False, False]])
        mask_out = mat._matrix_mask(data, mask_in)
        npt.assert_array_equal(mask_out, [[True, True], [False, False]])

    def test_cbar_ticks(self):

        f, (ax1, ax2) = plt.subplots(2)
        mat.heatmap(self.df_norm, ax=ax1, cbar_ax=ax2,
                    cbar_kws=dict(drawedges=True))
        assert len(ax2.collections) == 2


class TestDendrogram:
    rs = np.random.RandomState(sum(map(ord, "dendrogram")))

    x_norm = rs.randn(4, 8) + np.arange(8)
    x_norm = (x_norm.T + np.arange(4)).T
    letters = pd.Series(["A", "B", "C", "D", "E", "F", "G", "H"],
                        name="letters")

    df_norm = pd.DataFrame(x_norm, columns=letters)
    try:
        import fastcluster

        x_norm_linkage = fastcluster.linkage_vector(x_norm.T,
                                                    metric='euclidean',
                                                    method='single')
    except ImportError:
        x_norm_distances = distance.pdist(x_norm.T, metric='euclidean')
        x_norm_linkage = hierarchy.linkage(x_norm_distances, method='single')
    x_norm_dendrogram = hierarchy.dendrogram(x_norm_linkage, no_plot=True,
                                             color_threshold=-np.inf)
    x_norm_leaves = x_norm_dendrogram['leaves']
    df_norm_leaves = np.asarray(df_norm.columns[x_norm_leaves])

    default_kws = dict(linkage=None, metric='euclidean', method='single',
                       axis=1, label=True, rotate=False)

    def test_ndarray_input(self):
        p = mat._DendrogramPlotter(self.x_norm, **self.default_kws)
        npt.assert_array_equal(p.array.T, self.x_norm)
        pdt.assert_frame_equal(p.data.T, pd.DataFrame(self.x_norm))

        npt.assert_array_equal(p.linkage, self.x_norm_linkage)
        assert p.dendrogram == self.x_norm_dendrogram

        npt.assert_array_equal(p.reordered_ind, self.x_norm_leaves)

        npt.assert_array_equal(p.xticklabels, self.x_norm_leaves)
        npt.assert_array_equal(p.yticklabels, [])

        assert p.xlabel is None
        assert p.ylabel == ''

    def test_df_input(self):
        p = mat._DendrogramPlotter(self.df_norm, **self.default_kws)
        npt.assert_array_equal(p.array.T, np.asarray(self.df_norm))
        pdt.assert_frame_equal(p.data.T, self.df_norm)

        npt.assert_array_equal(p.linkage, self.x_norm_linkage)
        assert p.dendrogram == self.x_norm_dendrogram

        npt.assert_array_equal(p.xticklabels,
                               np.asarray(self.df_norm.columns)[
                                   self.x_norm_leaves])
        npt.assert_array_equal(p.yticklabels, [])

        assert p.xlabel == 'letters'
        assert p.ylabel == ''

    def test_df_multindex_input(self):

        df = self.df_norm.copy()
        index = pd.MultiIndex.from_tuples([("A", 1), ("B", 2),
                                           ("C", 3), ("D", 4)],
                                          names=["letter", "number"])
        index.name = "letter-number"
        df.index = index
        kws = self.default_kws.copy()
        kws['label'] = True

        p = mat._DendrogramPlotter(df.T, **kws)

        xticklabels = ["A-1", "B-2", "C-3", "D-4"]
        xticklabels = [xticklabels[i] for i in p.reordered_ind]
        npt.assert_array_equal(p.xticklabels, xticklabels)
        npt.assert_array_equal(p.yticklabels, [])
        assert p.xlabel == "letter-number"

    def test_axis0_input(self):
        kws = self.default_kws.copy()
        kws['axis'] = 0
        p = mat._DendrogramPlotter(self.df_norm.T, **kws)

        npt.assert_array_equal(p.array, np.asarray(self.df_norm.T))
        pdt.assert_frame_equal(p.data, self.df_norm.T)

        npt.assert_array_equal(p.linkage, self.x_norm_linkage)
        assert p.dendrogram == self.x_norm_dendrogram

        npt.assert_array_equal(p.xticklabels, self.df_norm_leaves)
        npt.assert_array_equal(p.yticklabels, [])

        assert p.xlabel == 'letters'
        assert p.ylabel == ''

    def test_rotate_input(self):
        kws = self.default_kws.copy()
        kws['rotate'] = True
        p = mat._DendrogramPlotter(self.df_norm, **kws)
        npt.assert_array_equal(p.array.T, np.asarray(self.df_norm))
        pdt.assert_frame_equal(p.data.T, self.df_norm)

        npt.assert_array_equal(p.xticklabels, [])
        npt.assert_array_equal(p.yticklabels, self.df_norm_leaves)

        assert p.xlabel == ''
        assert p.ylabel == 'letters'

    def test_rotate_axis0_input(self):
        kws = self.default_kws.copy()
        kws['rotate'] = True
        kws['axis'] = 0
        p = mat._DendrogramPlotter(self.df_norm.T, **kws)

        npt.assert_array_equal(p.reordered_ind, self.x_norm_leaves)

    def test_custom_linkage(self):
        kws = self.default_kws.copy()

        try:
            import fastcluster

            linkage = fastcluster.linkage_vector(self.x_norm, method='single',
                                                 metric='euclidean')
        except ImportError:
            d = distance.pdist(self.x_norm, metric='euclidean')
            linkage = hierarchy.linkage(d, method='single')
        dendrogram = hierarchy.dendrogram(linkage, no_plot=True,
                                          color_threshold=-np.inf)
        kws['linkage'] = linkage
        p = mat._DendrogramPlotter(self.df_norm, **kws)

        npt.assert_array_equal(p.linkage, linkage)
        assert p.dendrogram == dendrogram

    def test_label_false(self):
        kws = self.default_kws.copy()
        kws['label'] = False
        p = mat._DendrogramPlotter(self.df_norm, **kws)
        assert p.xticks == []
        assert p.yticks == []
        assert p.xticklabels == []
        assert p.yticklabels == []
        assert p.xlabel == ""
        assert p.ylabel == ""

    def test_linkage_scipy(self):
        p = mat._DendrogramPlotter(self.x_norm, **self.default_kws)

        scipy_linkage = p._calculate_linkage_scipy()

        from scipy.spatial import distance
        from scipy.cluster import hierarchy

        dists = distance.pdist(self.x_norm.T,
                               metric=self.default_kws['metric'])
        linkage = hierarchy.linkage(dists, method=self.default_kws['method'])

        npt.assert_array_equal(scipy_linkage, linkage)

    @pytest.mark.skipif(_no_fastcluster, reason="fastcluster not installed")
    def test_fastcluster_other_method(self):
        import fastcluster

        kws = self.default_kws.copy()
        kws['method'] = 'average'
        linkage = fastcluster.linkage(self.x_norm.T, method='average',
                                      metric='euclidean')
        p = mat._DendrogramPlotter(self.x_norm, **kws)
        npt.assert_array_equal(p.linkage, linkage)

    @pytest.mark.skipif(_no_fastcluster, reason="fastcluster not installed")
    def test_fastcluster_non_euclidean(self):
        import fastcluster

        kws = self.default_kws.copy()
        kws['metric'] = 'cosine'
        kws['method'] = 'average'
        linkage = fastcluster.linkage(self.x_norm.T, method=kws['method'],
                                      metric=kws['metric'])
        p = mat._DendrogramPlotter(self.x_norm, **kws)
        npt.assert_array_equal(p.linkage, linkage)

    def test_dendrogram_plot(self):
        d = mat.dendrogram(self.x_norm, **self.default_kws)

        ax = plt.gca()
        xlim = ax.get_xlim()
        # 10 comes from _plot_dendrogram in scipy.cluster.hierarchy
        xmax = len(d.reordered_ind) * 10

        assert xlim[0] == 0
        assert xlim[1] == xmax

        assert len(ax.collections[0].get_paths()) == len(d.dependent_coord)

    @pytest.mark.xfail(mpl.__version__ == "3.1.1",
                       reason="matplotlib 3.1.1 bug")
    def test_dendrogram_rotate(self):
        kws = self.default_kws.copy()
        kws['rotate'] = True

        d = mat.dendrogram(self.x_norm, **kws)

        ax = plt.gca()
        ylim = ax.get_ylim()

        # 10 comes from _plot_dendrogram in scipy.cluster.hierarchy
        ymax = len(d.reordered_ind) * 10

        # Since y axis is inverted, ylim is (80, 0)
        # and therefore not (0, 80) as usual:
        assert ylim[1] == 0
        assert ylim[0] == ymax

    def test_dendrogram_ticklabel_rotation(self):
        f, ax = plt.subplots(figsize=(2, 2))
        mat.dendrogram(self.df_norm, ax=ax)

        for t in ax.get_xticklabels():
            assert t.get_rotation() == 0

        plt.close(f)

        df = self.df_norm.copy()
        df.columns = [str(c) * 10 for c in df.columns]
        df.index = [i * 10 for i in df.index]

        f, ax = plt.subplots(figsize=(2, 2))
        mat.dendrogram(df, ax=ax)

        for t in ax.get_xticklabels():
            assert t.get_rotation() == 90

        plt.close(f)

        f, ax = plt.subplots(figsize=(2, 2))
        mat.dendrogram(df.T, axis=0, rotate=True)
        for t in ax.get_yticklabels():
            assert t.get_rotation() == 0
        plt.close(f)


class TestClustermap:
    rs = np.random.RandomState(sum(map(ord, "clustermap")))

    x_norm = rs.randn(4, 8) + np.arange(8)
    x_norm = (x_norm.T + np.arange(4)).T
    letters = pd.Series(["A", "B", "C", "D", "E", "F", "G", "H"],
                        name="letters")

    df_norm = pd.DataFrame(x_norm, columns=letters)
    try:
        import fastcluster

        x_norm_linkage = fastcluster.linkage_vector(x_norm.T,
                                                    metric='euclidean',
                                                    method='single')
    except ImportError:
        x_norm_distances = distance.pdist(x_norm.T, metric='euclidean')
        x_norm_linkage = hierarchy.linkage(x_norm_distances, method='single')
    x_norm_dendrogram = hierarchy.dendrogram(x_norm_linkage, no_plot=True,
                                             color_threshold=-np.inf)
    x_norm_leaves = x_norm_dendrogram['leaves']
    df_norm_leaves = np.asarray(df_norm.columns[x_norm_leaves])

    default_kws = dict(pivot_kws=None, z_score=None, standard_scale=None,
                       figsize=(10, 10), row_colors=None, col_colors=None,
                       dendrogram_ratio=.2, colors_ratio=.03,
                       cbar_pos=(0, .8, .05, .2))

    default_plot_kws = dict(metric='euclidean', method='average',
                            colorbar_kws=None,
                            row_cluster=True, col_cluster=True,
                            row_linkage=None, col_linkage=None,
                            tree_kws=None)

    row_colors = color_palette('Set2', df_norm.shape[0])
    col_colors = color_palette('Dark2', df_norm.shape[1])

    def test_ndarray_input(self):
        cg = mat.ClusterGrid(self.x_norm, **self.default_kws)
        pdt.assert_frame_equal(cg.data, pd.DataFrame(self.x_norm))
        assert len(cg.fig.axes) == 4
        assert cg.ax_row_colors is None
        assert cg.ax_col_colors is None

    def test_df_input(self):
        cg = mat.ClusterGrid(self.df_norm, **self.default_kws)
        pdt.assert_frame_equal(cg.data, self.df_norm)

    def test_corr_df_input(self):
        df = self.df_norm.corr()
        cg = mat.ClusterGrid(df, **self.default_kws)
        cg.plot(**self.default_plot_kws)
        diag = cg.data2d.values[np.diag_indices_from(cg.data2d)]
        npt.assert_array_equal(diag, np.ones(cg.data2d.shape[0]))

    def test_pivot_input(self):
        df_norm = self.df_norm.copy()
        df_norm.index.name = 'numbers'
        df_long = pd.melt(df_norm.reset_index(), var_name='letters',
                          id_vars='numbers')
        kws = self.default_kws.copy()
        kws['pivot_kws'] = dict(index='numbers', columns='letters',
                                values='value')
        cg = mat.ClusterGrid(df_long, **kws)

        pdt.assert_frame_equal(cg.data2d, df_norm)

    def test_colors_input(self):
        kws = self.default_kws.copy()

        kws['row_colors'] = self.row_colors
        kws['col_colors'] = self.col_colors

        cg = mat.ClusterGrid(self.df_norm, **kws)
        npt.assert_array_equal(cg.row_colors, self.row_colors)
        npt.assert_array_equal(cg.col_colors, self.col_colors)

        assert len(cg.fig.axes) == 6

    def test_categorical_colors_input(self):
        kws = self.default_kws.copy()

        row_colors = pd.Series(self.row_colors, dtype="category")
        col_colors = pd.Series(
            self.col_colors, dtype="category", index=self.df_norm.columns
        )

        kws['row_colors'] = row_colors
        kws['col_colors'] = col_colors

        exp_row_colors = list(map(mpl.colors.to_rgb, row_colors))
        exp_col_colors = list(map(mpl.colors.to_rgb, col_colors))

        cg = mat.ClusterGrid(self.df_norm, **kws)
        npt.assert_array_equal(cg.row_colors, exp_row_colors)
        npt.assert_array_equal(cg.col_colors, exp_col_colors)

        assert len(cg.fig.axes) == 6

    def test_nested_colors_input(self):
        kws = self.default_kws.copy()

        row_colors = [self.row_colors, self.row_colors]
        col_colors = [self.col_colors, self.col_colors]
        kws['row_colors'] = row_colors
        kws['col_colors'] = col_colors

        cm = mat.ClusterGrid(self.df_norm, **kws)
        npt.assert_array_equal(cm.row_colors, row_colors)
        npt.assert_array_equal(cm.col_colors, col_colors)

        assert len(cm.fig.axes) == 6

    def test_colors_input_custom_cmap(self):
        kws = self.default_kws.copy()

        kws['cmap'] = mpl.cm.PRGn
        kws['row_colors'] = self.row_colors
        kws['col_colors'] = self.col_colors

        cg = mat.clustermap(self.df_norm, **kws)
        npt.assert_array_equal(cg.row_colors, self.row_colors)
        npt.assert_array_equal(cg.col_colors, self.col_colors)

        assert len(cg.fig.axes) == 6

    def test_z_score(self):
        df = self.df_norm.copy()
        df = (df - df.mean()) / df.std()
        kws = self.default_kws.copy()
        kws['z_score'] = 1

        cg = mat.ClusterGrid(self.df_norm, **kws)
        pdt.assert_frame_equal(cg.data2d, df)

    def test_z_score_axis0(self):
        df = self.df_norm.copy()
        df = df.T
        df = (df - df.mean()) / df.std()
        df = df.T
        kws = self.default_kws.copy()
        kws['z_score'] = 0

        cg = mat.ClusterGrid(self.df_norm, **kws)
        pdt.assert_frame_equal(cg.data2d, df)

    def test_standard_scale(self):
        df = self.df_norm.copy()
        df = (df - df.min()) / (df.max() - df.min())
        kws = self.default_kws.copy()
        kws['standard_scale'] = 1

        cg = mat.ClusterGrid(self.df_norm, **kws)
        pdt.assert_frame_equal(cg.data2d, df)

    def test_standard_scale_axis0(self):
        df = self.df_norm.copy()
        df = df.T
        df = (df - df.min()) / (df.max() - df.min())
        df = df.T
        kws = self.default_kws.copy()
        kws['standard_scale'] = 0

        cg = mat.ClusterGrid(self.df_norm, **kws)
        pdt.assert_frame_equal(cg.data2d, df)

    def test_z_score_standard_scale(self):
        kws = self.default_kws.copy()
        kws['z_score'] = True
        kws['standard_scale'] = True
        with pytest.raises(ValueError):
            mat.ClusterGrid(self.df_norm, **kws)

    def test_color_list_to_matrix_and_cmap(self):
        # Note this uses the attribute named col_colors but tests row colors
        matrix, cmap = mat.ClusterGrid.color_list_to_matrix_and_cmap(
            self.col_colors, self.x_norm_leaves, axis=0)

        for i, leaf in enumerate(self.x_norm_leaves):
            color = self.col_colors[leaf]
            assert_colors_equal(cmap(matrix[i, 0]), color)

    def test_nested_color_list_to_matrix_and_cmap(self):
        # Note this uses the attribute named col_colors but tests row colors
        colors = [self.col_colors, self.col_colors[::-1]]
        matrix, cmap = mat.ClusterGrid.color_list_to_matrix_and_cmap(
            colors, self.x_norm_leaves, axis=0)

        for i, leaf in enumerate(self.x_norm_leaves):
            for j, color_row in enumerate(colors):
                color = color_row[leaf]
                assert_colors_equal(cmap(matrix[i, j]), color)

    def test_color_list_to_matrix_and_cmap_axis1(self):
        matrix, cmap = mat.ClusterGrid.color_list_to_matrix_and_cmap(
            self.col_colors, self.x_norm_leaves, axis=1)

        for j, leaf in enumerate(self.x_norm_leaves):
            color = self.col_colors[leaf]
            assert_colors_equal(cmap(matrix[0, j]), color)

    def test_color_list_to_matrix_and_cmap_different_sizes(self):
        colors = [self.col_colors, self.col_colors * 2]
        with pytest.raises(ValueError):
            matrix, cmap = mat.ClusterGrid.color_list_to_matrix_and_cmap(
                colors, self.x_norm_leaves, axis=1)

    def test_savefig(self):
        # Not sure if this is the right way to test....
        cg = mat.ClusterGrid(self.df_norm, **self.default_kws)
        cg.plot(**self.default_plot_kws)
        cg.savefig(tempfile.NamedTemporaryFile(), format='png')

    def test_plot_dendrograms(self):
        cm = mat.clustermap(self.df_norm, **self.default_kws)

        assert len(cm.ax_row_dendrogram.collections[0].get_paths()) == len(
            cm.dendrogram_row.independent_coord
        )
        assert len(cm.ax_col_dendrogram.collections[0].get_paths()) == len(
            cm.dendrogram_col.independent_coord
        )
        data2d = self.df_norm.iloc[cm.dendrogram_row.reordered_ind,
                                   cm.dendrogram_col.reordered_ind]
        pdt.assert_frame_equal(cm.data2d, data2d)

    def test_cluster_false(self):
        kws = self.default_kws.copy()
        kws['row_cluster'] = False
        kws['col_cluster'] = False

        cm = mat.clustermap(self.df_norm, **kws)
        assert len(cm.ax_row_dendrogram.lines) == 0
        assert len(cm.ax_col_dendrogram.lines) == 0

        assert len(cm.ax_row_dendrogram.get_xticks()) == 0
        assert len(cm.ax_row_dendrogram.get_yticks()) == 0
        assert len(cm.ax_col_dendrogram.get_xticks()) == 0
        assert len(cm.ax_col_dendrogram.get_yticks()) == 0

        pdt.assert_frame_equal(cm.data2d, self.df_norm)

    def test_row_col_colors(self):
        kws = self.default_kws.copy()
        kws['row_colors'] = self.row_colors
        kws['col_colors'] = self.col_colors

        cm = mat.clustermap(self.df_norm, **kws)

        assert len(cm.ax_row_colors.collections) == 1
        assert len(cm.ax_col_colors.collections) == 1

    def test_cluster_false_row_col_colors(self):
        kws = self.default_kws.copy()
        kws['row_cluster'] = False
        kws['col_cluster'] = False
        kws['row_colors'] = self.row_colors
        kws['col_colors'] = self.col_colors

        cm = mat.clustermap(self.df_norm, **kws)
        assert len(cm.ax_row_dendrogram.lines) == 0
        assert len(cm.ax_col_dendrogram.lines) == 0

        assert len(cm.ax_row_dendrogram.get_xticks()) == 0
        assert len(cm.ax_row_dendrogram.get_yticks()) == 0
        assert len(cm.ax_col_dendrogram.get_xticks()) == 0
        assert len(cm.ax_col_dendrogram.get_yticks()) == 0
        assert len(cm.ax_row_colors.collections) == 1
        assert len(cm.ax_col_colors.collections) == 1

        pdt.assert_frame_equal(cm.data2d, self.df_norm)

    def test_row_col_colors_df(self):
        kws = self.default_kws.copy()
        kws['row_colors'] = pd.DataFrame({'row_1': list(self.row_colors),
                                          'row_2': list(self.row_colors)},
                                         index=self.df_norm.index,
                                         columns=['row_1', 'row_2'])
        kws['col_colors'] = pd.DataFrame({'col_1': list(self.col_colors),
                                          'col_2': list(self.col_colors)},
                                         index=self.df_norm.columns,
                                         columns=['col_1', 'col_2'])

        cm = mat.clustermap(self.df_norm, **kws)

        row_labels = [l.get_text() for l in
                      cm.ax_row_colors.get_xticklabels()]
        assert cm.row_color_labels == ['row_1', 'row_2']
        assert row_labels == cm.row_color_labels

        col_labels = [l.get_text() for l in
                      cm.ax_col_colors.get_yticklabels()]
        assert cm.col_color_labels == ['col_1', 'col_2']
        assert col_labels == cm.col_color_labels

    def test_row_col_colors_df_shuffled(self):
        # Tests if colors are properly matched, even if given in wrong order

        m, n = self.df_norm.shape
        shuffled_inds = [self.df_norm.index[i] for i in
                         list(range(0, m, 2)) + list(range(1, m, 2))]
        shuffled_cols = [self.df_norm.columns[i] for i in
                         list(range(0, n, 2)) + list(range(1, n, 2))]

        kws = self.default_kws.copy()

        row_colors = pd.DataFrame({'row_annot': list(self.row_colors)},
                                  index=self.df_norm.index)
        kws['row_colors'] = row_colors.loc[shuffled_inds]

        col_colors = pd.DataFrame({'col_annot': list(self.col_colors)},
                                  index=self.df_norm.columns)
        kws['col_colors'] = col_colors.loc[shuffled_cols]

        cm = mat.clustermap(self.df_norm, **kws)
        assert list(cm.col_colors)[0] == list(self.col_colors)
        assert list(cm.row_colors)[0] == list(self.row_colors)

    def test_row_col_colors_df_missing(self):
        kws = self.default_kws.copy()
        row_colors = pd.DataFrame({'row_annot': list(self.row_colors)},
                                  index=self.df_norm.index)
        kws['row_colors'] = row_colors.drop(self.df_norm.index[0])

        col_colors = pd.DataFrame({'col_annot': list(self.col_colors)},
                                  index=self.df_norm.columns)
        kws['col_colors'] = col_colors.drop(self.df_norm.columns[0])

        cm = mat.clustermap(self.df_norm, **kws)

        assert list(cm.col_colors)[0] == [(1.0, 1.0, 1.0)] + list(self.col_colors[1:])
        assert list(cm.row_colors)[0] == [(1.0, 1.0, 1.0)] + list(self.row_colors[1:])

    def test_row_col_colors_df_one_axis(self):
        # Test case with only row annotation.
        kws1 = self.default_kws.copy()
        kws1['row_colors'] = pd.DataFrame({'row_1': list(self.row_colors),
                                           'row_2': list(self.row_colors)},
                                          index=self.df_norm.index,
                                          columns=['row_1', 'row_2'])

        cm1 = mat.clustermap(self.df_norm, **kws1)

        row_labels = [l.get_text() for l in
                      cm1.ax_row_colors.get_xticklabels()]
        assert cm1.row_color_labels == ['row_1', 'row_2']
        assert row_labels == cm1.row_color_labels

        # Test case with only col annotation.
        kws2 = self.default_kws.copy()
        kws2['col_colors'] = pd.DataFrame({'col_1': list(self.col_colors),
                                           'col_2': list(self.col_colors)},
                                          index=self.df_norm.columns,
                                          columns=['col_1', 'col_2'])

        cm2 = mat.clustermap(self.df_norm, **kws2)

        col_labels = [l.get_text() for l in
                      cm2.ax_col_colors.get_yticklabels()]
        assert cm2.col_color_labels == ['col_1', 'col_2']
        assert col_labels == cm2.col_color_labels

    def test_row_col_colors_series(self):
        kws = self.default_kws.copy()
        kws['row_colors'] = pd.Series(list(self.row_colors), name='row_annot',
                                      index=self.df_norm.index)
        kws['col_colors'] = pd.Series(list(self.col_colors), name='col_annot',
                                      index=self.df_norm.columns)

        cm = mat.clustermap(self.df_norm, **kws)

        row_labels = [l.get_text() for l in cm.ax_row_colors.get_xticklabels()]
        assert cm.row_color_labels == ['row_annot']
        assert row_labels == cm.row_color_labels

        col_labels = [l.get_text() for l in cm.ax_col_colors.get_yticklabels()]
        assert cm.col_color_labels == ['col_annot']
        assert col_labels == cm.col_color_labels

    def test_row_col_colors_series_shuffled(self):
        # Tests if colors are properly matched, even if given in wrong order

        m, n = self.df_norm.shape
        shuffled_inds = [self.df_norm.index[i] for i in
                         list(range(0, m, 2)) + list(range(1, m, 2))]
        shuffled_cols = [self.df_norm.columns[i] for i in
                         list(range(0, n, 2)) + list(range(1, n, 2))]

        kws = self.default_kws.copy()

        row_colors = pd.Series(list(self.row_colors), name='row_annot',
                               index=self.df_norm.index)
        kws['row_colors'] = row_colors.loc[shuffled_inds]

        col_colors = pd.Series(list(self.col_colors), name='col_annot',
                               index=self.df_norm.columns)
        kws['col_colors'] = col_colors.loc[shuffled_cols]

        cm = mat.clustermap(self.df_norm, **kws)

        assert list(cm.col_colors) == list(self.col_colors)
        assert list(cm.row_colors) == list(self.row_colors)

    def test_row_col_colors_series_missing(self):
        kws = self.default_kws.copy()
        row_colors = pd.Series(list(self.row_colors), name='row_annot',
                               index=self.df_norm.index)
        kws['row_colors'] = row_colors.drop(self.df_norm.index[0])

        col_colors = pd.Series(list(self.col_colors), name='col_annot',
                               index=self.df_norm.columns)
        kws['col_colors'] = col_colors.drop(self.df_norm.columns[0])

        cm = mat.clustermap(self.df_norm, **kws)
        assert list(cm.col_colors) == [(1.0, 1.0, 1.0)] + list(self.col_colors[1:])
        assert list(cm.row_colors) == [(1.0, 1.0, 1.0)] + list(self.row_colors[1:])

    def test_row_col_colors_ignore_heatmap_kwargs(self):

        g = mat.clustermap(self.rs.uniform(0, 200, self.df_norm.shape),
                           row_colors=self.row_colors,
                           col_colors=self.col_colors,
                           cmap="Spectral",
                           norm=mpl.colors.LogNorm(),
                           vmax=100)

        assert np.array_equal(
            np.array(self.row_colors)[g.dendrogram_row.reordered_ind],
            g.ax_row_colors.collections[0].get_facecolors()[:, :3]
        )

        assert np.array_equal(
            np.array(self.col_colors)[g.dendrogram_col.reordered_ind],
            g.ax_col_colors.collections[0].get_facecolors()[:, :3]
        )

    def test_row_col_colors_raise_on_mixed_index_types(self):

        row_colors = pd.Series(
            list(self.row_colors), name="row_annot", index=self.df_norm.index
        )

        col_colors = pd.Series(
            list(self.col_colors), name="col_annot", index=self.df_norm.columns
        )

        with pytest.raises(TypeError):
            mat.clustermap(self.x_norm, row_colors=row_colors)

        with pytest.raises(TypeError):
            mat.clustermap(self.x_norm, col_colors=col_colors)

    def test_mask_reorganization(self):

        kws = self.default_kws.copy()
        kws["mask"] = self.df_norm > 0

        g = mat.clustermap(self.df_norm, **kws)
        npt.assert_array_equal(g.data2d.index, g.mask.index)
        npt.assert_array_equal(g.data2d.columns, g.mask.columns)

        npt.assert_array_equal(g.mask.index,
                               self.df_norm.index[
                                   g.dendrogram_row.reordered_ind])
        npt.assert_array_equal(g.mask.columns,
                               self.df_norm.columns[
                                   g.dendrogram_col.reordered_ind])

    def test_ticklabel_reorganization(self):

        kws = self.default_kws.copy()
        xtl = np.arange(self.df_norm.shape[1])
        kws["xticklabels"] = list(xtl)
        ytl = self.letters.loc[:self.df_norm.shape[0]]
        kws["yticklabels"] = ytl

        g = mat.clustermap(self.df_norm, **kws)

        xtl_actual = [t.get_text() for t in g.ax_heatmap.get_xticklabels()]
        ytl_actual = [t.get_text() for t in g.ax_heatmap.get_yticklabels()]

        xtl_want = xtl[g.dendrogram_col.reordered_ind].astype("<U1")
        ytl_want = ytl[g.dendrogram_row.reordered_ind].astype("<U1")

        npt.assert_array_equal(xtl_actual, xtl_want)
        npt.assert_array_equal(ytl_actual, ytl_want)

    def test_noticklabels(self):

        kws = self.default_kws.copy()
        kws["xticklabels"] = False
        kws["yticklabels"] = False

        g = mat.clustermap(self.df_norm, **kws)

        xtl_actual = [t.get_text() for t in g.ax_heatmap.get_xticklabels()]
        ytl_actual = [t.get_text() for t in g.ax_heatmap.get_yticklabels()]
        assert xtl_actual == []
        assert ytl_actual == []

    def test_size_ratios(self):

        # The way that wspace/hspace work in GridSpec, the mapping from input
        # ratio to actual width/height of each axes is complicated, so this
        # test is just going to assert comparative relationships

        kws1 = self.default_kws.copy()
        kws1.update(dendrogram_ratio=.2, colors_ratio=.03,
                    col_colors=self.col_colors, row_colors=self.row_colors)

        kws2 = kws1.copy()
        kws2.update(dendrogram_ratio=.3, colors_ratio=.05)

        g1 = mat.clustermap(self.df_norm, **kws1)
        g2 = mat.clustermap(self.df_norm, **kws2)

        assert (g2.ax_col_dendrogram.get_position().height
                > g1.ax_col_dendrogram.get_position().height)

        assert (g2.ax_col_colors.get_position().height
                > g1.ax_col_colors.get_position().height)

        assert (g2.ax_heatmap.get_position().height
                < g1.ax_heatmap.get_position().height)

        assert (g2.ax_row_dendrogram.get_position().width
                > g1.ax_row_dendrogram.get_position().width)

        assert (g2.ax_row_colors.get_position().width
                > g1.ax_row_colors.get_position().width)

        assert (g2.ax_heatmap.get_position().width
                < g1.ax_heatmap.get_position().width)

        kws1 = self.default_kws.copy()
        kws1.update(col_colors=self.col_colors)
        kws2 = kws1.copy()
        kws2.update(col_colors=[self.col_colors, self.col_colors])

        g1 = mat.clustermap(self.df_norm, **kws1)
        g2 = mat.clustermap(self.df_norm, **kws2)

        assert (g2.ax_col_colors.get_position().height
                > g1.ax_col_colors.get_position().height)

        kws1 = self.default_kws.copy()
        kws1.update(dendrogram_ratio=(.2, .2))

        kws2 = kws1.copy()
        kws2.update(dendrogram_ratio=(.2, .3))

        g1 = mat.clustermap(self.df_norm, **kws1)
        g2 = mat.clustermap(self.df_norm, **kws2)

        assert (g2.ax_row_dendrogram.get_position().width
                == g1.ax_row_dendrogram.get_position().width)

        assert (g2.ax_col_dendrogram.get_position().height
                > g1.ax_col_dendrogram.get_position().height)

    def test_cbar_pos(self):

        kws = self.default_kws.copy()
        kws["cbar_pos"] = (.2, .1, .4, .3)

        g = mat.clustermap(self.df_norm, **kws)
        pos = g.ax_cbar.get_position()
        assert pytest.approx(tuple(pos.p0)) == kws["cbar_pos"][:2]
        assert pytest.approx(pos.width) == kws["cbar_pos"][2]
        assert pytest.approx(pos.height) == kws["cbar_pos"][3]

        kws["cbar_pos"] = None
        g = mat.clustermap(self.df_norm, **kws)
        assert g.ax_cbar is None

    def test_square_warning(self):

        kws = self.default_kws.copy()
        g1 = mat.clustermap(self.df_norm, **kws)

        with pytest.warns(UserWarning):
            kws["square"] = True
            g2 = mat.clustermap(self.df_norm, **kws)

        g1_shape = g1.ax_heatmap.get_position().get_points()
        g2_shape = g2.ax_heatmap.get_position().get_points()
        assert np.array_equal(g1_shape, g2_shape)

    def test_clustermap_annotation(self):

        g = mat.clustermap(self.df_norm, annot=True, fmt=".1f")
        for val, text in zip(np.asarray(g.data2d).flat, g.ax_heatmap.texts):
            assert text.get_text() == "{:.1f}".format(val)

        g = mat.clustermap(self.df_norm, annot=self.df_norm, fmt=".1f")
        for val, text in zip(np.asarray(g.data2d).flat, g.ax_heatmap.texts):
            assert text.get_text() == "{:.1f}".format(val)

    def test_tree_kws(self):

        rgb = (1, .5, .2)
        g = mat.clustermap(self.df_norm, tree_kws=dict(color=rgb))
        for ax in [g.ax_col_dendrogram, g.ax_row_dendrogram]:
            tree, = ax.collections
            assert tuple(tree.get_color().squeeze())[:3] == rgb

import numpy as np

from sklearn.utils._plot import plot_heatmap


def test_plot_heatmap(pyplot):
    """Check the overall heatmap plotting."""
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    fig, ax, im, text = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
    )

    assert ax.get_xlabel() == "xlabel"
    assert ax.get_ylabel() == "ylabel"
    for label in ax.get_xticklabels():
        assert label.get_text() in ["col_1", "col_2", "col_3"]
    for label in ax.get_yticklabels():
        assert label.get_text() in ["row_1", "row_2"]
    for annotation, expected_annotation in zip(text.ravel(), data.ravel()):
        assert annotation.get_text() == str(expected_annotation)
    assert im.get_clim() == (data.min(), data.max())
    assert im.get_cmap() == pyplot.cm.viridis
    assert len(fig.axes) == 2


def test_plot_heatmap_no_annotation(pyplot):
    """Check the heatmap plotting without annotations."""
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    _, _, _, text = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
        include_values=False,
    )

    assert text is None


def test_plot_heatmap_no_colorbar(pyplot):
    """Check the heatmap plotting without colorbar."""
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    fig, _, _, _ = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
        colorbar=False,
    )

    assert len(fig.axes) == 1


def test_plot_heatmap_cmap(pyplot):
    """Check that we set properly the heatmap when it is passed."""
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    _, _, im, _ = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
    )
    assert im.get_cmap() == pyplot.cm.viridis

    _, _, im, _ = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
        cmap="Reds",
    )
    assert im.get_cmap() == pyplot.cm.Reds


def test_plot_heatmap_values_format(pyplot):
    """Check that we can force the format of the printed values."""
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    _, ax, im, text = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
        values_format=".3f",
    )
    for annotation, expected_annotation in zip(text.ravel(), data.ravel()):
        assert annotation.get_text() == f"{expected_annotation:.3f}"


def test_plot_heatmap_axis(pyplot):
    """Check that we don't create a new axis if we pass one."""
    _, ax = pyplot.subplots()
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    _, ax_out, im, text = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
        values_format=".3f",
        ax=ax,
    )
    assert ax is ax_out


def test_plot_heatmap_im_kw(pyplot):
    """Check that we can pass extra keyword to the underlying `imshow` function."""
    data = np.array([[0.1, 0.3, 0.4], [0.2, 0.4, 0.5]])
    cmap = "RdBu"
    _, _, im, _ = plot_heatmap(
        data,
        ylabel="ylabel",
        xlabel="xlabel",
        yticklabels=["row_1", "row_2"],
        xticklabels=["col_1", "col_2", "col_3"],
        im_kw={"cmap": cmap},
    )
    assert im.cmap.name == cmap

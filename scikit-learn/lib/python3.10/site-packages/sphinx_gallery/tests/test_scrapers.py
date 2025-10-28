"""Testing image scrapers."""

import os

import pytest
from sphinx.errors import ConfigError, ExtensionError

import sphinx_gallery
from sphinx_gallery.gen_gallery import _fill_gallery_conf_defaults
from sphinx_gallery.scrapers import (
    _KNOWN_IMG_EXTS,
    SG_IMAGE,
    ImagePathIterator,
    _reset_matplotlib,
    figure_rst,
    matplotlib_scraper,
    save_figures,
)


@pytest.fixture(scope="function")
def make_gallery_conf(tmpdir):
    """Sets up a test sphinx-gallery configuration."""
    # Skip if numpy not installed
    pytest.importorskip("numpy")

    def make_gallery_conf(init=None):
        gallery_conf = _fill_gallery_conf_defaults(init or {})
        gallery_conf.update(
            src_dir=str(tmpdir), examples_dir=str(tmpdir), gallery_dir=str(tmpdir)
        )

        return gallery_conf

    return make_gallery_conf


class matplotlib_svg_scraper:
    """Test matplotlib svg scraper."""

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        """Call matplotlib scraper with 'svg' format."""
        return matplotlib_scraper(*args, format="svg", **kwargs)


@pytest.mark.parametrize("ext", ("png", "svg"))
def test_save_matplotlib_figures(make_gallery_conf, ext):
    """Test matplotlib figure save."""
    gallery_conf = make_gallery_conf(
        {"image_scrapers": (matplotlib_svg_scraper(),)} if ext == "svg" else {}
    )
    import matplotlib.pyplot as plt  # nest these so that Agg can be set

    plt.plot(1, 1)
    fname_template = os.path.join(gallery_conf["gallery_dir"], "image{0}.png")
    image_path_iterator = ImagePathIterator(fname_template)
    block = ("",) * 3
    block_vars = dict(image_path_iterator=image_path_iterator)
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(image_path_iterator) == 1
    fname = f"/image1.{ext}"
    assert fname in image_rst
    fname = gallery_conf["gallery_dir"] + fname
    assert os.path.isfile(fname)

    def _create_two_images():
        image_path_iterator.next()
        image_path_iterator.next()
        plt.plot(1, 1)
        plt.figure()
        plt.plot(1, 1)

    # Test capturing 2 images with shifted start number
    _create_two_images()
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(image_path_iterator) == 5
    for ii in range(4, 6):
        fname = f"/image{ii}.{ext}"
        assert fname in image_rst
        fname = gallery_conf["gallery_dir"] + fname
        assert os.path.isfile(fname)

    # Test `sphinx_gallery_multi_image(_block)` variables work; these variables prevent
    # images with `sphx-glr-single-img` classes from being converted to
    # `sphx-glr-multi-img` classes; requires > 1 image

    # Test file-wide `sphinx_gallery_multi_image` variable
    _create_two_images()
    block_vars["file_conf"] = {"multi_image": "single"}
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert "sphx-glr-single-img" in image_rst
    assert "sphx-glr-multi-img" not in image_rst

    # Test block-specific `sphinx_gallery_multi_image_block` variable
    # (test with default `sphinx_gallery_multi_image`, i.e. != "single")
    _create_two_images()
    block_vars["file_conf"] = {}
    block_vars["multi_image"] = "single"
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert "sphx-glr-single-img" in image_rst
    assert "sphx-glr-multi-img" not in image_rst
    # (test block-specific setting overrides file-wide setting)
    _create_two_images()
    block_vars["file_conf"] = {"multi_image": "single"}
    block_vars["multi_image"] = "multi"
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert "sphx-glr-single-img" not in image_rst
    assert "sphx-glr-multi-img" in image_rst


def test_image_srcset_config(make_gallery_conf):
    with pytest.raises(ConfigError, match="'image_srcset' config allowed"):
        make_gallery_conf({"image_srcset": "2x"})
    with pytest.raises(ConfigError, match="Invalid value for image_srcset parameter"):
        make_gallery_conf({"image_srcset": [False]})
    with pytest.raises(ConfigError, match="Invalid value for image_srcset parameter"):
        make_gallery_conf({"image_srcset": ["200"]})

    conf = make_gallery_conf({"image_srcset": ["2x"]})
    assert conf["image_srcset"] == [2.0]
    conf = make_gallery_conf({"image_srcset": ["1x", "2x"]})
    assert conf["image_srcset"] == [2.0]  # "1x" is implied.


def test_save_matplotlib_figures_hidpi(make_gallery_conf):
    """Test matplotlib hidpi figure save."""
    gallery_conf = make_gallery_conf({"image_srcset": ["2x"]})
    ext = "png"

    import matplotlib.pyplot as plt  # nest these so that Agg can be set

    plt.plot(1, 1)
    fname_template = os.path.join(gallery_conf["gallery_dir"], "image{0}.png")
    image_path_iterator = ImagePathIterator(fname_template)
    block = ("",) * 3
    block_vars = dict(image_path_iterator=image_path_iterator)
    image_rst = save_figures(block, block_vars, gallery_conf)

    fname = f"/image1.{ext}"
    assert fname in image_rst
    assert f"/image1_2_00x.{ext} 2.00x" in image_rst

    assert len(image_path_iterator) == 1
    fname = gallery_conf["gallery_dir"] + fname
    fnamehi = gallery_conf["gallery_dir"] + f"/image1_2_00x.{ext}"

    assert os.path.isfile(fname)
    assert os.path.isfile(fnamehi)

    # Test capturing 2 images with shifted start number
    image_path_iterator.next()
    image_path_iterator.next()
    plt.plot(1, 1)
    plt.figure()
    plt.plot(1, 1)
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(image_path_iterator) == 5
    for ii in range(4, 6):
        fname = f"/image{ii}.{ext}"
        assert fname in image_rst

        fname = gallery_conf["gallery_dir"] + fname
        assert os.path.isfile(fname)
        fname = f"/image{ii}_2_00x.{ext}"
        assert fname in image_rst
        fname = gallery_conf["gallery_dir"] + fname
        assert os.path.isfile(fname)


def _custom_func(x, y, z):
    return y["image_path_iterator"].next()


def test_custom_scraper(make_gallery_conf, monkeypatch):
    """Test custom scrapers."""
    # Test the API contract for custom scrapers
    with monkeypatch.context() as m:
        m.setattr(
            sphinx_gallery, "_get_sg_image_scraper", lambda: _custom_func, raising=False
        )
        for cust in (_custom_func, "sphinx_gallery"):
            # smoke test that it works
            make_gallery_conf({"image_scrapers": [cust]})
    # degenerate
    # without the monkey patch to add sphinx_gallery._get_sg_image_scraper,
    # we should get an error
    with pytest.raises(ConfigError, match="Unknown string option"):
        make_gallery_conf({"image_scrapers": ["sphinx_gallery"]})

    # other degenerate conditions
    with pytest.raises(ConfigError, match="Unknown string option for image_scraper"):
        make_gallery_conf({"image_scrapers": ["foo"]})
    for cust, msg in [
        (_custom_func, "did not produce expected image"),
        (lambda x, y, z: 1.0, "was not a string"),
    ]:
        conf = make_gallery_conf({"image_scrapers": [cust]})
        fname_template = os.path.join(conf["gallery_dir"], "image{0}.png")
        image_path_iterator = ImagePathIterator(fname_template)
        block = ("",) * 3
        block_vars = dict(image_path_iterator=image_path_iterator)
        with pytest.raises(ExtensionError, match=msg):
            save_figures(block, block_vars, conf)
    # degenerate string interface
    with monkeypatch.context() as m:
        m.setattr(sphinx_gallery, "_get_sg_image_scraper", "foo", raising=False)
        with pytest.raises(ConfigError, match="^Unknown string option for image_"):
            make_gallery_conf({"image_scrapers": ["sphinx_gallery"]})
    with monkeypatch.context() as m:
        m.setattr(sphinx_gallery, "_get_sg_image_scraper", lambda: "foo", raising=False)
        with pytest.raises(ConfigError, match="craper.*must be callable"):
            make_gallery_conf({"image_scrapers": ["sphinx_gallery"]})


@pytest.mark.parametrize("ext", _KNOWN_IMG_EXTS)
def test_figure_rst(ext):
    """Test reST generation of images."""
    figure_list = ["sphx_glr_plot_1." + ext]
    image_rst = figure_rst(figure_list, ".")
    single_image = f"""
.. image-sg:: /sphx_glr_plot_1.{ext}
   :alt: pl
   :srcset: /sphx_glr_plot_1.{ext}
   :class: sphx-glr-single-img
"""
    assert image_rst == single_image

    image_rst = figure_rst(figure_list + ["second." + ext], ".")

    image_list_rst = f"""
.. rst-class:: sphx-glr-horizontal


    *

      .. image-sg:: /sphx_glr_plot_1.{ext}
          :alt: pl
          :srcset: /sphx_glr_plot_1.{ext}
          :class: sphx-glr-multi-img

    *

      .. image-sg:: /second.{ext}
          :alt: pl
          :srcset: /second.{ext}
          :class: sphx-glr-multi-img
"""
    assert image_rst == image_list_rst


def test_figure_rst_path():
    """Test figure path correct in figure reSt."""
    # Tests issue #229
    local_img = [os.path.join(os.getcwd(), "third.png")]
    image_rst = figure_rst(local_img, ".")

    single_image = SG_IMAGE % ("third.png", "", "/third.png")
    assert image_rst == single_image


def test_figure_rst_srcset():
    """Test reST generation of images with srcset paths correct."""
    figure_list = ["sphx_glr_plot_1.png"]
    hipaths = [{0: "sphx_glr_plot_1.png", 2.0: "sphx_glr_plot_1_2_00.png"}]
    image_rst = figure_rst(figure_list, ".", srcsetpaths=hipaths)
    single_image = """
.. image-sg:: /sphx_glr_plot_1.png
   :alt: pl
   :srcset: /sphx_glr_plot_1.png, /sphx_glr_plot_1_2_00.png 2.00x
   :class: sphx-glr-single-img
"""
    assert image_rst == single_image

    hipaths += [{0: "second.png", 2.0: "second_2_00.png"}]
    image_rst = figure_rst(figure_list + ["second.png"], ".", srcsetpaths=hipaths)
    image_list_rst = """
.. rst-class:: sphx-glr-horizontal


    *

      .. image-sg:: /sphx_glr_plot_1.png
          :alt: pl
          :srcset: /sphx_glr_plot_1.png, /sphx_glr_plot_1_2_00.png 2.00x
          :class: sphx-glr-multi-img

    *

      .. image-sg:: /second.png
          :alt: pl
          :srcset: /second.png, /second_2_00.png 2.00x
          :class: sphx-glr-multi-img
"""
    assert image_rst == image_list_rst

    # test issue #229
    local_img = [os.path.join(os.getcwd(), "third.png")]
    image_rst = figure_rst(local_img, ".")

    single_image = SG_IMAGE % ("third.png", "", "/third.png")
    assert image_rst == single_image


def test_iterator():
    """Test ImagePathIterator."""
    ipi = ImagePathIterator("foo{0}")
    ipi._stop = 10
    with pytest.raises(ExtensionError, match="10 images"):
        for ii in ipi:
            pass


def test_reset_matplotlib(make_gallery_conf):
    """Test _reset_matplotlib."""
    import matplotlib

    matplotlib.rcParams["lines.linewidth"] = 42
    matplotlib.units.registry.clear()

    gallery_conf = make_gallery_conf()
    _reset_matplotlib(gallery_conf, "")

    assert matplotlib.rcParams["lines.linewidth"] != 42
    assert len(matplotlib.units.registry) > 0

import os

import pytest

from sphinx.errors import ConfigError, ExtensionError
import sphinx_gallery
from sphinx_gallery.gen_gallery import _complete_gallery_conf
from sphinx_gallery.scrapers import (figure_rst, mayavi_scraper, SG_IMAGE,
                                     matplotlib_scraper, ImagePathIterator,
                                     save_figures, _KNOWN_IMG_EXTS)
from sphinx_gallery.utils import _get_image


@pytest.fixture(scope='function')
def gallery_conf(tmpdir):
    """Sets up a test sphinx-gallery configuration"""
    # Skip if numpy not installed
    pytest.importorskip("numpy")

    gallery_conf = _complete_gallery_conf({}, str(tmpdir), True, False)
    gallery_conf.update(examples_dir=str(tmpdir), gallery_dir=str(tmpdir))
    return gallery_conf


class matplotlib_svg_scraper():

    def __repr__(self):
        return self.__class__.__name__

    def __call__(self, *args, **kwargs):
        return matplotlib_scraper(*args, format='svg', **kwargs)


@pytest.mark.parametrize('ext', ('png', 'svg'))
def test_save_matplotlib_figures(gallery_conf, ext):
    """Test matplotlib figure save."""
    if ext == 'svg':
        gallery_conf['image_scrapers'] = (matplotlib_svg_scraper(),)
    import matplotlib.pyplot as plt  # nest these so that Agg can be set
    plt.plot(1, 1)
    fname_template = os.path.join(gallery_conf['gallery_dir'], 'image{0}.png')
    image_path_iterator = ImagePathIterator(fname_template)
    block = ('',) * 3
    block_vars = dict(image_path_iterator=image_path_iterator)
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(image_path_iterator) == 1
    fname = '/image1.{0}'.format(ext)
    assert fname in image_rst
    fname = gallery_conf['gallery_dir'] + fname
    assert os.path.isfile(fname)

    # Test capturing 2 images with shifted start number
    image_path_iterator.next()
    image_path_iterator.next()
    plt.plot(1, 1)
    plt.figure()
    plt.plot(1, 1)
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(image_path_iterator) == 5
    for ii in range(4, 6):
        fname = '/image{0}.{1}'.format(ii, ext)
        assert fname in image_rst
        fname = gallery_conf['gallery_dir'] + fname
        assert os.path.isfile(fname)


def test_save_matplotlib_figures_hidpi(gallery_conf):
    """Test matplotlib hidpi figure save."""
    ext = 'png'
    gallery_conf['image_srcset'] = ["2x"]

    import matplotlib.pyplot as plt  # nest these so that Agg can be set
    plt.plot(1, 1)
    fname_template = os.path.join(gallery_conf['gallery_dir'], 'image{0}.png')
    image_path_iterator = ImagePathIterator(fname_template)
    block = ('',) * 3
    block_vars = dict(image_path_iterator=image_path_iterator)
    image_rst = save_figures(block, block_vars, gallery_conf)

    fname = f'/image1.{ext}'
    assert fname in image_rst
    assert f'/image1_2_0x.{ext} 2.0x' in image_rst

    assert len(image_path_iterator) == 1
    fname = gallery_conf['gallery_dir'] + fname
    fnamehi = gallery_conf['gallery_dir'] + f'/image1_2_0x.{ext}'

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
        fname = f'/image{ii}.{ext}'
        assert fname in image_rst

        fname = gallery_conf['gallery_dir'] + fname
        assert os.path.isfile(fname)
        fname = f'/image{ii}_2_0x.{ext}'
        assert fname in image_rst
        fname = gallery_conf['gallery_dir'] + fname
        assert os.path.isfile(fname)


def test_save_mayavi_figures(gallery_conf, req_mpl, req_pil):
    """Test file naming when saving figures. Requires mayavi."""
    import numpy as np
    Image = _get_image()
    try:
        from mayavi import mlab
    except ImportError:
        raise pytest.skip('Mayavi not installed')
    import matplotlib.pyplot as plt
    mlab.options.offscreen = True

    gallery_conf.update(
        image_scrapers=(matplotlib_scraper, mayavi_scraper))
    fname_template = os.path.join(gallery_conf['gallery_dir'], 'image{0}.png')
    image_path_iterator = ImagePathIterator(fname_template)
    block = ('',) * 3
    block_vars = dict(image_path_iterator=image_path_iterator)

    plt.axes([-0.1, -0.1, 1.2, 1.2])
    plt.pcolor([[0]], cmap='Greens')
    mlab.test_plot3d()
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(plt.get_fignums()) == 0
    assert len(image_path_iterator) == 2
    assert '/image0.png' not in image_rst
    assert '/image1.png' in image_rst
    assert '/image2.png' in image_rst
    assert '/image3.png' not in image_rst
    assert not os.path.isfile(fname_template.format(0))
    assert os.path.isfile(fname_template.format(1))
    assert os.path.isfile(fname_template.format(2))
    assert not os.path.isfile(fname_template.format(0))
    with Image.open(fname_template.format(1)) as img:
        pixels = np.asarray(img.convert("RGB"))
    assert (pixels == [247, 252, 245]).all()  # plt first

    # Test next-value handling, plus image_scrapers modification
    gallery_conf.update(image_scrapers=(matplotlib_scraper,))
    mlab.test_plot3d()
    plt.axes([-0.1, -0.1, 1.2, 1.2])
    plt.pcolor([[0]], cmap='Reds')
    image_rst = save_figures(block, block_vars, gallery_conf)
    assert len(plt.get_fignums()) == 0
    assert len(image_path_iterator) == 3
    assert '/image1.png' not in image_rst
    assert '/image2.png' not in image_rst
    assert '/image3.png' in image_rst
    assert '/image4.png' not in image_rst
    assert not os.path.isfile(fname_template.format(0))
    for ii in range(3):
        assert os.path.isfile(fname_template.format(ii + 1))
    assert not os.path.isfile(fname_template.format(4))
    with Image.open(fname_template.format(3)) as img:
        pixels = np.asarray(img.convert("RGB"))
    assert (pixels == [255, 245, 240]).all()


def _custom_func(x, y, z):
    return y['image_path_iterator'].next()


def test_custom_scraper(gallery_conf, monkeypatch):
    """Test custom scrapers."""
    # Test the API contract for custom scrapers
    complete_args = (gallery_conf, gallery_conf['gallery_dir'], True, False)
    with monkeypatch.context() as m:
        m.setattr(sphinx_gallery, '_get_sg_image_scraper',
                  lambda: _custom_func, raising=False)
        for cust in (_custom_func, 'sphinx_gallery'):
            gallery_conf.update(image_scrapers=[cust])
            # smoke test that it works
            _complete_gallery_conf(*complete_args, check_keys=False)
    # degenerate
    # without the monkey patch to add sphinx_gallery._get_sg_image_scraper,
    # we should get an error
    gallery_conf.update(image_scrapers=['sphinx_gallery'])
    with pytest.raises(ConfigError,
                       match="has no attribute '_get_sg_image_scraper'"):
        _complete_gallery_conf(*complete_args, check_keys=False)

    # other degenerate conditions
    gallery_conf.update(image_scrapers=['foo'])
    with pytest.raises(ConfigError, match='Unknown image scraper'):
        _complete_gallery_conf(*complete_args, check_keys=False)
    gallery_conf.update(image_scrapers=[_custom_func])
    fname_template = os.path.join(gallery_conf['gallery_dir'],
                                  'image{0}.png')
    image_path_iterator = ImagePathIterator(fname_template)
    block = ('',) * 3
    block_vars = dict(image_path_iterator=image_path_iterator)
    with pytest.raises(ExtensionError, match='did not produce expected image'):
        save_figures(block, block_vars, gallery_conf)
    gallery_conf.update(image_scrapers=[lambda x, y, z: 1.])
    with pytest.raises(ExtensionError, match='was not a string'):
        save_figures(block, block_vars, gallery_conf)
    # degenerate string interface
    gallery_conf.update(image_scrapers=['sphinx_gallery'])
    with monkeypatch.context() as m:
        m.setattr(sphinx_gallery, '_get_sg_image_scraper', 'foo',
                  raising=False)
        with pytest.raises(ConfigError, match='^Unknown image.*\n.*callable'):
            _complete_gallery_conf(*complete_args, check_keys=False)
    with monkeypatch.context() as m:
        m.setattr(sphinx_gallery, '_get_sg_image_scraper', lambda: 'foo',
                  raising=False)
        with pytest.raises(ConfigError, match='^Scraper.*was not callable'):
            _complete_gallery_conf(*complete_args, check_keys=False)


@pytest.mark.parametrize('ext', _KNOWN_IMG_EXTS)
def test_figure_rst(ext):
    """Testing rst of images"""
    figure_list = ['sphx_glr_plot_1.' + ext]
    image_rst = figure_rst(figure_list, '.')
    single_image = f"""
.. image-sg:: /sphx_glr_plot_1.{ext}
   :alt: pl
   :srcset: /sphx_glr_plot_1.{ext}
   :class: sphx-glr-single-img
"""
    assert image_rst == single_image

    image_rst = figure_rst(figure_list + ['second.' + ext], '.')

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

    # test issue #229
    local_img = [os.path.join(os.getcwd(), 'third.' + ext)]
    image_rst = figure_rst(local_img, '.')

    single_image = SG_IMAGE % ("third." + ext, '', "/third." + ext)
    assert image_rst == single_image


@pytest.mark.parametrize('ext', ['png'])
def test_figure_rst_srcset(ext):
    """Testing rst of images"""
    figure_list = ['sphx_glr_plot_1.' + ext]
    hipaths = [{0: 'sphx_glr_plot_1.png', 2.0: 'sphx_glr_plot_1_2_0.png'}]
    image_rst = figure_rst(figure_list, '.', srcsetpaths=hipaths)
    single_image = f"""
.. image-sg:: /sphx_glr_plot_1.{ext}
   :alt: pl
   :srcset: /sphx_glr_plot_1.{ext}, /sphx_glr_plot_1_2_0.{ext} 2.0x
   :class: sphx-glr-single-img
"""
    assert image_rst == single_image

    hipaths += [{0: 'second.png', 2.0: 'second_2_0.png'}]
    image_rst = figure_rst(figure_list + ['second.' + ext], '.',
                           srcsetpaths=hipaths+[])

    image_list_rst = f"""
.. rst-class:: sphx-glr-horizontal


    *

      .. image-sg:: /sphx_glr_plot_1.{ext}
          :alt: pl
          :srcset: /sphx_glr_plot_1.png, /sphx_glr_plot_1_2_0.png 2.0x
          :class: sphx-glr-multi-img

    *

      .. image-sg:: /second.{ext}
          :alt: pl
          :srcset: /second.{ext}, /second_2_0.{ext} 2.0x
          :class: sphx-glr-multi-img
"""
    assert image_rst == image_list_rst

    # test issue #229
    local_img = [os.path.join(os.getcwd(), 'third.' + ext)]
    image_rst = figure_rst(local_img, '.')

    single_image = SG_IMAGE % ("third." + ext, '', "/third." + ext)
    assert image_rst == single_image


def test_iterator():
    """Test ImagePathIterator."""
    ipi = ImagePathIterator('foo{0}')
    ipi._stop = 10
    with pytest.raises(ExtensionError, match='10 images'):
        for ii in ipi:
            pass

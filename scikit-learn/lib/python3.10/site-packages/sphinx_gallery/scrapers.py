# Author: Óscar Nájera
# License: 3-clause BSD
"""Scrapers for embedding images.

Collect images that have been produced by code blocks.

The only scraper we natively support is Matplotlib, others should
live in modules that will support them (e.g., PyVista, Plotly).  Scraped
images are injected as rst ``image-sg`` directives into the ``.rst``
file generated for each example script.
"""

import importlib
import inspect
import os
import re
import sys
from pathlib import PurePosixPath
from textwrap import indent
from warnings import filterwarnings

from sphinx.errors import ExtensionError

from .utils import optipng

__all__ = [
    "save_figures",
    "figure_rst",
    "ImagePathIterator",
    "clean_modules",
    "matplotlib_scraper",
]


###############################################################################
# Scrapers


def _import_matplotlib():
    """Import matplotlib safely."""
    # make sure that the Agg backend is set before importing any
    # matplotlib
    import matplotlib

    matplotlib.use("agg")
    matplotlib_backend = matplotlib.get_backend().lower()

    filterwarnings(
        "ignore",
        category=UserWarning,
        message="Matplotlib is currently using agg, which is a"
        " non-GUI backend, so cannot show the figure."
        "|(\n|.)*is non-interactive, and thus cannot be"
        " shown",
    )

    if matplotlib_backend != "agg":
        raise ExtensionError(
            "Sphinx-Gallery relies on the matplotlib 'agg' backend to "
            "render figures and write them to files. You are "
            "currently using the {} backend. Sphinx-Gallery will "
            "terminate the build now, because changing backends is "
            "not well supported by matplotlib. We advise you to move "
            "sphinx_gallery imports before any matplotlib-dependent "
            "import. Moving sphinx_gallery imports at the top of "
            "your conf.py file should fix this issue".format(matplotlib_backend)
        )

    import matplotlib.pyplot as plt

    return matplotlib, plt


def _matplotlib_fig_titles(fig):
    titles = []
    # get supertitle if exists
    suptitle = getattr(fig, "_suptitle", None)
    if suptitle is not None:
        titles.append(suptitle.get_text())
    # get titles from all axes, for all locs
    title_locs = ["left", "center", "right"]
    for ax in fig.axes:
        for loc in title_locs:
            text = ax.get_title(loc=loc)
            if text:
                titles.append(text)
    fig_titles = ", ".join(titles)
    return fig_titles


_ANIMATION_RST = """
.. container:: sphx-glr-animation

    .. raw:: html

        {0}
"""
_ANIMATION_VIDEO_RST = """
.. video:: {video}
   :class: sphx-glr-single-img
   :height: {height}
   :width: {width}
{options}
"""


def matplotlib_scraper(block, block_vars, gallery_conf, **kwargs):
    """Scrape Matplotlib images.

    Parameters
    ----------
    block : sphinx_gallery.py_source_parser.Block
        The code block to be executed. Format (label, content, lineno).
    block_vars : dict
        Dict of block variables.
    gallery_conf : dict
        Contains the configuration of Sphinx-Gallery
    **kwargs : dict
        Additional keyword arguments to pass to
        :meth:`~matplotlib.figure.Figure.savefig`, e.g. ``format='svg'``. The
        ``format`` keyword argument in particular is used to set the file
        extension of the output file (currently only 'png', 'jpg', 'svg',
        'gif', and 'webp' are supported).

        This is not used internally, but intended for use when overriding the scraper.

    Returns
    -------
    rst : str
        The reStructuredText that will be rendered to HTML containing
        the images. This is often produced by :func:`figure_rst`.
    """
    # Do not use _import_matplotlib() to avoid potentially changing the backend
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.animation import Animation

    image_path_iterator = block_vars["image_path_iterator"]
    image_rsts = []
    srcset = gallery_conf["image_srcset"]

    # Check for animations
    anims = {}
    if gallery_conf["matplotlib_animations"][0]:
        for ani in block_vars["example_globals"].values():
            if isinstance(ani, Animation):
                anims[ani._fig] = ani
    # Then standard images
    for fig_num, image_path in zip(plt.get_fignums(), image_path_iterator):
        image_path = PurePosixPath(image_path)
        if "format" in kwargs:
            image_path = image_path.with_suffix("." + kwargs["format"])
        # Convert figure number to Figure.
        fig = plt.figure(fig_num)
        # Deal with animations
        if anim := anims.get(fig):
            image_rsts.append(_anim_rst(anim, image_path, gallery_conf))
            continue
        # get fig titles
        fig_titles = _matplotlib_fig_titles(fig)
        to_rgba = matplotlib.colors.colorConverter.to_rgba
        # shallow copy should be fine here, just want to avoid changing
        # "kwargs" for subsequent figures processed by the loop
        these_kwargs = kwargs.copy()
        for attr in ["facecolor", "edgecolor"]:
            fig_attr = getattr(fig, "get_" + attr)()
            default_attr = matplotlib.rcParams["figure." + attr]
            if to_rgba(fig_attr) != to_rgba(default_attr) and attr not in kwargs:
                these_kwargs[attr] = fig_attr

        # save the figures, and populate the srcsetpaths
        try:
            fig.savefig(image_path, **these_kwargs)
            dpi0 = matplotlib.rcParams["savefig.dpi"]
            if dpi0 == "figure":
                dpi0 = fig.dpi
            dpi0 = these_kwargs.get("dpi", dpi0)
            srcsetpaths = {0: image_path}

            # save other srcset paths, keyed by multiplication factor:
            for mult in srcset:
                multst = f"{mult:.2f}".replace(".", "_")
                name = f"{image_path.stem}_{multst}x{image_path.suffix}"
                hipath = image_path.parent / PurePosixPath(name)
                hikwargs = {**these_kwargs, "dpi": mult * dpi0}
                fig.savefig(hipath, **hikwargs)
                srcsetpaths[mult] = hipath
            srcsetpaths = [srcsetpaths]
        except Exception:
            plt.close("all")
            raise

        if "images" in gallery_conf["compress_images"]:
            optipng(str(image_path), gallery_conf["compress_images_args"])
            for hipath in srcsetpaths[0].items():
                optipng(str(hipath), gallery_conf["compress_images_args"])

        image_rsts.append(
            figure_rst(
                [image_path],
                gallery_conf["src_dir"],
                fig_titles,
                srcsetpaths=srcsetpaths,
            )
        )

    # Determine whether single-img should be converted to multi-img
    convert_to_multi_image = True  # default is to convert
    if block_vars.get("multi_image") is not None:  # block setting takes precedence
        convert_to_multi_image = block_vars["multi_image"] != "single"
    elif block_vars.get("file_conf") is not None:  # then file setting
        convert_to_multi_image = block_vars["file_conf"].get("multi_image") != "single"

    plt.close("all")
    rst = ""
    if len(image_rsts) == 1:
        rst = image_rsts[0]
    elif len(image_rsts) > 1:
        if convert_to_multi_image:
            image_rsts = [
                re.sub(
                    r":class: sphx-glr-single-img", ":class: sphx-glr-multi-img", image
                )
                for image in image_rsts
            ]
        image_rsts = [
            HLIST_IMAGE_MATPLOTLIB + indent(image, " " * 6) for image in image_rsts
        ]
        rst = HLIST_HEADER + "".join(image_rsts)
    return rst


def _anim_rst(anim, image_path, gallery_conf):
    from matplotlib import rcParams
    from matplotlib.animation import FFMpegWriter, ImageMagickWriter

    # output the thumbnail as the image, as it will just be copied
    # if it's the file thumbnail
    fig = anim._fig
    image_path = image_path.with_suffix(".gif")
    fig_size = fig.get_size_inches()
    thumb_size = gallery_conf["thumbnail_size"]
    use_dpi = round(min(t_s / f_s for t_s, f_s in zip(thumb_size, fig_size)))
    if FFMpegWriter.isAvailable():
        writer = "ffmpeg"
    elif ImageMagickWriter.isAvailable():
        writer = "imagemagick"
    else:
        writer = None
    anim.save(str(image_path), writer=writer, dpi=use_dpi)

    _, fmt = gallery_conf["matplotlib_animations"]
    # Formats that are embedded in rst
    html = None
    if fmt is None:
        html = anim._repr_html_()
        if html is None:  # plt.rcParams['animation.html'] == 'none'
            html = anim.to_jshtml()
    elif fmt == "html5":
        html = anim.to_html5_video()
    elif fmt == "jshtml":
        html = anim.to_jshtml()
    if html is not None:
        html = indent(html, "     ")
        return _ANIMATION_RST.format(html)

    # Formats that are saved and use `video` directive
    video = image_path.with_suffix(f".{fmt}")
    anim.save(video)
    options = ["autoplay"]
    if getattr(anim, "_repeat", False):
        options.append("loop")
    dpi = rcParams["savefig.dpi"]
    if dpi == "figure":
        dpi = fig.dpi
    # relative_to doesn't work on windows
    # video_uri = video.relative_to(gallery_conf["src_dir"]).as_posix()
    video_uri = PurePosixPath(os.path.relpath(video, gallery_conf["src_dir"]))
    html = _ANIMATION_VIDEO_RST.format(
        video=f"/{video_uri}",
        width=int(fig_size[0] * dpi),
        height=int(fig_size[1] * dpi),
        options="".join(f"   :{opt}:\n" for opt in options),
    )
    return html


_scraper_dict = dict(
    matplotlib=matplotlib_scraper,
)


class ImagePathIterator:
    """Iterate over image paths for a given example.

    Parameters
    ----------
    image_path : str
        The template image path.
    """

    def __init__(self, image_path):
        self.image_path = image_path
        self.paths = list()
        self._stop = 1000000

    def __len__(self):
        """Return the number of image paths used.

        Returns
        -------
        n_paths : int
            The number of paths.
        """
        return len(self.paths)

    def __iter__(self):
        """Iterate over paths.

        Returns
        -------
        paths : iterable of str

        This enables the use of this Python pattern::

            >>> for epoch in epochs:  # doctest: +SKIP
            >>>     print(epoch)  # doctest: +SKIP

        Where ``epoch`` is given by successive outputs of
        :func:`mne.Epochs.next`.
        """
        # we should really never have 1e6, let's prevent some user pain
        for ii in range(self._stop):
            yield self.next()
        else:
            raise ExtensionError(f"Generated over {self._stop} images")

    def next(self):
        """Return the next image path, with numbering starting at 1."""
        return self.__next__()

    def __next__(self):
        # The +1 here is because we start image numbering at 1 in filenames
        path = self.image_path.format(len(self) + 1)
        self.paths.append(path)
        return path


# For now, these are what we support
# Update advanced.rst if this list is changed
_KNOWN_IMG_EXTS = ("png", "svg", "jpg", "gif", "webp")


def _find_image_ext(path):
    """Find an image, tolerant of different file extensions."""
    path = os.path.splitext(path)[0]
    for ext in _KNOWN_IMG_EXTS:
        this_path = f"{path}.{ext}"
        if os.path.isfile(this_path):
            break
    else:
        ext = "png"
    return (f"{path}.{ext}", ext)


def save_figures(block, block_vars, gallery_conf):
    """Save all open figures of the example code-block.

    Parameters
    ----------
    block : sphinx_gallery.py_source_parser.Block
        The code block to be executed. Format (label, content, lineno).
    block_vars : dict
        Dict of block variables.
    gallery_conf : dict
        Contains the configuration of Sphinx-Gallery

    Returns
    -------
    images_rst : str
        rst code to embed the images in the document.
    """
    from .gen_rst import _get_callables

    image_path_iterator = block_vars["image_path_iterator"]
    all_rst = ""
    prev_count = len(image_path_iterator)
    for scraper in _get_callables(gallery_conf, "image_scrapers"):
        rst = scraper(block, block_vars, gallery_conf)
        if not isinstance(rst, str):
            raise ExtensionError(
                f"rst from scraper {scraper!r} was not a "
                f"string, got type {type(rst)}:\n{rst!r}"
            )
        n_new = len(image_path_iterator) - prev_count
        for ii in range(n_new):
            current_path, _ = _find_image_ext(
                image_path_iterator.paths[prev_count + ii]
            )
            if not os.path.isfile(current_path):
                raise ExtensionError(
                    f"Scraper {scraper} did not produce expected image:\n{current_path}"
                )
        all_rst += rst
    return all_rst


def figure_rst(figure_list, sources_dir, fig_titles="", srcsetpaths=None):
    """Generate reST for a list of image filenames.

    Depending on whether we have one or more figures, we use a
    single rst call to 'image' or a horizontal list.

    Parameters
    ----------
    figure_list : list
        List of strings of the figures' absolute paths.
    sources_dir : str
        absolute path of Sphinx documentation sources
    fig_titles : str
        Titles of figures, empty string if no titles found. Currently
        only supported for matplotlib figures, default = ''.
    srcsetpaths : list or None
        List of dictionaries containing absolute paths.  If
        empty, then srcset field is populated with the figure path.
        (see ``image_srcset`` configuration option).  Otherwise,
        each dict is of the form
        {0: /images/image.png, 2.0: /images/image_2_00x.png}
        where the key is the multiplication factor and the contents
        the path to the image created above.

    Returns
    -------
    images_rst : str
        rst code to embed the images in the document

    The rst code will have a custom ``image-sg`` directive that allows
    multiple resolution images to be served e.g.:
    ``:srcset: /plot_types/imgs/img_001.png,
      /plot_types/imgs/img_2_00x.png 2.00x``

    """
    if srcsetpaths is None:
        # this should never happen, but figure_rst is public, so
        # this has to be a kwarg...
        srcsetpaths = [{0: fl} for fl in figure_list]

    figure_paths = [
        os.path.relpath(figure_path, sources_dir).replace(os.sep, "/").lstrip("/")
        for figure_path in figure_list
    ]

    # Get alt text
    alt = ""
    if fig_titles:
        alt = fig_titles
    elif figure_list:
        file_name = os.path.split(figure_list[0])[1]
        # remove ext & 'sphx_glr_' from start & n#'s from end
        file_name_noext = os.path.splitext(file_name)[0][9:-4]
        # replace - & _ with \s
        file_name_final = re.sub(r"[-,_]", " ", file_name_noext)
        alt = file_name_final
    alt = _single_line_sanitize(alt)

    images_rst = ""
    if len(figure_paths) == 1:
        figure_name = figure_paths[0]
        hinames = srcsetpaths[0]
        srcset = _get_srcset_st(sources_dir, hinames)
        images_rst = SG_IMAGE % (figure_name, alt, srcset)

    elif len(figure_paths) > 1:
        images_rst = HLIST_HEADER
        for nn, figure_name in enumerate(figure_paths):
            hinames = srcsetpaths[nn]
            srcset = _get_srcset_st(sources_dir, hinames)

            images_rst += HLIST_SG_TEMPLATE % (figure_name, alt, srcset)

    return images_rst


def _get_srcset_st(sources_dir, hinames):
    """Create the srcset string for including on the rst line.

    For example; `sources_dir` might be `/home/sample-proj/source`,
    hinames posix paths:
    0: /home/sample-proj/source/plot_types/images/img1.png,
    2.0: /home/sample-proj/source/plot_types/images/img1_2_00x.png,

    The result would be:
    '/plot_types/basic/images/sphx_glr_pie_001.png,
    /plot_types/basic/images/sphx_glr_pie_001_2_00x.png 2.00x'
    """
    srcst = ""
    for k in hinames.keys():
        path = os.path.relpath(hinames[k], sources_dir).replace(os.sep, "/").lstrip("/")
        srcst += "/" + path
        if k == 0:
            srcst += ", "
        else:
            srcst += f" {k:1.2f}x, "
    if srcst[-2:] == ", ":
        srcst = srcst[:-2]
    srcst += ""

    return srcst


def _single_line_sanitize(s):
    """Remove problematic newlines."""
    # For example, when setting a :alt: for an image, it shouldn't have \n
    # This is a function in case we end up finding other things to replace
    return s.replace("\n", " ")


# The following strings are used when we have several pictures: we use
# an html div tag that our CSS uses to turn the lists into horizontal
# lists.
HLIST_HEADER = """
.. rst-class:: sphx-glr-horizontal

"""

HLIST_IMAGE_MATPLOTLIB = """
    *
"""

HLIST_SG_TEMPLATE = """
    *

      .. image-sg:: /%s
          :alt: %s
          :srcset: %s
          :class: sphx-glr-multi-img
"""

SG_IMAGE = """
.. image-sg:: /%s
   :alt: %s
   :srcset: %s
   :class: sphx-glr-single-img
"""

# keep around for back-compat:
SINGLE_IMAGE = """
 .. image:: /%s
     :alt: %s
     :class: sphx-glr-single-img
"""

###############################################################################
# Module resetting


def _reset_matplotlib(gallery_conf, fname):
    """Reset matplotlib."""
    mpl, plt = _import_matplotlib()
    plt.rcdefaults()
    importlib.reload(mpl.units)
    importlib.reload(mpl.dates)
    importlib.reload(mpl.category)


def _reset_seaborn(gallery_conf, fname):
    """Reset seaborn."""
    seaborn_module = sys.modules.get("seaborn")
    if seaborn_module is not None:
        seaborn_module.reset_defaults()


_reset_dict = {
    "matplotlib": _reset_matplotlib,
    "seaborn": _reset_seaborn,
}


def clean_modules(gallery_conf, fname, when):
    """Remove, unload, or reset modules.

    After a script is executed it can load a variety of settings that one
    does not want to influence in other examples in the gallery.

    Parameters
    ----------
    gallery_conf : dict
        The gallery configuration.
    fname : str or None
        The example being run. Will be None when this is called entering
        a directory of examples to be built.
    when : str
        Whether this module is run before or after examples.

        This parameter is only forwarded when the callables accept 3
        parameters.
    """
    from .gen_rst import _get_callables

    for reset_module in _get_callables(gallery_conf, "reset_modules"):
        sig = inspect.signature(reset_module)
        if len(sig.parameters) == 3:
            third_param = list(sig.parameters.keys())[2]
            if third_param != "when":
                raise ValueError(
                    f"3rd parameter in {reset_module.__name__} "
                    "function signature must be 'when', "
                    f"got {third_param}"
                )
            reset_module(gallery_conf, fname, when=when)
        else:
            reset_module(gallery_conf, fname)

# This module defines an image scraper for sphinx-gallery
# https://sphinx-gallery.github.io/
# which can be used by projects using plotly in their documentation.
from glob import glob
import os
import shutil

import plotly

plotly.io.renderers.default = "sphinx_gallery_png"


def plotly_sg_scraper(block, block_vars, gallery_conf, **kwargs):
    """Scrape Plotly figures for galleries of examples using
    sphinx-gallery.

    Examples should use ``plotly.io.show()`` to display the figure with
    the custom sphinx_gallery renderer.

    Since the sphinx_gallery renderer generates both html and static png
    files, we simply crawl these files and give them the appropriate path.

    Parameters
    ----------
    block : tuple
        A tuple containing the (label, content, line_number) of the block.
    block_vars : dict
        Dict of block variables.
    gallery_conf : dict
        Contains the configuration of Sphinx-Gallery
    **kwargs : dict
        Additional keyword arguments to pass to
        :meth:`~matplotlib.figure.Figure.savefig`, e.g. ``format='svg'``.
        The ``format`` kwarg in particular is used to set the file extension
        of the output file (currently only 'png' and 'svg' are supported).

    Returns
    -------
    rst : str
        The ReSTructuredText that will be rendered to HTML containing
        the images.

    Notes
    -----
    Add this function to the image scrapers
    """
    examples_dir = os.path.dirname(block_vars["src_file"])
    pngs = sorted(glob(os.path.join(examples_dir, "*.png")))
    htmls = sorted(glob(os.path.join(examples_dir, "*.html")))
    image_path_iterator = block_vars["image_path_iterator"]
    image_names = list()
    seen = set()
    for html, png in zip(htmls, pngs):
        if png not in seen:
            seen |= set(png)
            this_image_path_png = next(image_path_iterator)
            this_image_path_html = os.path.splitext(this_image_path_png)[0] + ".html"
            image_names.append(this_image_path_html)
            shutil.move(png, this_image_path_png)
            shutil.move(html, this_image_path_html)
    # Use the `figure_rst` helper function to generate rST for image files
    return figure_rst(image_names, gallery_conf["src_dir"])


def figure_rst(figure_list, sources_dir):
    """Generate RST for a list of PNG filenames.

    Depending on whether we have one or more figures, we use a
    single rst call to 'image' or a horizontal list.

    Parameters
    ----------
    figure_list : list
        List of strings of the figures' absolute paths.
    sources_dir : str
        absolute path of Sphinx documentation sources

    Returns
    -------
    images_rst : str
        rst code to embed the images in the document
    """

    figure_paths = [
        os.path.relpath(figure_path, sources_dir).replace(os.sep, "/").lstrip("/")
        for figure_path in figure_list
    ]
    images_rst = ""
    if not figure_paths:
        return images_rst
    figure_name = figure_paths[0]
    figure_path = os.path.join("images", os.path.basename(figure_name))
    images_rst = SINGLE_HTML % figure_path
    return images_rst


SINGLE_HTML = """
.. raw:: html
    :file: %s
"""

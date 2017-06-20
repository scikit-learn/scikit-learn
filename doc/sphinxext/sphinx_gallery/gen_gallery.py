# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
========================
Sphinx-Gallery Generator
========================

Attaches Sphinx-Gallery to Sphinx in order to generate the galleries
when building the documentation.
"""


from __future__ import division, print_function, absolute_import
import re
import os
from . import glr_path_static
from .gen_rst import generate_dir_rst
from .docs_resolv import embed_code_links


def clean_gallery_out(build_dir):
    """Deletes images under the sphx_glr namespace in the build directory"""
    # Sphinx hack: sphinx copies generated images to the build directory
    #  each time the docs are made.  If the desired image name already
    #  exists, it appends a digit to prevent overwrites.  The problem is,
    #  the directory is never cleared.  This means that each time you build
    #  the docs, the number of images in the directory grows.
    #
    # This question has been asked on the sphinx development list, but there
    #  was no response: http://osdir.com/ml/sphinx-dev/2011-02/msg00123.html
    #
    # The following is a hack that prevents this behavior by clearing the
    #  image build directory from gallery images each time the docs are built.
    #  If sphinx changes their layout between versions, this will not
    #  work (though it should probably not cause a crash).
    # Tested successfully on Sphinx 1.0.7

    build_image_dir = os.path.join(build_dir, '_images')
    if os.path.exists(build_image_dir):
        filelist = os.listdir(build_image_dir)
        for filename in filelist:
            if filename.startswith('sphx_glr') and filename.endswith('png'):
                os.remove(os.path.join(build_image_dir, filename))


def generate_gallery_rst(app):
    """Generate the Main examples gallery reStructuredText

    Start the sphinx-gallery configuration and recursively scan the examples
    directories in order to populate the examples gallery
    """
    try:
        plot_gallery = eval(app.builder.config.plot_gallery)
    except TypeError:
        plot_gallery = bool(app.builder.config.plot_gallery)

    gallery_conf.update(app.config.sphinx_gallery_conf)
    gallery_conf.update(plot_gallery=plot_gallery)
    gallery_conf.update(abort_on_example_error=app.builder.config.abort_on_example_error)

    # this assures I can call the config in other places
    app.config.sphinx_gallery_conf = gallery_conf
    app.config.html_static_path.append(glr_path_static())

    clean_gallery_out(app.builder.outdir)

    examples_dirs = gallery_conf['examples_dirs']
    gallery_dirs = gallery_conf['gallery_dirs']

    if not isinstance(examples_dirs, list):
        examples_dirs = [examples_dirs]
    if not isinstance(gallery_dirs, list):
        gallery_dirs = [gallery_dirs]

    mod_examples_dir = os.path.relpath(gallery_conf['mod_example_dir'],
                                       app.builder.srcdir)
    seen_backrefs = set()

    for examples_dir, gallery_dir in zip(examples_dirs, gallery_dirs):
        examples_dir = os.path.relpath(examples_dir,
                                       app.builder.srcdir)
        gallery_dir = os.path.relpath(gallery_dir,
                                      app.builder.srcdir)

        for workdir in [examples_dir, gallery_dir, mod_examples_dir]:
            if not os.path.exists(workdir):
                os.makedirs(workdir)

        # we create an index.rst with all examples
        fhindex = open(os.path.join(gallery_dir, 'index.rst'), 'w')
        # Here we don't use an os.walk, but we recurse only twice: flat is
        # better than nested.
        fhindex.write(generate_dir_rst(examples_dir, gallery_dir, gallery_conf,
                                       seen_backrefs))
        for directory in sorted(os.listdir(examples_dir)):
            if os.path.isdir(os.path.join(examples_dir, directory)):
                src_dir = os.path.join(examples_dir, directory)
                target_dir = os.path.join(gallery_dir, directory)
                fhindex.write(generate_dir_rst(src_dir, target_dir,
                                               gallery_conf,
                                               seen_backrefs))
        fhindex.flush()


def touch_empty_backreferences(app, what, name, obj, options, lines):
    """Generate empty back-reference example files

    This avoids inclusion errors/warnings if there are no gallery
    examples for a class / module that is being parsed by autodoc"""

    examples_path = os.path.join(app.srcdir,
                                 app.config.sphinx_gallery_conf["mod_example_dir"],
                                 "%s.examples" % name)

    if not os.path.exists(examples_path):
        # touch file
        open(examples_path, 'w').close()


gallery_conf = {
    'filename_pattern': re.escape(os.sep) + 'plot',
    'examples_dirs': '../examples',
    'gallery_dirs': 'auto_examples',
    'mod_example_dir': os.path.join('modules', 'generated'),
    'doc_module': (),
    'reference_url': {},
}


def setup(app):
    """Setup sphinx-gallery sphinx extension"""
    app.add_config_value('plot_gallery', True, 'html')
    app.add_config_value('abort_on_example_error', False, 'html')
    app.add_config_value('sphinx_gallery_conf', gallery_conf, 'html')
    app.add_stylesheet('gallery.css')

    extensions_attr = '_extensions' if hasattr(app, '_extensions') else 'extensions'
    if 'sphinx.ext.autodoc' in getattr(app, extensions_attr):
        app.connect('autodoc-process-docstring', touch_empty_backreferences)

    app.connect('builder-inited', generate_gallery_rst)

    app.connect('build-finished', embed_code_links)


def setup_module():
    # HACK: Stop nosetests running setup() above
    pass

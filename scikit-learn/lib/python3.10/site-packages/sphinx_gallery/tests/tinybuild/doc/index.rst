Tiny test build
===============

The complete Sphinx-Gallery project is automatically documented for
every module. Examples `here <auto_examples/index.html>`_.

.. currentmodule:: sphinx_gallery

.. automodule:: sphinx_gallery
   :no-members:
   :no-inherited-members:

:py:mod:`sphinx_gallery`:

.. autosummary::
   :toctree: gen_modules/
   :template: module.rst

   backreferences
   _dummy
   _dummy.nested
   docs_resolv
   downloads
   gen_gallery
   gen_rst
   notebook
   py_source_parser
   scrapers
   sorting

Examples
--------

This tests that mini-gallery reference labels work:
:ref:`sphx_glr_backref_sphinx_gallery.backreferences.identify_names`.

.. toctree::
   :maxdepth: 2

   auto_examples/index.rst
   auto_examples_with_rst/index.rst
   auto_examples_rst_index/index.rst
   auto_examples_README_header/index.rst

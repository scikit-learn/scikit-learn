=====================================
Example with rst files passed through
=====================================

Sometimes it is desirable to mix rst files with python files in the same example
directory.  Here, the source directory has a mix of python examples (``plot_boo.py``)
and raw rst files (``rst_example1.rst``).  The rst files are passed directly from
the source directory  ``/examples_with_rst/rst_example1.rst`` to the
directory created by ``sphinx-gallery`` in
``/doc/examples_with_rst/rst_example1.rst`` without change.

Any rst files that are passed through this way must have a manual
``toctree`` entry somewhere, or you will get a warning that the file doesn't exist in
the toctree.  We add that here to the ``README.txt`` as::

    Rst files
    =========

    .. toctree::

      rst_example1
      rst_example2

    Sphinx-gallery files
    ====================


Note that the python example also shows up as a usual thumbnail below this table of contents.

Rst files
=========

.. toctree::

    rst_example1
    rst_example2

Sphinx-gallery files
====================

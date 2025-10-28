==========================================
Example gallery that has its own index.rst
==========================================

Usually made by sphinx-gallery from the ``README.txt``.
However, if ``index.rst`` is given, we use that instead and
ignore README.txt.

.. toctree::

    plot_examp

Subtopic one
============
.. toctree::

    examp_subdir1/index

Subtopic two
============

This subtopic's directory has its own index.rst.

.. toctree::

    examp_subdir2/index

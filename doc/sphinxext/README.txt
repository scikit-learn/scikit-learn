=====================================
numpydoc -- Numpy's Sphinx extensions
=====================================

Numpy's documentation uses several custom extensions to Sphinx.  These
are shipped in this ``numpydoc`` package, in case you want to make use
of them in third-party projects.

The following extensions are available:

  - ``numpydoc``: support for the Numpy docstring format in Sphinx, and add
    the code description directives ``np-function``, ``np-cfunction``, etc.
    that support the Numpy docstring syntax.

  - ``numpydoc.traitsdoc``: For gathering documentation about Traits attributes.

  - ``numpydoc.plot_directives``: Adaptation of Matplotlib's ``plot::``
    directive. Note that this implementation may still undergo severe
    changes or eventually be deprecated.

  - ``numpydoc.only_directives``: (DEPRECATED)

  - ``numpydoc.autosummary``: (DEPRECATED) An ``autosummary::`` directive.
    Available in Sphinx 0.6.2 and (to-be) 1.0 as ``sphinx.ext.autosummary``,
    and it the Sphinx 1.0 version is recommended over that included in
    Numpydoc.


numpydoc
========

Numpydoc inserts a hook into Sphinx's autodoc that converts docstrings
following the Numpy/Scipy format to a form palatable to Sphinx.

Options
-------

The following options can be set in conf.py:

- numpydoc_use_plots: bool

  Whether to produce ``plot::`` directives for Examples sections that
  contain ``import matplotlib``.

- numpydoc_show_class_members: bool

  Whether to show all members of a class in the Methods and Attributes
  sections automatically.

- numpydoc_edit_link: bool  (DEPRECATED -- edit your HTML template instead)

  Whether to insert an edit link after docstrings.

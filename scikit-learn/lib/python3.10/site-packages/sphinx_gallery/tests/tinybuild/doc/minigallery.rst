Test mini-galleries
===================

Test 1-N

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder

Test 1-D-D

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder
    :add-heading:

Test 1-D-C

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder
    :add-heading:
    :heading-level: -

Test 1-C-D

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder
    :add-heading: This is a custom heading

Test 2-N

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder sphinx_gallery.sorting.FileNameSortKey

Test 2-D-D

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder sphinx_gallery.sorting.FileNameSortKey
    :add-heading:

Test 2-C-C

.. minigallery:: sphinx_gallery.sorting.ExplicitOrder sphinx_gallery.sorting.FileNameSortKey
    :add-heading: This is a different custom heading
    :heading-level: =

Test 1-F

.. minigallery:: ../examples/plot_log.py

Test 2-F-G

.. minigallery:: ../examples/plot_matplotlib*.py


Test 3-F-G-B

.. minigallery:: ../examples/plot_log.py ../examples/*matplotlib*.py sphinx_gallery.sorting.ExplicitOrder sphinx_gallery.sorting.FileNameSortKey
    :add-heading: All the input types


Test 1-F-R

.. minigallery:: ../examples_with_rst/*.py


Test 1-S

.. minigallery:: ../examples_rst_index/*/*.py

Test 3-N

.. minigallery::

   ../examples/plot_log.py
   ../examples/*matplotlib*.py
   sphinx_gallery.sorting.ExplicitOrder
   sphinx_gallery.sorting.FileNameSortKey

Test duplicates

.. minigallery::

    sphinx_gallery.sorting.ExplicitOrder
    ../examples/plot_second_future_imports.py
    ../examples/plot_second_future_import*
    sphinx_gallery.py_source_parser.Block
    ../examples/plot_numpy_matplotib.py

.. _astro_tutorial_setup:

===============================
Tutorial Setup and Installation
===============================

.. topic:: Objectives

   At the end of this section, you will
  
   1. Have scikit-learn and all the prerequisites and dependencies for
      this tutorial installed on your machine.
   2. Download the source files and data required for this tutorial

Python Prerequisites
--------------------

This tutorial is based on scikit-learn, which has the following dependencies:

- `numpy <http://numpy.scipy.org>`_ : this is a python module which has powerful
  tools for the creation and manipulation of arrays.  It is the foundation of
  most scientific computing packages in python

- `scipy <http://www.scipy.org>`_ : this is a python module which builds on
  numpy and provides fast implementations of many basic scientific algorithms.

- `matplotlib <http://matplotlib.sourceforge.net/>`_ : this is a powerful
  package for generating plots, figures, and diagrams.  Our main form of
  visual interaction with data and results depends on matplotlib.

We will also make extensive use of `iPython <http://ipython.org>`_, an
interactive python interpreter.  In particular, much of the interactive
material will make use of the 
`ipython notebook`_ functionality, which was introduced in ipython version 0.12.

Installing scikit-learn and Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please refer to the :ref:`install page <installing_scikit_learn>` for
per-system instructions on installing scikit-learn.  In addition to
``numpy``, ``scipy``, and ``scikit-learn``, this tutorial will assume that
you have ``matplotlib`` and ``ipython`` installed as well.

  * Under **Debian or Ubuntu Linux** you should use::

      % sudo apt-get install build-essential python-dev python-numpy \
        python-numpy-dev python-scipy libatlas-dev g++ python-matplotlib \
        ipython

  * Under **MacOSX** you should probably use a scientific python distribution
    such as `Scipy Superpack`_

  * Under **Windows** the `Python(x,y)`_ is probably your best bet to get a
    working numpy / scipy environment up and running.

  * Power-users may wish to install bleeding edge versions of these
    packages from the source.  The source can be downloaded using
    ``git`` from their respective `GitHub`_ repositories.

Alternatively under Windows and MaxOSX you can use the EPD_ (Enthought
Python Distribution) which is a (non-open source) packaging of the
scientific python stack.

.. note::
   that to use `ipython notebook`_, you must install ``ipython`` version
   0.12 and several other dependencies.  Refer to the ipython documentation
   for details.

.. _`Scipy Superpack`: http://fonnesbeck.github.com/ScipySuperpack/
.. _`Python(x,y)`: http://www.pythonxy.com/
.. _EPD: https://www.enthought.com/products/epd.php
.. _GitHub: http://www.github.com
.. _`ipython notebook`: http://ipython.org/ipython-doc/stable/interactive/htmlnotebook.html


Tutorial Files
--------------
The source code for the example files in the following pages is best
accessed through cloning the scikit-learn repository using
`git <http://git-scm.com/>`_.  Once ``git`` is installed the
command to accomplish this is::

    % git clone https://github.com/jakevdp/scikit-learn/
    % cd scikit-learn
    % git checkout -b sklearn-tutorial origin/sklearn-tutorial

The last line switches to the branch of the repository which contains
the source of this tutorial.  In what follows, the path to the
``scikit-learn`` source folder will be named ``$SKL_HOME``.
Within this directory can be found the tutorial source::

    % cd doc/tutorial/astronomy

In what follows, this directory will be named ``$TUTORIAL_HOME``. It
should contain the following folders:

  * ``data`` - folder to put the datasets used during the tutorial

  * ``skeletons`` - sample incomplete scripts for the exercices

  * ``solutions`` - solutions of the exercices

  * ``notebooks`` - ipython notebooks which provide an interactive interface
    to parts of this tutorial.

You can aleardy copy the skeletons into a new folder named ``workspace``
where you will edit your own files for the exercices while keeping
the original skeletons intact::

    % cp -r skeletons workspace


Download the datasets
---------------------

Machine Learning algorithms need data. Go to each ``$TUTORIAL_HOME/data/``
sub-folder and run the ``fetch_data.py`` script from there (after
having read them first).  This will download a dataset to the current
directory.  This tutorial has three such datasets; they will be used
in the examples and exercises later on.

For example::

    % cd $TUTORIAL_HOME/data/sdss_colors
    % less fetch_data.py
    % python fetch_data.py

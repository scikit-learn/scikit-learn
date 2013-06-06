Tutorial setup
==============

To get started with this tutorial, you firstly must have the
*scikit-learn* and all of its requiered dependencies installed.
The source of this tutorial can be found within your
scikit-learn folder::

    scikit-learn/doc/tutorial/text_analytics/

The tutorial folder, should contain the following folders:

  * ``*.rst files`` - the source of the tutorial document written with sphinx

  * ``data`` - folder to put the datasets used during the tutorial

  * ``skeletons`` - sample incomplete scripts for the exercices

  * ``solutions`` - solutions of the exercices


You can aleardy copy the skeletons into a new folder somewhere
on your hard-drive named ``sklearn_tut_workspace`` where you
will edit your own files for the exercices while keeping
the original skeletons intact::

    % cp -r skeletons work_directory/sklear_tut_workspace


Install scikit-learn build dependencies
---------------------------------------

Please refer to the `scikit-learn install`_ page for per-system
instructions.

.. _`scikit-learn install`: http://scikit-learn.sourceforge.net/install.html

You must have ``numpy``, ``scipy``, ``matplotlib`` and ``ipython``
installed:

  * Under **Debian or Ubuntu Linux** you should use::

      % sudo apt-get install build-essential python-dev python-numpy \
        python-numpy-dev python-scipy libatlas-dev g++ python-matplotlib \
        ipython

  * Under **MacOSX** you should probably use a scientific python distribution
    such as `Scipy Superpack`_

  * Under **Windows** the `Python(x,y)`_ is probably your best bet to get a
    working numpy / scipy environment up and running.

Alternatively under Windows and MaxOSX you can use the EPD_ (Enthought
Python Distribution) which is a (non-open source) packaging of the
scientific python stack.

.. _`Scipy Superpack`: http://stronginference.com/scipy-superpack/
.. _`Python(x,y)`: http://www.pythonxy.com/
.. _EPD: https://www.enthought.com/products/epd.php


Build scikit-learn from source
------------------------------

Here are the instructions to install the current master from source
on a POSIX system (e.g. Linux and MacOSX). **In the folder next to
$TUTORIAL_HOME** do::

    % git clone https://github.com/scikit-learn/scikit-learn.git
    % cd scikit-learn

You can then build it locally and add it to your PYTHONPATH environment
variable::

    % python setup.py build_ext -i
    % export PYTHONPATH=`pwd`

Alternatively you can install the library globally::

    % python setup.py build
    % sudo python setup.py install

You should also be able to launch the tests from anywhere in the system
(if nose is installed) with the following::

    % nosetests sklearn

The output should end with ``OK`` as in::

    ----------------------------------------------------------------------
    Ran 589 tests in 36.876s

    OK (SKIP=2)

If this is not the case please send a mail to the `scikit-learn mailing list`_
including the error messages along with the version number of all the afore
mentioned dependencies and your operating system.

.. _`scikit-learn mailing list`: https://lists.sourceforge.net/lists/listinfo/scikit-learn-general

In the rest of the tutorial, the path to the ``scikit-learn`` source
folder will be named ``$SKL_HOME``.

As usual building from source under Windows is slightly more complicated.
Checkout the `build instructions`_ on the scikit-learn website.

.. _`build instructions`: http://scikit-learn.sourceforge.net/dev/install.html#building-on-windows


Download the datasets
---------------------

Machine Learning algorithms need data. Go to each ``$TUTORIAL_HOME/data``
sub-folder and run the ``fetch_data.py`` script from there (after
having read them first).

For instance::

    % cd $TUTORIAL_HOME/data/languages
    % less fetch_data.py
    % python fetch_data.py


Tutorial setup
==============

The following assumes you have extracted the source distribution
of this tutorial somewhere on your local disk. Alternatively you
can use git to clone this repo directly from github onto your
local disk.

In the following we will name this folder ``$TUTORIAL_HOME``. It
should contain the following folders:


  * ``tutorial`` - the source of the tutorial document written with sphinx

  * ``data`` - folder to put the datasets used during the tutorial

  * ``skeletons`` - sample incomplete scripts for the exercices

  * ``solutions`` - solutions of the exercices


You can aleardy copy the skeletons into a new folder named ``workspace``
where you will edit your own files for the exercices while keeping
the original skeletons intact::

    % cp -r skeletons workspace


Install scikit-learn 0.7
------------------------

Please refer to the `scikit-learn install`_ page for per-system instructions.

.. _`scikit-learn install`: http://scikit-learn.sourceforge.net/install.html

You must have ``numpy``, ``scipy`` and ``matplotlib`` installed first.

Here are the instructions to install the 0.7 release from source
on a POSIX system (e.g. Linux and MacOSX). First download the release
archive and extract it **in the folder next to $TUTORIAL_HOME**::

    % wget http://pypi.python.org/packages/source/s/scikits.learn/scikits.learn-0.7.tar.gz
    % tar zxvf scikits.learn-0.7.tar.gz
    % cd scikits.learn-0.7

You can then build it locally and add it to your PYTHONPATH environment
variable::

    % python setup.py build_ext -i
    % export PYTHONPATH=`pwd`

If you want to install the library globally, do the following instead::

    % python setup.py build
    % sudo python setup.py install

Whatever the installation procedure you should check that the '0.7' version is
active in your python path::

    % python -c "import scikits.learn; print scikits.learn.__version__"
    0.7

You should also be able to launch the tests from anywhere in the system
(if nose is installed) with the following::

    % python -c "import scikits.learn as skl; skl.test()"

The output should end with ``OK`` as in::

    ----------------------------------------------------------------------
    Ran 623 tests in 26.108s

    OK (SKIP=2)


In the rest of the tutorial, the path to the extracted archive folder
``scikits.learn-0.7`` will be named ``$SKL_HOME``.



Download the datasets
---------------------

Machine Learning algorithms need data. Go to each ``$TUTORIAL_HOME/data``
sub-folder and run the ``fetch_data.py`` script from there (after
having read them first).

For instance::

    % cd $TUTORIAL_HOME/data/languages
    % less fetch_data.py
    % python fetch_data.py


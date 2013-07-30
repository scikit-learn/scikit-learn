Tutorial setup
==============

To get started with this tutorial, you firstly must have the
*scikit-learn* and all of its requiered dependencies installed.

Please refer to the `scikit-learn install`_ page for more information
and for per-system instructions.

.. _`scikit-learn install`: http://scikit-learn.sourceforge.net/install.html

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

    % cp -r skeletons work_directory/sklearn_tut_workspace


Download the datasets
---------------------

Machine Learning algorithms need data. Go to each ``$TUTORIAL_HOME/data``
sub-folder and run the ``fetch_data.py`` script from there (after
having read them first).

For instance::

    % cd $TUTORIAL_HOME/data/languages
    % less fetch_data.py
    % python fetch_data.py


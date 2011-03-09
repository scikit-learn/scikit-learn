.. -*- mode: rst -*-

About
=====

``scikit-learn`` is a python module for machine learning built on
top of numpy / scipy.

The purpose of the ``scikit-learn-tutorial`` subproject is to learn
how to apply machine learning to practical situations using the
algorithms implemented in the ``scikit-learn`` library.

The target audience is experienced Python developers familiar with
numpy and scipy.


Downloading the PDF
===================

Prebuilt versions of this tutorial are available from the `github download
page`_.

While following the exercices you might find helpful to use the official
`scikit-learn user guide (PDF)`_ as a more comprehensive reference::

If you need a numpy refresher please first have a look at the
`Scientific Python lecture notes (PDF)`_, esp. chapter 4.

.. _`github download page`: https://github.com/scikit-learn/scikit-learn-tutorial/archives/master
.. _`scikit-learn User Guide (PDF)`: http://downloads.sourceforge.net/project/scikit-learn/documentation/user_guide-0.7.pdf
.. _`Scientific Python lecture notes (PDF)`: http://scipy-lectures.github.com/_downloads/PythonScientific.pdf


Source code of the tutorial
===========================

The project is hosted on github at https://github.com/scikit-learn/scikit-learn-tutorial


Building the tutorial
=====================

You can build the HTML and PDF (requires pdflatex) versions of this
tutorial by installing sphinx (1.0.0+)::

  $ sudo pip install -U sphinx

Then for the html variant::

  $ cd tutorial
  $ make html

The results is available in the ``_build/html/`` subdolder. Point your browser
to the ``index.html`` file for table of content.

To build the PDF variant::

  $ make latex
  $ cd _build/latex
  $ pdflatex scikit_learn_tutorial.tex

You should get a file named ``scikit_learn_tutorial.pdf`` as output.


Mailing list
============

If you have questions about this tutorial you can ask them on the
``scikit-learn`` mailing list on sourceforge:
https://lists.sourceforge.net/lists/listinfo/scikit-learn-general


IRC channel
===========

Some developers tend to hang around the channel ``#scikit-learn``
at ``irc.freenode.net``, especially during the week preparing a new
release. If nobody is available to answer your questions there don't
hesitate to ask it on the mailing list to reach a wider audience.


License
=======

This tutorial is distributed under the Creative Commons Attribution
3.0 license. The python source code and exercices solutions are
distributed under the same license as the ``scikit-learn`` project
(Simplidied BSD).


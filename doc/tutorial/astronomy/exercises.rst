.. _astronomy_exercises

===================================
Exercises: Taking it a step further
===================================

Here we describe three exercises that will walk you through more advanced
approaches to the classification, regression, and dimensionality reduction
tasks outlined in the previous sections.

Before beginning the exercises, be sure you have downloaded the tutorial
source as described in the `setup <setup.html>`_ section.
For each exercise, there are two files.  In ``$TUTORIAL_HOME/skeletons``,
there is a skeleton script which can be filled-in to complete the exercise.
In ``$TUTORIAL_HOME/solutions`` you can find a completed version of the
exercise files.

To work on the tutorials, begin by copying the content of the ``skeletons``
folder to your own workspace::

    % cp -r skeletons workspace

You can then edit the contents of ``workspace`` without losing the original
files.  You must also make sure you have run the ``fetch_data.py`` scripts
in the appropriate subdirectory of ``$TUTORIAL_HOME/data``.
Next start the ipython interpreter and run the script with::

    In [1]: %run workspace/exercise*.py datadir

If an exception is triggered, you can type ``%debug`` to fire-up an ipdb
session.  

Each exercise contains the necessary import statements, data loading, and code
to evaluate the predictive accuracy of the model being used.  The remaining
code must be filled-in by the user.  These places are marked by the comment
``#TODO``.  Each exercise also has in-line descriptions which go over the
details of each sub-task.


Exercise 1: Photometric Classification with GMM
-----------------------------------------------

In this exercise, you will improve on the results of the classification
example described in the `classification <classification.html>`_ section.
We previously used a simple Gaussian Naive Bayes classifier to distinguish
between stars and quasars.  Here we will use Gaussian Mixture Models
to recreate the Gaussian Naive Bayes classifier, and then tweak the
parameters of these mixture models, evaluating them using a cross-validation
set, and attempt to arrive at a better classifier.

This task is broken into several parts:

    1. fill-in the function ``print_results``, which calculates the precision
       and recall of a set of classifications.

    2. Create a pair of single-component gaussian mixture models,
       calculate the priors for each class, and use these to compute the
       likelihoods for the cross-validation set.  This can be compared
       directly with the Gaussian Naive Bayes results

    3. Change the number of components and the covariance type to optimize
       the classifier using the cross-validation set.

    4. Once you've converged on a good set of GMM parameters, predict the
       labels for the test set, and compare this with the labels from the
       literature.

ipython command is::

    In [1]: %run workspace/exercise_01.py data/sdss_colors

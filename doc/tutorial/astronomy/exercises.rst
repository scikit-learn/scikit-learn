.. _astro_exercises:

===================================
Exercises: Taking it a step further
===================================

Here we describe three exercises that will walk you through more advanced
approaches to the classification, regression, and dimensionality reduction
tasks outlined in the previous sections.

Before beginning the exercises, be sure you have downloaded the tutorial
source as described in the `setup <setup.html>`_ section.
For each exercise, there are three files.  In ``$TUTORIAL_HOME/skeletons``,
there is a skeleton script which can be filled-in to complete the exercise.
In ``$TUTORIAL_HOME/solutions`` you can find a completed version of the
exercise files.  ``$TUTORIAL_HOME/notebooks`` has `iPython notebook`_
files which can be used to interactively complete the tutorials.

Using the notebook files is the preferred way to complete the tutorial.
If ipython notebook is unavailable, the files in the ``skeletons`` and
``solutions`` directory can be used.  Begin by copying the
content of the ``skeletons`` folder to your own workspace::

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


.. _astro_exercise_1:

Exercise 1: Photometric Classification with GMM
-----------------------------------------------

.. note::
   The information in this section is available in an interactive notebook
   :download:`10_exercise01.ipynb <notebooks/10_exercise01.ipynb>`,
   which can be viewed using `iPython notebook`_.  This requires the
   sdss colors data, which can be downloaded as described in
   :ref:`astro_tutorial_setup`.

In this exercise, you will improve on the results of the classification
example described in the `classification <classification.html>`_ section.
We previously used a simple Gaussian Naive Bayes classifier to distinguish
between stars and quasars.  Here we will use
`Gaussian Mixture Models <http://scikit-learn.org/0.6/modules/mixture.html>`_
to recreate the Gaussian Naive Bayes classifier, and then tweak the
parameters of these mixture models, evaluating them using a cross-validation
set, and attempt to arrive at a better classifier.

The overall goal of the task is to improve the quasar classifier, measuring
its performance using the precision, recall, and F1-score.

This task is broken into several parts:

    1. Re-implement Gaussian Naive Bayes using Gaussian Mixture models.
       This involves training two single-component GMMs, calculating the
       priors for each class, and using these to compute likelihoods.
       This can be compared directly with the Gaussian Naive Bayes results
       to check that the implementation is correct.

    2. Experiment with various covariance types and numbers of components
       to optimize the classifier using the cross-validation set.

    3. Once you've converged on a good set of GMM parameters, predict the
       labels for the test set, and compare this with the labels from the
       literature.

the ipython command is::

    In [1]: %run workspace/exercise_01.py data/sdss_colors

Althernatively, the ipython notebook file is found in

    :download:`notebooks/10_exercise01.ipynb`

.. _astro_exercise_2:

Exercise 2: Photometric redshifts with Decision Trees
-----------------------------------------------------

.. note::
   The information in this section is available in an interactive notebook
   :download:`11_exercise02.ipynb <notebooks/11_exercise02.ipynb>`,
   which can be viewed using `iPython notebook`_.  This requires the
   sdss photoz data, which can be downloaded as described in
   :ref:`astro_tutorial_setup`.

In this exercise, you will seek to improve the results of the photometric
redshift regression problem described in the `regression <regression.html>`_
section.  This exercise will draw heavily on the concepts of bias and
variance, using the learning curve plots introduced in 
`section 3 <practical.html>`_ of this tutorial.

There are two goals of this exercise:

    1. Find the best possible decision tree classifier for the data

    2. Decide what the best use of future resources is.  Should the
       astronomers seek to observe more objects (i.e. increase the number of
       training samples), or record more observations for each object
       (i.e. increase the number of features)?

The exercise is broken into the following tasks:

    1. Compute the training error and cross-validation error as a function
       of the ``max_depth`` parameter used in the Decision Tree Classifier.

    2. Compute the training error and cross-validation error as a function
       of the number of training samples.

    3. Repeat these two tasks, recording the outlier rate rather than the
       rms error.

    4. Analyze these results: should future observations focus on increasing
       the number of samples, or increasing the number of features?  Does
       this answer change depending on whether the rms error or outlier
       rate is the metric used?

the ipython command is::

    In [1]: %run workspace/exercise_02.py data/sdss_photoz/

Althernatively, the ipython notebook file is found in

    :download:`notebooks/11_exercise02.ipynb`

.. _astro_exercise_3:

Exercise 3: Dimensionality Reduction of Spectra
-----------------------------------------------

In this exercise, you will use several dimensionality reduction techniques
to view low-dimensional projections of galaxy & quasar spectra from the
Sloan Digital Sky Survey.  This exercise is much less quantitative than the
previous ones: it mainly will help you to get a qualitative sense of the
characteristics of these learning methods.

There is a programming section, followed by an experimentation section.  The
skeleton is set up to use command-line options to compare different sets of
parameters

Programming
~~~~~~~~~~~
The file has several places with "TODO" marked.  In these, you will use the
specified unsupervised method to project the data ``X`` into the
lower-dimensional ``X_proj``.

   1. Use :class:`sklearn.decomposition.RandomizedPCA` to project the data

      the ipython command is::

      	  In [1]: %run workspace/exercise_03.py data/sdss_spectra/ -m pca

      Note the argument ``-m`` which specifies the method  to use.

   2. Use :class:`sklearn.manifold.LocallyLinearEmbedding` with
      ``method='standard'`` to project the data.

      the ipython command is::

      	  In [1]: %run workspace/exercise_03.py data/sdss_spectra/ -m lle

   3. Use :class:`sklearn.manifold.LocallyLinearEmbedding` with
      ``method='standard'`` to project the data.

      the ipython command is::

      	  In [1]: %run workspace/exercise_03.py data/sdss_spectra/ -m mlle

   4. Use :class:`sklearn.manifold.Isomap` to project the data.

      the ipython command is::

      	  In [1]: %run workspace/exercise_03.py data/sdss_spectra/ -m isomap

Experimentation
~~~~~~~~~~~~~~~
Your goal is to find a projection that does a good job of separating the
various classes of spectra, and lays them out in a way that might allow
intuitive evaluation of the relationships between points.  The script is
set-up as a command-line interface.  You should address the following
questions:

   1. How sensitive is PCA to the set of data used?  To the number of 
      training points?  You can test this out as follows::

          In [1]: %run workspace/exercise_03.py data/sdss_spectra -m pca -n 1000 -s

      This will perform PCA on a subset of 1000 points.  ``-s`` indicates that
      the data should be shuffled, so that the set of points is different every
      time.  How stable is the projection between different subsets of the
      data?  How does the projection change as the number of points is
      increased?

   2. Address the same questions with LLE, MLLE, and Isomap.  Which of these
      manifold methods appears to give the most stable results?

   3. Now we can vary the number of neighbors used with LLE, MLLE, and Isomap.
      This is accomplished as follows::

          In [1]: %run workspace/exercise_03.py data/sdss_spectra -m lle -k 20

      This call will execute LLE with 20 neighbors.  Try this for several
      values of `k`.  How does the number of
      neighbors change the projection?  Among LLE, MLLE, and Isomap, which
      produces the most stable results as the number of neighbors are changed?

   4. Finally, we'll test the effects of normalization.  This can be done
      as follows::

          In [1]: %run workspace/exercise_03.py data/sdss_spectra -N l2

      this will perform PCA with L2-normalization.  The other options are
      ``-N l1`` for L1-normalization, and ``-N none`` for no normalization.
      Normalization has the effect of bringing all the spectra closer
      together: unnormalized spectra may be very bright (for nearby objects)
      or very dim (for far away objects).  Normalization corrects for this
      source of variance in the data.  How do the projected results change
      as you vary the normalization?

   5. By now, you should have an idea of which method and which combination of
      parameters give the best qualitative separation between the points.
      Re-run this method using the full `n`=4000 dataset now::

          In [1]: %run python workspace/exercise_03.py data/sdss_spectra -n 4000 -m [method] [other options]

      This should give you a projection of the data that gives a good
      visualization of the relationship between points.  An astronomer may
      go further and try to develop rough cut-offs that would give a broad
      classification to an unlabeled test point.  This sort of procedure could
      be used as the first step of a physically-motivated classification
      pipeline, or to flag potentially interesting objects for quick
      followup.


.. _`iPython notebook`: http://ipython.org/ipython-doc/stable/interactive/htmlnotebook.html

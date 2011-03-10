Exercises
=========

To do the exercises, copy the content of the 'skeletons' folder as
a new folder named 'workspaces'.

Then fire an ipython shell and run the work-in-progress script with::

  >>> %run workspace/exercise_XX_script.py arg1 arg2 arg3

If an exception is triggered, use the ``%debug`` to fire-up a post
mortem ipdb session.

Refine the implementation and iterate until the exercise is solved.


Exercise 1: Sentiment Analysis on movie reviews
-----------------------------------------------

- Write a text classification pipeline to classify movie reviews as either
  positive or negative.

- Find a good set of parameters using grid search.

- Evaluate the performance on a held out test set.


Exercise 2: Language identification
-----------------------------------

- Write a text classification pipeline using a custom preprocessor and
  CharNGramAnalyzer using data from Wikipedia articles as training set.

- Evaluate the performance on some held out test set.


Exercise 3: CLI text classification utility
-------------------------------------------

Using the results of the previous exercises and the ``cPickle``
module of the standard library, write a command line utility that
detect the language of some text provided on ``stdin`` and estimate
the polarity (positive or negative) if the text is written in
English.

Bonus point if the utility is able to give a confidence level for its
predictions.


Exercise 4: Face recognition
----------------------------

Build a classifier that recognize person on faces pictures from the
Labeled Faces in the Wild dataset.


Exercises
=========

To do the exercises, copy the content of the 'skeletons' folder as
a new folder named 'workspace'::

  % cp -r skeletons workspace

You can then edit the content of the workspace without fear of loosing
the original exercise instructions.

Then fire an ipython shell and run the work-in-progress script with::

  [1] %run workspace/exercise_XX_script.py arg1 arg2 arg3

If an exception is triggered, use ``%debug`` to fire-up a post
mortem ipdb session.

Refine the implementation and iterate until the exercise is solved.

**For each exercise, the skeleton file provides all the necessary import
statements, boilerplate code to load the data and sample code to evaluate
the predictive accurracy of the model.**


Exercise 1: Language identification
-----------------------------------

- Write a text classification pipeline using a custom preprocessor and
  ``CharNGramAnalyzer`` using data from Wikipedia articles as training set.

- Evaluate the performance on some held out test set.

ipython command line::

  %run workspace/exercise_01_language_train_model.py data/languages/paragraphs/


Exercise 2: Sentiment Analysis on movie reviews
-----------------------------------------------

- Write a text classification pipeline to classify movie reviews as either
  positive or negative.

- Find a good set of parameters using grid search.

- Evaluate the performance on a held out test set.

ipython command line::

  %run workspace/exercise_02_sentiment.py data/movie_reviews/txt_sentoken/


Exercise 3: CLI text classification utility
-------------------------------------------

Using the results of the previous exercises and the ``cPickle``
module of the standard library, write a command line utility that
detects the language of some text provided on ``stdin`` and estimate
the polarity (positive or negative) if the text is written in
English.

Bonus point if the utility is able to give a confidence level for its
predictions.


Exercise 4: Face recognition
----------------------------

Build a classifier that recognizes persons on face pictures from the
Labeled Faces in the Wild dataset.

ipython command line::

  %run workspace/exercise_04_face_recognition.py data/data/labeled_faces_wild/lfw_preprocessed/

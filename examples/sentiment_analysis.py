"""
Train a text classification model using a TF-IDF vectorizer and linear SVM 
classifier.

The model is trained on a subset of the '20 newsgroups' dataset filtered 
by the following categories:
1. 'alt.atheism'
2. 'soc.religion.christian'
3. 'comp.graphics'
4. 'sci.med'

It then evaluates the model using precision, recall, and F1-score.

Steps:
1. Load dataset.
2. Split into training and testing sets.
3. Create a pipeline for text vectorization and classification.
4. Train the model on the training data.
5. Make predictions on the test data.
6. Evaluate the model's performance.

This is the best model
"""



import numpy as np

from sklearn.feature_extraction.text import TfidfVectorizer

from sklearn.svm import SVC

from sklearn.pipeline import make_pipeline

from sklearn.datasets import fetch_20newsgroups

from sklearn.model_selection import train_test_split

from sklearn import metrics

 

# Load dataset

categories = ['alt.atheism', 'soc.religion.christian', 'comp.graphics', 
'sci.med']

newsgroups_train = fetch_20newsgroups(subset='train', 
categories=categories)

newsgroups_test = fetch_20newsgroups(subset='test', categories=categories)

 

# Split data

X_train, X_test, y_train, y_test = train_test_split(newsgroups_train.data, 
newsgroups_train.target, test_size=0.3, random_state=42)

 

# Create a pipeline

text_clf = make_pipeline(TfidfVectorizer(), SVC(kernel='linear'))

 

# Train the model

text_clf.fit(X_train, y_train)

 

# Predict

predicted = text_clf.predict(X_test)

 

# Evaluate

print(metrics.classification_report(y_test, predicted))



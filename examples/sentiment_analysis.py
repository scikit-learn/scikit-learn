import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline
from sklearn.datasets import fetch_20newsgroups
from sklearn.model_selection import train_test_split
from sklearn import metrics

# Load dataset
categories = ['alt.atheism', 'soc.religion.christian', 'comp.graphics', 'sci.med']
newsgroups_data = fetch_20newsgroups(subset='train', categories=categories)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(newsgroups_data.data, newsgroups_data.target, test_size=0.2, random_state=42)

# Create a pipeline with TF-IDF Vectorizer and SVM
model = make_pipeline(TfidfVectorizer(), SVC())

# Train the model
model.fit(X_train, y_train)

# Predict on the test data
predicted = model.predict(X_test)

# Print the classification report
print(metrics.classification_report(y_test, predicted, target_names=newsgroups_data.target_names))

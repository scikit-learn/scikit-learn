=============================================
A Beginner's guide to Comparing ML Algorithms
=============================================

As a new machine learning engineer, some of the algorithms may seem really complex. What's more, deciding which algorithm to use (even when they're available) can be somewhat of a challenge. This guide walks you through the basic machine learning algorithms and when to use them.

.. _Linear_Regression:

Linear Regression
=================

- Linear Regression is widely considered one of the most basic machine learning algorithms.  It is used to find the relationship between dependent and independent variables. Linear Regression is often used for the following problems in Machine Learning:

- Evaluate business trends and make estimates/forecasts. (Predict advertising and marketing profits, housing prices, etc)

- Evaluate the relationship between age,gender,height (independent) and weight (dependent) to predict the risk of obesity.

And any other problem, where it appears that there is a direct relationship between a set of features and an output variable. 

- Link to our linear regression model documentation:
 <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html>`_

.. _Logistic_Regression:

Logistic Regression
===================

- Logistic Regression is a basic machine learning algorithm which is used for classification problems,based on the concept of probability. It is often used when the dependent variable is categorical, i.e the decision to be made is a binary (Yes/No) one.

- Logistic Regression is often used for the following problems:

- Spam classification: Determining whether an email is spam or not.

- Tumor Malignancy prediction: Predicting whether a tumor is malignant or not.

- Speech Classification: Determining whether a certain dialouge or comment is sexist/racist or offensive to any community in general.

- Predicting Customer loyalty: Determing whether a company may be able to retain a customer, depending on the kind of services they have provided.

and any other problem, where the answer seems to be a simple 'yes' or 'no'.

- Link to our logistic regression model documentation:
 <https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html>`_

.. _Support_Vector_Machine:

Support Vector Machine
======================

- The objective of the support vector machine algorithm is to find a hyperplane in an N-dimensional space(N — the number of features) that distinctly classifies the data points. SVMs are also often used for classification problems. 

- The main difference between Logistic Regression and SVM is that while the former fits the data points as if they are along a continuous function, the latter fits a function (hyperplane) that attempts to separate two classes of data that could be of multiple dimensions.

- SVMs are used for the following purposes:

- Face Detection (Image detection in general)

- Handwriting recognition

- Gene classification 

- Text Classification 

or any other problems, where it seems necessary to be able to distinguish between multiple classes in order to classify the independent variable.

- Link to our Supoort Vector Machine model documentation:
 <https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html>`_





"""
===============================================
Interpretation of coefficients in linear models
===============================================

Linear models describe situations in which the target value is expected to be
a linear combination of the features (see the :ref:`linear_model` User Guide
section for a description of a set of linear model methods available in
scikit-learn).
It is important to emphasize that coefficients in multiple linear models
represent the relationship between the given feature and the target
assuming that other features remain constant.

This example will provide some hints in interpreting coefficient in linear
models, using data from the "Current Population Survey" from 1985.
We will be interested in the prediction of the wage as a function of various
features such as experience, age, or education.

A description of the dataset follows.
"""

print(__doc__)

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#############################################################################
# Determinants of Wages from the 1985 Current Population Survey
# -------------------------------------------------------------
#
# The dataset
# ...........
#
# We fetch the data from `OpenML <http://openml.org/>`_.
# Note that setting the parameter `as_frame` to True will retrieve the data
# as a pandas dataframe.

from sklearn.datasets import fetch_openml

survey = fetch_openml(data_id=534, as_frame=True)

##############################################################################
# Then, we identify features `X` and targets `y`: the column WAGE is our
# target variable (i.e., the variable which we want to predict).
#
X = survey.data[survey.feature_names]
X.describe(include="all")

##############################################################################
# Notice that the dataset contains categorical and numerical variables.
# Some of the categorical variables are binary variables.
# About the numerical ones we can observe that AGE and EXPERIENCE have similar
# distributions while the EDUCATION distribution is narrower.
# This will give us directions on how to preprocess the data thereafter.

X.head()

##############################################################################
y = survey.target.values.ravel()
survey.target.head()

###############################################################################
# We split the sample in a train and a test dataset,
# Only the train dataset will be used in the following exploratory analysis.
# This is a way to emulate a real situation where predictions are performed on
# an unknown target.

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, random_state=42
)

##############################################################################
# First, let's get some insights by looking at the marginal links between the
# different variables. Only numerical variables will be used.

train_dataset = X_train.copy()
train_dataset.insert(0, "WAGE", y_train)
sns.pairplot(train_dataset, diag_kind='kde')
plt.show()

##############################################################################
# Looking closely at the WAGE distribution it could be noticed that it has a
# long tail and we could take its logarithm
# to simplify our problem and approximate a normal distribution.
# The WAGE is increasing when EDUCATION is increasing.
# Also, the EXPERIENCE and AGE are linearly correlated.
#
# The pipeline
# ............
#
# To design our machine-learning pipeline, we manually
# check the type of data that we are dealing with:

survey.data.info()

#############################################################################
# As seen previously, the dataset contains columns with different data types
# and we need to apply a specific preprocessing for each data types.
# In particular categorical variables cannot be included in linear model if not
# coded as integers first. In addition, to avoid categorical features to be
# treated as ordered values, we need to one-hot-encode them.
# Our pre-processor will
#
# - one-hot encode (i.e., generate a column by category) the categorical
#   columns;
# - replace by 0 and 1 the categories of binary columns;
# - as a first approach (we will see after how the normalisation of numerical
#   values will affect our discussion), keep numerical values as they are.

from sklearn.compose import make_column_transformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder

categorical_columns = ['RACE', 'OCCUPATION', 'SECTOR']
binary_columns = ['MARR', 'UNION', 'SEX', 'SOUTH']
numerical_columns = ['EDUCATION', 'EXPERIENCE', 'AGE']

preprocessor = make_column_transformer(
    (OneHotEncoder(), categorical_columns),
    (OrdinalEncoder(), binary_columns),
    remainder='passthrough'
)

##############################################################################
# To describe the dataset as a linear model we choose to use a ridge regressor
# and to model the logarithm of the WAGE.
# We sample the complexity parameter space between 1.e-10 and 1.e10.

from sklearn.pipeline import make_pipeline
from sklearn.linear_model import RidgeCV
from sklearn.compose import TransformedTargetRegressor

model = make_pipeline(
    preprocessor,
    TransformedTargetRegressor(
        regressor=RidgeCV(alphas=np.logspace(-10,10,21)),
        func=np.log10,
        inverse_func=sp.special.exp10
    )
)

##############################################################################
# Processing the dataset
# ......................
#
# First, we fit the model and we verify which value for :math:`\alpha` has been
# selected.

model.fit(X_train, y_train)
model[-1].regressor_.alpha_

##############################################################################
# Once verified that the :math:`\alpha` parameter is not at the boundary of
# the sampled parameter space, we can check the performance of the computed
# model using, for example, the median absolute error of the model and the R
# squared coefficient.

from sklearn.metrics import median_absolute_error

y_pred = model.predict(X_train)
mae = median_absolute_error(y_train, y_pred)
string_score = 'MAE on training set: {0:.2f} $/hour'.format(mae)
y_pred = model.predict(X_test)
mae = median_absolute_error(y_test, y_pred)
r2score = model.score(X_test,y_test)

string_score += '\nMAE on testing set: {0:.2f} $/hour'.format(mae)
string_score += '\nR2 score: {0:.4f}'.format(r2score)
fig, ax = plt.subplots(figsize=(6, 6))
sns.regplot(y_test, y_pred)

plt.text(3, 20, string_score)

plt.ylabel('Model predictions')
plt.xlabel('Truths')
plt.xlim([0, 27])
plt.ylim([0, 27])

##############################################################################
# The model learnt is far from being a good model making accurate predictions:
# the R squared score is very low.
# As interpretation tools characterize model rather than the generative process
# of the data itself, it needs to be emphasized that interpretations are
# correct if the model is correct as well.
# In this case, we are more interested in providing a methodology than in
# having a good description of the data: a bad example illustrates the
# importance of cross checking the results.
#
# Interpreting coefficients
# .........................
#
# First of all, we can plot the values of the coefficients of the regressor we
# have fitted.

feature_names = (model.named_steps['columntransformer']
                      .named_transformers_['onehotencoder']
                      .get_feature_names(input_features=categorical_columns))
feature_names = np.concatenate(
    [feature_names, binary_columns, numerical_columns])

coefs = pd.DataFrame(
    model.named_steps['transformedtargetregressor'].regressor_.coef_,
    columns=['Coefficients'], index=feature_names
)
coefs.plot(kind='barh', figsize=(9, 7))
plt.axvline(x=0, color='.5')
plt.subplots_adjust(left=.3)

###############################################################################
# Soon we realize that we cannot compare different coefficients since we did
# not scale the data before the fit, features having different value ranges.
# For instance, the AGE coefficient is expressed in $/hours/leaving years
# while the EDUCATION is expressed in $/hours/years of education.
# This is evident if we compare feature standard deviations.

X_train_preprocessed = pd.DataFrame(
    model.named_steps['columntransformer'].transform(X_train),
    columns=feature_names
)
X_train_preprocessed.std().plot(kind='barh', figsize=(9, 7))
plt.title('Features std. dev.')
plt.subplots_adjust(left=.3)

###############################################################################
# We should then normalize the coefficients by the standard deviation and we
# will be able to compare them.

coefs = pd.DataFrame(
    model.named_steps['transformedtargetregressor'].regressor_.coef_ *
    X_train_preprocessed.std(),
    columns=['Coefficients'], index=feature_names
)
coefs.plot(kind='barh', figsize=(9, 7))
plt.axvline(x=0, color='.5')
plt.subplots_adjust(left=.3)

###############################################################################
# The plot above tells us that an increase of the AGE will induce a decrease
# of the WAGE when all other features remain constant. Also an increase of
# the EXPERIENCE will induce an increase of the WAGE when all other
# features remain constant.
#
# The first interpretation might look counter-intuitive at first, if one
# relates the relationship between AGE and WAGE as a marginal link.
# However, as previously mentioned, a linear model computes a conditional
# link between AGE and WAGE given all other features.
# Therefore, one should also interpret that for a given experience and all
# other features constant as well, a younger person would have a higher
# wage.
#
# Checking the coefficient stability
# ..................................
#
# The stability of the coefficients is a guarantee of the robustness of the
# model. We can check the coefficient stability through cross-validation.

from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedKFold

cv_model = cross_validate(
    model, X, y, cv=RepeatedKFold(n_splits=5, n_repeats=5),
    return_estimator=True, n_jobs=-1
)
coefs = pd.DataFrame(
    [est.named_steps['transformedtargetregressor'].regressor_.coef_ *
     X_train_preprocessed.std()
     for est in cv_model['estimator']],
    columns=feature_names
)
plt.figure(figsize=(9, 7))
sns.swarmplot(data=coefs, orient='h', color='k', alpha=0.5)
sns.boxplot(data=coefs, orient='h', color='cyan')
plt.axvline(x=0, color='.5')
plt.title('Stability of coefficients')
plt.subplots_adjust(left=.3)

###############################################################################
# The AGE and EXPERIENCE coefficients are highly unstable which might be
# due to the collinearity between the 2 features.
#
# In order to verify our interpretation we remove one of the 2 features and
# check what is the impact on the model stability.

column_to_drop = ['AGE']

cv_model = cross_validate(
    model, X.drop(columns=column_to_drop), y,
    cv=RepeatedKFold(n_splits=5, n_repeats=5),
    return_estimator=True, n_jobs=-1
)
coefs = pd.DataFrame(
    [est.named_steps['transformedtargetregressor'].regressor_.coef_ *
     X_train_preprocessed.drop(columns=column_to_drop).std()
     for est in cv_model['estimator']],
    columns=feature_names[:-1]
)
plt.figure(figsize=(9, 7))
sns.swarmplot(data=coefs, orient='h', color='k', alpha=0.5)
sns.boxplot(data=coefs, orient='h', color='cyan')
plt.axvline(x=0, color='.5')
plt.title('Stability of coefficients')
plt.subplots_adjust(left=.3)

###############################################################################
# The estimation of the EXPERIENCE coefficient is now more stable and
# remain important for all predictors trained during cross-validation.

"""
===============================================
Interpretation of coefficients in linear models
===============================================

Linear models describe situations in which the target value is expected to be
a linear combination of the features (see the :ref:`linear_model` User guide
section for a description of a set of linear model methods available in
scikit-learn).
It is important to emphasize that linear models compute conditional links.
The interpretation of the coefficient gives the relationship between the
feature and the target given that other features remain constant.

This example will show some hints in interpreting coefficient in linear models, 
using data from the "Current Population Survey" from 1985.
"""

print(__doc__)

from time import time
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import RidgeCV
from sklearn.compose import TransformedTargetRegressor
from sklearn.metrics import median_absolute_error
from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedKFold

#############################################################################
# Determinants of Wages from the 1985 Current Population Survey
# -------------------------------------------------------------
#
# First of all we fetch the data from `OpenML <http://openml.org/>`_.
# Note that setting the parameter `as_frame` to `True` will retrieve the data
# as a pandas dataframe.
# Then, we identify features (`X`) and targets (`y`): the column 'WAGE' is our
# target variable (i.e., the variable which we want to predict).
 
survey = fetch_openml(data_id=534, as_frame=True)

X = survey.data[survey.feature_names]
y = survey.target.values.ravel()
ax = sns.kdeplot(y, shade=True, color="r")
plt.xlabel(survey.target_names)
plt.show()

##############################################################################
# Note that the "WAGE" distribution has a long tail and we could take its log
# to simplify our problem getting closer to a normal distribution.
#
# The dataset is composed by columns with different data types and we need to
# apply a specific preprocessing for each data types.
# Our pre-processor will
# 
# - one-hot encode (i.e., generate a column by category) the categorical
#   columns;
# - replace by 0 and 1 the categories of binary columns;
# - keep numerical values as they are.

categorical_columns = ['RACE', 'OCCUPATION', 'SECTOR']
binary_columns = ['MARR', 'UNION', 'SEX', 'SOUTH']
numerical_columns = ['EDUCATION', 'EXPERIENCE', 'AGE']

preprocessor = make_column_transformer(
    (OneHotEncoder(), categorical_columns),
    (OrdinalEncoder(), binary_columns),
    remainder='passthrough'
)

##############################################################################
# Modeling the data
# .................
#
# We will fit a ridge regressor and transform the target before the fit using
# a log transform.
# But before computing the model we split the sample in a train and a test
# dataset.

model = make_pipeline(
    preprocessor,
    TransformedTargetRegressor(
        regressor=RidgeCV(),
        func=np.log10,
        inverse_func=sp.special.exp10
    )
)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, random_state=42
)

model.fit(X_train, y_train);


##############################################################################
# Scoring the model
# .................
#
# We can check the performance of the computed model using, for example, the
# median absolute error of the model.  

def mae_scorer(model, X_train, X_test, y_train, y_test):
    y_pred = model.predict(X_train)
    string_score = f'MAE on training set: {median_absolute_error(y_train, y_pred):.2f} $/hour'
    y_pred = model.predict(X_test)
    string_score += f'\nMAE on testing set: {median_absolute_error(y_test, y_pred):.2f} $/hour'
    return string_score

fig, ax = plt.subplots(figsize=(6, 6))
y_pred = model.predict(X_test)
sns.regplot(y_test, y_pred)

plt.text(3, 20, mae_scorer(model, X_train, X_test, y_train, y_test))

plt.ylabel('Model predictions')
plt.xlabel('Truths')
plt.xlim([0, 27])
plt.ylim([0, 27]);

##############################################################################
# The model learnt is far to be a good model making accurate prediction.
# As interpretation tools characterize model rather than the generative process
# of the data itself, it needs to be emphasized that interpretations are correct
# if the model is correct as well.

##############################################################################
# Interpreting coefficients
# .........................
#
# First of all, we can plot the values of the coefficients of the regressor we
# have fitted.

feature_names = (model.named_steps['columntransformer']
                      .named_transformers_['onehotencoder']
                      .get_feature_names(input_features=categorical_columns))
feature_names = np.concatenate([feature_names, binary_columns, numerical_columns])

coefs = pd.DataFrame(
    model.named_steps['transformedtargetregressor'].regressor_.coef_,
    columns=['Coefficients'], index=feature_names
)
coefs.plot(kind='barh', figsize=(9, 7))
plt.axvline(x=0, color='.5');

###############################################################################
# Soon we realize that we cannot compare different coefficients since we did
# not scale the data before the fit features having different value ranges.
# For instance, the "AGE" coefficient is expressed in $/hours/leaving years
# while the "EDUCATION" is expressed in $/hours/years of education.
# This is evident if we compare feature standard deviations.

X_train_preprocessed = pd.DataFrame(
    model.named_steps['columntransformer'].transform(X_train),
    columns=feature_names
)
X_train_preprocessed.std().plot(kind='barh', figsize=(9, 7))
plt.title('Features std. dev.');

###############################################################################
# We can then normalize the coefficients by the standard deviation and we will
# be able to compare them.

coefs = pd.DataFrame(
    model.named_steps['transformedtargetregressor'].regressor_.coef_ *
    X_train_preprocessed.std(),
    columns=['Coefficients'], index=feature_names
)
coefs.plot(kind='barh', figsize=(9, 7))
plt.axvline(x=0, color='.5');

###############################################################################
# The plot above tells us that an increase of the "AGE" will induce a decrease
# of the "WAGE" when all other features remain constant. Also an increase of
# the "EXPERIENCE" will induce an increase of the "WAGE" when all other
# features remain constant.
#
# The first interpretation might look counter-intuitive at first, if one relates
# the relationship between "AGE" and "WAGE" as a marginal link.
# However, as previously mentioned, a linear model computes a conditional
# link between "AGE" and "WAGE" given all other features.
# Therefore, one should also interpret that for a given experience (and all other
# features constant as well ...), a younger person would have an higher wage.
#
# Checking the coefficient stability
# ..................................
#
# The stability of the coefficients is a guarantee of the robustness of the
# model. We can check the coefficient stability through cross-validation.

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
sns.boxenplot(data=coefs, orient='h', color='C0')
plt.axvline(x=0, color='.5')
plt.title('Stability of coefficients');

###############################################################################
# The "AGE" and "EXPERIENCE" coefficients are highly instable which might be
# due to the collinearity between the 2 features.

age = survey.data['AGE'].values
experience = survey.data['EXPERIENCE'].values
sns.regplot(age,experience,scatter_kws={"color": "black", "alpha": 0.2, "s": 30},
    line_kws={"color": "red"})

##############################################################################
# We can remove one of the 2 features and check what is the impact on the
# features stability.


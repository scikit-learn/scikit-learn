"""
==================================================
Partial Dependence Plots with Categorical Features
==================================================

In this example, we plot partial dependence of a trained regression model on
various categorical features using the Ames housing dataset. Partial dependence
plots show the dependence between the target function [1]_ and a set of
features of interest, marginalizing over the values of all other features
(the complement features). In the case of categorical features, the partial
dependence plots are bar charts with individual bars showing the partial
dependence of each category.

.. [1] For classification you can think of it as the regression score before
       the link function.

.. note::
    See also
    :ref:`sphx_glr_auto_examples_inspection_plot_partial_dependence.py`
"""
print(__doc__)
import matplotlib.pyplot as plt
import numpy as np

from sklearn.compose import make_column_transformer
from sklearn.datasets import fetch_openml
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import plot_partial_dependence
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle


def load_ames_housing():
    df = fetch_openml(name="house_prices", as_frame=True)
    X = df.data
    y = df.target

    features = ['YrSold', 'HeatingQC', 'Street', 'YearRemodAdd', 'Heating',
                'MasVnrType', 'BsmtUnfSF', 'Foundation', 'MasVnrArea',
                'MSSubClass', 'ExterQual', 'Condition2', 'GarageCars',
                'GarageType', 'OverallQual', 'TotalBsmtSF', 'BsmtFinSF1',
                'HouseStyle', 'MiscFeature', 'MoSold']

    X = X[features]
    X, y = shuffle(X, y, random_state=0)

    X = X[:600]
    y = y[:600]
    return X, np.log(y)


X, y = load_ames_housing()


# %%
# Preprocess the data and build a regression model
##############################################################################
#
# In preparing the dataset, we will impute the missing values (with new
# category 'missing' for categorical columns and mean value for numerical
# columns) and encode the categories with
# :class:`sklearn.preprocessing.OneHotEncoder
# <sklearn.preprocessing.OneHotEncoder>`

cat_cols = X.columns[X.dtypes == 'O']
num_cols = X.columns[X.dtypes == 'float64']

categories = [X[column].unique() for column in X[cat_cols]]

for cat in categories:
    cat[cat == None] = 'missing'  # noqa

cat_proc_nlin = make_pipeline(
    SimpleImputer(missing_values=None, strategy='constant',
                  fill_value='missing'),
    OneHotEncoder(categories=categories)
)

num_proc_nlin = make_pipeline(
    SimpleImputer(strategy='mean'),
    StandardScaler()
)

processor_nlin = make_column_transformer(
    (cat_proc_nlin, cat_cols),
    (num_proc_nlin, num_cols),
    remainder='passthrough'
)

# Fit a :class:`~sklearn.ensemble.RandomForestRegressor`.

rf_pipeline = make_pipeline(processor_nlin, RandomForestRegressor())
rf_pipeline.fit(X, y)


# %%
# Single-variable partial dependence plots
##############################################################################
#
# Let's compute single-variable partial dependence plots for some of the
# categorical variables.

features_to_plot = ['ExterQual', 'HeatingQC', 'Street', 'HouseStyle']
categorical = [X[f].dtype == 'O' for f in features_to_plot]
plot_partial_dependence(rf_pipeline, X, features=features_to_plot,
                        is_categorical=categorical, n_cols=2)
fig = plt.gcf()
plt.subplots_adjust(hspace=0.5)
plt.show()

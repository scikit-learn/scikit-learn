"""
=================================
Combine predictors using stacking
=================================

Stacking refers to a method to blend estimators. In this strategy, some
estimators are individually fitted on some training data while a final
estimator is trained using the stacked predictions of these base estimators.

In this example, we illustrate the use case in which different regressors are
stacked together and a final linear penalized regressor is used to output the
prediction. We compare the performance of each individual regressor with the
stacking strategy. Stacking slightly improves the overall performance.

"""
print(__doc__)

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
#          Maria Telenczuk    <https://github.com/maikia>
# License: BSD 3 clause

###############################################################################
# The function ``plot_regression_results`` is used to plot the predicted and
# true targets.

import matplotlib.pyplot as plt


def plot_regression_results(ax, y_true, y_pred, title, scores, elapsed_time):
    """Scatter plot of the predicted vs true targets."""
    ax.plot([y_true.min(), y_true.max()],
            [y_true.min(), y_true.max()],
            '--r', linewidth=2)
    ax.scatter(y_true, y_pred, alpha=0.2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.set_xlim([y_true.min(), y_true.max()])
    ax.set_ylim([y_true.min(), y_true.max()])
    ax.set_xlabel('Measured')
    ax.set_ylabel('Predicted')
    extra = plt.Rectangle((0, 0), 0, 0, fc="w", fill=False,
                          edgecolor='none', linewidth=0)
    ax.legend([extra], [scores], loc='upper left')
    title = title + '\n Evaluation in {:.2f} seconds'.format(elapsed_time)
    ax.set_title(title)


###############################################################################
# Download and prepare the dataset
###############################################################################
#
# We will use Ames Housing dataset which is a set of 1460 residential homes
# in Ames, Iowa. Each of them is described by 80 features and the task is to
# predict the final price of the houses. Here, we select only 20 most
# intersting features chosen using GradientBoostingRegressor
# (TODO: add link, check correct
# https://github.com/scikit-learn/scikit-learn/pull/16400).
# https://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_regression.html
#
# Ames Housing dataset is not a part of Sklearn and we will download it from OpenML (TODO:
# link: https://www.openml.org/d/42165).


from sklearn.datasets import fetch_openml


def load_ames_housing():
    df = fetch_openml(data_id=42165, as_frame=True)
    X = df.data
    y = df.target

    features = ['YrSold', 'HeatingQC', 'Street', 'YearRemodAdd', 'Heating',
       'MasVnrType', 'BsmtUnfSF', 'Foundation', 'MasVnrArea', 'MSSubClass',
       'ExterQual', 'Condition2', 'GarageCars', 'GarageType', 'OverallQual',
       'TotalBsmtSF', 'BsmtFinSF1', 'HouseStyle', 'MiscFeature', 'MoSold']

    X = X[features]
    return X, y

X, y = load_ames_housing()

from sklearn.pipeline import make_pipeline
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.compose import make_column_transformer


###############################################################################
# Stack of predictors on a single data set
###############################################################################
#
# It is sometimes tedious to find the model which will best perform on a given
# dataset. Stacking provide an alternative by combining the outputs of several
# learners, without the need to choose a model specifically. The performance of
# stacking is usually close to the best model and sometimes it can outperform
# the prediction performance of each individual model.
#
# Here, we combine 3 learners (linear and non-linear) and use a ridge regressor
# to combine their outputs together.

from sklearn.ensemble import StackingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LassoCV
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from sklearn.compose import make_column_transformer

#TODO: add ridge encoder

# The data we are going to use is going to be downloaded from the OpenMl, and
# therefore it still need to be preprocessed for our use. First, we have some
# categorical columns, second there are plenty of missing values in the
# dataset. Therefore we wil start from making a pipelines, one for our linear
# estimators and the other one for decision trees





#X = X[:500]
#y = y[:500]

cat_cols = X.columns[X.dtypes == 'O']
num_cols = X.columns[X.dtypes == 'float64']

#
categories = [
    X[column].unique() for column in X[cat_cols]]

for cat in categories:
    cat[cat == None] = 'missing'

categorical_proc_tree = Pipeline(steps=[
    ('imputer_none', SimpleImputer(missing_values=None, strategy='constant', fill_value='missing')),
    ('encoder', OrdinalEncoder(categories=categories))
    ])

# only remove the missing values
numerical_proc_tree = Pipeline(steps=[
    ('imputer_nan', SimpleImputer(strategy='mean'))])

#categorical_proc_tree.fit(X[cat_cols],y)
#categorical_proc_tree.transform(X[cat_cols])
#numerical_proc_tree.fit(X[num_cols],y)
#numerical_proc_tree.transform(X[num_cols])

categorical_proc_lin = make_pipeline(
    SimpleImputer(missing_values=None, strategy='constant', fill_value='missing'),
    OneHotEncoder(categories=categories)
)

numerical_proc_lin = make_pipeline(
    StandardScaler(),
    SimpleImputer(strategy='mean'),
)

#import pdb; pdb.set_trace()
#estimators = [
#    ('Random Forest', RandomForestRegressor(random_state=42)),
#    ('Lasso', LassoCV()),
#    ('Gradient Boosting', HistGradientBoostingRegressor(random_state=0))
#]

# transformation to use for decision tree based estimators
processor_tree = make_column_transformer(
                              (categorical_proc_tree, cat_cols),
                              (numerical_proc_tree, num_cols),
                              remainder = 'passthrough')

# transformation to use for linear estimators
processor_lin = make_column_transformer(
                              (categorical_proc_lin, cat_cols),
                              (numerical_proc_lin, num_cols),
                              remainder = 'passthrough')
#import pdb; pdb.set_trace()
lasso_pip = make_pipeline(processor_lin,
                            LassoCV())

rf_pip = make_pipeline(processor_tree,
                        RandomForestRegressor(random_state=42))

ridge_pip = make_pipeline(processor_lin,
                          RidgeCV())

gradient_pip = make_pipeline(processor_tree,
                            HistGradientBoostingRegressor(random_state=0))
#stacking_regressor = StackingRegressor(
#    estimators=estimators, final_estimator=RidgeCV()
#)
estimators = [('Random Forest', rf_pip),
              ('Lasso', lasso_pip),
              ('Gradient Boosting', gradient_pip)]
stacking_regressor = StackingRegressor(estimators = estimators,
                                       final_estimator=RidgeCV())

#from sklearn.linear_model import LogisticRegression

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

#ridge_pip.fit(X, y)
#lasso_pip.fit(X,y)
#import pdb;pdb.set_trace()

###############################################################################
# Load dataset
###############################################################################
#
# We will be using Ames housing dataset. The objective is to predict the price
# of the residential homes in Ames, Iowa. `This dataset`_ was first compiled by
# Dean De Cock and became better known after it was used for the `Kaggle
# challenge`_.
#
# Ames housing dataset is not part of the Sklearn and therefore we will fetch
# it from `OpenML`_. It consists of multiple features characterizing the
# houses. Some of those features are not numeric values and therefore we will
# first replace them with numbers. Next, we need to impute all the missing
# values (here we will use 'most_frequent' strategy).
#
# .. _`This dataset`: http://jse.amstat.org/v19n3/decock.pdf
# .. _`Kaggle challenge`:
# https://www.kaggle.com/c/house-prices-advanced-regression-techniques
# .. _`OpenML`: https://www.openml.org/d/42165


'''
    o_columns = X.columns[X.dtypes == 'O']
    for column in o_columns:
        string_value = X[column].unique()
        string_dict = dict(zip(string_value, range(len(string_value))))
        X = X.replace({column: string_dict})

    imp = SimpleImputer(strategy="most_frequent")
    X = imp.fit_transform(X)

    X = np.asarray(X)
    y = np.asarray(y)

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp.fit(X)

    return X, y
    '''


###############################################################################
# Now we can use Ames dataset to make the predictions. We check the performance
# of each individual predictor as well as the stack of the regressors.
# To speed up the calculations we will use only part of the dataset, however
# feel free to take all the data-points and see if the result changes.

import time
import numpy as np
from sklearn.model_selection import cross_validate, cross_val_predict




fig, axs = plt.subplots(2, 2, figsize=(9, 7))
axs = np.ravel(axs)

for ax, (name, est) in zip(axs, estimators + [('Stacking Regressor',
                                               stacking_regressor)]):
    start_time = time.time()

    print(name)
    score = cross_validate(est, X, y,
                           scoring=['r2', 'neg_mean_absolute_error'],
                           n_jobs=-1, verbose=0)
    elapsed_time = time.time() - start_time

    y_pred = cross_val_predict(est, X, y, n_jobs=-1, verbose=0)

    plot_regression_results(
        ax, y, y_pred,
        name,
        (r'$R^2={:.2f} \pm {:.2f}$' + '\n' + r'$MAE={:.2f} \pm {:.2f}$')
        .format(np.mean(score['test_r2']),
                np.std(score['test_r2']),
                -np.mean(score['test_neg_mean_absolute_error']),
                np.std(score['test_neg_mean_absolute_error'])),
        elapsed_time)


plt.suptitle('Single predictors versus stacked predictors')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.show()

###############################################################################
# The stacked regressor will combine the strengths of the different regressors.
# However, we also see that training the stacked regressor is much more
# computationally expensive.

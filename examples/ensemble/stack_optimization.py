"""
=============================================
Stack ensemble tuning on a regression problem
=============================================

This example shows the flexiblity provided when using `FeatureUnion`/`Pipeline`
APIs together with the stack's meta estimator.

We are also using a regression dataset, showing the stacking API is really
flexible.

"""
import numpy as np

from sklearn import datasets
from sklearn.ensemble import StackMetaEstimator
from sklearn.pipeline import make_union, make_pipeline
from sklearn.metrics import r2_score

# base estimators
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Ridge, Lasso

# dataset
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split, ParameterGrid

# parameters
RANDOM_SEED = 17583

X, y = load_diabetes(return_X_y=True)

Xtrain, Xtest, ytrain, ytest = train_test_split(X, y, random_state=RANDOM_SEED)


# Let's create some regressors using parameter grids
configs = [{'constructor': SVR,
            'parameter_grid': ParameterGrid({'C': [1, 1e3],
                                             'kernel': ['rbf', 'linear']})},
           {'constructor': Ridge,
            'parameter_grid': ParameterGrid({'alpha': [1., .1, .01],
                                             'normalize': [True, False],
                                             'random_state': [RANDOM_SEED]})},
           {'constructor': Lasso,
            'parameter_grid': ParameterGrid({'alpha': [1., .1, .01],
                                             'normalize': [True, False],
                                             'random_state': [RANDOM_SEED]})}]

base_regressors = []
for config in configs:
    base_regressors.extend([config['constructor'](**params)
                            for params in list(config['parameter_grid'])])
    print("Generated %d regressors of type %s"
          % (len(config['parameter_grid']), config['constructor'].__name__))

print "Total number of base regressors: %d" % len(base_regressors)

# For each regressor, starting from a model made of all regressors, let's try
# removing one by one and checking if the ensemble performs better. The models
# will be combined with a LinearRegression. We'll fit them beforehand, so we
# speed things up.

combiner = LinearRegression()

layer0_metas = [StackMetaEstimator(x) for x in base_regressors]
initial_layer0 = make_union(*layer0_metas)
initial_Xtrain2 = initial_layer0.fit_transform(Xtrain, ytrain)
initial_Xtest2 = initial_layer0.transform(Xtest)

# calculate baseline score to compare
combiner.fit(initial_Xtrain2, ytrain)

ypred = combiner.predict(initial_Xtest2)
baseline = r2_score(ytest, ypred)

print "Baseline is %.7f" % baseline

removed_metas = []  # we'll store blacklisted models here

for i in range(len(layer0_metas)):
    Xtrain2 = np.delete(initial_Xtrain2, i, axis=1)
    Xtest2 = np.delete(initial_Xtest2, i, axis=1)
    combiner.fit(Xtrain2, ytrain)
    ypred = combiner.predict(Xtest2)
    score = r2_score(ytest, ypred)
    if score > baseline:
        print "Removing model #%d improves score by %.7f" % (i, score-baseline)
        removed_metas.append(i)

layer0_metas = [x for i, x in enumerate(layer0_metas)
                if i not in removed_metas]


# build model with chosen metamodels
Xtrain2 = np.delete(initial_Xtrain2, removed_metas, axis=1)
combiner.fit(Xtrain2, ytrain)
layer0 = make_union(*layer0_metas)
final_model = make_pipeline(layer0, combiner)

ypred = final_model.predict(Xtest)
print "Final score: %.7f" % r2_score(ytest, ypred)

"""
======================================
Poisson regression and non-normal loss
======================================

This example illustrates the use of linear Poisson regression
on the French Motor Third-Party Liability Claims dataset [1] and compares
it with models learned with least squared error.

We start by defining a few helper functions for loading the data and
visualizing results.


.. [1]  A. Noll, R. Salzmann and M.V. Wuthrich, Case Study: French Motor
    Third-Party Liability Claims (November 8, 2018).
    `doi:10.2139/ssrn.3164764 <http://dx.doi.org/10.2139/ssrn.3164764>`_

"""
print(__doc__)

# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#          Roman Yurchak <rth.yurchak@gmail.com>
# License: BSD 3 clause
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.datasets import fetch_openml
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import PoissonRegressor, LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder
from sklearn.preprocessing import StandardScaler, KBinsDiscretizer
from sklearn.ensemble import GradientBoostingRegressor

from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import mean_poisson_deviance


def load_mtpl2(n_samples=100000):
    """Fetcher for French Motor Third-Party Liability Claims dataset

    Parameters
    ----------
    n_samples: int, default=100000
      number of samples to select (for faster run time).
    """

    # freMTPL2freq dataset from https://www.openml.org/d/41214
    df_freq = fetch_openml(data_id=41214, as_frame=True)['data']
    df_freq['IDpol'] = df_freq['IDpol'].astype(np.int)
    df_freq.set_index('IDpol', inplace=True)

    # freMTPL2sev dataset from https://www.openml.org/d/41215
    df_sev = fetch_openml(data_id=41215, as_frame=True)['data']

    # sum ClaimAmount over identical IDs
    df_sev = df_sev.groupby('IDpol').sum()

    df = df_freq.join(df_sev, how="left")
    df["ClaimAmount"].fillna(0, inplace=True)

    # unquote string fields
    for column_name in df.columns[df.dtypes.values == np.object]:
        df[column_name] = df[column_name].str.strip("'")
    return df.iloc[:n_samples]


##############################################################################
#
# 1. Loading datasets and pre-processing
# --------------------------------------
#
# We construct the freMTPL2 dataset by joining the  freMTPL2freq table,
# containing the number of claims (``ClaimNb``) with the freMTPL2sev table
# containing the claim amount (``ClaimAmount``) for the same user ids.

df = load_mtpl2(n_samples=100000)

# Note: filter out claims with zero amount, as the severity model
# requires a strictly positive target values.
df.loc[(df.ClaimAmount == 0) & (df.ClaimNb >= 1), "ClaimNb"] = 0

# correct for unreasonable observations (that might be data error)
df["ClaimNb"] = df["ClaimNb"].clip(upper=4)
df["Exposure"] = df["Exposure"].clip(upper=1)

column_trans = ColumnTransformer(
    [
        ("Veh_Driv_Age", KBinsDiscretizer(n_bins=10), ["VehAge", "DrivAge"]),
        (
            "Veh_Brand_Gas_Region",
            OneHotEncoder(),
            ["VehBrand", "VehPower", "VehGas", "Region", "Area"],
        ),
        ("BonusMalus", "passthrough", ["BonusMalus"]),
        (
            "Density_log",
            make_pipeline(
                FunctionTransformer(np.log, validate=False), StandardScaler()
            ),
            ["Density"],
        ),
    ],
    remainder="drop",
)
X = column_trans.fit_transform(df)

##############################################################################
#
# The number of claims (``ClaimNb``) is a positive integer that can be modeled
# as a Poisson distribution. It is then assumed to be the number of discrete
# events occurring with a constant rate in a given time interval
# (``Exposure``). Here we model the frequency ``y = ClaimNb / Exposure``,
# which is still a (scaled) Poisson distribution.
#
# A very important property of the Poisson distribution is its mean-variance
# relation: The variance is proportional to the mean.

df["Frequency"] = df.ClaimNb / df.Exposure

print(
   pd.cut(df.Frequency, [-1e-6, 1e-6, 1, 2, 3, 4, 5]).value_counts()
)

##############################################################################
#
# It worth noting that 96 % of users have 0 claims, and if we were to convert
# this problem into a binary classification task, it would be significantly
# imbalanced.
#
# To evaluate the pertinence of the used metrics, we will consider as a
# baseline an estimator that returns 0 for any input.

df_train, df_test, X_train, X_test = train_test_split(df, X, random_state=2)


eps = 1e-5
print("MSE: %.3f" % mean_squared_error(
        df_test.Frequency.values, np.zeros(len(df_test)),
        df_test.Exposure.values))
print("MAE: %.3f" % mean_absolute_error(
        df_test.Frequency.values, np.zeros(len(df_test)),
        df_test.Exposure.values))
print("mean Poisson deviance: %.3f" % mean_poisson_deviance(
        df_test.Frequency.values, eps + np.zeros(len(df_test)),
        df_test.Exposure.values))


##############################################################################
#
# We start by modeling the target variable with the least squares linear
# regression model,


linregr = LinearRegression()
linregr.fit(X_train, df_train.Frequency, sample_weight=df_train.Exposure)

print("LinearRegression")
print("MSE: %.3f" % mean_squared_error(
          df_test.Frequency.values, linregr.predict(X_test),
          df_test.Exposure.values))
print("MSE: %.3f" % mean_absolute_error(
          df_test.Frequency.values, linregr.predict(X_test),
          df_test.Exposure.values))
print("mean Poisson deviance: %.3f" % mean_poisson_deviance(
        df_test.Frequency.values, np.fmax(linregr.predict(X_test), eps),
        df_test.Exposure.values))

##############################################################################
#
# The Poisson deviance cannot be computed because negative values are
# predicted by the model,

print('Number Negatives: %s / total: %s' % (
      (linregr.predict(X_test) < 0).sum(), X_test.shape[0]))

##############################################################################
#
# Next we fit the Poisson regressor on the target variable,

glm_freq = PoissonRegressor(alpha=0, max_iter=1000)
glm_freq.fit(X_train, df_train.Frequency, sample_weight=df_train.Exposure)

print("PoissonRegressor")
print("MSE: %.3f" % mean_squared_error(
        df_test.Frequency.values, glm_freq.predict(X_test),
        df_test.Exposure.values))
print("MAE: %.3f" % mean_absolute_error(
        df_test.Frequency.values, glm_freq.predict(X_test),
        df_test.Exposure.values))
print("mean Poisson deviance: %.3f" % mean_poisson_deviance(
        df_test.Frequency.values, glm_freq.predict(X_test),
        df_test.Exposure.values))

##############################################################################
#
# Finally we will consider a non linear model  with Gradient boosting that
# still minimizes the least square error.


gbr = GradientBoostingRegressor(max_depth=3)
gbr.fit(X_train, df_train.Frequency.values,
        sample_weight=df_train.Exposure.values)


print("GradientBoostingRegressor")
print("MSE: %.3f" % mean_squared_error(
      df_test.Frequency.values, gbr.predict(X_test), df_test.Exposure.values))
print("MAE: %.3f" % mean_absolute_error(
      df_test.Frequency.values, gbr.predict(X_test), df_test.Exposure.values))
print("mean Poisson deviance: %.3f" % mean_poisson_deviance(
      df_test.Frequency.values, np.fmax(gbr.predict(X_test), eps),
      df_test.Exposure.values))

##############################################################################
#
# In this example, although Gradient boosting minimizes the least square error,
# because of a higher predictive power it also results in a smaller Poisson
# deviance than the Poisson regression model.
#
# Evaluating models with a single train / test split is prone to numerical
# errors, we can verify that we would also get equivalent resuts with the
# cross-validation score.
#
# The difference between these models can also be visualized by comparing the
# histogram of observed target values with that of predicted values,


fig, axes = plt.subplots(1, 4, figsize=(16, 3))

df_train.Frequency.hist(bins=np.linspace(-1, 10, 50), ax=axes[0])

axes[0].set_title('Experimental data')

for idx, model in enumerate([linregr, glm_freq, gbr]):
    y_pred = model.predict(X_train)

    pd.Series(y_pred).hist(bins=np.linspace(-1, 8, 50), ax=axes[idx+1])
    axes[idx + 1].set_title(model.__class__.__name__)

for axi in axes:
    axi.set(
        yscale='log',
        xlabel="y (Frequency)"
    )

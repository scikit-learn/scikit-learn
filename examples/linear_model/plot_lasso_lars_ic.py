"""
==============================================
Lasso model selection via information criteria
==============================================

This example reproduces the example of Fig. 2 of [ZHT2007]_. A
:class:`~sklearn.linear_model.LassoLarsIC` estimator is fit on a
diabetes dataset and the AIC and the BIC criteria are used to select
the best model.

.. note::
    It is important to note that the optimization to find `alpha` with
    :class:`~sklearn.linear_model.LassoLarsIC` relies on the AIC or BIC
    criteria that are computed in-sample, thus on the training set directly.
    This approach differs from the cross-validation procedure. For a comparison
    of the two approaches, you can refer to the following example:
    :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_model_selection.py`.

.. topic:: References

    .. [ZHT2007] :arxiv:`Zou, Hui, Trevor Hastie, and Robert Tibshirani.
       "On the degrees of freedom of the lasso."
       The Annals of Statistics 35.5 (2007): 2173-2192.
       <0712.0881>`
"""

# Author: Alexandre Gramfort
#         Guillaume Lemaitre
# License: BSD 3 clause

# %%
# We will use the diabetes dataset.
from sklearn.datasets import load_diabetes

X, y = load_diabetes(return_X_y=True, as_frame=True)
n_samples = X.shape[0]
X.head()

# %%
# Scikit-learn provides an estimator called
# :class:`~sklearn.linear_model.LinearLarsIC` that uses either Akaike's
# information criterion (AIC) or the Bayesian information criterion (BIC) to
# select the best model. Before fitting
# this model, we will scale the dataset.
#
# In the following, we are going to fit two models to compare the values
# reported by AIC and BIC.
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoLarsIC
from sklearn.pipeline import make_pipeline

lasso_lars_ic = make_pipeline(StandardScaler(), LassoLarsIC(criterion="aic")).fit(X, y)


# %%
# To be in line with the definition in [ZHT2007]_, we need to rescale the
# AIC and the BIC. Indeed, Zou et al. are ignoring some constant terms
# compared to the original definition of AIC derived from the maximum
# log-likelihood of a linear model. You can refer to
# :ref:`mathematical detail section for the User Guide <lasso_lars_ic>`.
def zou_et_al_criterion_rescaling(criterion, n_samples, noise_variance):
    """Rescale the information criterion to follow the definition of Zou et al."""
    return criterion - n_samples * np.log(2 * np.pi * noise_variance) - n_samples


# %%
import numpy as np

aic_criterion = zou_et_al_criterion_rescaling(
    lasso_lars_ic[-1].criterion_,
    n_samples,
    lasso_lars_ic[-1].noise_variance_,
)

index_alpha_path_aic = np.flatnonzero(
    lasso_lars_ic[-1].alphas_ == lasso_lars_ic[-1].alpha_
)[0]

# %%
lasso_lars_ic.set_params(lassolarsic__criterion="bic").fit(X, y)

bic_criterion = zou_et_al_criterion_rescaling(
    lasso_lars_ic[-1].criterion_,
    n_samples,
    lasso_lars_ic[-1].noise_variance_,
)

index_alpha_path_bic = np.flatnonzero(
    lasso_lars_ic[-1].alphas_ == lasso_lars_ic[-1].alpha_
)[0]

# %%
# Now that we collected the AIC and BIC, we can as well check that the minima
# of both criteria happen at the same alpha. Then, we can simplify the
# following plot.
index_alpha_path_aic == index_alpha_path_bic

# %%
# Finally, we can plot the AIC and BIC criterion and the subsequent selected
# regularization parameter.
import matplotlib.pyplot as plt

plt.plot(aic_criterion, color="tab:blue", marker="o", label="AIC criterion")
plt.plot(bic_criterion, color="tab:orange", marker="o", label="BIC criterion")
plt.vlines(
    index_alpha_path_bic,
    aic_criterion.min(),
    aic_criterion.max(),
    color="black",
    linestyle="--",
    label="Selected alpha",
)
plt.legend()
plt.ylabel("Information criterion")
plt.xlabel("Lasso model sequence")
_ = plt.title("Lasso model selection via AIC and BIC")

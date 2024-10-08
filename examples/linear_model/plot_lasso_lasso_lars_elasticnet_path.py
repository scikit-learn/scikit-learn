'''
# Regularization Paths for Lasso, Lasso-LARS, and ElasticNet

In this example, we will explore and compare the regularization paths of three important 
linear models used for regularization: 
- :func:`~sklearn.linear_model.Lasso`
- :func:`~sklearn.linear_model.LassoLars`
- :func:`~sklearn.linear_model.ElasticNet`

## What is a Regularization Path?

Regularization path is a plot between model coefficients  and the regularization parameter (alpha).
For models like Lasso and ElasticNet, the path shows how coefficients 
are shrunk towards zero as regularization becomes stronger. This helps in feature selection 
and model interpretability.

We will dive into comparing Lasso vs ElasticNet, 
and Lasso vs Lasso-LARS, focusing on their regularization paths.
'''

import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.datasets import load_diabetes
from sklearn.linear_model import enet_path, lasso_path, lars_path

# Load the dataset
X, y = load_diabetes(return_X_y=True)
X /= X.std(axis=0)  # Standardize data (this ensures the features have mean 0 and variance 1)

'''
### 1. Lasso and ElasticNet: A Comparison of Regularization path

Lasso (Least Absolute Shrinkage and Selection Operator) uses L1 regularization, meaning it 
penalizes the absolute value of the coefficients. As a result, Lasso tends to produce sparse 
models, where some coefficients are exactly zero.

ElasticNet, on the other hand, is a combination of L1 and L2 regularization. The L2 penalty 
helps overcome some limitations of Lasso, particularly when features are highly correlated.

$$ \text{ElasticNet Loss} = \frac{1}{2n_{\text{samples}}} \|y - Xw\|^2_2 + \alpha \rho \|w\|_1 + \alpha (1 - \rho) \|w\|_2^2 $$

where $\rho$ is the mix ratio between Lasso (L1) and Ridge (L2) penalties.
'''

eps = 5e-3  # A smaller value of eps leads to a longer regularization path

# Compute the regularization path for Lasso
alphas_lasso, coefs_lasso, _ = lasso_path(X, y, eps=eps)

# Compute the regularization path for ElasticNet with l1_ratio=0.8
l1_ratio = 0.8  # This controls the mix between L1 and L2 (Lasso and Ridge)
alphas_enet, coefs_enet, _ = enet_path(X, y, eps=eps, l1_ratio=l1_ratio)

# Plot the paths for Lasso and ElasticNet
plt.figure(figsize=(10, 6))
colors = cycle(["b", "r", "g", "c", "k"])

for coef_l, coef_e, c in zip(coefs_lasso, coefs_enet, colors):
    plt.semilogx(alphas_lasso, coef_l, c=c)
    plt.semilogx(alphas_enet, coef_e, linestyle="--", c=c)

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("Lasso vs ElasticNet Regularization Path")
plt.legend(["Lasso", "ElasticNet (L1 ratio = 0.8)"], loc="upper left")
plt.axis("tight")
plt.show()

'''
We can see in the plot that as alpha increases (more regularization), both Lasso and ElasticNet drive coefficients towards 
zero. However, ElasticNet's combination of L1 and L2 regularization causes coefficients to 
shrink more smoothly as compared to Lasso. This allows ElasticNet to handle correlated features 
better, whereas Lasso might arbitrarily select one of the correlated features and set the rest to zero.
'''

'''
### 2. Lasso vs Lasso-LARS: Regularization Path

The main difference between Lasso and Lasso-LARS is the method it uses
to minimize loss.
Lasso uses cordinate descent to minimize the loss function which is an computationally
expensive method but Lasso-LARS (Least Angle Regression) is a more efficient 
algorithm. It finds the minimum solution by moving in a path of most correlated 
features. The regularization path of lasso and lasso-lars would similar, but 
lasso-lars would be must faster when there are many correlated features.

Let compute and compare their regularization paths.
'''

# Compute the regularization path for Lasso-LARS
alphas_lars, _, coefs_lars = lars_path(X, y, method="lasso")

# Plot the paths for Lasso and Lasso-LARS
plt.figure(figsize=(10, 6))
colors = cycle(["b", "r", "g", "c", "k"])

for coef_lasso, coef_lars, c in zip(coefs_lasso, coefs_lars, colors):
    plt.semilogx(alphas_lasso, coef_lasso, c=c)
    plt.semilogx(alphas_lars, coef_lars, linestyle="--", c=c)

plt.xlabel("alpha")
plt.ylabel("coefficients")
plt.title("Lasso vs Lasso-LARS Regularization Path")
plt.legend(["Lasso", "Lasso-LARS"], loc="upper right")
plt.axis("tight")
plt.show()

'''
As We can see the paths for Lasso and Lasso-LARS are close to each other. But lasso-LARS has a more direct 
path instead of a curve smooth path, that is because of its method of implementation.
Both methods set some coefficients to exactly zero, but the LARS algorithm moves in the 
direction of the strongest feature correlation, making it particularly suited for sparse models.
'''

'''
### 3. Positive Constraints

Both Lasso and ElasticNet can also enforce positive constraints on the coefficients by 
specifying `positive=True`.

Lets see how positive constraints impact the regularization paths for Lasso.
'''

alphas_positive_lasso, coefs_positive_lasso, _ = lasso_path(X, y, eps=eps, positive=True)
alphas_positive_enet, coefs_positive_enet, _ = enet_path(X, y, eps=eps, l1_ratio=l1_ratio, positive=True)
alphas_positive_lars, _, coefs_positive_lars = lars_path(X, y, method="lasso", positive=True)

# Plot all three subplots in one row
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

colors = cycle(["b", "r", "g", "c", "k"])

# First plot: Lasso vs Positive Lasso
for coef_lasso, coef_positive_lasso, c in zip(coefs_lasso, coefs_positive_lasso, colors):
    axes[0].semilogx(alphas_lasso, coef_lasso, c=c)
    axes[0].semilogx(alphas_positive_lasso, coef_positive_lasso, linestyle="--", c=c)

axes[0].set_xlabel("alpha")
axes[0].set_ylabel("coefficients")
axes[0].set_title("Lasso vs Positive Lasso")
axes[0].legend(["Lasso", "Positive Lasso"], loc="upper right")
axes[0].axis("tight")

# Second plot:  ElasticNet vs Positive ElasticNet
for coef_e, coef_pe, c in zip(coefs_enet, coefs_positive_enet, colors):
    axes[1].semilogx(alphas_enet, coef_e, c=c)
    axes[1].semilogx(alphas_positive_enet, coef_pe, linestyle="--", c=c)

axes[1].set_xlabel("alpha")
axes[1].set_ylabel("coefficients")
axes[1].set_title(" ElasticNet vs Positive ElasticNet")
axes[1].legend(["ElasticNet", "Positive ElasticNet (L1 ratio = 0.8)"], loc="upper left")
axes[1].axis("tight")

# Third plot: Lasso-LARS vs Positive Lasso-LARS
for coef_lars, coef_positive_lars, c in zip(coefs_lars, coefs_positive_lars, colors):
    axes[2].semilogx(alphas_lars, coef_lars, c=c)
    axes[2].semilogx(alphas_positive_lars, coef_positive_lars, linestyle="--", c=c)

axes[2].set_xlabel("alpha")
axes[2].set_ylabel("coefficients")
axes[2].set_title("Lasso-LARS vs Positive Lasso-LARS")
axes[2].legend(["Lasso-LARS", "Positive Lasso-LARS"], loc="upper right")
axes[2].axis("tight")

# Display the plots
plt.tight_layout()
plt.show()

'''
When we enforce positive constraints on Lasso, the regularization path differs, as coefficients 
are restricted to positive values only. This constraint leads to a different path, particularly 
for coefficients that would have otherwise become negative.
'''

'''
## Conclusion:

This example illustrates how the choice of regularization method and solver impacts the 
regularization path. Lasso and ElasticNet differ in their penalties (L1 vs a mix of L1 and L2), 
while Lasso and Lasso-LARS differ in their solvers, with LARS being more efficient for 
high-dimensional problems. Additionally, positive constraints can lead to different paths, 
forcing non-negative coefficients in models like Lasso and ElasticNet.
'''

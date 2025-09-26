"""
============================================================
Univariate Feature Selection vs. SVM Weights Comparison
============================================================

This example demonstrates the effectiveness of univariate feature selection
in identifying informative features and improving a subsequent model.

We start with the Iris dataset, which has 4 informative features. A large
number of noisy, non-informative features are then added.

We compare three things:
1. The scores from a univariate feature selection method (F-test p-values).
2. The feature weights from an SVM trained on all features (informative + noisy).
3. The feature weights from an SVM trained only on the features selected by the
   univariate method.

The plot shows that univariate selection correctly identifies the original 4
features as most significant. While a standard SVM assigns non-zero weights
to many noisy features, using a selection step beforehand forces the SVM to
focus only on the important features, leading to a more robust and
interpretable model.
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets, svm
from sklearn.feature_selection import SelectFpr, f_classif
from sklearn.pipeline import make_pipeline

def main():
    """
    Main function to run the feature selection and SVM comparison.
    """
    print(__doc__)

    # --------------------------------------------------------------------------
    # 1. Data Preparation
    # --------------------------------------------------------------------------
    # Load the IRIS dataset, which has 4 informative features
    iris = datasets.load_iris()
    X, y = iris.data, iris.target
    n_original_features = X.shape[1]

    # Add noisy features that are not correlated with the target
    np.random.seed(42)  # for reproducibility
    E = np.random.normal(size=(len(X), 35))
    X_noisy = np.hstack((X, E))
    n_total_features = X_noisy.shape[1]

    # --------------------------------------------------------------------------
    # 2. Feature Scoring and Weight Calculation
    # --------------------------------------------------------------------------
    # (A) Get scores from Univariate Feature Selection (F-test)
    # We select features based on the false positive rate test (alpha=0.1)
    selector_fpr = SelectFpr(f_classif, alpha=0.1)
    selector_fpr.fit(X_noisy, y)
    # The scores are the negative log of the p-values
    univariate_scores = -np.log10(selector_fpr.pvalues_)
    univariate_scores /= univariate_scores.max()

    # (B) Get weights from an SVM trained on ALL features
    clf = svm.SVC(kernel='linear')
    clf.fit(X_noisy, y)
    # Calculate weights as the squared sum of coefficients
    svm_weights_all_features = (clf.coef_**2).sum(axis=0)
    svm_weights_all_features /= svm_weights_all_features.max()

    # (C) Use a pipeline to perform selection and then train an SVM
    # This is the recommended scikit-learn approach
    clf_pipeline = make_pipeline(
        SelectFpr(f_classif, alpha=0.1),
        svm.SVC(kernel='linear')
    )
    clf_pipeline.fit(X_noisy, y)
    # Get the SVM model from the pipeline
    svm_model_after_selection = clf_pipeline.named_steps['svc']
    # Get the feature selector from the pipeline
    selector_in_pipeline = clf_pipeline.named_steps['selectfpr']
    # Get the weights from the SVM trained on selected features
    svm_weights_after_selection = (svm_model_after_selection.coef_**2).sum(axis=0)
    svm_weights_after_selection /= svm_weights_after_selection.max()
    
    # Map the pipeline SVM weights back to the original feature indices
    full_svm_weights_pipeline = np.zeros(n_total_features)
    selected_features_mask = selector_in_pipeline.get_support()
    full_svm_weights_pipeline[selected_features_mask] = svm_weights_after_selection

    # --------------------------------------------------------------------------
    # 3. Plotting the Results
    # --------------------------------------------------------------------------
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))

    x_indices = np.arange(n_total_features)
    bar_width = 0.25

    # Plot the scores and weights as bars
    ax.bar(
        x_indices - bar_width,
        univariate_scores,
        width=bar_width,
        label=r'Univariate Score ($-Log(p_{value})$)',
        color='darkgreen',
        edgecolor='black'
    )
    ax.bar(
        x_indices,
        svm_weights_all_features,
        width=bar_width,
        label='SVM Weight (All Features)',
        color='firebrick',
        edgecolor='black'
    )
    ax.bar(
        x_indices + bar_width,
        full_svm_weights_pipeline,
        width=bar_width,
        label='SVM Weight (After Selection)',
        color='cornflowerblue',
        edgecolor='black'
    )

    # Add a vertical line to separate original and noisy features
    ax.axvline(
        n_original_features - 0.5,
        ls='--',
        color='gray',
        label='Original/Noisy Feature Boundary'
    )

    # Final plot styling
    ax.set_title("Comparing Univariate Feature Selection with SVM Weights", fontsize=16)
    ax.set_xlabel('Feature Number', fontsize=12)
    ax.set_ylabel('Normalized Score / Weight', fontsize=12)
    ax.set_xticks(x_indices)
    ax.set_xticklabels(x_indices + 1) # Label features starting from 1
    ax.axis('tight')
    ax.legend(loc='upper right', fontsize=10)
    
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
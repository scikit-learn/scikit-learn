"""
================================================================================
Linear and Quadratic Discriminant Analysis with scikit-learn and custom datasets
================================================================================


This example provides quick workflow for Discriminant Analysis with 
classification report, and plots the decision boundary and confusion matrix learned
by :class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis` (LDA) and
:class:`~sklearn.discriminant_analysis.QuadraticDiscriminantAnalysis` (QDA) with 
scikit-learn and custom datasets. This example mainly focuses on class separation 
with LDA and QDA using train/test split, and visualisation of class separation and 
confusion matrix on test dataset.
"""


# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


import numpy as np
import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Union
from sklearn.pipeline import Pipeline
from sklearn.discriminant_analysis import *
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import RobustScaler
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_breast_cancer, make_blobs
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, classification_report


# A wrapper function for working with scikit-learn dataset to be used with the main function named discriminant_analysis_train_test
def fit_model_and_visualise_with_sklearn_dataset(
    X_train, 
    X_test, 
    y_train, 
    y_test, 
    model,
    dataset,
    model_name,
    figsize: Optional[tuple] = None,
    fig_title_1: Optional[str] = None
):
    
    model.fit(X_train, y_train)
            
    xx_min, xx_max= X_test[:,0].min() - 1, X_test[:,0].max() + 1
    yy_min, yy_max= X_test[:,1].min() - 1, X_test[:,1].max() + 1
            
    xx, yy= np.meshgrid(
        np.linspace(xx_min,xx_max,500),
        np.linspace(yy_min,yy_max,500)
    )
            
    grid= np.c_[xx.ravel(),yy.ravel()]
    Z= model.predict(grid)
    Z= Z.reshape(xx.shape)

    if figsize is None:
        figsize=(8,6)
    
    fig, ax= plt.subplots(1,2,figsize=figsize)
    ax[0].contourf(
        xx,
        yy,
        Z,
        alpha=0.25,
        cmap="coolwarm"
    )
    for i, class_name in enumerate(dataset.target_names):
        ax[0].scatter(
            X_test[y_test==i,0],
            X_test[y_test==i,1],
            label=class_name,
            edgecolors="black",
            alpha=0.75
            )
            
    ax[0].set_xlabel(dataset.feature_names[0])
    ax[0].set_ylabel(dataset.feature_names[1])
    ax[0].legend(loc="best")
    if model_name.lower() == "lda":
        ax[0].set_title("Linear Discriminant Analysis")
    else:
        ax[0].set_title("Quadratic Discriminant Analysis")
    
    preds= model.predict(X_test)
    cm= confusion_matrix(y_test,preds)
    disp= ConfusionMatrixDisplay(confusion_matrix=cm,display_labels=dataset.target_names)
            
    disp.plot(cmap="Blues",ax=ax[1],colorbar=False)
    if model_name.lower() == "lda":
        if fig_title_1 is None:
            fig.suptitle("LDA Classification")
        else:
            fig.suptitle(f"LDA: {fig_title_1}", fontsize=16,fontweight="bold")
        ax[1].set_title("LDA - Confusion Matrix")
        
    else:
        if fig_title_1 is None:
            fig.suptitle("QDA Classification")
        else:
            fig.suptitle(f"QDA: {fig_title_1}", fontsize=16,fontweight="bold")
        ax[1].set_title("QDA - Confusion Matrix")
    plt.tight_layout()
    plt.show()
            
    print("\nClassication Report:\n\n", classification_report(y_test, preds, target_names=dataset.target_names))
            
    print(f"\n\nTrain Accuracy: {model.score(X_train, y_train):.4f}")
    print(f"Test Accuracy: {model.score(X_test, y_test):.4f}")


# A wrapper function for working with custom dataset to be used with the main function named discriminant_analysis_train_test
def fit_model_and_visualise_with_custom_dataset(
    X,
    y,
    X_train, 
    X_test, 
    y_train, 
    y_test, 
    model,
    model_name,
    feature_names: Optional[Union[str, list, dict]] = None,
    target_names: Optional[Union[list,dict]] = None,
    figsize: Optional[tuple]= None,
    fig_title_1: Optional[str] = None
):
    X_train = np.asarray(X_train)
    X_test = np.asarray(X_test)
    y_train = np.asarray(y_train)
    y_test = np.asarray(y_test)

    model.fit(X_train, y_train)
    
    if feature_names is None:
        if hasattr(X,"columns"):
            feature_names= X.columns 
        else:
            feature_names= ["Feature 1", "Feature 2"]
    elif isinstance(feature_names, str):
        feature_names= list(X.columns)
    
    classes= np.unique(y)
    if isinstance(target_names, dict):
        target_names= [target_names[label] for label in classes]
    elif target_names is None:
        target_names= [str(label) for label in classes]
        
        
    xx_min, xx_max= X_test[:,0].min() - 1, X_test[:,0].max() + 1
    yy_min, yy_max= X_test[:,1].min() - 1, X_test[:,1].max() + 1
            
    xx, yy= np.meshgrid(
        np.linspace(xx_min,xx_max,500),
        np.linspace(yy_min,yy_max,500)
    )
            
    grid= np.c_[xx.ravel(),yy.ravel()]
    Z= model.predict(grid)
    Z= Z.reshape(xx.shape)

    if figsize is None:
        figsize=(8,6)
    
    fig, ax= plt.subplots(1,2,figsize=figsize)
    ax[0].contourf(
        xx,
        yy,
        Z,
        alpha=0.25,
        cmap="coolwarm"
    )
    for label, class_name in zip(classes,target_names):
        ax[0].scatter(
            X_test[y_test==label,0],
            X_test[y_test==label,1],
            label=class_name,
            edgecolors="black",
            alpha=0.75
            )
            
    ax[0].set_xlabel(feature_names[0])
    ax[0].set_ylabel(feature_names[1])
    ax[0].legend(loc="best")
    if model_name.lower() == "lda":
        ax[0].set_title("Linear Discriminant Analysis")
    else:
        ax[0].set_title("Quadratic Discriminant Analysis")
    
    preds= model.predict(X_test)
    cm= confusion_matrix(y_test,preds)
    disp= ConfusionMatrixDisplay(confusion_matrix=cm,display_labels=target_names)
            
    disp.plot(cmap="Blues",ax=ax[1],colorbar=False)
    if model_name.lower() == "lda":
        if fig_title_1 is None:
            fig.suptitle("LDA Classification")
        else:
            fig.suptitle(f"LDA: {fig_title_1}", fontsize=16,fontweight="bold")
        ax[1].set_title("LDA - Confusion Matrix")
        
    else:
        if fig_title_1 is None:
            fig.suptitle("QDA Classification")
        else:
            fig.suptitle(f"QDA: {fig_title_1}", fontsize=16,fontweight="bold")
        ax[1].set_title("QDA - Confusion Matrix")

    plt.tight_layout()
    plt.show()
            
    print("\nClassication Report:\n\n", classification_report(y_test, preds, labels=classes, target_names=target_names))
            
    print(f"\n\nTrain Accuracy: {model.score(X_train, y_train):.4f}")
    print(f"Test Accuracy: {model.score(X_test, y_test):.4f}")


def discriminant_analysis_train_test(
    dataset: object = None,      # Provide ready-made scikit-learn dataset such breast cancer, iris or any other dataset. In case dataset parameter is being passed then X and y should not be provided and vice versa.
    *,
    X:  Optional[Union[pd.DataFrame, pl.DataFrame, np.ndarray]] = None,  # Provide feature variables from your own custom dataset either generated via scikit-learn modules or your own personal dataset. If X is being passed then y too must be passed.
    y: Optional[Union[pd.Series, pl.Series, np.ndarray]] = None, # Provide target variable from your dataset. If y is being passed then X too must be passed.
    test_size: float = 0.25,  # Provide your desired test size, if not provided then it defaults to 0.25 .i.e 25% data is used for train and remaining data for training purposes
    model_name: str = "lda", # Specify your desired model between lda and qda, if not provided then it defaults to qda
    feature_names: Optional[Union[str, list]] = None, # Specify feature names which would be used in the plot displaying decision boundary and class separation, if using custom or your personal dataset .i.e when passing both X and y parameters.
    target_names: Optional[dict] = None, # Specify feature names which would be used as labels in the decision boundary and confusion matrix, if required.
    figsize: Optional[tuple] = (8,6), # Specify your preferred figsize in tuples or else it defaults to (8,6)
    fig_title_1: Optional[str] = None # Specify the title you would like to provide for the decision boundary plot.
):
    
    
    try:
        if model_name.lower() == "lda":
            model= LinearDiscriminantAnalysis()
        elif model_name.lower() == "qda":
            model= QuadraticDiscriminantAnalysis()
        else:
            raise ValueError("Value must be either qda or lda")
        
        if dataset is not None and X is None and y is None:
            X,y= dataset.data, dataset.target
            X= X[:,:2]
            X_train, X_test, y_train, y_test= train_test_split(X,y,test_size=test_size,shuffle=True,stratify=y, random_state=42)
            
            fit_model_and_visualise_with_sklearn_dataset(
                X_train=X_train, 
                X_test=X_test, 
                y_train=y_train, 
                y_test=y_test, 
                model=model, 
                dataset=dataset, 
                model_name=model_name,
                figsize=figsize,
                fig_title_1=fig_title_1
            )
        
        elif X is not None and y is not None:
            if isinstance(X, (pd.DataFrame, pl.DataFrame)):
                X = X.to_numpy()
                X = np.asarray(X)
                if X.ndim != 2 or X.shape[1] < 2:
                    raise ValueError("X must have at least two feature columns")
            elif isinstance(X, (pd.Series, pl.Series)):
                raise ValueError("X must be a 2D DataFrame or array, not a Series")
                
            X= X[:,:2]
            X_train, X_test, y_train, y_test= train_test_split(X,y,test_size=test_size,shuffle=True,stratify=y, random_state=42)
            fit_model_and_visualise_with_custom_dataset(
                X=X,
                y=y,
                X_train=X_train, 
                X_test=X_test, 
                y_train=y_train, 
                y_test=y_test, 
                model=model, 
                model_name=model_name,
                figsize=figsize,
                feature_names=feature_names,
                target_names=target_names,
                fig_title_1=fig_title_1
                )
        else:
            raise ValueError("Provide either sklearn dataset or both X and y if using your own dataset")
    except Exception as e:
        print(f"An Error Occurred: {e}")



# %%
# Using workflow with scikit-learn dataset
discriminant_analysis_train_test(load_breast_cancer(), fig_title_1="Breast Cancer Binary classification", model_name="qda")

# %%
# Generating custom dataset with scikit-learn's dataset

X, y = make_blobs(
    n_samples=500,
    centers=[(-1.5, 0), (1.5, 0)],
    cluster_std=1.5,  # Increase for more overlap, decrease for cleaner separation
    random_state=42
)

# Optional: Add small independent measurement noise
rng = np.random.default_rng(42)
X += rng.normal(loc=0, scale=0.25, size=X.shape)

# %%
# Using workflow with custom dataset dataset
discriminant_analysis_train_test(
    X=X,
    y=y, 
    target_names={0:"False",1:"True"}, # We can pass multiple labels in this parameter if it's a multiclass classification.
    feature_names=["1st Feature", "2nd Feature"], 
    fig_title_1="Binary Classification", 
    figsize=(12,10), 
    model_name="lda", 
    test_size=0.20
)


# %%
# The major difference between LDA and QDA is that LDA assumes that the covariance matrix of 
# each class is equal, while QDA estimates a covariance matrix per class. If the data requires
# linear separation then LDA should be chosen, whereas QDA should be used when class needs to be 
# separated based on their respective covariance matrix.

# While selecting the features for Discriminant Analysis and / or this workflow it's better to 
# make sure that those features aren't highly correlated to each other .i.e the problem of Multicollinearity
# or else they won't be able to provide the meaningful results. This workflow needs not more or less than 2 feature variables
# so it's important to have basic checks and selecting the features carefully. While dealing with categorical features
# we must first encode them using appropriate encoders. 

# It's also important to note that while LDA can also be used in supervised dimensionality reduction however, the main
# purpose of this workflow was quick and automated Discriminant Analysis with train/test split/
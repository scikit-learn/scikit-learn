"""
==========================
Pipeline Anova SVM
==========================

Simple usages of pipeline:
- ANOVA SVM-C
"""

from scikits.learn import svm
from scikits.learn.datasets import samples_generator
from scikits.learn.feature_selection.univariate_selection import SelectKBest,f_regression
from scikits.learn.pipeline import Pipeline

# import some data to play with
X, y = samples_generator.test_dataset_classif(k=5)


# ANOVA SVM-C
# 1) anova filter, take 5 best ranked features 
anova_filter = SelectKBest(f_regression, k=5)
# 2) svm
clf = svm.SVC(kernel='linear')

anova_svm = Pipeline([anova_filter], clf)
anova_svm.fit(X, y)
anova_svm.predict(X)


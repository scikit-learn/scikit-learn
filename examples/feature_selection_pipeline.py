"""
==========================
Pipeline Anova SVM
==========================

Simple usage of Pipeline made ANOVA | SVM-C
"""
#import numpy as np
#import pylab as pl
from  scikits.learn import svm, datasets
from  scikits.learn.datasets import samples_generator
import scikits.learn.feature_selection.new_univariate_selection as ufs
import scikits.learn.feature_selection.univ_scores as ufs_scores
from scikits.learn import pipeline



# import some data to play with
X,y = samples_generator.test_dataset_classif(k=5)

# 1) transformer 
univariate_filter = ufs.UnivariateFilter(ufs.SelectKBest(k=5),
    ufs_scores.f_regression)
# 2) predictor
clf = svm.SVC(kernel='linear')

pipe = pipeline.Pipeline([univariate_filter],clf)
pipe.fit(X,y)
pipe.predict(X)



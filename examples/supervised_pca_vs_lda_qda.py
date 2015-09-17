#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=====================
Supervised PCA against LDA and QDA
=====================

A comparison of supervised PCA against LDA and QDA.
"""
print(__doc__)


# Code source: Gaël Varoquaux
#              Andreas Müller
# Modified for documentation by Jaques Grobler
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_classification
from sklearn.lda import LDA
from sklearn.qda import QDA
from supervised_pca import SupervisedPCAClassifier

total_range=100
performances={}

names = ["LDA", "QDA","SuperPCA thres=0","SuperPCA thres=0.3","SuperPCA thres=0.7"]
ncomponents={names[2]:[],names[3]:[],names[4]:[]}

classifiers = [
    LDA(),
    QDA(),
    SupervisedPCAClassifier(threshold=0),
    SupervisedPCAClassifier(threshold=0.3),
    SupervisedPCAClassifier(threshold=0.7)
    ]  

for name in names:
    performances[name]=[]

    # iterate over classifiers

    
for i in range(1,total_range):
    X, y = make_classification(n_features=i*10, n_redundant=i*5, n_informative=i,
                               random_state=1, n_clusters_per_class=1)
    

    X = StandardScaler().fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.4)
                      
        
    for name, clf in zip(names, classifiers):
        
        clf.fit(X_train, y_train)
        pred=clf.predict(X_test)
        score = sum(y_test==pred)/(1.0*len(pred))
       
        performances[name].append(score)
        try:
            ncomponents[name].append(clf.get_n_components())
        except:
            pass
    
 
x=[k*10 for k in range(1,total_range)]

plt.figure()
plt.subplot(311)
plt.title("Score against number of features")  
for name in names:    
    plt.plot(x,performances[name])

plt.legend(labels=names,loc="best")

plt.subplot(312)

plt.title("Score boxplots for each classifier")

dummy=[]
for name in names:
    dummy.append(performances[name])
plt.boxplot(dummy,labels=names)

plt.subplot(313)

plt.title("Number of components against features")
plotcomponentsSPCA0=plt.plot(x,ncomponents[names[2]])
plotComponentsSPCA01=plt.plot(x,ncomponents[names[3]])
plotComponentsSPCA07=plt.plot(x,ncomponents[names[4]])
plt.legend([names[2],names[3],names[4]],loc="best")
        
        

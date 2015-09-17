#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=====================
Supervised PCA against LDA and QDA
=====================

A comparison of supervised PCA against elastic net.

The data are artificially created and the variables can be described by a 
lower dimensional space.

Supervised PCA shows better results until about 100 features. The performance of all models
is similar after that point. A supervisedPCA model with a threshold=0.7 works with a smaller
number of components, while having similar performance to the rest of the models. This can be particularly
useful in situations where interpretability is important.
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
from sklearn.datasets import make_regression
from supervised_pca import SupervisedPCARegressor
from sklearn.linear_model import ElasticNet

total_range=50
performances={}

names = ["ElasticNet","SuperPCA thres=0","SuperPCA thres=0.1","SuperPCA thres=0.7"]
ncomponents={names[1]:[],names[2]:[],names[3]:[]}



for name in names:
    performances[name]=[]

    # iterate over classifiers

    
for i in range(1,total_range):
    print(i)
    
    classifiers = [
        ElasticNet(),
        SupervisedPCARegressor(threshold=0),
        SupervisedPCARegressor(threshold=0.1),
        SupervisedPCARegressor(threshold=0.7)
        ]      
    
    X, y = make_regression(n_features=i*5, n_informative=i*4,
                               random_state=1,effective_rank=i)
    

    X = StandardScaler().fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.4)
                      
        
    for name, clf in zip(names, classifiers):
        
        clf.fit(X_train, y_train)
        pred=clf.predict(X_test)
        score = np.mean(abs(y_test-pred))
       
        performances[name].append(score)
        try:
            ncomponents[name].append(clf.get_n_components())
        except:
            pass
    
 
x=[k*5 for k in range(1,total_range)]

plt.figure()
plt.subplot(311)
plt.title("MAE against number of features")  
for name in names:    
    plt.plot(x,performances[name])

plt.legend(labels=names,loc="best")

plt.subplot(312)

plt.title("MAE boxplots for each classifier")

dummy=[]
for name in names:
    dummy.append(performances[name])
plt.boxplot(dummy,labels=names)

plt.subplot(313)

plt.title("Number of components against features")
plotcomponentsSPCA0=plt.plot(x,ncomponents[names[1]])
plotComponentsSPCA01=plt.plot(x,ncomponents[names[2]])
plotComponentsSPCA07=plt.plot(x,ncomponents[names[3]])
plt.legend([names[1],names[2],names[3]],loc="best")
        
        

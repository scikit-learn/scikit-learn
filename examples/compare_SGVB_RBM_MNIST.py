"""Example SGVB 
"""
# Authors:
# Otto Fabius
# Joost van Amersfoort

"""Demonstration comparing SGVB to RBM for feature extraction and classification"""

import numpy as np
from sklearn import neural_network, linear_model
from sklearn.datasets import fetch_mldata

print "Loading data"
mnist = fetch_mldata('MNIST original')

#Normalize to 0-1
mnist["data"] = mnist["data"]/255.

#Take the standard MNIST training set size for the RBM and auto-encoder
x_train, t_train = mnist["data"][:50000], mnist["target"][:50000]

#Reduced training sample for logistic regression to show effect of RBM and the variational auto-encoder on this problem.
data_train, data_test = mnist["data"][50000:51000], mnist["data"][51000:52000]
t_train, t_test = mnist["target"][50000:51000], mnist["target"][51000:52000]

print 'extracting features with SGVB auto-encoder, default is 40 iterations...'
encoder = neural_network.SGVB(verbose=True)
encoder.fit(x_train)
SGVB_train_features = encoder.transform(data_train)
SGVB_test_features  = encoder.transform(data_test)
print 'done'

print 'extracting features with RBM...'
n_components=200
learning_rate = 0.03
batch_size=100
rbm = neural_network.BernoulliRBM(n_components,learning_rate,batch_size, verbose=True)
rbm.fit(x_train)
RBM_train_features = rbm.transform(data_train)
RBM_test_features  = rbm.transform(data_test)
print 'done'

print 'performing logistic regression on raw data...'
LogReg_raw = linear_model.LogisticRegression()
LogReg_raw.fit(data_train,t_train)
raw_score = LogReg_raw.score(data_test, t_test)
print 'Test score on raw data = ', raw_score

print 'performing logistic regression on RBM features...'
LogReg_RBM = linear_model.LogisticRegression()
LogReg_RBM.fit(RBM_train_features,t_train)
RBM_score = LogReg_RBM.score(RBM_test_features, t_test)
print 'Test score on RBM features = ', RBM_score

print 'performing logistic regression on SGVB features...'
LogReg_SGVB = linear_model.LogisticRegression()
LogReg_SGVB.fit(SGVB_train_features,t_train)
SGVB_score = LogReg_SGVB.score(SGVB_test_features, t_test)
print 'Test score on SGVB features = ', SGVB_score

# libraries
import numpy as np

# scikit-learn base libraries
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin

# scikit-learn modules
from sklearn.cross_validation import cross_val_predict
from sklearn.ensemble import VotingClassifier


class StackingClassifier(BaseEstimator,ClassifierMixin):
    '''
    stacking ensemble classifier based on scikit-learn
    '''
    def __init__(self,stage_one_clfs,stage_two_clfs,weights=None, n_runs=10):
        '''
        
        weights: weights of the stage_two_clfs
        n_runs: train stage_two_clfs n_runs times and average them (only for probabilistic output)
        '''
        self.stage_one_clfs = stage_one_clfs
        self.stage_two_clfs = stage_two_clfs
        if weights == None:
            self.weights = [1] * len(stage_two_clfs)
        else:
            self.weights = weights
        self.n_runs = n_runs
        
    def fit(self,X,y):
        '''
        fit the model
        '''
        self.__X = X
        self.__y = y
        
        # fit the first stage models
        for clf in self.stage_one_clfs:
            y_pred = cross_val_predict(clf, X, y, cv=5, n_jobs=1)
            clf.fit(X,y)
            y_pred  = np.reshape(y_pred,(len(y_pred),1))
            self.__X = np.hstack((self.__X,y_pred)) 
        
        # fit the second stage models
        for clf in self.stage_two_clfs:
            clf.fit(self.__X,self.__y)        
    
    def predict(self,X_test):
        '''
        predict the class for each sample
        '''
        self.__X_test = X_test
        
        # first stage
        for clf in self.stage_one_clfs:
            y_pred = clf.predict(X_test)
            y_pred  = np.reshape(y_pred,(len(y_pred),1))
            self.__X_test = np.hstack((self.__X_test,y_pred)) 
            
        # second stage
        preds = []
        for i in range(self.n_runs):
            j = 0
            for clf in self.stage_two_clfs:
                y_pred = clf.predict(self.__X_test)  
                preds.append(self.weights[j] * y_pred)
                j += 1
        # average predictions
        y_final = preds.pop(0)
        for pred in preds:
            y_final += pred
        y_out = y_final/(np.array(self.weights).sum() * self.n_runs)
        
        return y_out
        
    
    def predict_proba(self,X_test):
        '''
        predict the probability for each class for each sample
        '''
        self.__X_test = X_test
        
        # first stage
        for clf in self.stage_one_clfs:
            y_pred = clf.predict(X_test)
            y_pred  = np.reshape(y_pred,(len(y_pred),1))
            self.__X_test = np.hstack((self.__X_test,y_pred)) 
            
        # second stage
        majority_voting = VotingClassifier(estimators=self.stage_two_clfs, voting="hard", weights=self.weights)
        y_out = majority_voting.predict(self.__X_test)
        
        return y_out
    
    def run_gridsearch(self,X,y,params):
        '''
        
        '''
        print("Not implemented yet")
        
        
class StackingRegressor(BaseEstimator,RegressorMixin):
    '''
    stacking ensemble regressor based on scikit-learn
    '''
    def __init__(self,stage_one_clfs,stage_two_clfs,weights=None, n_runs=10):
        '''
        
        weights: weights of the stage_two_clfs
        n_runs: train stage_two_clfs n_runs times and average them
        '''
        self.stage_one_clfs = stage_one_clfs
        self.stage_two_clfs = stage_two_clfs
        if weights == None:
            self.weights = [1] * len(stage_two_clfs)
        else:
            self.weights = weights
        self.n_runs = n_runs
        
    def fit(self,X,y):
        '''
        fit the model
        '''
        self.__X = X
        self.__y = y
        
        # fit the first stage models
        for clf in self.stage_one_clfs:
            y_pred = cross_val_predict(clf, X, y, cv=5, n_jobs=1)
            clf.fit(X,y)
            y_pred  = np.reshape(y_pred,(len(y_pred),1))
            self.__X = np.hstack((self.__X,y_pred)) 
        
        # fit the second stage models
        for clf in self.stage_two_clfs:
            clf.fit(self.__X,self.__y)        
    
    def predict(self,X_test):
        '''
        Predict the value for each sample
        '''
        self.__X_test = X_test
        
        # first stage
        for clf in self.stage_one_clfs:
            y_pred = clf.predict(X_test)
            y_pred  = np.reshape(y_pred,(len(y_pred),1))
            self.X_test = np.hstack((self.__X_test,y_pred)) 
            
        # second stage
        preds = []
        for i in range(self.n_runs):
            j = 0
            for clf in self.stage_two_clfs:
                y_pred = clf.predict(self.__X_test)  
                preds.append(self.weights[j] * y_pred)
                j += 1
        # average predictions
        y_final = preds.pop(0)
        for pred in preds:
            y_final += pred
        y_final = y_final/(np.array(self.weights).sum() * self.n_runs)
        
        return y_final
    
    def run_gridsearch(self,X,y,params):
        '''
        
        '''
        print("Not implemented yet")
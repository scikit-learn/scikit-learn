"""
Operation Suggester 
--------------------

This module focuses on the people who are new to the data science field. 
This module will suggest the user about the operations that they can perform
with respect to dataset. You must call the package method to execute the entire module.

Example : 

Df --> Dataframe --> type = Data frame
target_feature --> target label to predict --> type = String

object = Operation_Suggester(Df , target_feature).package() --> Code to execute the entire module. 



Developed by : Santosh Saxena
Date : 18/8/2021
"""


# Module Declaration
#-----------------------------------------------------------------------#
import numpy as np
from sklearn.svm import SVC , SVR
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression , LogisticRegression
from sklearn.tree import DecisionTreeClassifier , DecisionTreeRegressor
#-----------------------------------------------------------------------#


class Operation_Suggester:

    """
    This is an operation suggester module. This module will give suggesstion when user 
    will give dataframe as an input.

    """

    def __init__(self, df , target_feature):
        """
        This is an initializaion function. This methdod is used to load
        the necessary inputs.

        Input : 
        1) Dataframe as df with dataframe data type.
        2) Target Feature as target_feature with string data type.

        Output : 
        NaN

        """
        try:
            self.df = df
            self.y = target_feature
            self.df = self.df[self.df[self.y].isnull() == False]
            self.y_data_type = self.df[self.y].dtypes
            
        except Exception as e:
            print("pass")

    def missing_values(self):

        """
        Missing values method is supposed to process the data with respect to NUll values.

        If the amount of missing values in a particular column is more than 50 % then that column will be
        removed.

        If the data type is object then the missing values will be removed
        If the data type is not object then the missing values is filled with mean

        Input : This mehtod works on class variables hence NaN.
        Output : Information regarding missing values and suggesstions.

        """

        try:
            print("\n\n")
            print("=="*5 , end = "")
            print("\nMissing Values")
            print("=="*5 , end ="")
            print("\n\n")
            print("These features is having missing values : ",self.df.columns[self.df.isnull().sum() != 0])
            print("\n")
            print("The feature in which missing values are too high are : ", self.df.columns[self.df.isnull().sum() > 0.5*len(self.df)])
            print(self.df.columns[self.df.isnull().sum() > 0.5*len(self.df)], " features is not considered")
            print("We suggest that you should eliminate these features", self.df.columns[self.df.isnull().sum() > 0.5*len(self.df)])
            
            self.df.drop(self.df.columns[self.df.isnull().sum() > 0.5*len(self.df)], axis = 1, inplace = True) 
            print("\nProcessing others missing features\n")
            for i in self.df.columns[self.df.isnull().sum() != 0]:
                if(self.df[i].dtypes == "object"):
                    print(i, " Feature should be dropped because of missing values")
                    self.df = self.df[self.df[i].isnull() == False]
            
                if(self.df[i].dtypes != "object"):
                    print(i, " Feature, missing values should be replaced with mean")
                    self.df[i]= self.df[i].fillna(self.df[i].mean())
            print("=="*5 , end = "")
        
        except Exception as e:
            print(e)

    def label_encoder(self):

        """
        This module will encode the labels 
        Input: Nan
        Ouptut : Information regarding categorical features and suggesstions. 

        working : class variables will be encoded. 
        
        """
        try:
            print("\n\n")
            print("=="*5 , end = "")
            print("\nCategorical Feature processing")
            print("=="*5 , end = "")
            print("\n\n")
            for i in self.df.columns[self.df.dtypes == "object"]:
                print(str(i) + " should be encoded")
                self.encoder = LabelEncoder()
                self.df[i] = self.encoder.fit_transform(self.df[i])
            print("=="*5 , end = "")
            print()
        except Exception as e:
            print(e)


    def pca(self):

        """
        This method is supposed to reduce the dimension. 
        because of this method multi-collinearity will not occur. 

        Input : Nan
        Output : Information and suggesstion regarding dimensions.

        """

        try:
            print("\n\n")
            print("=="*5 , end = "")
            print("\nDimensionalty Reduction")
            print("=="*5 , end = "")
            print("\n\n")
            self.X = self.df.drop(self.y , axis = 1)
            self.y = self.df[self.y]
            print("The shape of X is : " , self.X.shape)
            pca = PCA()
            pca.fit(self.X)
            self.no_of_dimensions = 0
            for i in np.cumsum(pca.explained_variance_ratio_):
                if(i < 0.99):
                    self.no_of_dimensions += 1

            print("The dimensions of the dataset should be reduced to : ",self.no_of_dimensions)

            pca = PCA(self.no_of_dimensions)
            self.X = pca.fit_transform(self.X)
            print("The shape of X after dimensionality reduction : ", self.X.shape)
            print("=="*5 , end = "")
            print()

        except Exception as e:
            print(e)

    def data_preprocessing(self):

        """
        This methods is combining all the preprocessing mehtod to make pipeline.

        Input : Nan
        Output : Nan

        """

        try:
            self.missing_values()
            self.label_encoder()
            self.pca()
        except Exception as e:
            print(e)


    def find_categorical_algo(self):

        """
        This method will suggest the base algorithm that can be best with respect to your model.
        Input : Nan
        Output : suggesstion of algorithm 
    
        """

        try:
            print("\n\n")
            print("=="*5, end = "")
            print("\nSuggesting Categorical algorithm")
            print("=="*5 , end = "")
            print("\n\n")
            self.accuracy = []
            self.models = [LogisticRegression() , SVC() , DecisionTreeClassifier()]
            x_train, x_test, y_train, y_test = train_test_split(self.X , self.y)
            for i in self.models:
                model = i
                model.fit(x_train, y_train)
                self.accuracy.append(model.score(x_test , y_test))
        
            for i in range(len(self.accuracy)):
                if(self.accuracy[i] == max(self.accuracy)):
                    print(str(self.models[i]) + " is a good model for this dataset")
                    print("You can add ensemble techniques with the base model as " + str(self.models[i]))
            print("=="*5 , end = "")
            print()

        except Exception as e:
            print(e)

    def find_regression_algo(self):
        """
        This method will suggest the base algorithm that can be best with respect to your model.
        Input : Nan
        Output : suggesstion of algorithm 
    
        """
        try:
            print("\n\n")
            print("=="*5, end = "")
            print("\nSuggesting Regression algorithm")
            print("=="*5 , end = "")
            print("\n\n")
            self.accuracy = []
            self.models = [LinearRegression() , SVR() , DecisionTreeRegressor()]
            x_train, x_test, y_train, y_test = train_test_split(self.X , self.y)
            for i in self.models:
                model = i
                model.fit(x_train, y_train)
                self.accuracy.append(model.score(x_test , y_test))
        
            for i in range(len(self.accuracy)):
                if(self.accuracy[i] == max(self.accuracy)):
                    print(str(self.models[i]) + " is a good model for this dataset")
                    print("You can add ensemble techniques with the base model as " + str(self.models[i]))
        
            print("=="*5 , end = "")
            print()
        except Exception as e:
            print(e)

    def find_best_algo(self):
        """
        This method is combining all the package related to algorithm selection. 

        Input : Nan
        Output : Nan

        """

        try:
            if(self.y_data_type == "object"):
                self.find_categorical_algo()
            elif(self.y_data_type != "object"):
                self.find_regression_algo()
            else:
                raise("Some issue with find best algo in data type")
        
        except Exception as e:
            print(e)

    def package(self):
        """
        This is the package function . User needs to call this function to execute each and every task. 

        Input : NaN
        Output : Nan.

        """
        try:
            self.data_preprocessing()        
            self.find_best_algo()
        except Exception as e:
            print(e)


# Similarly we can further add ensemble learning and Deep learning suggesstions. 
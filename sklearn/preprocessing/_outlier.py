'''
Author 1: Sahib Nanda <sabbykabby12@gmail.com>
LinkedIn: https://www.linkedin.com/in/sahib-nanda-44b2bb264/
GitHub: https://github.com/thesahibnanda
'''

#Importing Required Libraries
import warnings

import pandas as pd
import numpy as np
from scipy import stats

__all__ = ["FixOutliers"]

class FixOutliers:
    """ 
    Used To Remove Or Impute Outliers According To The Approaches (To Detect Outliers) Defined By The User
    Remove Means Deleting The Datapoints Or Rows That Have Outliers
    Impute Means Replacing Values Of Outliers With Some Other Value Like Mean, Median, Mode.

    Approaches To Detect Outliers: 

    Approach 1 - On The Basis Of Z Score
    Z = (Datapoint - Mean) / Standard Deviation
    Find The Data Points Whose Z-Scores Are Higher Than The Specified Threshold Value (Default Threshold Value Is 3)

    Approach 2 - On The Basis Of InterQuartile Range
    Q1 Is The Point That Accounts 25 Percentile Of All Data
    Q2 Is The Point That Accounts 50 Percentile Of All Data
    Q3 Is The Point That Accounts 75 Percentile Of All Data
    Interquartile Range = Q3 - Q1
    Minimum of Quartile = Q1 - (1.5 X Interquartile Range)
    Maximum of Quartile = Q3 + (1.5 X Interquartile Range)
    Data Points Not In Range Of Minimum of Quartile And Maximum of Quartile Are Outliers 

    Approach 3 - On The Basis Of Both Approaches Mentioned Above
    First Outliers Will Detected And Taken Care Of Using 1st Approach Then 2nd Approach

    Parameters
    ----------
    approach: 'str', Optional (Default = 'interquartile_range')
    User Specifies The Approach He/She Wants To Be Executed, Either It Can Be interquartile_range ,z_score or every.
    """

    #Parameters And Their Description Is Provided In Respective Function

    def __init__(
        self,
        approach='interquartile_range'
        ):

        if approach not in ['z_score', 'interquartile_range', 'every']:
            raise ValueError('Invalid Approach Parameter')
        
        self.approach = approach
    
    def fit_transform(self, df, col_name, treatment = 'impute', imputation = 'mean', consider_outliers = False, threshold = 3):
        """
        Parameters
        ----------
        df: 'pandas DataFrame'
        DataFrame That Needs To Undergo Outlier Treatment

        col_name: 'str'
        Name Of Column Of The DataFrame That Needs To Undergo Outlier Treatment

        treatment: 'str', Optional (Default = 'impute')
        Specifies What Treatment Wil Outliers Undergo (Possible Values = ['impute','remove'])

        imputation: 'str', Optional (Default = 'mean')
        Specifies Which Value Will Be Used For Imputation (Possible Values = ['mean', 'median', 'mode', 'std'])
        Only Useful If treatment = 'impute'

        consider_outliers: 'bool', Optional (Default = False)
        Specifies Should Machine Consider Outliers While Computing Imputation Value
        Only Useful If treatment = 'impute'

        threshold: 'int', Optional (Default = 3)
        Specifies Threshold For Z-Score Outlier Computation
        Only Useful If self.approach = 'z_score'

        Return
        ------
        df: 'pandas DataFrame'
        DataFrame After Outlier Treatment
        """
        if self.approach == 'interquartile_range' and threshold != 3:
            raise AttributeError("InterQuartile Range Does Not Require Threshold")
        
        if treatment not in ['remove', 'impute']:
            raise ValueError("Unknown Value Passed As treatment")
        
        if imputation not in ['mean', 'median', 'mode', 'std']:
            raise ValueError("Unknown Value Passed As imputation")
        
        if consider_outliers not in [True, False]:
            raise ValueError("consider_outliers Needs To <'bool'>")

        if treatment == 'remove' and (imputation != 'mean' or consider_outliers != False):
            raise AttributeError("No Need Of imputation or consider_outliers when treatment = 'remove' ")
        
        #Now Functions Defined Out Of This Class Will Be Called

        if self.approach == 'z_score':

            if treatment == 'remove':
                output = Outliers_Removal_Z_Score(df, col_name, threshold) 
                #Removal Of Outliers From Z-Score
                return output
            
            elif treatment == 'impute':

                output = Impute_Z_Score(df, col_name, method=imputation, consider_outliers=consider_outliers, threshold=threshold)
                #Imputation Of Outliers From Z-Score
                return output
            
        elif self.approach == 'interquartile_range':

            if treatment == 'remove':
                output = Outliers_Removal_IQR(df, col_name)
                #Removal Of Outliers From InterQuartile Range
                return output
            
            elif treatment == 'impute':
                output = Impute_IQR(df, col_name, method=imputation, consider_outliers=consider_outliers)
                #Imputation Of Outliers From InterQuartile Range
                return output
        
        elif self.approach == 'every':

            if treatment == 'remove':

                df = Outliers_Removal_Z_Score(df, col_name, threshold) #Z-Score Removal
                output = Outliers_Removal_IQR(df, col_name) #InterQuartile Range Removal

                return output
            
            elif treatment == 'impute':
                
                #Z-Score Imputation
                df = Impute_Z_Score(df, col_name, method=imputation, consider_outliers=consider_outliers, threshold=threshold)

                #InterQuartile Range Imputation
                output = Impute_IQR(df, col_name, method=imputation, consider_outliers=consider_outliers)

                return output



            


#Below Defined Functions Will Be Used In Class FixOutliers

#Outlier Removal From Z Score
def Outliers_Removal_Z_Score(df, col_name, threshold):
    """
        Parameters
        ----------
        df: 'pandas DataFrame'
        DataFrame From Which Outliers Need To Be Detected Is Passed

        col_name: 'str'
        Name Of The Column That Needs To Undergo Outlier Detection And Removal Is Passed

        threshold: 'int'
        Threshold Value For Outlier Computation Is Passed

        Returns
        -------
        df: 'pandas DataFrame'
        DataFrame After Removal Of Outliers Is Returned
    """

    #Calculate Z - Score
    z_scores = np.abs((df[col_name] - df[col_name].mean()) / df[col_name].std())

    #Identify Outliers
    outliers = z_scores > threshold

    #Removal
    df = df[~outliers]

    return df


#Outlier Removal From InterQuartile Range
def Outliers_Removal_IQR(df, col_name):
    """
        Parameters
        ----------
        df: 'pandas DataFrame'
        DataFrame From Which Outliers Need To Be Detected Is Passed

        col_name: 'str'
        Name Of The Column That Needs To Undergo Outlier Detection And Removal Is Passed

        Returns
        -------
        df: 'pandas DataFrame'
        DataFrame After Removal Of Outliers Is Returned
    """

    #Calculate First(25 Percentile) And Third(75 Percentile) Quartile
    Q1 = df[col_name].quantile(0.25)
    Q3 = df[col_name].quantile(0.75)

    #Calculating InterQuartile Range
    IQR = Q3 - Q1

    #Defining Upper Bound And Lower Bound
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    #Removing Outliers
    df = df[(df[col_name] >= lower_bound) & (df[col_name] <= upper_bound)]

    return df


#Outliers Imputation From Z - Score
def Impute_Z_Score(df, col_name, method, consider_outliers, threshold):
    """
        Parameters
        ----------
        df: 'pandas DataFrame'
        DataFrame From Which Outliers Need To Be Detected Is Passed

        col_name: 'str'
        Name Of The Column That Needs To Undergo Outlier Detection And Removal Is Passed

        method: 'str'
        What Value Should Be Used For Imputation - Mean, Median, Mode Or Std. Deviation

        consider_outliers: 'bool'
        Should Outliers Be Considered While Calculating Imputation Values

        threshold: 'int'
        Threshold Value For Outlier Computation Is Passed

        Returns
        -------
        df: 'pandas DataFrame'
        DataFrame After Imputation
    """
    
    #Convert Method In Lower Case
    method = method.lower()
        
    #Finding Outliers In The DataFrame
    
    #Calculate Z - Score
    z_scores = np.abs((df[col_name] - df[col_name].mean()) / df[col_name].std())
    
    #Identify Outliers
    outliers = z_scores > threshold
    
    #Considering Outliers
    if consider_outliers:
        if method == 'mean':
            impute = df[col_name].mean()
            df.loc[outliers, col_name] = impute
        elif method == 'median':
            impute = df[col_name].median()
            df.loc[outliers, col_name] = impute
        elif method == 'mode':
            impute = stats.mode(df[col_name] , keepdims = True)[0][0]
            df.loc[outliers, col_name] = impute
        elif method == 'std':
            impute = df[col_name].std()
            df.loc[outliers, col_name] = impute
        else:
            raise ValueError('Invalid Imputation, Valid Options: Mean, Median, Mode Or Std')
    else:
        if method == 'mean':
            impute = df.loc[~outliers, col_name].mean()
            df.loc[outliers, col_name] = impute
        elif method == 'median':
            impute = df.loc[~outliers, col_name].median()
            df.loc[outliers, col_name] = impute
        elif method == 'mode':
            impute = stats.mode(df.loc[~outliers, col_name], keepdims = True)[0][0]
            df.loc[outliers, col_name] = impute
        elif method == 'std':
            impute = df.loc[~outliers, col_name].std()
            df.loc[outliers, col_name] = impute
        else:
            raise ValueError('Invalid Imputation, Valid Options: Mean, Median, Mode Or Std')
            
    return df


#Outliers Imputation From InterQuartile Range
def Impute_IQR(df, col_name, method, consider_outliers):
    """
        Parameters
        ----------
        df: 'pandas DataFrame'
        DataFrame From Which Outliers Need To Be Detected Is Passed

        col_name: 'str'
        Name Of The Column That Needs To Undergo Outlier Detection And Removal Is Passed

        method: 'str'
        What Value Should Be Used For Imputation - Mean, Median, Mode Or Std. Deviation

        consider_outliers: 'bool'
        Should Outliers Be Considered While Calculating Imputation Values

        Returns
        -------
        df: 'pandas DataFrame'
        DataFrame After Imputation
    """
    
    #Convert Method In Lower Case
    method = method.lower()
        
    #Finding Outliers In The DataFrame
    
    #Calculate IQR
    q1 = df[col_name].quantile(0.25)
    q3 = df[col_name].quantile(0.75)
    iqr = q3 - q1
    
    #Identify Outliers
    outliers = (df[col_name] < (q1 - 1.5*iqr)) | (df[col_name] > (q3 + 1.5*iqr))
    
    #Considering Outliers
    if consider_outliers:
        if method == 'mean':
            impute = df[col_name].mean()
            df.loc[outliers, col_name] = impute
        elif method == 'median':
            impute = df[col_name].median()
            df.loc[outliers, col_name] = impute
        elif method == 'mode':
            impute = stats.mode(df[col_name], keepdims = True)[0][0]
            df.loc[outliers, col_name] = impute
        elif method == 'std':
            impute = df[col_name].std()
            df.loc[outliers, col_name] = impute
        else:
            raise ValueError('Invalid Imputation, Valid Options: Mean, Median, Mode Or Std')
    else:
        if method == 'mean':
            impute = df.loc[~outliers, col_name].mean()
            df.loc[outliers, col_name] = impute
        elif method == 'median':
            impute = df.loc[~outliers, col_name].median()
            df.loc[outliers, col_name] = impute
        elif method == 'mode':
            impute = stats.mode(df.loc[~outliers, col_name], keepdims = True)[0][0]
            df.loc[outliers, col_name] = impute
        elif method == 'std':
            impute = df.loc[~outliers, col_name].std()
            df.loc[outliers, col_name] = impute
        else:
            raise ValueError('Invalid Imputation, Valid Options: Mean, Median, Mode Or Std')
            
    return df
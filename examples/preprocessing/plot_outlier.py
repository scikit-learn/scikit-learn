import matplotlib.pyplot as plt
from sklearn.preprocessing import FixOutliers
from sklearn.datasets import load_iris
import pandas as pd
import numpy as np
import seaborn as sns

#Load Dataset
iris = load_iris()

# Create a DataFrame from the data and target variables
df1 = pd.DataFrame(data=iris.data, columns=iris.feature_names)
df1['target'] = iris.target



#Fixing Outliers
df2 = df1.copy()
clf = FixOutliers(approach = 'every')
for i in iris.feature_names:
    df2 = clf.fit_transform(df2,i, treatment = 'remove')
df2 = pd.DataFrame(df2)


#Plot The Correlation Matrix Of df1
plt.figure(figsize=(10, 8))
corr1 = df1.corr()
sns.heatmap(corr1, annot=True, cmap='coolwarm')
plt.title('Correlation Matrix of df1')

#Plot The Correlation Matrix Of df2
plt.figure(figsize=(10, 8))
corr2 = df2.corr()
sns.heatmap(corr2, annot=True, cmap='coolwarm')
plt.title('Correlation Matrix of df2')

plt.show()
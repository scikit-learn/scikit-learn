from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import numpy as np
import pandas as pd
from StringIO import StringIO
import gzip
from urllib import urlopen
import time

# calculates the precision at the top k examples of the positive class
def precision_k(y_true, y_score, k):
    ranks = y_score.argsort()
    top_k = ranks[-k:]
    return y_true[top_k].sum()*1.0/k

# read a gzipped csv from a url into a pandas dataframe
def csv_from_gzip_url(url):
    f = StringIO(urlopen(url).read())
    s = gzip.GzipFile(fileobj=f, mode='rb')
    df = pd.read_csv(s)
    return df

# binarize all columns with object dtype
def binarize(df):
    categorical_columns = df.dtypes[df.dtypes == object].index
    for column in categorical_columns:
        categories = df[column].unique()
        for category in categories:
            df[category] = (df[column] == category)
        df.drop(column, axis=1, inplace=True)

# code the specified column as an integer
def code(df, column):
    categories = df[column].unique()
    for i, category in enumerate(categories):
        df.loc[df[column]==category, [column]] = i
    df[column] = df[column].astype(int)

kddtrain = csv_from_gzip_url('http://kdd.ics.uci.edu/databases/kddcup99/kddcup.data_10_percent.gz')
kddtest = csv_from_gzip_url('http://kdd.ics.uci.edu/databases/kddcup99/corrected.gz')

# rename columns because the csvs don't have headers
kddtrain.columns = range(42)
kddtest.columns = range(42)

kdd = pd.concat((kddtrain,kddtest))
code(kdd, 41)
binarize(kdd)

X = kdd.drop(41, axis=1).values
y = (kdd[41].values > 5)
X_train = X[0:len(kddtrain),:]
y_train = y[0:len(kddtrain)]
X_test = X[-len(kddtest):,:]
y_test = y[-len(kddtest):]

print 'baseline: {}'.format(y_train.sum()*1.0 / len(y_train)) # the minority class makes up 1.7% of the training set
print ''

common_params={'n_estimators':100, 'criterion':'entropy', 'n_jobs':-1}
params = [{}, {'class_weight':'auto'}, {'class_weight':'balanced_subsample'}, {'balanced':True}] # default, weighted random forest, balanced subsample, balanced random forest
k = y_test.sum()
for p in params:
    print 'forest parameters: {}'.format(p)
    p.update(common_params)
    clf = RandomForestClassifier(**p)

    start = time.clock()
    clf.fit(X_train,y_train)
    print 'time elapsed: {}'.format(time.clock() - start)

    y_score = clf.predict_proba(X_test)[:,1]
    y_predict = clf.predict(X_test)

    print 'precision at {}: {}'.format(k, precision_k(y_test, y_score, k))
    print 'auc: {}'.format(roc_auc_score(y_test, y_score))
    print ''

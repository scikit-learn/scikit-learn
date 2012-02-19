"""
Perform Naive Bayes classification of SDSS Quasar colors

Usage:

% python classification_naive_bayes.py data/sdss_colors/

First argument is the data directory where the colors data has been
downloaded.  Note that you must first run the 'fetch_data.py' script
in data/sdss_colors/
"""
import sys, os
import numpy as np

from sklearn import naive_bayes

try:
    data_dir = sys.argv[1]
    train_file = os.path.join(data_dir, 'sdssdr6_colors_class_train.npy')
    train_data = np.load(train_file)

    test_file = os.path.join(data_dir, 'sdssdr6_colors_class.200000.npy')
    test_data = np.load(test_file)
except:
    print __doc__
    sys.exit()

X_train = np.vstack([train_data['u-g'],
                     train_data['g-r'],
                     train_data['r-i'],
                     train_data['i-z']]).T

y_train = (train_data['redshift'] > 0).astype(int)
# y=0 -> quasar
# y=1 -> star

X_test = np.vstack([test_data['u-g'],
                    test_data['g-r'],
                    test_data['r-i'],
                    test_data['i-z']]).T
y_test = (test_data['label'] == 0).astype(int)
# y=0 -> quasar
# y=1 -> star


# Train the Naive Bayes Classifier
gnb = naive_bayes.GaussianNB()
gnb.fit(X_train, y_train)
y_pred = gnb.predict(X_test)


accuracy = float(np.sum(y_pred == y_test)) / len(y_pred)
print accuracy

print np.sum(y_test == 0)
print np.sum(y_test == 1)


# since there are many fewer quasars than stars, we'll call quasars positive
TP = np.sum((y_pred == 1) & (y_test == 1))
FP = np.sum((y_pred == 0) & (y_test == 1))
TN = np.sum((y_pred == 0) & (y_test == 0))
FN = np.sum((y_pred == 1) & (y_test == 0))

precision = TP / float(TP + FP)
recall = TP / float(TP + FN)

print precision
print recall

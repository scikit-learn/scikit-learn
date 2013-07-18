import numpy as np

from sklearn.cross_validation import cross_val_score
from sklearn.grid_search import GridSearchCV
from sklearn.datasets import load_digits
from sklearn.datasets import make_multilabel_classification

from sklearn.utils.validation import check_random_state

from sklearn.dummy import DummyClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier


from sklearn.metrics import accuracy_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics import f1_score

ESTIMATORS = {
    "dummy": DummyClassifier(),
    "1-nn": KNeighborsClassifier(n_neighbors=5),
    "5-nn": KNeighborsClassifier(n_neighbors=5),
    "extra trees": ExtraTreesClassifier(random_state=0),
    # "random forest": RandomForestClassifier(random_state=0),
    "decision tree": DecisionTreeClassifier(random_state=0),
}

rng = check_random_state(0)

X, y = make_multilabel_classification(n_samples=500,
                                      n_features=50,
                                      n_classes=10,
                                      return_indicator=True,
                                      allow_unlabeled=False,
                                      random_state=0)

# TODO compare hamming loss, jaccard and subset accuracy
#      as a function of the number of label

# TODO make histogram

for score_func in [hamming_loss, accuracy_score]:
    print score_func

    for name, est in ESTIMATORS.items():
        score = cross_val_score(est, X, y, score_func=hamming_loss)
        print "{0} = {1:0.4f}".format(name, np.mean(score))

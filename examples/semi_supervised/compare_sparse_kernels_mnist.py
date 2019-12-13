"""
=================================================
Label Propagation MNIST: Comparing Sparse Kernels
=================================================

This example compares the runtime and performance of two sparse kernels for
semisupervised learning on the MNIST digit dataset.

The MNIST dataset consists of 28x28 pixel grayscale images. Here, we will use a
subset of 10K images, reserving a fraction of these for testing.  We will
compare the performance and runtime of two sparse kernels, across a range of
low-supervision scenarios.

In each scenario, we will run each model multiple times, to increase our
confidence in the comparison between kernels.

The models will be evaluated for their accuracy at spreading labels during
training ("transductive learning"), as well as spreading labels to unseen
points at test time ("inductive learning").

The first kernel option produces a binary k-Nearest Neighbors adjacency matrix.
The second produces a kernel which is also k-sparse, but contains the same
weights as used in an RBF kernel.  

Notice that the performance of the sparse-RBF kernel is very sensitive to
parameters; the parameters used here were found by a quick manual search, so
the model can likely be improved with further optimization, and using this
kernel effectively on a new dataset requires hyperparameter tuning.
"""
import numpy as np
from sklearn.datasets import fetch_openml
from sklearn.semi_supervised import LabelSpreading
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.metrics import make_scorer
import time

Xorig, Yorig = fetch_openml('mnist_784', version=1, return_X_y=True)
Yorig = Y.astype(int)

# For a quick demonstration, use only a subset of the data
n_total = 10000
X = Xorig[:n_total, :]
Y = Yorig[:n_total]

# Save test set for inductive learning
test_fraction = 0.333
Xtrain, Xtest, Ytrain, Ytest = train_test_split(X, Y, test_size=test_fraction,
                                                random_state=0)

# Mask subset of train data for transductive learning
n_train = len(Ytrain)
#kwargs = {'gamma': 1e-9, 'n_neighbors': 50, 'n_jobs': -1, 'max_iter': 100}

#models = [LabelSpreading(kernel='knn', **kwargs),
#          LabelSpreading(kernel='sparse-rbf', **kwargs)]

#supervision_fractions = [0.001, 0.005, 0.01, 0.05, 0.1]

# First, we perform a grid search to optimize parameters for sparse-rbf kernel.
# For this purpose, we use a smaller subset of the data.
# Notice also that we 

class WrapLabelSpreading(LabelSpreading):
    """
    In order to perform a grid search over this semi-supervised model,
    we need to provide a thin wrapper that masks a subset of the data before 
    `fit` is called.
    """
    def __init__(self, supervision_fraction, kernel='sparse-rbf', gamma=20,
            n_neighbors=7, alpha=0.2, max_iter=30, tol=1e-3, n_jobs=None):

        self.supervision_fraction = supervision_fraction

        super().__init__(kernel=kernel, gamma=gamma,
                         n_neighbors=n_neighbors, alpha=alpha,
                         max_iter=max_iter, tol=tol, n_jobs=n_jobs)

    def fit(self, X, y):
        # mask a random subset of labels, based on self.supervision_fraction
        n_total = len(y)
        n_labeled = self.supervision_fraction * n_total

        indices = np.arange(n_total)
        np.random.seed(0)
        np.random.shuffle(indices)
        unlabeled_subset = indices[n_labeled:]

        y[unlabeled_subset] = -1

        super().fit(X,y)
        return self


# In all cases, we simply use max_iter=100
sparse_rbf_model = GridSearchCV(WrapLabelSpreading(kernel='sparse-rbf'),
                                param_grid= {
                                    'gamma': np.logspace(-8, 1, 10),
                                    'alpha': np.linspace(0, 1, 10),
                                     'n_neighbors': list(range(5,55,5))})

knn_model = GridSearchCV(WrapLabelSpreading(kernel='knn'),
                         param_grid= {
                             'n_neighbors': list(range(5,55,5))},
                             'alpha': np.linspace(0, 1, 10),
                            )


# Then, we compare the performance of optimized sparse-rbf kernel to knn kernel
supervision_fractions = [0.05, 0.1]
accuracies = {
        'transduction': { 'knn':[], 'sparse-rbf':[] }, 
        'induction': { 'knn':[], 'sparse-rbf':[] }
}
for supervision_fraction in supervision_fractions:
    supervision_fraction = 0.05
    n_labeled = int(supervision_fraction * n_train)
    indices = np.arange(n_train)
    unlabeled_set = indices[n_labeled:]

    Ymasked = np.copy(Ytrain)
    Ymasked[unlabeled_set] = -1

    for kernel_name, model in zip(['knn', 'sparse-rbf'], 
                                  [knn_model, sparse_rbf_model]):
        knn_acc_trans = []
        knn_acc_ind = []
        sparse_rbf_acc_trans = []
        sparse_rbf_acc_ind = []
        # Repeat each scenario 5 times to collect rough statistics
    #    for _ in range(5):
        print("="*80)
        t0 = time.time()
        print(f"MODEL: {model}")
        model.fit(Xtrain, Ymasked)
        t1 = time.time()

        predicted_labels = model.transduction_[unlabeled_set]
        true_labels = Ytrain[unlabeled_set]
        acc = np.sum(predicted_labels == true_labels) / len(unlabeled_set)
        print(f"accuracy: {acc}")



        print("-"*80)
        print(f"TRANSDUCTION: {n_labeled} labeled and " +
              f"{n_train - n_labeled} unlabeled points ({n_train} total)")
        print("-"*80)
        print("Confusion Matrix:")
        print(confusion_matrix(true_labels, predicted_labels,
                               labels=model.classes_))
        print("-"*80)
        print("Classification Report:")
        print(classification_report(true_labels, predicted_labels))
        print("-"*80)

        predicted_labels = model.predict(Xtest)
        t2 = time.time()

        print("-"*80)
        print(f"INDUCTION: {int(test_fraction * n_total)} test points")
        print("-"*80)
        print("Confusion Matrix:")
        print(confusion_matrix(Ytest, predicted_labels, labels=model.classes_))
        print("-"*80)
        print("Classification Report:")
        print(classification_report(Ytest, predicted_labels))
        print("-"*80)
        print(f"Runtimes: Transduction: {t1 - t0:.2f}s. Induction: {t2 - t1:.2f}s")

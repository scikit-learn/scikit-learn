import pstats
import cProfile

from sklearn.mlp.classes import MLPClassifier
from sklearn import datasets, preprocessing

digits = datasets.fetch_mldata("MNIST Original")
n_samples = len(digits.data)
data = digits.data.reshape((n_samples, -1))

data = data / 255.

#data = preprocessing.scale(data)
classifier = MLPClassifier(n_hidden=64, lr=0.3, lr_moment=0.3, batch_size=1)
cProfile.runctx("""classifier.fit(data[:n_samples / 2], digits.target[:n_samples / 2],
    max_epochs=10)""", globals(), locals(), "Profile.prof")
s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()

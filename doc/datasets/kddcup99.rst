
.. _kddcup99:

Kddcup 99 dataset
=================

The KDD Cup '99 dataset was created by processing the tcpdump portions
of the 1998 DARPA Intrusion Detection System (IDS) Evaluation dataset,
created by MIT Lincoln Lab. The artificial data (described on the `dataset's
homepage <http://kdd.ics.uci.edu/databases/kddcup99/kddcup99.html>`_) was
generated using a closed network and hand-injected attacks to produce a
large number of different types of attack with normal activity in the
background. As the initial goal was to produce a large training set for
supervised learning algorithms, there is a large proportion (80.1%) of
abnormal data which is unrealistic in real world, and inapropriate for
unsupervised anomaly detection which aims at detecting 'abnormal' data, ie
1) qualitatively different from normal data
2) in large minority among the observations.
We thus transform the KDD Data set into two differents data set: SA and SF.

-SA is obtained by simply selecting all the normal data, and a small
proportion of abnormal data to gives an anomaly proportion of 1%.

-SF is obtained as in [2]
by simply picking up the data whose attribute logged_in is positive, thus
focusing on the intrusion attack, which gives a proportion of 0.3% of
attack.

-http and smtp are two subsets of SF corresponding with third feature
equal to 'http' (resp. to 'smtp')

:func:`sklearn.datasets.fetch_kddcup99` will load the kddcup99 dataset;
it returns a dictionary-like object
with the feature matrix in the ``data`` member
and the target values in ``target``.
The dataset will be downloaded from the web if necessary.

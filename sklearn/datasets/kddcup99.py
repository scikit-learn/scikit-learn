"""KDDCUP 99 dataset.

A classic dataset for anomaly detection.

The dataset page is available from UCI Machine Learning Repository

https://archive.ics.uci.edu/ml/machine-learning-databases/kddcup99-mld/kddcup.data.gz

"""

import sys
import errno
from gzip import GzipFile
from io import BytesIO
import logging
import os
from os.path import exists, join
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

import numpy as np

from .base import get_data_home
from .base import Bunch
from ..externals import joblib, six
from ..utils import check_random_state
from ..utils import shuffle as shuffle_method


URL10 = ('http://archive.ics.uci.edu/ml/'
         'machine-learning-databases/kddcup99-mld/kddcup.data_10_percent.gz')

URL = ('http://archive.ics.uci.edu/ml/'
       'machine-learning-databases/kddcup99-mld/kddcup.data.gz')


logger = logging.getLogger()


def fetch_kddcup99(subset=None, shuffle=False, random_state=None,
                   percent10=True, download_if_missing=True):
    """Load and return the kddcup 99 dataset (classification).

    The KDD Cup '99 dataset was created by processing the tcpdump portions
    of the 1998 DARPA Intrusion Detection System (IDS) Evaluation dataset,
    created by MIT Lincoln Lab [1] . The artificial data was generated using
    a closed network and hand-injected attacks to produce a large number of
    different types of attack with normal activity in the background.
    As the initial goal was to produce a large training set for supervised
    learning algorithms, there is a large proportion (80.1%) of abnormal
    data which is unrealistic in real world, and inappropriate for unsupervised
    anomaly detection which aims at detecting 'abnormal' data, ie

    1) qualitatively different from normal data.

    2) in large minority among the observations.

    We thus transform the KDD Data set into two different data sets: SA and SF.

    - SA is obtained by simply selecting all the normal data, and a small
      proportion of abnormal data to gives an anomaly proportion of 1%.

    - SF is obtained as in [2]
      by simply picking up the data whose attribute logged_in is positive, thus
      focusing on the intrusion attack, which gives a proportion of 0.3% of
      attack.

    - http and smtp are two subsets of SF corresponding with third feature
      equal to 'http' (resp. to 'smtp')


    General KDD structure :

    ================      ==========================================
    Samples total         4898431
    Dimensionality        41
    Features              discrete (int) or continuous (float)
    Targets               str, 'normal.' or name of the anomaly type
    ================      ==========================================

    SA structure :

    ================      ==========================================
    Samples total         976158
    Dimensionality        41
    Features              discrete (int) or continuous (float)
    Targets               str, 'normal.' or name of the anomaly type
    ================      ==========================================

    SF structure :

    ================      ==========================================
    Samples total         699691
    Dimensionality        4
    Features              discrete (int) or continuous (float)
    Targets               str, 'normal.' or name of the anomaly type
    ================      ==========================================

    http structure :

    ================      ==========================================
    Samples total         619052
    Dimensionality        3
    Features              discrete (int) or continuous (float)
    Targets               str, 'normal.' or name of the anomaly type
    ================      ==========================================

    smtp structure :

    ================      ==========================================
    Samples total         95373
    Dimensionality        3
    Features              discrete (int) or continuous (float)
    Targets               str, 'normal.' or name of the anomaly type
    ================      ==========================================

    .. versionadded:: 0.18

    Parameters
    ----------
    subset : None, 'SA', 'SF', 'http', 'smtp'
        To return the corresponding classical subsets of kddcup 99.
        If None, return the entire kddcup 99 dataset.

    random_state : int, RandomState instance or None, optional (default=None)
        Random state for shuffling the dataset.
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    shuffle : bool, default=False
        Whether to shuffle dataset.

    percent10 : bool, default=False
        Whether to load only 10 percent of the data.

    download_if_missing : bool, default=True
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn and 'target', the regression target for each
        sample.


    References
    ----------
    .. [1] Analysis and Results of the 1999 DARPA Off-Line Intrusion
           Detection Evaluation Richard Lippmann, Joshua W. Haines,
           David J. Fried, Jonathan Korba, Kumar Das

    .. [2] A Geometric Framework for Unsupervised Anomaly Detection: Detecting
           Intrusions in Unlabeled Data (2002) by Eleazar Eskin, Andrew Arnold,
           Michael Prerau, Leonid Portnoy, Sal Stolfo

    """
    kddcup99 = _fetch_brute_kddcup99(shuffle=shuffle, percent10=percent10,
                                     download_if_missing=download_if_missing)

    data = kddcup99.data
    target = kddcup99.target

    if subset == 'SA':
        s = target == b'normal.'
        t = np.logical_not(s)
        normal_samples = data[s, :]
        normal_targets = target[s]
        abnormal_samples = data[t, :]
        abnormal_targets = target[t]

        n_samples_abnormal = abnormal_samples.shape[0]
        # selected abnormal samples:
        random_state = check_random_state(random_state)
        r = random_state.randint(0, n_samples_abnormal, 3377)
        abnormal_samples = abnormal_samples[r]
        abnormal_targets = abnormal_targets[r]

        data = np.r_[normal_samples, abnormal_samples]
        target = np.r_[normal_targets, abnormal_targets]

    if subset == 'SF' or subset == 'http' or subset == 'smtp':
        # select all samples with positive logged_in attribute:
        s = data[:, 11] == 1
        data = np.c_[data[s, :11], data[s, 12:]]
        target = target[s]

        data[:, 0] = np.log((data[:, 0] + 0.1).astype(float))
        data[:, 4] = np.log((data[:, 4] + 0.1).astype(float))
        data[:, 5] = np.log((data[:, 5] + 0.1).astype(float))

        if subset == 'http':
            s = data[:, 2] == b'http'
            data = data[s]
            target = target[s]
            data = np.c_[data[:, 0], data[:, 4], data[:, 5]]

        if subset == 'smtp':
            s = data[:, 2] == b'smtp'
            data = data[s]
            target = target[s]
            data = np.c_[data[:, 0], data[:, 4], data[:, 5]]

        if subset == 'SF':
            data = np.c_[data[:, 0], data[:, 2], data[:, 4], data[:, 5]]

    return Bunch(data=data, target=target)


def _fetch_brute_kddcup99(subset=None, data_home=None,
                          download_if_missing=True, random_state=None,
                          shuffle=False, percent10=False):

    """Load the kddcup99 dataset, downloading it if necessary.

    Parameters
    ----------
    subset : None, 'SA', 'SF', 'http', 'smtp'
        To return the corresponding classical subsets of kddcup 99.
        If None, return the entire kddcup 99 dataset.

    data_home : string, optional
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing : boolean, default=True
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    random_state : int, RandomState instance or None, optional (default=None)
        Random state for shuffling the dataset.
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    shuffle : bool, default=False
        Whether to shuffle dataset.

    percent10 : bool, default=False
        Whether to load only 10 percent of the data.

    Returns
    -------
    dataset : dict-like object with the following attributes:
        dataset.data : numpy array of shape (494021, 41)
            Each row corresponds to the 41 features in the dataset.
        dataset.target : numpy array of shape (494021,)
            Each value corresponds to one of the 21 attack types or to the
            label 'normal.'.
        dataset.DESCR : string
            Description of the kddcup99 dataset.

    """

    data_home = get_data_home(data_home=data_home)
    if sys.version_info[0] == 3:
        # The zlib compression format use by joblib is not compatible when
        # switching from Python 2 to Python 3, let us use a separate folder
        # under Python 3:
        dir_suffix = "-py3"
    else:
        # Backward compat for Python 2 users
        dir_suffix = ""
    if percent10:
        kddcup_dir = join(data_home, "kddcup99_10" + dir_suffix)
    else:
        kddcup_dir = join(data_home, "kddcup99" + dir_suffix)
    samples_path = join(kddcup_dir, "samples")
    targets_path = join(kddcup_dir, "targets")
    available = exists(samples_path)

    if download_if_missing and not available:
        _mkdirp(kddcup_dir)
        URL_ = URL10 if percent10 else URL
        logger.warning("Downloading %s" % URL_)
        f = BytesIO(urlopen(URL_).read())

        dt = [('duration', int),
              ('protocol_type', 'S4'),
              ('service', 'S11'),
              ('flag', 'S6'),
              ('src_bytes', int),
              ('dst_bytes', int),
              ('land', int),
              ('wrong_fragment', int),
              ('urgent', int),
              ('hot', int),
              ('num_failed_logins', int),
              ('logged_in', int),
              ('num_compromised', int),
              ('root_shell', int),
              ('su_attempted', int),
              ('num_root', int),
              ('num_file_creations', int),
              ('num_shells', int),
              ('num_access_files', int),
              ('num_outbound_cmds', int),
              ('is_host_login', int),
              ('is_guest_login', int),
              ('count', int),
              ('srv_count', int),
              ('serror_rate', float),
              ('srv_serror_rate', float),
              ('rerror_rate', float),
              ('srv_rerror_rate', float),
              ('same_srv_rate', float),
              ('diff_srv_rate', float),
              ('srv_diff_host_rate', float),
              ('dst_host_count', int),
              ('dst_host_srv_count', int),
              ('dst_host_same_srv_rate', float),
              ('dst_host_diff_srv_rate', float),
              ('dst_host_same_src_port_rate', float),
              ('dst_host_srv_diff_host_rate', float),
              ('dst_host_serror_rate', float),
              ('dst_host_srv_serror_rate', float),
              ('dst_host_rerror_rate', float),
              ('dst_host_srv_rerror_rate', float),
              ('labels', 'S16')]
        DT = np.dtype(dt)

        file_ = GzipFile(fileobj=f, mode='r')
        Xy = []
        for line in file_.readlines():
            if six.PY3:
                line = line.decode()
            Xy.append(line.replace('\n', '').split(','))
        file_.close()
        print('extraction done')
        Xy = np.asarray(Xy, dtype=object)
        for j in range(42):
            Xy[:, j] = Xy[:, j].astype(DT[j])

        X = Xy[:, :-1]
        y = Xy[:, -1]
        # XXX bug when compress!=0:
        # (error: 'Incorrect data length while decompressing[...] the file
        #  could be corrupted.')

        joblib.dump(X, samples_path, compress=0)
        joblib.dump(y, targets_path, compress=0)

    try:
        X, y
    except NameError:
        X = joblib.load(samples_path)
        y = joblib.load(targets_path)

    if shuffle:
        X, y = shuffle_method(X, y, random_state=random_state)

    return Bunch(data=X, target=y, DESCR=__doc__)


def _mkdirp(d):
    """Ensure directory d exists (like mkdir -p on Unix)
    No guarantee that the directory is writable.
    """
    try:
        os.makedirs(d)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

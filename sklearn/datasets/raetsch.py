import os
from scipy.io import loadmat
from sklearn.datasets.base import (_fetch_remote, get_data_home, Bunch,
                                   RemoteFileMetadata)


DATASETS = {'banana', 'breast_cancer', 'diabetis', 'flare_solar', 'german',
            'heart', 'image', 'ringnorm', 'splice', 'thyroid', 'titanic',
            'twonorm', 'waveform'}
ARCHIVE = RemoteFileMetadata(filename='benchmarks.mat',
                             url='https://github.com/tdiethe/gunnar_raetsch_benchmark_datasets/raw/master/benchmarks.mat',
                             checksum=('47c19e4bc4716edc4077cfa5ea61edf4d02af4ec51a0ecfe035626ae8b561c75'))


def fetch_raetsch(name, data_home=None):
    """Fetch Gunnar Raetsch's dataset.

    Fetch a Gunnar Raetsch's benchmark dataset by name. Availabe datasets are
    'banana', 'breast_cancer', 'diabetis', 'flare_solar', 'german', 'heart',
    'image', 'ringnorm', 'splice', 'thyroid', 'titanic', 'twonorm' and
    'waveform'. More info at
    https://github.com/tdiethe/gunnar_raetsch_benchmark_datasets.

    Parameters
    ----------
    name : string
        Dataset name.
    data_home : string or None, default None
        Specify another download and cache folder for the data sets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    Returns
    -------
    data : Bunch
        Dictionary-like object with all the data and metadata.

    """
    if name not in DATASETS:
        raise Exception('Avaliable datasets are ' + str(list(DATASETS)))
    dirname = os.path.join(get_data_home(data_home=data_home), 'raetsch')
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    filename = _fetch_remote(ARCHIVE, dirname=dirname)
    X, y, train_splits, test_splits = loadmat(filename)[name][0][0]
    inner_cv = zip(train_splits[:5] - 1, test_splits[:5] - 1)
    outer_cv = zip(train_splits - 1, test_splits - 1)
    return Bunch(data=X, target=y, inner_cv=inner_cv, outer_cv=outer_cv,
                 DESCR=name)

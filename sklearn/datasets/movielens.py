# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>,
#          Artem Sobolev
# License: BSD 3 clause

from io import BytesIO
from os.path import join, exists
from zipfile import ZipFile

try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen

import numpy as np
import scipy

from .base import get_data_home, Bunch
from ..externals import joblib

DATASETS = {
    '100k' : {
        'url' : 'http://www.grouplens.org/system/files/ml-100k.zip',
        'data_file' : 'ml-100k/u.data',
        'shape' : (943, 1682),
        'separator' : '\t',
        'nnz' : 100000,
    },
    '1m' : {
        'url' : 'http://www.grouplens.org/system/files/ml-1m.zip',
        'data_file' : 'ml-1m/ratings.dat',
        'shape' : (6040, 3952),
        'separator' : '::',
        'nnz' : 1000209,
    },
    '10m' : {
        'url' : 'http://www.grouplens.org/sites/www.grouplens.org/external_files/data/ml-10m.zip',
        'data_file' : 'ml-10M100K/ratings.dat',
        'shape' : (71567, 65133),
        'separator' : '::',
        'nnz' : 10000054,
    },
}

def fetch_movielens(dataset='100k', data_home=None, download_if_missing=True):

    data_home = get_data_home(data_home=data_home)
    zip_filename = join(data_home, 'ml-{}.zip'.format(dataset))
    target_filename = join(data_home,'ml-{}.pkz'.format(dataset))

    dataset_meta = DATASETS[dataset]
    data_url = dataset_meta['url']
    file_in_zip = dataset_meta['data_file']
    shape = dataset_meta['shape']
    delimiter = dataset_meta['separator']
    nnz = dataset_meta['nnz']

    if not exists(zip_filename):
        print('downloading movielens from {} to {}'.format(data_url, data_home))
        file_content = urlopen(data_url).read()
        with open(zip_filename, 'w') as f:
            f.write(file_content)

    if not exists(target_filename):

        with open(zip_filename, 'r') as f:
            buf = BytesIO(f.read())

        zip_file = ZipFile(buf)
        try:
            data_fd = zip_file.open(file_in_zip.format(dataset), 'r')

            # Extract the ratings
            timestamps = np.empty(nnz)
            data = np.empty((nnz, 3))
            for k in range(nnz):
                line = data_fd.readline()
                line = [float(el) for el in line.split(delimiter) if el]
                data[k, :] = line[:3]
                timestamps[k] = line[3]

            # Save them to the disk
            X = scipy.sparse.coo_matrix(
                (data[:, 2],
                 (data[:, 0]-1,
                  data[:, 1]-1)),
                shape=shape, dtype=np.float64_t)
            bunch = Bunch(data=X, timestamps=timestamps)
            joblib.dump(bunch, join(data_home, target_filename), compress=6)

            return bunch
        finally:
            zip_file.close()
    else:
        return joblib.load(join(data_home, target_filename))

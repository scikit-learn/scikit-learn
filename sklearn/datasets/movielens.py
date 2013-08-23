
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
    '100k' : ('http://www.grouplens.org/system/files/ml-100k.zip',
              'ml-100k/u.data',
              (943, 1682),
              '\t',
              100000,),
    '1m' : ('http://www.grouplens.org/system/files/ml-1m.zip',
            'ml-1m/ratings.dat',
            (6040, 3952),
            '::',
            1000209,),
    '10m' : ('http://www.grouplens.org/sites/www.grouplens.org/external_files/data/ml-10m.zip',
            'ml-10M100K/ratings.dat',
            (71567, 65133),
            '::',
            10000054,),
}

def fetch_movielens(dataset='100k', data_home=None, download_if_missing=True):

    data_home = get_data_home(data_home=data_home)
    zip_filename = join(data_home, 'ml-{}.zip'.format(dataset))
    target_filename = join(data_home,'ml-{}.pkz'.format(dataset))
    data_url = DATASETS[dataset][0]
    file_in_zip = DATASETS[dataset][1]
    shape = DATASETS[dataset][2]
    delimiter = DATASETS[dataset][3]
    nnz = DATASETS[dataset][4]


    if not exists(zip_filename):
        print('downloading movielens from %s to %s' % (data_url, data_home))
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
                shape=shape,)
            bunch = Bunch(data=X)
            joblib.dump(bunch, join(data_home, target_filename), compress=6)

            return bunch
        finally:
            zip_file.close()
    else:
        return joblib.load(join(data_home, target_filename))

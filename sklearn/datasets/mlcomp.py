# Copyright (c) 2010 Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause
"""Glue code to load http://mlcomp.org data as a scikit.learn dataset"""

import os
import numbers
from sklearn.datasets.base import load_files


def _load_document_classification(dataset_path, metadata, set_=None, **kwargs):
    if set_ is not None:
        dataset_path = os.path.join(dataset_path, set_)
    return load_files(dataset_path, metadata.get('description'), **kwargs)


LOADERS = {
    'DocumentClassification': _load_document_classification,
    # TODO: implement the remaining domain formats
}


def load_mlcomp(name_or_id, set_="raw", mlcomp_root=None, **kwargs):
    """Load a datasets as downloaded from http://mlcomp.org

    Parameters
    ----------

    name_or_id : the integer id or the string name metadata of the MLComp
                 dataset to load

    set_ : select the portion to load: 'train', 'test' or 'raw'

    mlcomp_root : the filesystem path to the root folder where MLComp datasets
                  are stored, if mlcomp_root is None, the MLCOMP_DATASETS_HOME
                  environment variable is looked up instead.

    **kwargs : domain specific kwargs to be passed to the dataset loader.

    Returns
    -------

    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'filenames', the files holding the raw to learn, 'target', the
        classification labels (integer index), 'target_names',
        the meaning of the labels, and 'DESCR', the full description of the
        dataset.

    Note on the lookup process: depending on the type of name_or_id,
    will choose between integer id lookup or metadata name lookup by
    looking at the unzipped archives and metadata file.

    TODO: implement zip dataset loading too
    """

    if mlcomp_root is None:
        try:
            mlcomp_root = os.environ['MLCOMP_DATASETS_HOME']
        except KeyError:
            raise ValueError("MLCOMP_DATASETS_HOME env variable is undefined")

    mlcomp_root = os.path.expanduser(mlcomp_root)
    mlcomp_root = os.path.abspath(mlcomp_root)
    mlcomp_root = os.path.normpath(mlcomp_root)

    if not os.path.exists(mlcomp_root):
        raise ValueError("Could not find folder: " + mlcomp_root)

    # dataset lookup
    if isinstance(name_or_id, numbers.Integral):
        # id lookup
        dataset_path = os.path.join(mlcomp_root, str(name_or_id))
    else:
        # assume name based lookup
        dataset_path = None
        expected_name_line = "name: " + name_or_id
        for dataset in os.listdir(mlcomp_root):
            metadata_file = os.path.join(mlcomp_root, dataset, 'metadata')
            if not os.path.exists(metadata_file):
                continue
            with open(metadata_file) as f:
                for line in f:
                    if line.strip() == expected_name_line:
                        dataset_path = os.path.join(mlcomp_root, dataset)
                        break
        if dataset_path is None:
            raise ValueError("Could not find dataset with metadata line: " +
                             expected_name_line)

    # loading the dataset metadata
    metadata = dict()
    metadata_file = os.path.join(dataset_path, 'metadata')
    if not os.path.exists(metadata_file):
        raise ValueError(dataset_path + ' is not a valid MLComp dataset')
    with open(metadata_file) as f:
        for line in f:
            if ":" in line:
                key, value = line.split(":", 1)
                metadata[key.strip()] = value.strip()

    format = metadata.get('format', 'unknow')
    loader = LOADERS.get(format)
    if loader is None:
        raise ValueError("No loader implemented for format: " + format)
    return loader(dataset_path, metadata, set_=set_, **kwargs)

# Copyright (c) 2010 Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD
"""Glue code to load http://mlcomp.org data as a scikit.learn dataset"""

import os
import numpy as np
from scikits.learn.datasets.base import Bunch
from scikits.learn.feature_extraction.text import HashingVectorizer
from scikits.learn.feature_extraction.text import SparseHashingVectorizer


def _load_document_classification(dataset_path, metadata, set_, sparse, **kw):
    """Loader implementation for the DocumentClassification format"""
    target = []
    target_names = []
    filenames = []
    vectorizer = kw.get('vectorizer')
    if vectorizer is None:
        if sparse:
            vectorizer = SparseHashingVectorizer()
        else:
            vectorizer = HashingVectorizer()

    # TODO: make it possible to plug a several pass system to filter-out tokens
    # that occur in more than 30% of the documents for instance.

    # TODO: use joblib.Parallel or multiprocessing to parallelize the following
    # (provided this is not IO bound)

    dataset_path = os.path.join(dataset_path, set_)
    folders = [f for f in sorted(os.listdir(dataset_path))
               if os.path.isdir(os.path.join(dataset_path, f))]
    for label, folder in enumerate(folders):
        target_names.append(folder)
        folder_path = os.path.join(dataset_path, folder)
        documents = [os.path.join(folder_path, d)
                     for d in sorted(os.listdir(folder_path))]
        vectorizer.vectorize_files(documents)
        target.extend(len(documents) * [label])
        filenames.extend(documents)

    return Bunch(data=vectorizer.get_vectors(), target=np.array(target),
                 target_names=target_names, filenames=filenames,
                 DESCR=metadata.get('description'))


LOADERS = {
    'DocumentClassification': _load_document_classification,
    # TODO: implement the remaining domain formats
}


def load_mlcomp(name_or_id, set_="raw", mlcomp_root=None, sparse=False,
                **kwargs):
    """Load a datasets as downloaded from http://mlcomp.org

    Parameters
    ----------

    name_or_id : the integer id or the string name metadata of the MLComp
                 dataset to load

    set_ : select the portion to load: 'train', 'test' or 'raw'

    mlcomp_root : the filesystem path to the root folder where MLComp datasets
                  are stored, if mlcomp_root is None, the MLCOMP_DATASETS_HOME
                  environment variable is looked up instead.

    sparse : boolean if True then use a scipy.sparse matrix for the data field,
             False by default

    **kwargs : domain specific kwargs to be passed to the dataset loader.

    Returns
    -------

    data : Bunch
        Dictionnary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels,
        'target_names', the meaning of the labels, and 'DESCR', the
        full description of the dataset.

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
    if isinstance(name_or_id, int):
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
            for line in file(metadata_file):
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
    for line in file(metadata_file):
        if ":" in line:
            key, value = line.split(":", 1)
            metadata[key.strip()] = value.strip()

    format = metadata.get('format', 'unknow')
    loader = LOADERS.get(format)
    if loader is None:
        raise ValueError("No loader implemented for format: " + format)
    return loader(dataset_path, metadata, set_=set_, sparse=sparse, **kwargs)



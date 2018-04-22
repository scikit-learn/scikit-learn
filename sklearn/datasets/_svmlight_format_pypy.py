# Optimized inner loop of load_svmlight_file.
#
# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck
#          Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

import array

import numpy as np
from ..externals import six

COMMA = u','.encode('ascii')
COLON = u':'.encode('ascii')


def _load_svmlight_file(f, dtype, multilabel, zero_based, query_id,
                        offset, length):
    qid_prefix = six.b('qid')
    offset_max = offset + length if length > 0 else -1

    # Special-case float32 but use float64 for everything else;
    # the Python code will do further conversions.
    if dtype != np.float32:
        dtype = np.float64
    data = []
    indices = array.array("i")
    indptr = [0]
    query = np.arange(0, dtype=np.int64)

    labels = []

    if offset > 0:
        f.seek(offset)
        # drop the current line that might be truncated and is to be
        # fetched by another call
        f.readline()

    for line in f:
        # skip comments
        hash_pos = line.find(b'#')
        if hash_pos >= 0:
            line = line[:hash_pos]

        line_parts = line.split()
        if len(line_parts) == 0:
            continue

        target, features = line_parts[0], line_parts[1:]
        if multilabel:
            if COLON in target:
                target, features = [], line_parts[0:]
            else:
                target = [float(y) for y in target.split(COMMA)]
            target.sort()
            labels.append(tuple(target))
        else:
            labels.append(float(target))

        prev_idx = -1
        n_features = len(features)
        if n_features and features[0].startswith(qid_prefix):
            _, value = features[0].split(COLON, 1)
            if query_id:
                query = np.resize(query, len(query) + 1)
                query[len(query) - 1] = np.int64(value)
            features.pop(0)
            n_features -= 1

        for i in range(n_features):
            idx_s, value = features[i].split(COLON, 1)
            idx = int(idx_s)
            if idx < 0 or not zero_based and idx == 0:
                raise ValueError(
                    "Invalid index %d in SVMlight/LibSVM data file." % idx)
            if idx <= prev_idx:
                raise ValueError("Feature indices in SVMlight/LibSVM data "
                                 "file should be sorted and unique.")

            indices.append(idx)

            data.append(float(value))

            prev_idx = idx

        indptr.append(len(data))

        if offset_max != -1 and f.tell() > offset_max:
            # Stop here and let another call deal with the following.
            break

    if not multilabel:
        labels = array.array('d', labels)
    indices = array.array('i', indices)
    if dtype == np.float32:
        data = array.array("f", data)
    else:
        data = array.array("d", data)
    indptr = array.array("i", indptr)

    return (dtype, data, indices, indptr, labels, query)

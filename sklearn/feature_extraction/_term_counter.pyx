import scipy.sparse as sp

def term_counts(vocabulary, term_count_dicts, dtype):
    i_indices = []
    j_indices = []
    values = []
    vocabulary = vocabulary

    for i, term_count_dict in enumerate(term_count_dicts):
        for term, count in term_count_dict.iteritems():
            j = vocabulary.get(term)
            if j is not None:
                i_indices.append(i)
                j_indices.append(j)
                values.append(count)
        # free memory as we go
        term_count_dict.clear()

    shape = (len(term_count_dicts), max(vocabulary.itervalues()) + 1)
    spmatrix = sp.coo_matrix((values, (i_indices, j_indices)),
                             shape=shape, dtype=dtype)

    return spmatrix

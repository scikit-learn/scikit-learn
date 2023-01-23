from ._tree cimport Tree


cdef class QuantileTree(Tree):
    def __cinit__(self):
        self.leaf_node_ids = vector[SIZE_t](self.capacity)
        self.leaf_samples = vector[DOUBLE_t](self.capacity)
         
    cdef int _set_leaf_node(
        self,
        SplitRecord* split_node,
        Node* node
    ) nogil except -1:
        """Set leaf node in a quantile tree.
        
        Leaf nodes are comprised of not only the traditional data
        stored in a random forest, but also every observation that falls
        into the leaf node. So we want ``y[samples[start:pos], :]`` or
        ``y[samples[pos:end], :]`` depending on which way a leaf node went.
        """
        cdef SIZE_t node_id = self.node_count

        # set relevant leaf node data
        node.left_child = _TREE_LEAF
        node.right_child = _TREE_LEAF
        node.feature = _TREE_UNDEFINED
        node.threshold = _TREE_UNDEFINED
        node.observations

        # now set leaf node data
        self.leaf_node_ids.push_back(node_id)
        
        # somehow need to pass in the samples 'start' and 'end' point that is in the
        # leaf node and add them to the data structure
        # XXX: currently this pseudocode assumes ``y`` is a a 1D vector, which it may not
        # be. So what container for ``leaf_samples`` would support ``y`` being 2D?
        for s_idx in range(start, end):
            self.leaf_samples[self.node_count].push_back(y[s_idx])
        return 1

    cpdef cnp.ndarray transform(self, object X):
        """Return an array of output values from the leaves X falls into.
        
        Parameters
        ----------
        X : cnp.ndarray of shape (n_samples, n_features)
            The input data array that we want to transform each sample.
        
        Returns
        -------
        out : cnp.ndarray of shape (n_samples, n_leaf_samples)
            For each sample, based on the leaf that it falls in, return the stored
            leaf node samples from the training dataset.
        """
        # get the indices of the X leaves - array of shape (n_samples, n_classes)
        X_leaves = self.apply(X)

        # extract array of values in leaves
        out = _extract_values(X)
        return out

    cpdef cnp.ndarray predict(self, object X, double quantiles=0.5):
        """Predict target for X is now using quantiles."""
        # get the indices of the X leaves - array of shape (n_samples, n_classes)
        X_leaves = self.apply(X)

        # [Not Implemented Function] compute the quantile
        out = _calculate_quantile(X_leaves, quantiles)

        if self.n_outputs == 1:
            out = out.reshape(X.shape[0], self.max_n_classes)
        return out
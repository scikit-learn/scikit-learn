# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
#
# Author: Thomas Moreau <thomas.moreau.2010@gmail.com>
# Author: Olivier Grisel <olivier.grisel@ensta.fr>


from cpython cimport Py_INCREF, PyObject, PyTypeObject

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from libc.stdio cimport printf
from libc.stdint cimport SIZE_MAX

from ..tree._utils cimport safe_realloc, sizet_ptr_to_ndarray
from ..utils import check_array

import numpy as np
cimport numpy as np
np.import_array()

cdef extern from "math.h":
    float fabsf(float x) nogil

cdef extern from "numpy/arrayobject.h":
    object PyArray_NewFromDescr(PyTypeObject* subtype, np.dtype descr,
                                int nd, np.npy_intp* dims,
                                np.npy_intp* strides,
                                void* data, int flags, object obj)
    int PyArray_SetBaseObject(np.ndarray arr, PyObject* obj)

# Build the corresponding numpy dtype for Cell.
# This works by casting `dummy` to an array of Cell of length 1, which numpy
# can construct a `dtype`-object for. See https://stackoverflow.com/q/62448946
# for a more detailed explanation.
cdef Cell dummy;
CELL_DTYPE = np.asarray(<Cell[:1]>(&dummy)).dtype

assert CELL_DTYPE.itemsize == sizeof(Cell)


cdef class _QuadTree:
    """Array-based representation of a QuadTree.

    This class is currently working for indexing 2D data (regular QuadTree) and
    for indexing 3D data (OcTree). It is planned to split the 2 implementations
    using `Cython.Tempita` to save some memory for QuadTree.

    Note that this code is currently internally used only by the Barnes-Hut
    method in `sklearn.manifold.TSNE`. It is planned to be refactored and
    generalized in the future to be compatible with nearest neighbors API of
    `sklearn.neighbors` with 2D and 3D data.
    """
    def __cinit__(self, int n_dimensions, int verbose):
        """Constructor."""
        # Parameters of the tree
        self.n_dimensions = n_dimensions
        self.verbose = verbose
        self.n_cells_per_cell = 2 ** self.n_dimensions

        # Inner structures
        self.max_depth = 0
        self.cell_count = 0
        self.capacity = 0
        self.n_points = 0
        self.cells = NULL

    def __dealloc__(self):
        """Destructor."""
        # Free all inner structures
        free(self.cells)

    property cumulative_size:
        def __get__(self):
            return self._get_cell_ndarray()['cumulative_size'][:self.cell_count]

    property leafs:
        def __get__(self):
            return self._get_cell_ndarray()['is_leaf'][:self.cell_count]

    def build_tree(self, X):
        """Build a tree from an array of points X."""
        cdef:
            int i
            DTYPE_t[3] pt
            DTYPE_t[3] min_bounds, max_bounds

        # validate X and prepare for query
        # X = check_array(X, dtype=DTYPE_t, order='C')
        n_samples = X.shape[0]

        capacity = 100
        self._resize(capacity)
        m = np.min(X, axis=0)
        M = np.max(X, axis=0)
        # Scale the maximum to get all points strictly in the tree bounding box
        # The 3 bounds are for positive, negative and small values
        M = np.maximum(M * (1. + 1e-3 * np.sign(M)), M + 1e-3)
        for i in range(self.n_dimensions):
            min_bounds[i] = m[i]
            max_bounds[i] = M[i]

            if self.verbose > 10:
                printf("[QuadTree] bounding box axis %i : [%f, %f]\n",
                       i, min_bounds[i], max_bounds[i])

        # Create the initial node with boundaries from the dataset
        self._init_root(min_bounds, max_bounds)

        for i in range(n_samples):
            for j in range(self.n_dimensions):
                pt[j] = X[i, j]
            self.insert_point(pt, i)

        # Shrink the cells array to reduce memory usage
        self._resize(capacity=self.cell_count)

    cdef int insert_point(self, DTYPE_t[3] point, SIZE_t point_index,
                          SIZE_t cell_id=0) nogil except -1:
        """Insert a point in the QuadTree."""
        cdef int ax
        cdef DTYPE_t n_frac
        cdef SIZE_t selected_child
        cdef Cell* cell = &self.cells[cell_id]
        cdef SIZE_t n_point = cell.cumulative_size

        if self.verbose > 10:
            printf("[QuadTree] Inserting depth %li\n", cell.depth)

        # Assert that the point is in the right range
        if DEBUGFLAG:
            self._check_point_in_cell(point, cell)

        # If the cell is an empty leaf, insert the point in it
        if cell.cumulative_size == 0:
            cell.cumulative_size = 1
            self.n_points += 1
            for i in range(self.n_dimensions):
                cell.barycenter[i] = point[i]
            cell.point_index = point_index
            if self.verbose > 10:
                printf("[QuadTree] inserted point %li in cell %li\n",
                       point_index, cell_id)
            return cell_id

        # If the cell is not a leaf, update cell internals and
        # recurse in selected child
        if not cell.is_leaf:
            for ax in range(self.n_dimensions):
                # barycenter update using a weighted mean
                cell.barycenter[ax] = (
                    n_point * cell.barycenter[ax] + point[ax]) / (n_point + 1)

            # Increase the size of the subtree starting from this cell
            cell.cumulative_size += 1

            # Insert child in the correct subtree
            selected_child = self._select_child(point, cell)
            if self.verbose > 49:
                printf("[QuadTree] selected child %li\n", selected_child)
            if selected_child == -1:
                self.n_points += 1
                return self._insert_point_in_new_child(point, cell, point_index)
            return self.insert_point(point, point_index, selected_child)

        # Finally, if the cell is a leaf with a point already inserted,
        # split the cell in n_cells_per_cell if the point is not a duplicate.
        # If it is a duplicate, increase the size of the leaf and return.
        if self._is_duplicate(point, cell.barycenter):
            if self.verbose > 10:
                printf("[QuadTree] found a duplicate!\n")
            cell.cumulative_size += 1
            self.n_points += 1
            return cell_id

        # In a leaf, the barycenter correspond to the only point included
        # in it.
        self._insert_point_in_new_child(cell.barycenter, cell, cell.point_index,
                                        cell.cumulative_size)
        return self.insert_point(point, point_index, cell_id)

    # XXX: This operation is not Thread safe
    cdef SIZE_t _insert_point_in_new_child(self, DTYPE_t[3] point, Cell* cell,
                                          SIZE_t point_index, SIZE_t size=1
                                          ) nogil:
        """Create a child of cell which will contain point."""

        # Local variable definition
        cdef:
            SIZE_t cell_id, cell_child_id, parent_id
            DTYPE_t[3] save_point
            DTYPE_t width
            Cell* child
            int i

        # If the maximal capacity of the Tree have been reached, double the capacity
        # We need to save the current cell id and the current point to retrieve them
        # in case the reallocation
        if self.cell_count + 1 > self.capacity:
            parent_id = cell.cell_id
            for i in range(self.n_dimensions):
                save_point[i] = point[i]
            self._resize(SIZE_MAX)
            cell = &self.cells[parent_id]
            point = save_point

        # Get an empty cell and initialize it
        cell_id = self.cell_count
        self.cell_count += 1
        child  = &self.cells[cell_id]

        self._init_cell(child, cell.cell_id, cell.depth + 1)
        child.cell_id = cell_id

        # Set the cell as an inner cell of the Tree
        cell.is_leaf = False
        cell.point_index = -1

        # Set the correct boundary for the cell, store the point in the cell
        # and compute its index in the children array.
        cell_child_id = 0
        for i in range(self.n_dimensions):
            cell_child_id *= 2
            if point[i] >= cell.center[i]:
                cell_child_id += 1
                child.min_bounds[i] = cell.center[i]
                child.max_bounds[i] = cell.max_bounds[i]
            else:
                child.min_bounds[i] = cell.min_bounds[i]
                child.max_bounds[i] = cell.center[i]
            child.center[i] = (child.min_bounds[i] + child.max_bounds[i]) / 2.
            width = child.max_bounds[i] - child.min_bounds[i]

            child.barycenter[i] = point[i]
            child.squared_max_width = max(child.squared_max_width, width*width)

        # Store the point info and the size to account for duplicated points
        child.point_index = point_index
        child.cumulative_size = size

        # Store the child cell in the correct place in children
        cell.children[cell_child_id] = child.cell_id

        if DEBUGFLAG:
            # Assert that the point is in the right range
            self._check_point_in_cell(point, child)
        if self.verbose > 10:
            printf("[QuadTree] inserted point %li in new child %li\n",
                   point_index, cell_id)

        return cell_id


    cdef bint _is_duplicate(self, DTYPE_t[3] point1, DTYPE_t[3] point2) nogil:
        """Check if the two given points are equals."""
        cdef int i
        cdef bint res = True
        for i in range(self.n_dimensions):
            # Use EPSILON to avoid numerical error that would overgrow the tree
            res &= fabsf(point1[i] - point2[i]) <= EPSILON
        return res


    cdef SIZE_t _select_child(self, DTYPE_t[3] point, Cell* cell) nogil:
        """Select the child of cell which contains the given query point."""
        cdef:
            int i
            SIZE_t selected_child = 0

        for i in range(self.n_dimensions):
            # Select the correct child cell to insert the point by comparing
            # it to the borders of the cells using precomputed center.
            selected_child *= 2
            if point[i] >= cell.center[i]:
                selected_child += 1
        return cell.children[selected_child]

    cdef void _init_cell(self, Cell* cell, SIZE_t parent, SIZE_t depth) nogil:
        """Initialize a cell structure with some constants."""
        cell.parent = parent
        cell.is_leaf = True
        cell.depth = depth
        cell.squared_max_width = 0
        cell.cumulative_size = 0
        for i in range(self.n_cells_per_cell):
            cell.children[i] = SIZE_MAX

    cdef void _init_root(self, DTYPE_t[3] min_bounds, DTYPE_t[3] max_bounds
                         ) nogil:
        """Initialize the root node with the given space boundaries"""
        cdef:
            int i
            DTYPE_t width
            Cell* root = &self.cells[0]

        self._init_cell(root, -1, 0)
        for i in range(self.n_dimensions):
            root.min_bounds[i] = min_bounds[i]
            root.max_bounds[i] = max_bounds[i]
            root.center[i] = (max_bounds[i] + min_bounds[i]) / 2.
            width = max_bounds[i] - min_bounds[i]
            root.squared_max_width = max(root.squared_max_width, width*width)
        root.cell_id = 0

        self.cell_count += 1

    cdef int _check_point_in_cell(self, DTYPE_t[3] point, Cell* cell
                                  ) nogil except -1:
        """Check that the given point is in the cell boundaries."""

        if self.verbose >= 50:
            if self.n_dimensions == 3:
                printf("[QuadTree] Checking point (%f, %f, %f) in cell %li "
                        "([%f/%f, %f/%f, %f/%f], size %li)\n",
                        point[0], point[1], point[2], cell.cell_id,
                        cell.min_bounds[0], cell.max_bounds[0], cell.min_bounds[1],
                        cell.max_bounds[1], cell.min_bounds[2], cell.max_bounds[2],
                        cell.cumulative_size)
            else:
                printf("[QuadTree] Checking point (%f, %f) in cell %li "
                        "([%f/%f, %f/%f], size %li)\n",
                        point[0], point[1],cell.cell_id, cell.min_bounds[0],
                        cell.max_bounds[0], cell.min_bounds[1],
                        cell.max_bounds[1], cell.cumulative_size)

        for i in range(self.n_dimensions):
            if (cell.min_bounds[i] > point[i] or
                    cell.max_bounds[i] <= point[i]):
                with gil:
                    msg = "[QuadTree] InsertionError: point out of cell "
                    msg += "boundary.\nAxis %li: cell [%f, %f]; point %f\n"

                    msg %= i, cell.min_bounds[i],  cell.max_bounds[i], point[i]
                    raise ValueError(msg)

    def _check_coherence(self):
        """Check the coherence of the cells of the tree.

        Check that the info stored in each cell is compatible with the info
        stored in descendent and sibling cells. Raise a ValueError if this
        fails.
        """
        for cell in self.cells[:self.cell_count]:
            # Check that the barycenter of inserted point is within the cell
            # boundaries
            self._check_point_in_cell(cell.barycenter, &cell)

            if not cell.is_leaf:
                # Compute the number of point in children and compare with
                # its cummulative_size.
                n_points = 0
                for idx in range(self.n_cells_per_cell):
                    child_id = cell.children[idx]
                    if child_id != -1:
                        child = self.cells[child_id]
                        n_points += child.cumulative_size
                        assert child.cell_id == child_id, (
                            "Cell id not correctly initiliazed.")
                if n_points != cell.cumulative_size:
                    raise ValueError(
                        "Cell {} is incoherent. Size={} but found {} points "
                        "in children. ({})"
                        .format(cell.cell_id, cell.cumulative_size,
                                n_points, cell.children))

        # Make sure that the number of point in the tree correspond to the
        # cummulative size in root cell.
        if self.n_points != self.cells[0].cumulative_size:
            raise ValueError(
                "QuadTree is incoherent. Size={} but found {} points "
                "in children."
                .format(self.n_points, self.cells[0].cumulative_size))

    cdef long summarize(self, DTYPE_t[3] point, DTYPE_t* results,
                        float squared_theta=.5, SIZE_t cell_id=0, long idx=0
                        ) nogil:
        """Summarize the tree compared to a query point.

        Input arguments
        ---------------
        point : array (n_dimensions)
             query point to construct the summary.
        cell_id : integer, optional (default: 0)
            current cell of the tree summarized. This should be set to 0 for
            external calls.
        idx : integer, optional (default: 0)
            current index in the result array. This should be set to 0 for
            external calls
        squared_theta: float, optional (default: .5)
            threshold to decide whether the node is sufficiently far
            from the query point to be a good summary. The formula is such that
            the node is a summary if
                node_width^2 / dist_node_point^2 < squared_theta.
            Note that the argument should be passed as theta^2 to avoid
            computing square roots of the distances.

        Output arguments
        ----------------
        results : array (n_samples * (n_dimensions+2))
            result will contain a summary of the tree information compared to
            the query point:
            - results[idx:idx+n_dimensions] contains the coordinate-wise
                difference between the query point and the summary cell idx.
                This is useful in t-SNE to compute the negative forces.
            - result[idx+n_dimensions+1] contains the squared euclidean
                distance to the summary cell idx.
            - result[idx+n_dimensions+2] contains the number of point of the
                tree contained in the summary cell idx.

        Return
        ------
        idx : integer
            number of elements in the results array.
        """
        cdef:
            int i, idx_d = idx + self.n_dimensions
            bint duplicate = True
            Cell* cell = &self.cells[cell_id]

        results[idx_d] = 0.
        for i in range(self.n_dimensions):
            results[idx + i] = point[i] - cell.barycenter[i]
            results[idx_d] += results[idx + i] * results[idx + i]
            duplicate &= fabsf(results[idx + i]) <= EPSILON

        # Do not compute self interactions
        if duplicate and cell.is_leaf:
            return idx

        # Check whether we can use this node as a summary
        # It's a summary node if the angular size as measured from the point
        # is relatively small (w.r.t. to theta) or if it is a leaf node.
        # If it can be summarized, we use the cell center of mass
        # Otherwise, we go a higher level of resolution and into the leaves.
        if cell.is_leaf or (
                (cell.squared_max_width / results[idx_d]) < squared_theta):
            results[idx_d + 1] = <DTYPE_t> cell.cumulative_size
            return idx + self.n_dimensions + 2

        else:
            # Recursively compute the summary in nodes
            for c in range(self.n_cells_per_cell):
                child_id = cell.children[c]
                if child_id != -1:
                    idx = self.summarize(point, results, squared_theta,
                                         child_id, idx)

        return idx

    def get_cell(self, point):
        """return the id of the cell containing the query point or raise
        ValueError if the point is not in the tree
        """
        cdef DTYPE_t[3] query_pt
        cdef int i

        assert len(point) == self.n_dimensions, (
            "Query point should be a point in dimension {}."
            .format(self.n_dimensions))

        for i in range(self.n_dimensions):
            query_pt[i] = point[i]

        return self._get_cell(query_pt, 0)

    cdef int _get_cell(self, DTYPE_t[3] point, SIZE_t cell_id=0
                       ) nogil except -1:
        """guts of get_cell.

        Return the id of the cell containing the query point or raise ValueError
        if the point is not in the tree"""
        cdef:
            SIZE_t selected_child
            Cell* cell = &self.cells[cell_id]

        if cell.is_leaf:
            if self._is_duplicate(cell.barycenter, point):
                if self.verbose > 99:
                    printf("[QuadTree] Found point in cell: %li\n",
                           cell.cell_id)
                return cell_id
            with gil:
                raise ValueError("Query point not in the Tree.")

        selected_child = self._select_child(point, cell)
        if selected_child > 0:
            if self.verbose > 99:
                printf("[QuadTree] Selected_child: %li\n", selected_child)
            return self._get_cell(point, selected_child)
        with gil:
            raise ValueError("Query point not in the Tree.")

    # Pickling primitives

    def __reduce__(self):
        """Reduce re-implementation, for pickling."""
        return (_QuadTree, (self.n_dimensions, self.verbose),
                           self.__getstate__())

    def __getstate__(self):
        """Getstate re-implementation, for pickling."""
        d = {}
        # capacity is inferred during the __setstate__ using nodes
        d["max_depth"] = self.max_depth
        d["cell_count"] = self.cell_count
        d["capacity"] = self.capacity
        d["n_points"] = self.n_points
        d["cells"] = self._get_cell_ndarray()
        return d

    def __setstate__(self, d):
        """Setstate re-implementation, for unpickling."""
        self.max_depth = d["max_depth"]
        self.cell_count = d["cell_count"]
        self.capacity = d["capacity"]
        self.n_points = d["n_points"]

        if 'cells' not in d:
            raise ValueError('You have loaded Tree version which '
                             'cannot be imported')

        cell_ndarray = d['cells']

        if (cell_ndarray.ndim != 1 or
                cell_ndarray.dtype != CELL_DTYPE or
                not cell_ndarray.flags.c_contiguous):
            raise ValueError('Did not recognise loaded array layout')

        self.capacity = cell_ndarray.shape[0]
        if self._resize_c(self.capacity) != 0:
            raise MemoryError("resizing tree to %d" % self.capacity)

        cells = memcpy(self.cells, (<np.ndarray> cell_ndarray).data,
                       self.capacity * sizeof(Cell))


    # Array manipulation methods, to convert it to numpy or to resize
    # self.cells array

    cdef np.ndarray _get_cell_ndarray(self):
        """Wraps nodes as a NumPy struct array.

        The array keeps a reference to this Tree, which manages the underlying
        memory. Individual fields are publicly accessible as properties of the
        Tree.
        """
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.cell_count
        cdef np.npy_intp strides[1]
        strides[0] = sizeof(Cell)
        cdef np.ndarray arr
        Py_INCREF(CELL_DTYPE)
        arr = PyArray_NewFromDescr(<PyTypeObject *> np.ndarray,
                                   CELL_DTYPE, 1, shape,
                                   strides, <void*> self.cells,
                                   np.NPY_DEFAULT, None)
        Py_INCREF(self)
        if PyArray_SetBaseObject(arr, <PyObject*> self) < 0:
            raise ValueError("Can't intialize array!")
        return arr

    cdef int _resize(self, SIZE_t capacity) nogil except -1:
        """Resize all inner arrays to `capacity`, if `capacity` == -1, then
           double the size of the inner arrays.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if self._resize_c(capacity) != 0:
            # Acquire gil only if we need to raise
            with gil:
                raise MemoryError()

    cdef int _resize_c(self, SIZE_t capacity=SIZE_MAX) nogil except -1:
        """Guts of _resize

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if capacity == self.capacity and self.cells != NULL:
            return 0

        if capacity == SIZE_MAX:
            if self.capacity == 0:
                capacity = 9  # default initial value to min
            else:
                capacity = 2 * self.capacity

        safe_realloc(&self.cells, capacity)

        # if capacity smaller than cell_count, adjust the counter
        if capacity < self.cell_count:
            self.cell_count = capacity

        self.capacity = capacity
        return 0

    def _py_summarize(self, DTYPE_t[:] query_pt, DTYPE_t[:, :] X, float angle):
        # Used for testing summarize
        cdef:
            DTYPE_t[:] summary
            int n_samples, n_dimensions

        n_samples = X.shape[0]
        n_dimensions = X.shape[1]
        summary = np.empty(4 * n_samples, dtype=np.float32)

        idx = self.summarize(&query_pt[0], &summary[0], angle * angle)
        return idx, summary

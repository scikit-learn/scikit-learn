# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# Author: Christopher Moody <chrisemoody@gmail.com>
# Author: Nick Travers <nickt@squareup.com>
# Implementation by Chris Moody & Nick Travers
# See http://homepage.tudelft.nl/19j49/t-SNE.html for reference
# implementations and papers describing the technique


from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from libc.math cimport sqrt, log
cimport numpy as np
import numpy as np

cdef char* EMPTY_STRING = ""

cdef extern from "math.h":
    float fabsf(float x) nogil

# Round points differing by less than this amount
# effectively ignoring differences near the 32bit 
# floating point precision
cdef float EPSILON = 1e-6

# This is effectively an ifdef statement in Cython
# It allows us to write printf debugging lines
# and remove them at compile time
cdef enum:
    DEBUGFLAG = 0

cdef extern from "time.h":
    # Declare only what is necessary from `tm` structure.
    ctypedef long clock_t
    clock_t clock() nogil
    double CLOCKS_PER_SEC


cdef extern from "cblas.h":
    float snrm2 "cblas_snrm2"(int N, float *X, int incX) nogil


cdef struct Node:
    # Keep track of the center of mass
    float* barycenter
    # If this is a leaf, the position of the point within this leaf 
    float* leaf_point_position
    # The number of points including all 
    # nodes below this one
    long cumulative_size
    # Number of points at this node
    long size
    # Index of the point at this node
    # Only defined for non-empty leaf nodes
    long point_index
    # level = 0 is the root node
    # And each subdivision adds 1 to the level
    long level
    # Left edge of this node
    float* left_edge
    # The center of this node, equal to le + w/2.0
    float* center
    # The width of this node -- used to calculate the opening
    # angle. Equal to width = re - le
    float* width
    # The value of the maximum width w
    float max_width

    # Does this node have children?
    # Default to leaf until we add points
    int is_leaf
    # Array of pointers to pointers of children
    Node **children
    # Keep a pointer to the parent
    Node *parent
    # Pointer to the tree this node belongs too
    Tree* tree

cdef struct Tree:
    # Holds a pointer to the root node
    Node* root_node 
    # Number of dimensions in the ouput
    int n_dimensions
    # Total number of cells
    long n_cells
    # Total number of points
    long n_points
    # Spit out diagnostic information?
    int verbose
    # How many cells per node? Should go as 2 ** n_dimensionss
    int n_cell_per_node

cdef Tree* init_tree(float[:] left_edge, float[:] width, int n_dimensions, 
                     int verbose) nogil:
    # tree is freed by free_tree
    cdef Tree* tree = <Tree*> malloc(sizeof(Tree))
    tree.n_dimensions = n_dimensions
    tree.n_cells = 0
    tree.n_points = 0
    tree.verbose = verbose
    tree.root_node = create_root(left_edge, width, n_dimensions)
    tree.root_node.tree = tree
    tree.n_cells += 1
    tree.n_cell_per_node = 2 ** n_dimensions
    if DEBUGFLAG:
        printf("[t-SNE] Tree initialised. Left_edge = (%1.9e, %1.9e, %1.9e)\n",
               left_edge[0], left_edge[1], left_edge[2])
        printf("[t-SNE] Tree initialised. Width = (%1.9e, %1.9e, %1.9e)\n",
                width[0], width[1], width[2])
    return tree

cdef Node* create_root(float[:] left_edge, float[:] width, int n_dimensions) nogil:
    # Create a default root node
    cdef int ax
    cdef int n_cell_per_node = 2 ** n_dimensions
    # root is freed by free_tree
    root = <Node*> malloc(sizeof(Node))
    root.is_leaf = 1
    root.parent = NULL
    root.level = 0
    root.cumulative_size = 0
    root.size = 0
    root.point_index = -1
    root.max_width = 0.0
    root.width = <float*> malloc(sizeof(float) * n_dimensions)
    root.left_edge = <float*> malloc(sizeof(float) * n_dimensions)
    root.center = <float*> malloc(sizeof(float) * n_dimensions)
    root.barycenter = <float*> malloc(sizeof(float) * n_dimensions)
    root.leaf_point_position= <float*> malloc(sizeof(float) * n_dimensions)
    root.children = NULL
    for ax in range(n_dimensions):
        root.width[ax] = width[ax]
        root.left_edge[ax] = left_edge[ax]
        root.center[ax] = 0.0
        root.barycenter[ax] = 0.
        root.leaf_point_position[ax] = -1
    for ax in range(n_dimensions):
        root.max_width = max(root.max_width, root.width[ax])
    if DEBUGFLAG:
        printf("[t-SNE] Created root node %p\n", root)
    return root

cdef Node* create_child(Node *parent, int[3] offset) nogil:
    # Create a new child node with default parameters
    cdef int ax
    # these children are freed by free_recursive
    child = <Node *> malloc(sizeof(Node))
    child.is_leaf = 1
    child.parent = parent
    child.level = parent.level + 1
    child.size = 0
    child.cumulative_size = 0
    child.point_index = -1
    child.tree = parent.tree
    child.max_width = 0.0
    child.width = <float*> malloc(sizeof(float) * parent.tree.n_dimensions)
    child.left_edge = <float*> malloc(sizeof(float) * parent.tree.n_dimensions)
    child.center = <float*> malloc(sizeof(float) * parent.tree.n_dimensions)
    child.barycenter = <float*> malloc(sizeof(float) * parent.tree.n_dimensions)
    child.leaf_point_position = <float*> malloc(sizeof(float) * parent.tree.n_dimensions)
    child.children = NULL
    for ax in range(parent.tree.n_dimensions):
        child.width[ax] = parent.width[ax] / 2.0
        child.left_edge[ax] = parent.left_edge[ax] + offset[ax] * parent.width[ax] / 2.0
        child.center[ax] = child.left_edge[ax] + child.width[ax] / 2.0
        child.barycenter[ax] = 0.
        child.leaf_point_position[ax] = -1.
    for ax in range(parent.tree.n_dimensions):
        child.max_width = max(child.max_width, child.width[ax])
    child.tree.n_cells += 1
    return child

cdef Node* select_child(Node *node, float[3] pos, long index) nogil:
    # Find which sub-node a position should go into
    # And return the appropriate node
    cdef int* offset = <int*> malloc(sizeof(int) * node.tree.n_dimensions)
    cdef int ax, idx
    cdef Node* child
    cdef int error
    for ax in range(node.tree.n_dimensions):
        offset[ax] = (pos[ax] - (node.left_edge[ax] + node.width[ax] / 2.0)) > 0.
    idx = offset2index(offset, node.tree.n_dimensions)
    child = node.children[idx]
    if DEBUGFLAG:
        printf("[t-SNE] Offset [%i, %i] with LE [%f, %f]\n",
               offset[0], offset[1], child.left_edge[0], child.left_edge[1])
    free(offset)
    return child


cdef inline void index2offset(int* offset, int index, int n_dimensions) nogil:
    # Convert a 1D index into N-D index; useful for indexing
    # children of a quadtree, octree, N-tree
    # Quite likely there's a fancy bitshift way of doing this
    # since the offset is equivalent to the binary representation
    # of the integer index
    # We read the offset array left-to-right
    # such that the least significat bit is on the right
    cdef int rem, k, shift
    for k in range(n_dimensions):
        shift = n_dimensions -k -1
        rem = ((index >> shift) << shift)
        offset[k] = rem > 0
        if DEBUGFLAG:
            printf("i2o index %i k %i rem %i offset", index, k, rem)
            for j in range(n_dimensions):
                printf(" %i", offset[j])
            printf(" n_dimensions %i\n", n_dimensions)
        index -= rem


cdef inline int offset2index(int* offset, int n_dimensions) nogil:
    # Calculate the 1:1 index for a given offset array
    # We read the offset array right-to-left
    # such that the least significat bit is on the right
    cdef int dim
    cdef int index = 0
    for dim in range(n_dimensions):
        index += (2 ** dim) * offset[n_dimensions - dim - 1]
        if DEBUGFLAG:
            printf("o2i index %i dim %i            offset", index, dim)
            for j in range(n_dimensions):
                printf(" %i", offset[j])
            printf(" n_dimensions %i\n", n_dimensions)
    return index


cdef void subdivide(Node* node) nogil:
    # This instantiates 2**n_dimensions = n_cell_per_node nodes for the current node
    cdef int idx = 0
    cdef int* offset = <int*> malloc(sizeof(int) * node.tree.n_dimensions)
    node.is_leaf = False
    node.children = <Node**> malloc(sizeof(Node*) * node.tree.n_cell_per_node)
    for idx in range(node.tree.n_cell_per_node):
        index2offset(offset, idx, node.tree.n_dimensions)
        node.children[idx] = create_child(node, offset)
    free(offset)


cdef int insert(Node *root, float pos[3], long point_index, long depth, long
        duplicate_count) nogil:
    # Introduce a new point into the tree
    # by recursively inserting it and subdividng as necessary
    # Carefully treat the case of identical points at the same node
    # by increasing the root.size and tracking duplicate_count
    cdef Node *child
    cdef long i
    cdef int ax
    cdef int not_identical = 1
    cdef int n_dimensions = root.tree.n_dimensions
    if DEBUGFLAG:
        printf("[t-SNE] [d=%i] Inserting pos %i [%f, %f] duplicate_count=%i "
                "into child %p\n", depth, point_index, pos[0], pos[1],
                duplicate_count, root)    
    # Increment the total number points including this
    # node and below it
    root.cumulative_size += duplicate_count
    # Evaluate the new center of mass, weighting the previous
    # center of mass against the new point data
    cdef double frac_seen = <double>(root.cumulative_size - 1) / (<double>
            root.cumulative_size)
    cdef double frac_new  = 1.0 / <double> root.cumulative_size
    # Assert that duplicate_count > 0
    if duplicate_count < 1:
        return -1
    # Assert that the point is inside the left & right edges
    for ax in range(n_dimensions):
        root.barycenter[ax] *= frac_seen
        if (pos[ax] > (root.left_edge[ax] + root.width[ax] + EPSILON)):
            printf("[t-SNE] Error: point (%1.9e) is above right edge of node "
                    "(%1.9e)\n", pos[ax], root.left_edge[ax] + root.width[ax])
            return -1
        if (pos[ax] < root.left_edge[ax] - EPSILON):
            printf("[t-SNE] Error: point (%1.9e) is below left edge of node "
                   "(%1.9e)\n", pos[ax], root.left_edge[ax])
            return -1
    for ax in range(n_dimensions):
        root.barycenter[ax] += pos[ax] * frac_new

    # If this node is unoccupied, fill it.
    # Otherwise, we need to insert recursively.
    # Two insertion scenarios: 
    # 1) Insert into this node if it is a leaf and empty
    # 2) Subdivide this node if it is currently occupied
    if (root.size == 0) & root.is_leaf:
        # Root node is empty and a leaf
        if DEBUGFLAG:
            printf("[t-SNE] [d=%i] Inserting [%f, %f] into blank cell\n", depth,
                   pos[0], pos[1])
        for ax in range(n_dimensions):
            root.leaf_point_position[ax] = pos[ax]
        root.point_index = point_index
        root.size = duplicate_count
        return 0
    else:
        # Root node is occupied or not a leaf
        if DEBUGFLAG:
            printf("[t-SNE] [d=%i] Node %p is occupied or is a leaf.\n", depth,
                    root)
            printf("[t-SNE] [d=%i] Node %p leaf = %i. Size %i\n", depth, root,
                    root.is_leaf, root.size)
        if root.is_leaf & (root.size > 0):
            # is a leaf node and is occupied
            for ax in range(n_dimensions):
                not_identical &= (fabsf(pos[ax] - root.leaf_point_position[ax]) < EPSILON)
                not_identical &= (root.point_index != point_index)
            if not_identical == 1:
                root.size += duplicate_count
                if DEBUGFLAG:
                    printf("[t-SNE] Warning: [d=%i] Detected identical "
                            "points. Returning. Leaf now has size %i\n",
                            depth, root.size)
                return 0
        # If necessary, subdivide this node before
        # descending
        if root.is_leaf:
            if DEBUGFLAG:
                printf("[t-SNE] [d=%i] Subdividing this leaf node %p\n", depth,
                        root)
            subdivide(root)
        # We have two points to relocate: the one previously
        # at this node, and the new one we're attempting
        # to insert
        if root.size > 0:
            child = select_child(root, root.leaf_point_position, root.point_index)
            if DEBUGFLAG:
                printf("[t-SNE] [d=%i] Relocating old point to node %p\n",
                        depth, child)
            insert(child, root.leaf_point_position, root.point_index, depth + 1, root.size)
        # Insert the new point
        if DEBUGFLAG:
            printf("[t-SNE] [d=%i] Selecting node for new point\n", depth)
        child = select_child(root, pos, point_index)
        if root.size > 0:
            # Remove the point from this node
            for ax in range(n_dimensions):
                root.leaf_point_position[ax] = -1            
            root.size = 0
            root.point_index = -1            
        return insert(child, pos, point_index, depth + 1, 1)

cdef int insert_many(Tree* tree, float[:,:] pos_array) nogil:
    # Insert each data point into the tree one at a time
    cdef long nrows = pos_array.shape[0]
    cdef long i
    cdef int ax
    cdef float row[3]
    cdef long err = 0
    for i in range(nrows):
        for ax in range(tree.n_dimensions):
            row[ax] = pos_array[i, ax]
        if DEBUGFLAG:
            printf("[t-SNE] inserting point %i: [%f, %f]\n", i, row[0], row[1])
        err = insert(tree.root_node, row, i, 0, 1)
        if err != 0:
            printf("[t-SNE] ERROR\n%s", EMPTY_STRING)
            return err
        tree.n_points += 1
    return err

cdef int free_tree(Tree* tree) nogil:
    cdef int check
    cdef long* cnt = <long*> malloc(sizeof(long) * 3)
    for i in range(3):
        cnt[i] = 0
    free_recursive(tree, tree.root_node, cnt)
    check = cnt[0] == tree.n_cells
    check &= cnt[2] == tree.n_points
    free(tree)
    free(cnt)
    return check

cdef void free_post_children(Node *node) nogil:
    free(node.width)
    free(node.left_edge)
    free(node.center)
    free(node.barycenter)
    free(node.leaf_point_position)
    free(node)

cdef void free_recursive(Tree* tree, Node *root, long* counts) nogil:
    # Free up all of the tree nodes recursively
    # while counting the number of nodes visited
    # and total number of data points removed
    cdef int idx
    cdef Node* child
    if not root.is_leaf:
        for idx in range(tree.n_cell_per_node):
            child = root.children[idx]
            free_recursive(tree, child, counts)
            counts[0] += 1
            if child.is_leaf:
                counts[1] += 1
                if child.size > 0:
                    counts[2] +=1
            else:
                free(child.children)

            free_post_children(child)

    if root == tree.root_node:
        if not root.is_leaf:
            free(root.children)

        free_post_children(root)

cdef long count_points(Node* root, long count) nogil:
    # Walk through the whole tree and count the number 
    # of points at the leaf nodes
    if DEBUGFLAG:
        printf("[t-SNE] Counting nodes at root node %p\n", root)
    cdef Node* child
    cdef int idx
    if root.is_leaf:
        count += root.size
        if DEBUGFLAG : 
            printf("[t-SNE] %p is a leaf node, no children\n", root)
            printf("[t-SNE] %i points in node %p\n", count, root)
        return count
    # Otherwise, get the children
    for idx in range(root.tree.n_cell_per_node):
        child = root.children[idx]
        if DEBUGFLAG:
            printf("[t-SNE] Counting points for child %p\n", child)
        if child.is_leaf and child.size > 0:
            if DEBUGFLAG:
                printf("[t-SNE] Child has size %d\n", child.size)
            count += child.size
        elif not child.is_leaf:
            if DEBUGFLAG:
                printf("[t-SNE] Child is not a leaf. Descending\n%s", EMPTY_STRING)
            count = count_points(child, count)
        # else case is we have an empty leaf node
        # which happens when we create a quadtree for
        # one point, and then the other neighboring cells
        # don't get filled in
    if DEBUGFLAG:
        printf("[t-SNE] %i points in this node\n", count)
    return count


cdef float compute_gradient(float[:,:] val_P,
                            float[:,:] pos_reference,
                            np.int64_t[:,:] neighbors,
                            float[:,:] tot_force,
                            Node* root_node,
                            float theta,
                            float dof,
                            long start,
                            long stop) nogil:
    # Having created the tree, calculate the gradient
    # in two components, the positive and negative forces
    cdef long i, coord
    cdef int ax
    cdef long n = pos_reference.shape[0]
    cdef int n_dimensions = root_node.tree.n_dimensions
    if root_node.tree.verbose > 11:
        printf("[t-SNE] Allocating %i elements in force arrays\n",
                n * n_dimensions * 2)
    cdef float* sum_Q = <float*> malloc(sizeof(float))
    cdef float* neg_f = <float*> malloc(sizeof(float) * n * n_dimensions)
    cdef float* neg_f_fast = <float*> malloc(sizeof(float) * n * n_dimensions)
    cdef float* pos_f = <float*> malloc(sizeof(float) * n * n_dimensions)
    cdef clock_t t1, t2
    cdef float sQ, error

    sum_Q[0] = 0.0
    t1 = clock()
    compute_gradient_negative(val_P, pos_reference, neg_f, root_node, sum_Q,
                              dof, theta, start, stop)
    t2 = clock()
    if root_node.tree.verbose > 15:
        printf("[t-SNE] Computing negative gradient: %e ticks\n", ((float) (t2 - t1)))
    sQ = sum_Q[0]
    t1 = clock()
    error = compute_gradient_positive(val_P, pos_reference, neighbors, pos_f,
                              n_dimensions, dof, sQ, start, root_node.tree.verbose)
    t2 = clock()
    if root_node.tree.verbose > 15:
        printf("[t-SNE] Computing positive gradient: %e ticks\n", ((float) (t2 - t1)))
    for i in range(start, n):
        for ax in range(n_dimensions):
            coord = i * n_dimensions + ax
            tot_force[i, ax] = pos_f[coord] - (neg_f[coord] / sum_Q[0])
    free(sum_Q)
    free(neg_f)
    free(neg_f_fast)
    free(pos_f)
    return sQ


cdef float compute_gradient_positive(float[:,:] val_P,
                                     float[:,:] pos_reference,
                                     np.int64_t[:,:] neighbors,
                                     float* pos_f,
                                     int n_dimensions,
                                     float dof,
                                     float sum_Q,
                                     np.int64_t start,
                                     int verbose) nogil:
    # Sum over the following expression for i not equal to j
    # grad_i = p_ij (1 + ||y_i - y_j||^2)^-1 (y_i - y_j)
    # This is equivalent to compute_edge_forces in the authors' code
    # It just goes over the nearest neighbors instead of all the data points
    # (unlike the non-nearest neighbors version of `compute_gradient_positive')
    cdef:
        int ax
        long i, j, k
        long K = neighbors.shape[1]
        long n = val_P.shape[0]
        float[3] buff
        float D, Q, pij
        float C = 0.0
        float exponent = (dof + 1.0) / -2.0
    cdef clock_t t1, t2
    t1 = clock()
    for i in range(start, n):
        for ax in range(n_dimensions):
            pos_f[i * n_dimensions + ax] = 0.0
        for k in range(K):
            j = neighbors[i, k]
            # we don't need to exclude the i==j case since we've 
            # already thrown it out from the list of neighbors
            D = 0.0
            Q = 0.0
            pij = val_P[i, j]
            for ax in range(n_dimensions):
                buff[ax] = pos_reference[i, ax] - pos_reference[j, ax]
                D += buff[ax] ** 2.0  
            Q = (((1.0 + D) / dof) ** exponent)
            D = pij * Q
            Q /= sum_Q
            C += pij * log((pij + EPSILON) / (Q + EPSILON))
            for ax in range(n_dimensions):
                pos_f[i * n_dimensions + ax] += D * buff[ax]
    t2 = clock()
    dt = ((float) (t2 - t1))
    if verbose > 10:
        printf("[t-SNE] Computed error=%1.4f in %1.1e ticks\n", C, dt)
    return C



cdef void compute_gradient_negative(float[:,:] val_P, 
                                    float[:,:] pos_reference,
                                    float* neg_f,
                                    Node *root_node,
                                    float* sum_Q,
                                    float dof,
                                    float theta, 
                                    long start, 
                                    long stop) nogil:
    if stop == -1:
        stop = pos_reference.shape[0] 
    cdef:
        int ax
        long i, j
        long n = stop - start
        float* force
        float* iQ 
        float* pos
        float* dist2s
        long* sizes
        float* deltas
        long* l
        int n_dimensions = root_node.tree.n_dimensions
        float qijZ, mult
        long idx, 
        long dta = 0
        long dtb = 0
        clock_t t1, t2, t3
        float* neg_force

    iQ = <float*> malloc(sizeof(float))
    force = <float*> malloc(sizeof(float) * n_dimensions)
    pos = <float*> malloc(sizeof(float) * n_dimensions)
    dist2s = <float*> malloc(sizeof(float) * n)
    sizes = <long*> malloc(sizeof(long) * n)
    deltas = <float*> malloc(sizeof(float) * n * n_dimensions)
    l = <long*> malloc(sizeof(long))
    neg_force= <float*> malloc(sizeof(float) * n_dimensions)

    for i in range(start, stop):
        # Clear the arrays
        for ax in range(n_dimensions):
            force[ax] = 0.0
            neg_force[ax] = 0.0
            pos[ax] = pos_reference[i, ax]
        iQ[0] = 0.0
        l[0] = 0
        # Find which nodes are summarizing and collect their centers of mass
        # deltas, and sizes, into vectorized arrays
        t1 = clock()
        compute_non_edge_forces(root_node, theta, i, pos, force, dist2s,
                                     sizes, deltas, l)
        t2 = clock()
        # Compute the t-SNE negative force
        # for the digits dataset, walking the tree
        # is about 10-15x more expensive than the 
        # following for loop
        exponent = (dof + 1.0) / -2.0
        for j in range(l[0]):
            qijZ = ((1.0 + dist2s[j]) / dof) ** exponent
            sum_Q[0] += sizes[j] * qijZ
            mult = sizes[j] * qijZ * qijZ
            for ax in range(n_dimensions):
                idx = j * n_dimensions + ax
                neg_force[ax] += mult * deltas[idx]
        t3 = clock()
        for ax in range(n_dimensions):
            neg_f[i * n_dimensions + ax] = neg_force[ax]
        dta += t2 - t1
        dtb += t3 - t2
    if root_node.tree.verbose > 20:
        printf("[t-SNE] Tree: %i clock ticks | ", dta)
        printf("Force computation: %i clock ticks\n", dtb)
    free(iQ)
    free(force)
    free(pos)
    free(dist2s)
    free(sizes)
    free(deltas)
    free(l)
    free(neg_force)


cdef void compute_non_edge_forces(Node* node, 
                                  float theta,
                                  long point_index,
                                  float* pos,
                                  float* force,
                                  float* dist2s,
                                  long* sizes,
                                  float* deltas,
                                  long* l) nogil:
    # Compute the t-SNE force on the point in pos given by point_index
    cdef:
        Node* child
        int i, j
        int n_dimensions = node.tree.n_dimensions
        long idx, idx1
        float dist_check
    
    # There are no points below this node if cumulative_size == 0
    # so do not bother to calculate any force contributions
    # Also do not compute self-interactions
    if node.cumulative_size > 0 and not (node.is_leaf and (node.point_index ==
        point_index)):
        # Compute distance between node center of mass and the reference point
        # I've tried rewriting this in terms of BLAS functions, but it's about
        # 1.5x worse when we do so, probbaly because the vectors are small
        idx1 = l[0] * n_dimensions
        deltas[idx1] = pos[0] - node.barycenter[0]
        idx = idx1
        for i in range(1, n_dimensions):
            idx += 1
            deltas[idx] = pos[i] - node.barycenter[i] 
        # do np.sqrt(np.sum(deltas**2.0))
        dist2s[l[0]] = snrm2(n_dimensions, &deltas[idx1], 1)
        # Check whether we can use this node as a summary
        # It's a summary node if the angular size as measured from the point
        # is relatively small (w.r.t. to theta) or if it is a leaf node.
        # If it can be summarized, we use the cell center of mass 
        # Otherwise, we go a higher level of resolution and into the leaves.
        if node.is_leaf or ((node.max_width / dist2s[l[0]]) < theta):
            # Compute the t-SNE force between the reference point and the
            # current node
            sizes[l[0]] = node.cumulative_size
            dist2s[l[0]] = dist2s[l[0]] * dist2s[l[0]]
            l[0] += 1
        else:
            # Recursively apply Barnes-Hut to child nodes
            for idx in range(node.tree.n_cell_per_node):
                child = node.children[idx]
                if child.cumulative_size == 0: 
                    continue
                compute_non_edge_forces(child, theta,
                        point_index, pos, force, dist2s, sizes, deltas,
                        l)


cdef float compute_error(float[:, :] val_P,
                        float[:, :] pos_reference,
                        np.int64_t[:,:] neighbors,
                        float sum_Q,
                        int n_dimensions,
                        int verbose) nogil:
    cdef int i, j, ax
    cdef int I = neighbors.shape[0]
    cdef int K = neighbors.shape[1]
    cdef float pij, Q
    cdef float C = 0.0
    cdef clock_t t1, t2
    cdef float dt, delta
    t1 = clock()
    for i in range(I):
        for k in range(K):
            j = neighbors[i, k]
            pij = val_P[i, j]
            Q = 0.0
            for ax in range(n_dimensions):
                delta = (pos_reference[i, ax] - pos_reference[j, ax])
                Q += delta * delta
            Q = (1.0 / (sum_Q + Q * sum_Q))
            C += pij * log((pij + EPSILON) / (Q + EPSILON))
    t2 = clock()
    dt = ((float) (t2 - t1))
    if verbose > 10:
        printf("[t-SNE] Computed error=%1.4f in %1.1e ticks\n", C, dt)
    return C


def calculate_edge(pos_output):
    # Make the boundaries slightly outside of the data
    # to avoid floating point error near the edge
    left_edge = np.min(pos_output, axis=0)
    right_edge = np.max(pos_output, axis=0) 
    center = (right_edge + left_edge) * 0.5
    width = np.maximum(np.subtract(right_edge, left_edge), EPSILON)
    # Exagerate width to avoid boundary edge
    width = width.astype(np.float32) * 1.001
    left_edge = center - width / 2.0
    right_edge = center + width / 2.0
    return left_edge, right_edge, width

def gradient(float[:,:] pij_input, 
             float[:,:] pos_output, 
             np.int64_t[:,:] neighbors, 
             float[:,:] forces, 
             float theta,
             int n_dimensions,
             int verbose,
             float dof = 1.0,
             long skip_num_points=0):
    # This function is designed to be called from external Python
    # it passes the 'forces' array by reference and fills thats array
    # up in-place
    cdef float C
    n = pos_output.shape[0]
    left_edge, right_edge, width = calculate_edge(pos_output)
    assert width.itemsize == 4
    assert pij_input.itemsize == 4
    assert pos_output.itemsize == 4
    assert forces.itemsize == 4
    m = "Number of neighbors must be < # of points - 1"
    assert n - 1 >= neighbors.shape[1], m
    m = "neighbors array and pos_output shapes are incompatible"
    assert n == neighbors.shape[0], m
    m = "Forces array and pos_output shapes are incompatible"
    assert n == forces.shape[0], m
    m = "Pij and pos_output shapes are incompatible"
    assert n == pij_input.shape[0], m
    m = "Pij and pos_output shapes are incompatible"
    assert n == pij_input.shape[1], m
    if verbose > 10:
        printf("[t-SNE] Initializing tree of n_dimensions %i\n", n_dimensions)
    cdef Tree* qt = init_tree(left_edge, width, n_dimensions, verbose)
    if verbose > 10:
        printf("[t-SNE] Inserting %i points\n", pos_output.shape[0])
    err = insert_many(qt, pos_output)
    assert err == 0, "[t-SNE] Insertion failed"
    if verbose > 10:
        # XXX: format hack to workaround lack of `const char *` type
        # in the generated C code that triggers error with gcc 4.9
        # and -Werror=format-security
        printf("[t-SNE] Computing gradient\n%s", EMPTY_STRING)
    sum_Q = compute_gradient(pij_input, pos_output, neighbors, forces,
                             qt.root_node, theta, dof, skip_num_points, -1)
    C = compute_error(pij_input, pos_output, neighbors, sum_Q, n_dimensions,
                      verbose)
    if verbose > 10:
        # XXX: format hack to workaround lack of `const char *` type
        # in the generated C code
        # and -Werror=format-security
        printf("[t-SNE] Checking tree consistency\n%s", EMPTY_STRING)
    cdef long count = count_points(qt.root_node, 0)
    m = ("Tree consistency failed: unexpected number of points=%i "
         "at root node=%i" % (count, qt.root_node.cumulative_size))
    assert count == qt.root_node.cumulative_size, m 
    m = "Tree consistency failed: unexpected number of points on the tree"
    assert count == qt.n_points, m
    free_tree(qt)
    return C


# Helper functions
def check_quadtree(X, np.int64_t[:] counts):
    """
    Helper function to access quadtree functions for testing
    """
    
    X = X.astype(np.float32)
    left_edge, right_edge, width = calculate_edge(X)
    # Initialise a tree
    qt = init_tree(left_edge, width, 2, 2)
    # Insert data into the tree
    insert_many(qt, X)

    cdef long count = count_points(qt.root_node, 0)
    counts[0] = count
    counts[1] = qt.root_node.cumulative_size
    counts[2] = qt.n_points
    free_tree(qt)
    return counts


cdef int helper_test_index2offset(int* check, int index, int n_dimensions):
    cdef int* offset = <int*> malloc(sizeof(int) * n_dimensions)
    cdef int error_check = 1
    for i in range(n_dimensions):
        offset[i] = 0
    index2offset(offset, index, n_dimensions)
    for i in range(n_dimensions):
        error_check &= offset[i] == check[i]
    free(offset)
    return error_check


def test_index2offset():
    ret = 1
    ret &= helper_test_index2offset([1, 0, 1], 5, 3) == 1
    ret &= helper_test_index2offset([0, 0, 0], 0, 3) == 1
    ret &= helper_test_index2offset([0, 0, 1], 1, 3) == 1
    ret &= helper_test_index2offset([0, 1, 0], 2, 3) == 1
    ret &= helper_test_index2offset([0, 1, 1], 3, 3) == 1
    ret &= helper_test_index2offset([1, 0, 0], 4, 3) == 1
    return ret


def test_index_offset():
    cdef int n_dimensions, idx, tidx, k
    cdef int error_check = 1
    cdef int* offset 
    for n_dimensions in range(2, 10):
        offset = <int*> malloc(sizeof(int) * n_dimensions)
        for k in range(n_dimensions):
            offset[k] = 0
        for idx in range(2 ** n_dimensions):
            index2offset(offset, idx, n_dimensions)
            tidx = offset2index(offset, n_dimensions)
            error_check &= tidx == idx
            assert error_check == 1
        free(offset)
    return error_check

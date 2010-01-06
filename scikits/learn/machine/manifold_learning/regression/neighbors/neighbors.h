/**
 * \file neighbors.h
 * Contains the declarations of the wrapper for the neighbor tree
 */

#ifndef NEIGHBOORS
#define NEIGHBOORS

# ifdef __cplusplus
extern "C"
{
# endif

#ifdef __GNUC__
#define neighboors_dll
#else
#ifdef neighboors_EXPORTS
#define neighboors_dll _declspec(dllexport)
#else
#define neighboors_dll _declspec(dllimport)
#endif
#endif

typedef unsigned long (*numpy_callback)(double, unsigned long);

/**
 * Allocates the tree and populates it
 * @param data is the data that populates the space to analyze
 * @param dim1 is the first dimension
 * @param dim2 is the second dimension
 * @param levels is the number of levels to use
 * @return a pointer to the new tree
 */
neighboors_dll Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* allocate_neighborer(double* data, unsigned long dim1, unsigned long dim2, unsigned long levels);

/**
 * Deletes the tree
 * @param tree is the tree to destroy
 * @return 0
 */
neighboors_dll int delete_neighborer(Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* tree);


/**
 * Finds the k-neighboors of a point
 * @param tree is the tree to use
 * @param data is the new point
 * @param dim1 is the size of the new point
 * @param neighbors is the number of neighboors to use
 * @param callback is a callback to the Python list of tuples containing the results
 */
neighboors_dll unsigned long find_kneighbors(Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* tree, double* data, unsigned long dim1, unsigned long neighboors, numpy_callback callback);


/**
 * Finds the neighboors of a point in a Window
 * @param tree is the tree to use
 * @param data is the new point
 * @param dim1 is the size of the new point
 * @param windowSize is the size of the window
 * @param callback is a callback to the Python list of tuples containing the results
 */
neighboors_dll unsigned long find_parzen(Neighbor::NeighborsFromMatrix<Matrix::Matrix<double> >* tree, double* data, unsigned long dim1, double windowSize, numpy_callback callback);

# ifdef __cplusplus
}
# endif

#endif

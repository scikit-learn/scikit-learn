/**
 * \file cost_function.h
 * Wrapper for ctypes for the modified cost function
 */

#ifndef COSTFUNCTION
#define COSTFUNCTION

# ifdef __cplusplus
extern "C"
{
# endif

#ifdef __GNUC__
#define cost_function_dll
#else
#ifdef cost_function_EXPORTS
#define cost_function_dll _declspec(dllexport)
#else
#define cost_function_dll _declspec(dllimport)
#endif
#endif

typedef double* (*numpy_allocator)(unsigned long, unsigned long*);
/**
 * Wrapper for allocating the instance
 * @param data is a pointer to the distance matrix
 * @param dim1 is the first dimension of the matrix
 * @param dim2 is the second dimension
 * @param nbCoords is the number of coordinates to get
 * @param epsilon is the epsilon factor added so that the function is derivable
 * @param sigma is the weight factor for noise filtering
 * @param x1 is the limit when difference between distances is to big and when qudratic cost is used
 * @return a new instance of the cost function
 */
cost_function_dll Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* allocate_cost_function(double* data, unsigned long dim1, unsigned long dim2, unsigned long nbCoords, double epsilon, double sigma, double x1);
/**
 * Wrapper for deleting the instance
 * @param function is the instance to delete
 */
cost_function_dll int delete_cost_function(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* function);
/**
 * Wrapper for calling the instance
 * @param function is the instance to call
 * @param data is the paramaters
 * @param dim1 is the dimension of the vector of parameters
 */
cost_function_dll double call_cost_function(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* function, double* data, unsigned long dim1);
/**
 * Wrapper for calling the gradient
 * @param function is the instance to call
 * @param data is the paramaters
 * @param dim1 is the dimension of the vector of parameters
 * @param allocator allocates and returns the new numpy array
 */
cost_function_dll void gradient_cost_function(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* function, double* data, unsigned long dim1, numpy_allocator allocator);

# ifdef __cplusplus
}
# endif

#endif

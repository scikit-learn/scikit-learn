#include <iostream>
#include "modifiedCompression.h"
#include "cost_function.h"

Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* allocate_cost_function(double* data, unsigned long dim1, unsigned long dim2, unsigned long nbCoords, double epsilon, double sigma, double x1)
{
  return new Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >(Matrix::Matrix<double>(dim1, dim2, data), nbCoords, epsilon, sigma, x1);
}

int delete_cost_function(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* function)
{
  delete function;
  return 0;
}

double call_cost_function(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* function, double* data, unsigned long dim1)
{
  return function->value(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >::ParameterType(dim1, 1, data));
}

void gradient_cost_function(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >* function, double* data, unsigned long dim1, numpy_allocator allocator)
{
  double* array = allocator(1, &dim1);
  Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >::ParameterType gradient = function->gradient(Reconstruction::EnhancedCompressionFunction<Matrix::Matrix<double> >::ParameterType(dim1, 1, data));
  std::copy(gradient.begin(), gradient.end(), array);
}

/**
 * \file modifiedCompression.h
 * Our own compression method, a little bit modified
 */

// Matthieu Brucher
// Last Change : 2008-02-21 15:29

#ifndef MODIFIEDCOMPRESSIONMETHOD
#define MODIFIEDCOMPRESSIONMETHOD

#include <matrix/matrix_lib.h>
#include <matrix/sub_matrix_lib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>

namespace Reconstruction
{
  /// The enhanced compression cost function, derives from a stochastic function in case of someone needs it :)
  /**
   * This function provides a "static" filter for the distances, because the 2 last terms are based on the real distance instead of the other
   */
  template<class MatrixType>
  struct EnhancedCompressionFunction
  {
    typedef typename MatrixType::Data_Type DataType;
    enum{Size = MatrixType::staticHeight};
    /// Typedef of the parameters
    typedef typename Matrix::ColumnVector<DataType, Size>::Type ParameterType;
    /// Typedef of the parameters
    typedef const typename Matrix::ColumnVector<DataType, Size>::Type ConstParameterType;
    /**
     * The construction of the cost function
     * @param realDistance are the real distances that must be kept
     * @param dimension is the height of the compressed coordinates, the size of the hidden space of the points
     * @param epsilon is a small value for not having discontinuities at the origin
     * @param sigma is the variance of the underlying noise
     * @param x1 is an indicator to use quadratic convergence if the estimated distance is too different from the real distance
     */
    EnhancedCompressionFunction(const MatrixType& realDistance, unsigned long dimension, DataType epsilon, DataType sigma, DataType x1)
      :parameterSize(dimension), realDistance(realDistance), epsilon(epsilon), sigma(sigma), x1(x1)
    {
      std::cout << x1 << std::endl;
    }

    /**
     * Calculates the value of the function at a point
     * @param parameters are the paramaters at which the function is to be calculated
     * @return the value of the cost function
     */
    DataType value(const ParameterType& parameters)
    {
      this->parameters = parameters;

      DataType cost = DataTypeTraits<DataType>::zero(parameters(0,0));
//#pragma omp parallel for schedule(dynamic, 500) reduction(+:cost)
      for(long i = 0; i < realDistance.width(); ++i)
      {
        DataType tempCost = DataTypeTraits<DataType>::zero(parameters(0,0));
        for(long j = 0; j < i; ++j)
        {
          DataType redDistance(reductedDistance(i, j));
          DataType reaDistance(realDistance(i, j));
          DataType squareDifference = DataTypeTraits<DataType>::norm2(redDistance - reaDistance);

          if(reaDistance > 0)
            tempCost += DataTypeTraits<DataType>::sqrt(epsilon + squareDifference) *
              DataTypeTraits<DataType>::sqrt((x1 + squareDifference) * DataTypeTraits<DataType>::inverse(x1)) *
              reaDistance * DataTypeTraits<DataType>::inverse(sigma + reaDistance);
        }
        cost += tempCost;
      }
      return cost;
    }

    /**
     * Calculates the gradient of the function at a point for a piece of parameters
     * @param parameters are the paramaters at which the gradient of the function is to be calculated
     * @param chunkToUpdate is the piece number of parameters to update
     * @return the gradient of the cost function
     */
    const typename Matrix::ColumnVector<DataType, 0U>::Type gradient(const ParameterType& parameters, unsigned long chunkToUpdate)
    {
      this->parameters = parameters;
      return partialGrad(parameters, chunkToUpdate);
    }
    /**
     * Calculates the gradient of the function at a point
     * @param parameters are the paramaters at which the gradient of the function is to be calculated
     * @return the gradient of the cost function
     */
    const ParameterType gradient(const ParameterType& parameters)
    {
      this->parameters = parameters;
      ParameterType grad(parameters.height(), 1U);
      for(unsigned long chunkToUpdate = 0; chunkToUpdate < parameters.height() / parameterSize; ++chunkToUpdate)
      {
        Matrix::SubMatrix<ParameterType> innerGrad(grad, 0, parameterSize * chunkToUpdate, 1, parameterSize);
        innerGrad = partialGrad(parameters, chunkToUpdate);
      }
      return grad;
    }
    private:
    /**
     * Calculates the distance with a norm2 between two points in the reducted space
     * The coordinates of a point in the reducted space in considered to be the number of parameters the stochastic function changes at each update
     * @param xi is the first point
     * @param xj is the second point
     * @return the distance between the two points
     */
    DataType reductedDistance(unsigned long xi, unsigned long xj) const
    {
      return DataTypeTraits<DataType>::sqrt(norm2(Matrix::SubMatrix<const ParameterType>(parameters, 0, xi * parameterSize, 1, parameterSize) - Matrix::SubMatrix<const ParameterType>(parameters, 0, xj * parameterSize, 1, parameterSize)));
    }

    /**
     * Calculates the gradient for a point
     * @param parameters is the vector containing the coordinates of all the points
     * @param chunkToUpdate is the point whose gradient will be upgraded
     */
    const ParameterType partialGrad(const ParameterType& parameters, unsigned long chunkToUpdate)
    {
      ParameterType grad(parameterSize, 1U);
      Matrix::SubMatrix<ConstParameterType> innerParam(parameters, 0, chunkToUpdate * parameterSize, 1, parameterSize);

      ParameterType xk(innerParam);
      for(unsigned long i = 0; i < realDistance.width(); ++i)
      {
        if(i == chunkToUpdate)
          continue;

        DataType cost=0;
        DataType redDistance(reductedDistance(i, chunkToUpdate));
        DataType reaDistance(realDistance(i, chunkToUpdate));
        DataType squareDifference = DataTypeTraits<DataType>::norm2(redDistance - reaDistance);

        if(reaDistance > 0)
        {
          cost = (redDistance - reaDistance) * DataTypeTraits<DataType>::inverse(DataTypeTraits<DataType>::sqrt(epsilon + squareDifference) * redDistance) *
            DataTypeTraits<DataType>::sqrt((x1 + squareDifference) * DataTypeTraits<DataType>::inverse(x1)) *
            reaDistance * DataTypeTraits<DataType>::inverse(sigma + reaDistance);

          cost += DataTypeTraits<DataType>::sqrt(epsilon + squareDifference) *
            (redDistance - reaDistance) * DataTypeTraits<DataType>::inverse(DataTypeTraits<DataType>::sqrt((x1 + squareDifference)) * DataTypeTraits<DataType>::inverse(x1) * redDistance) *
            reaDistance * DataTypeTraits<DataType>::inverse(sigma + reaDistance);
        }

        innerParam.setColumnOffset(i * parameterSize);
        ParameterType xi(innerParam);
        grad += (xk - xi) * cost;
      }

      return grad;
    }

    /// Size of a point in the reducted space
    unsigned long parameterSize;
    /// The parameters that are used
    ParameterType parameters;
    /// The real points that must be compressed
    typename MatrixType::Result realDistance;
    /// The number used to assess non nullity of th square root and its derivate
    DataType epsilon;
    /// Kind of variance of the noise, under this number, a distance is not really taken in account
    DataType sigma;
    /// Kind of horizon, over this number, a distance is not really taken in account
    DataType x1;
  };
}
#endif

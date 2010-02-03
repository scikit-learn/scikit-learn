/**
 * \file matrix_traits.h
 * This files contains the definitions of zero, one for each datatype
 */

#ifndef MATRIX_TRAITS
#define MATRIX_TRAITS

#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <matrix/type_traits.h>

namespace Matrix
{
  template<class DataType, unsigned long Height = 0U, unsigned long Width = 0U>
  struct Matrix;
  template<class MatrixType>
  struct SubMatrix;
  /// Workaround for template typedef of ColumVectors :( waiting for Cxx :(
  template<class DataType, unsigned long Size = 0U>
  struct ColumnVector
  {
    typedef Matrix<DataType, Size, 1U> Type;
  };
  /// Workaround for template typedef of LineVectors :( waiting for Cxx :(
  template<class DataType, unsigned long Size = 0U>
  struct LineVector
  {
    typedef Matrix<DataType, 1U, Size> Type;
  };

  template<class MatrixType, template <class> class SVDFunc>
  struct SVD;
  template<class MatrixType>
  struct SVDHouseolderFunctions;
  template<class MatrixType>
  struct SVDHermitianFunctions;
}

/// This structure contains the value for zero, one and some other operation for an Matrix::Matrix type
template <class DataType, unsigned long Height, unsigned long Width>
struct DataTypeTraits<Matrix::Matrix<DataType, Height, Width> >
{
  /// The difference between 1 and the least value greater than 1 that is representable
  static const Matrix::Matrix<DataType, Height, Width> epsilon(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    return Matrix::Matrix<DataType, Height, Width>::DataTraits::epsilon(data(0, 0)) * Matrix::Matrix<DataType, Height, Width>(data.height(), data.width(), Matrix::Matrix<DataType, Height, Width>::DataTraits::one(data(0, 0)));
  }

  static Matrix::Matrix<DataType, Height, Width> zero(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    return Matrix::Matrix<DataType, Height, Width>(data.height(), data.height(), Matrix::Matrix<DataType, Height, Width>::DataTraits::zero(data(0, 0)));
  }

  static Matrix::Matrix<DataType, Height, Width> one(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    Matrix::Matrix<DataType, Height, Width> temp_mat(data.height(), data.width(), Matrix::Matrix<DataType, Height, Width>::DataTraits::zero(data(0, 0)));
    unsigned long size(std::min(data.height(), data.width()));
    for (unsigned long i = 0; i < size; i++)
      temp_mat(i,i) = Matrix::Matrix<DataType, Height, Width>::DataTraits::one(data(0, 0));
    return temp_mat;
  }

  static const unsigned long height(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    return data.height();
  }

  static const unsigned long width(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    return data.width();
  }

  static void clone(const Matrix::Matrix<DataType, Height, Width>& value, Matrix::Matrix<DataType, Height, Width>& resultat)
  {
    resultat.clone(value);
  }

  /**
   * This function calculates the absolute value of a Matrix::Matrix
   * @param data is the piece of data
   * @return the absolute value of the piece of data
   */
  static const Matrix::Matrix<DataType, Height, Width> absolute(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    unsigned long height = data.height();
    unsigned long width = data.width();
    Matrix::Matrix<DataType, Height, Width> abso(height, width);

    for(unsigned long i = 0; i < height; ++i)
    {
      for(unsigned long j = 0; j < width; ++j)
      {
        abso(i, j) = Matrix::Matrix<DataType, Height, Width>::DataTraits::absolute(data(i, j));
      }
    }
    return abso;
  }
  
  /**
   * This function conjugates a piece of data
   * @param data is the piece of data to be conjugated
   * @return the conjugated piece of data
   */
  static const Matrix::Matrix<DataType, Height, Width> conjuge(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    return transpose(data);
  }

  /**
   * This function inverses a piece of data
   * @param data is the piece of data to be conjugated
   * @return the inversed piece of data
   */
  static const Matrix::Matrix<DataType, Height, Width> inverse(const Matrix::Matrix<DataType, Height, Width>& data)
  {
    return inversion(data);
  }

  /**
   * This function returns the data needed to transpose the original data
   * @param data is the piece of data to be transposed
   * @return the anti-diagonal system
   */
  static Matrix::Matrix<DataType, Width, Height> create_transpose(const Matrix::Matrix<DataType, Height, Width>& data);
};

#endif

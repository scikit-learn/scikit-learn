/**
 * \file inversion.h
 * This file is there to include the correct headers for determinant and inversion calculus
 */

#ifndef INVERSION
#define INVERSION

#include <stdexcept>
#include <matrix/matrix_svd_lib.h>
#include <matrix/svd_houseolder.h>

namespace Matrix
{
  /**
   * Calculates the determinant of the square matrix
   * @return the determinant of the square matrix
   * @pre the matrix is square
   * \todo change the SVD part by a LU decomposition
   */
  template<class MatrixType>
  typename MatrixType::Data_Type determinant(const MatrixType& matrix)
  {
    assert(matrix.width() == matrix.height());
    typename MatrixType::Data_Type deter;
    switch(matrix.height())
    {
      case 1:
      {
        deter = matrix(0,0);
        break;
      }
      case 2:
      {
        deter = matrix(0,0) * matrix(1,1) - matrix(1,0) * matrix(0,1);
        break;
      }
      case 3:
      {
        deter = matrix(0,0) * (matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1))
            + matrix(1,0) * (matrix(2,1) * matrix(0,2) - matrix(0,1) * matrix(2,2))
            + matrix(2,0) * (matrix(0,1) * matrix(1,2) - matrix(1,1) * matrix(0,2));
        break;
      }
      default:
      {
        SVD<MatrixType, SVDHouseolderFunctions> det(matrix);
        det.biorthogonalization();
        deter = det.det();
      }
    }
    return (deter);
  }

  /**
   * Inverts a matrix <br>
   * The method is described <a href=doi:10.1016/S0377-0427(00)00613-0>here</a>
   * @return the inverse of the square matrix
   * @pre the matrix is square
   */
  template<class MatrixType>
  typename MatrixType::Self inversion(const MatrixType& matrix)
  {
    typedef typename MatrixType::DataTraits DataTraits;

    assert(matrix.width() == matrix.height());
    switch(matrix.width())
    {
      case 1:
      {
        return MatrixType(1, 1, DataTraits::inverse(matrix(0,0)));
      }
      case 2:
      {
        MatrixType temp_square(2, 2, DataTraits::zero(matrix(0,0)));
        typename MatrixType::Data_Type deter = DataTraits::inverse(determinant(matrix));
        temp_square(0, 0) = matrix(1, 1) * deter;
        temp_square(1, 1) = matrix(0, 0) * deter;
        temp_square(0, 1) = - matrix(0, 1) * deter;
        temp_square(1, 0) = - matrix(1, 0) * deter;
        return temp_square;
      }
      case 3:
      {
        MatrixType temp_square(3, 3, DataTraits::zero(matrix(0,0)));
        typename MatrixType::Data_Type deter = DataTraits::inverse(determinant(matrix));
        temp_square(0, 0) = (matrix(1, 1) * matrix(2,2) - matrix(1,2) * matrix(2,1)) * deter;
        temp_square(1, 0) =-(matrix(1, 0) * matrix(2,2) - matrix(2,0) * matrix(1,2)) * deter;
        temp_square(2, 0) = (matrix(1, 0) * matrix(2,1) - matrix(2,0) * matrix(1,1)) * deter;
        temp_square(0, 1) =-(matrix(0, 1) * matrix(2,2) - matrix(0,2) * matrix(2,1)) * deter;
        temp_square(1, 1) = (matrix(0, 0) * matrix(2,2) - matrix(0,2) * matrix(2,0)) * deter;
        temp_square(2, 1) =-(matrix(0, 0) * matrix(2,1) - matrix(0,1) * matrix(2,0)) * deter;
        temp_square(0, 2) = (matrix(0, 1) * matrix(1,2) - matrix(1,1) * matrix(0,2)) * deter;
        temp_square(1, 2) =-(matrix(0, 0) * matrix(1,2) - matrix(1,0) * matrix(0,2)) * deter;
        temp_square(2, 2) = (matrix(0, 0) * matrix(1,1) - matrix(0,1) * matrix(1,0)) * deter;
        return temp_square;
      }
      default:
      {
        Matrix<typename MatrixType::Data_Type, 0U, 0U> temp_square(prepare_inversion(matrix));

        for(unsigned long i = 0; i < matrix.height(); ++i)
        {
          unsigned long abs_pivot, ord_pivot;
          temp_square.max(matrix.width(), matrix.height(), (matrix.width() * 2 - i), (matrix.height() * 2 - i), &abs_pivot, &ord_pivot);
          if(temp_square(ord_pivot, abs_pivot) == DataTraits::zero(temp_square(ord_pivot, abs_pivot)))
            throw std::runtime_error("Division by zero error");

          typename MatrixType::Data_Type inverse(DataTraits::inverse(temp_square(ord_pivot, abs_pivot)));

          temp_square = Matrix<typename MatrixType::Data_Type>::leftPivotForSquareMatrix(temp_square, matrix.height(), matrix.width(), abs_pivot, ord_pivot, inverse) * temp_square * Matrix<typename MatrixType::Data_Type>::rightPivotForSquareMatrix(temp_square, matrix.height(), matrix.width(), abs_pivot, ord_pivot, inverse);
        }
        return MatrixType(temp_square);
      }
    }
  }

  template<class MatrixType>
  static Matrix<typename MatrixType::Data_Type, MatrixType::staticHeight * 2, MatrixType::staticWidth * 2> prepare_inversion(const MatrixType& matrix)
  {
    unsigned long width = matrix.width();
    unsigned long height = matrix.height();
    Matrix<typename MatrixType::Data_Type, MatrixType::staticHeight * 2, MatrixType::staticWidth * 2> m_inversion(height * 2, height * 2, matrix(0, 0));

    /* putting zero in the first quadrant */
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        m_inversion(i,j) = MatrixType::DataTraits::zero(matrix(0, 0));
      }
    }
    /* putting identity in the second and third quadrants */
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        m_inversion(i + height, j) =  MatrixType::DataTraits::zero(matrix(0, 0));
        m_inversion(i, j + width) =  MatrixType::DataTraits::zero(matrix(0, 0));
      }
      m_inversion(j + height, j) =  MatrixType::DataTraits::one(matrix(0, 0));
      m_inversion(j, j + width) =  MatrixType::DataTraits::one(matrix(0, 0));
    }
    /* putting matrix in the fourth quadrant */
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        m_inversion((i + height),(j + width)) = - matrix(i, j);
      }
    }
    return m_inversion;
  }

  /**
   * Calculates the pseudo inverse of a hermitian matrix
   * @param matrix is the matrix to inverse
   * @param epsilon is the precision of the pseudo inverse
   * @return the pseudoinverse
   * \pre matrix must be hermitian
   */
  template<class MatrixType>
  const typename MatrixType::Self pseudoinverse(const MatrixType& matrix, typename MatrixType::Data_Type epsilon = 0.)
  {
    typedef typename MatrixType::DataTraits DataTraits;
    if(epsilon == 0.)
      epsilon = DataTraits::epsilon(matrix(0,0));

    SVD<MatrixType, SVDHermitianFunctions> svd(matrix);
    svd.biorthogonalization();
    svd.diagonalization();
    svd.S(0, 0) = DataTraits::inverse(svd.S(0, 0));
    for(unsigned long i = 1; i < matrix.height(); ++i)
    {
      if(DataTraits::absolute(svd.S(i, i)) > svd.S(0, 0) * epsilon)
        svd.S(i, i) = DataTraits::inverse(svd.S(i, i));
      else
        svd.S(i, i) = DataTraits::zero(svd.S(i, i));
    }
    return typename MatrixType::Self(svd.U * svd.S * svd.V);
  }
  /**
   * Calculates the pseudo inverse of a matrix with the innerproduct
   * @param matrix is the matrix to inverse
   * @param epsilon is the precision of the pseudo inverse
   * @return the pseudoinverse
   */
  template<class MatrixType>
  const typename MatrixType::Transpose pseudoinverseI(const MatrixType& matrix, typename MatrixType::Data_Type epsilon = 0.)
  {
    typedef typename MatrixType::DataTraits DataTraits;
    if(epsilon == 0.)
      epsilon = DataTraits::epsilon(matrix(0,0));

    Matrix<typename MatrixType::Data_Type, MatrixType::staticWidth, MatrixType::staticHeight> innerProduct = inner(matrix);
    return pseudoinverse(innerProduct) * transpose(matrix);
  }

  /**
   * Calculates the pseudo inverse of a matrix with the outerproduct
   * @param matrix is the matrix to inverse
   * @param epsilon is the precision of the pseudo inverse
   * @return the pseudoinverse
   */
  template<class MatrixType>
  const typename MatrixType::Transpose pseudoinverseO(const MatrixType& matrix, typename MatrixType::Data_Type epsilon = 0.)
  {
    typedef typename MatrixType::DataTraits DataTraits;
    if(epsilon == 0.)
      epsilon = DataTraits::epsilon(matrix(0,0));

    Matrix<typename MatrixType::Data_Type, MatrixType::staticWidth, MatrixType::staticHeight> outerProduct = outer(matrix);
    return transpose(matrix) * pseudoinverse(outerProduct);
  }
}

#endif

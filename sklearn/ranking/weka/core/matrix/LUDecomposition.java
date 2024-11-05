/*
 * This software is a cooperative product of The MathWorks and the National
 * Institute of Standards and Technology (NIST) which has been released to the
 * public domain. Neither The MathWorks nor NIST assumes any responsibility
 * whatsoever for its use by other parties, and makes no guarantees, expressed
 * or implied, about its quality, reliability, or any other characteristic.
 */

/*
 * LUDecomposition.java
 * Copyright (C) 1999 The Mathworks and NIST
 *
 */

package weka.core.matrix;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/** 
 * LU Decomposition.
 * <P>
 * For an m-by-n matrix A with m &gt;= n, the LU decomposition is an m-by-n
 * unit lower triangular matrix L, an n-by-n upper triangular matrix U, and a
 * permutation vector piv of length m so that A(piv,:) = L*U.  If m &lt; n,
 * then L is m-by-m and U is m-by-n.
 * <P>
 * The LU decompostion with pivoting always exists, even if the matrix is
 * singular, so the constructor will never fail.  The primary use of the LU
 * decomposition is in the solution of square systems of simultaneous linear
 * equations.  This will fail if isNonsingular() returns false.
 * <p/>
 * Adapted from the <a href="http://math.nist.gov/javanumerics/jama/" target="_blank">JAMA</a> package.
 *
 * @author The Mathworks and NIST 
 * @author Fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class LUDecomposition 
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -2731022568037808629L;

  /** 
   * Array for internal storage of decomposition.
   * @serial internal array storage.
   */
  private double[][] LU;

  /** 
   * Row and column dimensions, and pivot sign.
   * @serial column dimension.
   * @serial row dimension.
   * @serial pivot sign.
   */
  private int m, n, pivsign; 

  /** 
   * Internal storage of pivot vector.
   * @serial pivot vector.
   */
  private int[] piv;

  /** 
   * LU Decomposition
   * @param  A   Rectangular matrix
   */
  public LUDecomposition(Matrix A) {

    // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

    LU = A.getArrayCopy();
    m = A.getRowDimension();
    n = A.getColumnDimension();
    piv = new int[m];
    for (int i = 0; i < m; i++) {
      piv[i] = i;
    }
    pivsign = 1;
    double[] LUrowi;
    double[] LUcolj = new double[m];

    // Outer loop.

    for (int j = 0; j < n; j++) {

      // Make a copy of the j-th column to localize references.

      for (int i = 0; i < m; i++) {
        LUcolj[i] = LU[i][j];
      }

      // Apply previous transformations.

      for (int i = 0; i < m; i++) {
        LUrowi = LU[i];

        // Most of the time is spent in the following dot product.

        int kmax = Math.min(i,j);
        double s = 0.0;
        for (int k = 0; k < kmax; k++) {
          s += LUrowi[k]*LUcolj[k];
        }

        LUrowi[j] = LUcolj[i] -= s;
      }

      // Find pivot and exchange if necessary.

      int p = j;
      for (int i = j+1; i < m; i++) {
        if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
          p = i;
        }
      }
      if (p != j) {
        for (int k = 0; k < n; k++) {
          double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
        }
        int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
        pivsign = -pivsign;
      }

      // Compute multipliers.

      if (j < m & LU[j][j] != 0.0) {
        for (int i = j+1; i < m; i++) {
          LU[i][j] /= LU[j][j];
        }
      }
    }
  }

  /** 
   * Is the matrix nonsingular?
   * @return     true if U, and hence A, is nonsingular.
   */
  public boolean isNonsingular() {
    for (int j = 0; j < n; j++) {
      if (LU[j][j] == 0)
        return false;
    }
    return true;
  }

  /** 
   * Return lower triangular factor
   * @return     L
   */
  public Matrix getL() {
    Matrix X = new Matrix(m,n);
    double[][] L = X.getArray();
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (i > j) {
          L[i][j] = LU[i][j];
        } else if (i == j) {
          L[i][j] = 1.0;
        } else {
          L[i][j] = 0.0;
        }
      }
    }
    return X;
  }

  /** 
   * Return upper triangular factor
   * @return     U
   */
  public Matrix getU() {
    Matrix X = new Matrix(n,n);
    double[][] U = X.getArray();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i <= j) {
          U[i][j] = LU[i][j];
        } else {
          U[i][j] = 0.0;
        }
      }
    }
    return X;
  }

  /** 
   * Return pivot permutation vector
   * @return     piv
   */
  public int[] getPivot() {
    int[] p = new int[m];
    for (int i = 0; i < m; i++) {
      p[i] = piv[i];
    }
    return p;
  }

  /** 
   * Return pivot permutation vector as a one-dimensional double array
   * @return     (double) piv
   */
  public double[] getDoublePivot() {
    double[] vals = new double[m];
    for (int i = 0; i < m; i++) {
      vals[i] = (double) piv[i];
    }
    return vals;
  }

  /** 
   * Determinant
   * @return     det(A)
   * @exception  IllegalArgumentException  Matrix must be square
   */
  public double det() {
    if (m != n) {
      throw new IllegalArgumentException("Matrix must be square.");
    }
    double d = (double) pivsign;
    for (int j = 0; j < n; j++) {
      d *= LU[j][j];
    }
    return d;
  }

  /** 
   * Solve A*X = B
   * @param  B   A Matrix with as many rows as A and any number of columns.
   * @return     X so that L*U*X = B(piv,:)
   * @exception  IllegalArgumentException Matrix row dimensions must agree.
   * @exception  RuntimeException  Matrix is singular.
   */
  public Matrix solve(Matrix B) {
    if (B.getRowDimension() != m) {
      throw new IllegalArgumentException("Matrix row dimensions must agree.");
    }
    if (!this.isNonsingular()) {
      throw new RuntimeException("Matrix is singular.");
    }

    // Copy right hand side with pivoting
    int nx = B.getColumnDimension();
    Matrix Xmat = B.getMatrix(piv,0,nx-1);
    double[][] X = Xmat.getArray();

    // Solve L*Y = B(piv,:)
    for (int k = 0; k < n; k++) {
      for (int i = k+1; i < n; i++) {
        for (int j = 0; j < nx; j++) {
          X[i][j] -= X[k][j]*LU[i][k];
        }
      }
    }
    // Solve U*X = Y;
    for (int k = n-1; k >= 0; k--) {
      for (int j = 0; j < nx; j++) {
        X[k][j] /= LU[k][k];
      }
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < nx; j++) {
          X[i][j] -= X[k][j]*LU[i][k];
        }
      }
    }
    return Xmat;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
}

/*
 * This software is a cooperative product of The MathWorks and the National
 * Institute of Standards and Technology (NIST) which has been released to the
 * public domain. Neither The MathWorks nor NIST assumes any responsibility
 * whatsoever for its use by other parties, and makes no guarantees, expressed
 * or implied, about its quality, reliability, or any other characteristic.
 */

/*
 * QRDecomposition.java
 * Copyright (C) 1999 The Mathworks and NIST
 *
 */

package weka.core.matrix;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/** 
 * QR Decomposition.
 * <P>
 * For an m-by-n matrix A with m &gt;= n, the QR decomposition is an m-by-n
 * orthogonal matrix Q and an n-by-n upper triangular matrix R so that A = Q*R.
 * <P>
 * The QR decompostion always exists, even if the matrix does not have full
 * rank, so the constructor will never fail.  The primary use of the QR
 * decomposition is in the least squares solution of nonsquare systems of
 * simultaneous linear equations.  This will fail if isFullRank() returns false.
 * <p/>
 * Adapted from the <a href="http://math.nist.gov/javanumerics/jama/" target="_blank">JAMA</a> package.
 *
 * @author The Mathworks and NIST 
 * @author Fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class QRDecomposition 
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -5013090736132211418L;

  /** 
   * Array for internal storage of decomposition.
   *    @serial internal array storage.
   */
  private double[][] QR;

  /** 
   * Row and column dimensions.
   *    @serial column dimension.
   *    @serial row dimension.
   */
  private int m, n;

  /** 
   * Array for internal storage of diagonal of R.
   *    @serial diagonal of R.
   */
  private double[] Rdiag;

  /** 
   * QR Decomposition, computed by Householder reflections.
   * @param A    Rectangular matrix
   */
  public QRDecomposition(Matrix A) {
    // Initialize.
    QR = A.getArrayCopy();
    m = A.getRowDimension();
    n = A.getColumnDimension();
    Rdiag = new double[n];

    // Main loop.
    for (int k = 0; k < n; k++) {
      // Compute 2-norm of k-th column without under/overflow.
      double nrm = 0;
      for (int i = k; i < m; i++) {
        nrm = Maths.hypot(nrm,QR[i][k]);
      }

      if (nrm != 0.0) {
        // Form k-th Householder vector.
        if (QR[k][k] < 0) {
          nrm = -nrm;
        }
        for (int i = k; i < m; i++) {
          QR[i][k] /= nrm;
        }
        QR[k][k] += 1.0;

        // Apply transformation to remaining columns.
        for (int j = k+1; j < n; j++) {
          double s = 0.0; 
          for (int i = k; i < m; i++) {
            s += QR[i][k]*QR[i][j];
          }
          s = -s/QR[k][k];
          for (int i = k; i < m; i++) {
            QR[i][j] += s*QR[i][k];
          }
        }
      }
      Rdiag[k] = -nrm;
    }
  }

  /** 
   * Is the matrix full rank?
   * @return     true if R, and hence A, has full rank.
   */
  public boolean isFullRank() {
    for (int j = 0; j < n; j++) {
      if (Rdiag[j] == 0)
        return false;
    }
    return true;
  }

  /** 
   * Return the Householder vectors
   * @return     Lower trapezoidal matrix whose columns define the reflections
   */
  public Matrix getH() {
    Matrix X = new Matrix(m,n);
    double[][] H = X.getArray();
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (i >= j) {
          H[i][j] = QR[i][j];
        } else {
          H[i][j] = 0.0;
        }
      }
    }
    return X;
  }

  /** 
   * Return the upper triangular factor
   * @return     R
   */
  public Matrix getR() {
    Matrix X = new Matrix(n,n);
    double[][] R = X.getArray();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i < j) {
          R[i][j] = QR[i][j];
        } else if (i == j) {
          R[i][j] = Rdiag[i];
        } else {
          R[i][j] = 0.0;
        }
      }
    }
    return X;
  }

  /** 
   * Generate and return the (economy-sized) orthogonal factor
   * @return     Q
   */
  public Matrix getQ() {
    Matrix X = new Matrix(m,n);
    double[][] Q = X.getArray();
    for (int k = n-1; k >= 0; k--) {
      for (int i = 0; i < m; i++) {
        Q[i][k] = 0.0;
      }
      Q[k][k] = 1.0;
      for (int j = k; j < n; j++) {
        if (QR[k][k] != 0) {
          double s = 0.0;
          for (int i = k; i < m; i++) {
            s += QR[i][k]*Q[i][j];
          }
          s = -s/QR[k][k];
          for (int i = k; i < m; i++) {
            Q[i][j] += s*QR[i][k];
          }
        }
      }
    }
    return X;
  }

  /** 
   * Least squares solution of A*X = B
   * @param B    A Matrix with as many rows as A and any number of columns.
   * @return     X that minimizes the two norm of Q*R*X-B.
   * @exception  IllegalArgumentException  Matrix row dimensions must agree.
   * @exception  RuntimeException  Matrix is rank deficient.
   */
  public Matrix solve(Matrix B) {
    if (B.getRowDimension() != m) {
      throw new IllegalArgumentException("Matrix row dimensions must agree.");
    }
    if (!this.isFullRank()) {
      throw new RuntimeException("Matrix is rank deficient.");
    }

    // Copy right hand side
    int nx = B.getColumnDimension();
    double[][] X = B.getArrayCopy();

    // Compute Y = transpose(Q)*B
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < nx; j++) {
        double s = 0.0; 
        for (int i = k; i < m; i++) {
          s += QR[i][k]*X[i][j];
        }
        s = -s/QR[k][k];
        for (int i = k; i < m; i++) {
          X[i][j] += s*QR[i][k];
        }
      }
    }
    // Solve R*X = Y;
    for (int k = n-1; k >= 0; k--) {
      for (int j = 0; j < nx; j++) {
        X[k][j] /= Rdiag[k];
      }
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < nx; j++) {
          X[i][j] -= X[k][j]*QR[i][k];
        }
      }
    }
    return (new Matrix(X,n,nx).getMatrix(0,n-1,0,nx-1));
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

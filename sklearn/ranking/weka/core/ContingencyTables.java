/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    ContingencyTables.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

/**
 * Class implementing some statistical routines for contingency tables.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class ContingencyTables
  implements RevisionHandler {

  /** The natural logarithm of 2 */
  private static double log2 = Math.log(2);

  /**
   * Returns chi-squared probability for a given matrix.
   *
   * @param matrix the contigency table
   * @param yates is Yates' correction to be used?
   * @return the chi-squared probability
   */

  public static double chiSquared(double [][] matrix, boolean yates) {

    int df = (matrix.length - 1) * (matrix[0].length - 1);

    return Statistics.chiSquaredProbability(chiVal(matrix, yates), df);
  }

  /**
   * Computes chi-squared statistic for a contingency table.
   *
   * @param matrix the contigency table
   * @param useYates is Yates' correction to be used?
   * @return the value of the chi-squared statistic
   */
  public static double chiVal(double [][] matrix, boolean useYates) {
    
    int df, nrows, ncols, row, col;
    double[] rtotal, ctotal;
    double expect = 0, chival = 0, n = 0;
    boolean yates = true;
    
    nrows = matrix.length;
    ncols = matrix[0].length;
    rtotal = new double [nrows];
    ctotal = new double [ncols];
    for (row = 0; row < nrows; row++) {
      for (col = 0; col < ncols; col++) {
	rtotal[row] += matrix[row][col];
	ctotal[col] += matrix[row][col];
	n += matrix[row][col];
      }
    }
    df = (nrows - 1)*(ncols - 1);
    if ((df > 1) || (!useYates)) {
      yates = false;
    } else if (df <= 0) {
      return 0;
    }
    chival = 0.0;
    for (row = 0; row < nrows; row++) {
      if (Utils.gr(rtotal[row], 0)) {
	for (col = 0; col < ncols; col++) {
	  if (Utils.gr(ctotal[col], 0)) {
	    expect = (ctotal[col] * rtotal[row]) / n;
	    chival += chiCell (matrix[row][col], expect, yates);
	  }
	}
      }
    }
    return chival;
  }

  /**
   * Tests if Cochran's criterion is fullfilled for the given
   * contingency table. Rows and columns with all zeros are not considered
   * relevant.
   *
   * @param matrix the contigency table to be tested
   * @return true if contingency table is ok, false if not
   */
  public static boolean cochransCriterion(double[][] matrix) {

    double[] rtotal, ctotal;
    double n = 0, expect, smallfreq = 5;
    int smallcount = 0, nonZeroRows = 0, nonZeroColumns = 0, nrows, ncols, 
      row, col;

    nrows = matrix.length;
    ncols = matrix[0].length;

    rtotal = new double [nrows];
    ctotal = new double [ncols];
    for (row = 0; row < nrows; row++) {
      for (col = 0; col < ncols; col++) {
	rtotal[row] += matrix[row][col];
	ctotal[col] += matrix[row][col];
	n += matrix[row][col];
      }
    }
    for (row = 0; row < nrows; row++) {
      if (Utils.gr(rtotal[row], 0)) {
	nonZeroRows++;
      }
    }
    for (col = 0; col < ncols; col++) {
      if (Utils.gr(ctotal[col], 0)) {
	nonZeroColumns++;
      }
    }
    for (row = 0; row < nrows; row++) {
      if (Utils.gr(rtotal[row], 0)) {
	for (col = 0; col < ncols; col++) {
	  if (Utils.gr(ctotal[col], 0)) {
	    expect = (ctotal[col] * rtotal[row]) / n;
	    if (Utils.sm(expect, smallfreq)) {
	      if (Utils.sm(expect, 1)) {
		return false;
	      } else {
		smallcount++;
		if (smallcount > (nonZeroRows * nonZeroColumns) / smallfreq) {
		  return false;
		}
	      }
	    }
	  }
	}
      }
    }
    return true;
  }

  /**
   * Computes Cramer's V for a contingency table.
   *
   * @param matrix the contingency table
   * @return Cramer's V
   */
  public static double CramersV(double [][] matrix) {

    int row, col, nrows,ncols, min;
    double n = 0;
    
    nrows = matrix.length;
    ncols = matrix[0].length;
    for (row = 0; row < nrows; row++) {
      for (col = 0; col < ncols; col++) {
	n += matrix[row][col];
      }
    }
    min = nrows < ncols ? nrows-1 : ncols-1;
    if ((min == 0) || Utils.eq(n, 0))
      return 0;
    return Math.sqrt(chiVal(matrix, false) / (n * (double)min)); 
  } 

  /**
   * Computes the entropy of the given array.
   *
   * @param array the array
   * @return the entropy
   */
  public static double entropy(double[] array) {

    double returnValue = 0, sum = 0;

    for (int i = 0; i < array.length; i++) {
      returnValue -= lnFunc(array[i]);
      sum += array[i];
    }
    if (Utils.eq(sum, 0)) {
      return 0;
    } else {
      return (returnValue + lnFunc(sum)) / (sum * log2);
    }
  }

  /**
   * Computes conditional entropy of the rows given
   * the columns.
   *
   * @param matrix the contingency table
   * @return the conditional entropy of the rows given the columns
   */
  public static double entropyConditionedOnColumns(double[][] matrix) {
    
    double returnValue = 0, sumForColumn, total = 0;

    for (int j = 0; j < matrix[0].length; j++) {
      sumForColumn = 0;
      for (int i = 0; i < matrix.length; i++) {
	returnValue = returnValue + lnFunc(matrix[i][j]);
	sumForColumn += matrix[i][j];
      }
      returnValue = returnValue - lnFunc(sumForColumn);
      total += sumForColumn;
    }
    if (Utils.eq(total, 0)) {
      return 0;
    }
    return -returnValue / (total * log2);
  }

  /**
   * Computes conditional entropy of the columns given
   * the rows.
   *
   * @param matrix the contingency table
   * @return the conditional entropy of the columns given the rows
   */
  public static double entropyConditionedOnRows(double[][] matrix) {
    
    double returnValue = 0, sumForRow, total = 0;

    for (int i = 0; i < matrix.length; i++) {
      sumForRow = 0;
      for (int j = 0; j < matrix[0].length; j++) {
	returnValue = returnValue + lnFunc(matrix[i][j]);
	sumForRow += matrix[i][j];
      }
      returnValue = returnValue - lnFunc(sumForRow);
      total += sumForRow;
    }
    if (Utils.eq(total, 0)) {
      return 0;
    }
    return -returnValue / (total * log2);
  }

  /**
   * Computes conditional entropy of the columns given the rows
   * of the test matrix with respect to the train matrix. Uses a
   * Laplace prior. Does NOT normalize the entropy.
   *
   * @param train the train matrix 
   * @param test the test matrix
   * @param numClasses the number of symbols for Laplace
   * @return the entropy
   */
  public static double entropyConditionedOnRows(double[][] train, 
						double[][] test,
						double numClasses) {
    
    double returnValue = 0, trainSumForRow, testSumForRow, testSum = 0;

    for (int i = 0; i < test.length; i++) {
      trainSumForRow = 0;
      testSumForRow = 0;
      for (int j = 0; j < test[0].length; j++) {
	returnValue -= test[i][j] * Math.log(train[i][j] + 1);
	trainSumForRow += train[i][j];
	testSumForRow += test[i][j];
      }
      testSum = testSumForRow;
      returnValue += testSumForRow * Math.log(trainSumForRow + 
					     numClasses);
    }
    return returnValue / (testSum * log2);
  }

  /**
   * Computes the rows' entropy for the given contingency table.
   *
   * @param matrix the contingency table
   * @return the rows' entropy
   */
  public static double entropyOverRows(double[][] matrix) {
    
    double returnValue = 0, sumForRow, total = 0;

    for (int i = 0; i < matrix.length; i++) {
      sumForRow = 0;
      for (int j = 0; j < matrix[0].length; j++) {
	sumForRow += matrix[i][j];
      }
      returnValue = returnValue - lnFunc(sumForRow);
      total += sumForRow;
    }
    if (Utils.eq(total, 0)) {
      return 0;
    }
    return (returnValue + lnFunc(total)) / (total * log2);
  }

  /**
   * Computes the columns' entropy for the given contingency table.
   *
   * @param matrix the contingency table
   * @return the columns' entropy
   */
  public static double entropyOverColumns(double[][] matrix){
    
    double returnValue = 0, sumForColumn, total = 0;

    for (int j = 0; j < matrix[0].length; j++){
      sumForColumn = 0;
      for (int i = 0; i < matrix.length; i++) {
	sumForColumn += matrix[i][j];
      }
      returnValue = returnValue - lnFunc(sumForColumn);
      total += sumForColumn;
    }
    if (Utils.eq(total, 0)) {
      return 0;
    }
    return (returnValue + lnFunc(total)) / (total * log2);
  }

  /**
   * Computes gain ratio for contingency table (split on rows).
   * Returns Double.MAX_VALUE if the split entropy is 0.
   *
   * @param matrix the contingency table
   * @return the gain ratio
   */
  public static double gainRatio(double[][] matrix){
    
    double preSplit = 0, postSplit = 0, splitEnt = 0,
      sumForRow, sumForColumn, total = 0, infoGain;

    // Compute entropy before split
    for (int i = 0; i < matrix[0].length; i++) {
      sumForColumn = 0;
      for (int j = 0; j < matrix.length; j++) 
	sumForColumn += matrix[j][i];
      preSplit += lnFunc(sumForColumn);
      total += sumForColumn;
    }
    preSplit -= lnFunc(total);

    // Compute entropy after split and split entropy
    for (int i = 0; i < matrix.length; i++) {
      sumForRow = 0;
      for (int j = 0; j < matrix[0].length; j++) {
	postSplit += lnFunc(matrix[i][j]);
	sumForRow += matrix[i][j];
      }
      splitEnt += lnFunc(sumForRow);
    }
    postSplit -= splitEnt;
    splitEnt -= lnFunc(total);

    infoGain = preSplit - postSplit;
    if (Utils.eq(splitEnt, 0))
      return 0;
    return infoGain / splitEnt;
  }

  /**
   * Returns negative base 2 logarithm of multiple hypergeometric
   * probability for a contingency table.
   *
   * @param matrix the contingency table
   * @return the log of the hypergeometric probability of the contingency table 
   */
  public static double log2MultipleHypergeometric(double[][] matrix) {

    double sum = 0, sumForRow, sumForColumn, total = 0;

    for (int i = 0; i < matrix.length; i++) {
      sumForRow = 0;
      for (int j = 0; j < matrix[i].length; j++) {
	sumForRow += matrix[i][j];
      }
      sum += SpecialFunctions.lnFactorial(sumForRow);
      total += sumForRow;
    }
    for (int j = 0; j < matrix[0].length; j++) {
      sumForColumn = 0;
      for (int i = 0; i < matrix.length; i++) {
	sumForColumn += matrix [i][j];
      }
      sum += SpecialFunctions.lnFactorial(sumForColumn);
    }
    for (int i = 0; i < matrix.length; i++) {
      for (int j = 0; j < matrix[i].length; j++) {
	sum -= SpecialFunctions.lnFactorial(matrix[i][j]);
      }
    }
    sum -= SpecialFunctions.lnFactorial(total);
    return -sum / log2;
  }

  /**
   * Reduces a matrix by deleting all zero rows and columns.
   *
   * @param matrix the matrix to be reduced
   * @return the matrix with all zero rows and columns deleted
   */
  public static double[][] reduceMatrix(double[][] matrix) {

    int row, col, currCol, currRow, nrows, ncols, 
      nonZeroRows = 0, nonZeroColumns = 0;
    double[] rtotal, ctotal;
    double[][] newMatrix;

    nrows = matrix.length;
    ncols = matrix[0].length;
    rtotal = new double [nrows];
    ctotal = new double [ncols];
    for (row = 0; row < nrows; row++) {
      for (col = 0; col < ncols; col++) {
	rtotal[row] += matrix[row][col];
	ctotal[col] += matrix[row][col];
      }
    }
    for (row = 0; row < nrows; row++) {
      if (Utils.gr(rtotal[row],0)) {
	nonZeroRows++;
      }
    }
    for (col = 0; col < ncols; col++) {
      if (Utils.gr(ctotal[col],0)) {
	nonZeroColumns++;
      }
    }
    newMatrix = new double[nonZeroRows][nonZeroColumns];
    currRow = 0;
    for (row = 0; row < nrows; row++) {
      if (Utils.gr(rtotal[row],0)) {
	currCol = 0;
	for (col = 0; col < ncols; col++) {
	  if (Utils.gr(ctotal[col],0)) {
	    newMatrix[currRow][currCol] = matrix[row][col];
	    currCol++;
	  }
	}
	currRow++;
      }
    }
    return newMatrix;
  }
    
   /**
    * Calculates the symmetrical uncertainty for base 2.
    *
    * @param matrix the contingency table
    * @return the calculated symmetrical uncertainty
    *
    */
  public static double symmetricalUncertainty(double matrix[][]) {

    double sumForColumn, sumForRow, total = 0, columnEntropy = 0, 
      rowEntropy = 0, entropyConditionedOnRows = 0, infoGain = 0;

    // Compute entropy for columns
    for (int i = 0; i < matrix[0].length; i++) {
      sumForColumn = 0;
      for (int j = 0; j < matrix.length; j++) {
	sumForColumn += matrix[j][i];
      }
      columnEntropy += lnFunc(sumForColumn);
      total += sumForColumn;
    }
    columnEntropy -= lnFunc(total);

    // Compute entropy for rows and conditional entropy
    for (int i = 0; i < matrix.length; i++) {
      sumForRow = 0;
      for (int j = 0; j < matrix[0].length; j++) { 
	sumForRow += matrix[i][j];
	entropyConditionedOnRows += lnFunc(matrix[i][j]);
      }
      rowEntropy += lnFunc(sumForRow);
    }
    entropyConditionedOnRows -= rowEntropy;
    rowEntropy -= lnFunc(total);
    infoGain = columnEntropy - entropyConditionedOnRows;
    if (Utils.eq(columnEntropy, 0) || Utils.eq(rowEntropy, 0))
      return 0;
    return 2.0 * (infoGain / (columnEntropy + rowEntropy));
  }

  /**
   * Computes Goodman and Kruskal's tau-value for a contingency table.
   *
   * @param matrix the contingency table
   * @return Goodman and Kruskal's tau-value
   */
  public static double tauVal(double[][] matrix) {
    
    int nrows, ncols, row, col;
    double [] ctotal;
    double maxcol = 0, max, maxtotal = 0, n = 0;
    
    nrows = matrix.length;
    ncols = matrix[0].length;
    ctotal = new double [ncols];
    for (row = 0; row < nrows; row++) {
      max = 0;
      for (col = 0; col < ncols; col++) {
	if (Utils.gr(matrix[row][col], max)) 
	  max = matrix[row][col];
	ctotal[col] += matrix[row][col];
	n += matrix[row][col];
      }
      maxtotal += max;
    }
    if (Utils.eq(n, 0)) {
      return 0;
    }
    maxcol = ctotal[Utils.maxIndex(ctotal)];
    return (maxtotal - maxcol)/(n - maxcol);
  }

  /**
   * Help method for computing entropy.
   */
  private static double lnFunc(double num){
    
    // Constant hard coded for efficiency reasons
    if (num < 1e-6) {
      return 0;
    } else {
      return num * Math.log(num);
    }
  }

  /**
   * Computes chi-value for one cell in a contingency table.
   *
   * @param freq the observed frequency in the cell
   * @param expected the expected frequency in the cell
   * @return the chi-value for that cell; 0 if the expected value is
   * too close to zero 
   */
  private static double chiCell(double freq, double expected, 
                                boolean yates){

    // Cell in empty row and column?
    if (Utils.smOrEq(expected, 0)) {
      return 0;
    }

    // Compute difference between observed and expected value
    double diff = Math.abs(freq - expected);
    if (yates) {

      // Apply Yates' correction if wanted
      diff -= 0.5;

      // The difference should never be negative
      if (diff < 0) {
        diff = 0;
      }
    }

    // Return chi-value for the cell
    return (diff * diff / expected);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }

  /**
   * Main method for testing this class.
   */
  public static void main(String[] ops) {

    double[] firstRow = {10, 5, 20};
    double[] secondRow = {2, 10, 6};
    double[] thirdRow = {5, 10, 10};
    double[][] matrix = new double[3][0];

    matrix[0] = firstRow; matrix[1] = secondRow; matrix[2] = thirdRow;
    for (int i = 0; i < matrix.length; i++) {
      for (int j = 0; j < matrix[i].length; j++) {
	System.out.print(matrix[i][j] + " ");
      }
      System.out.println();
    }
    System.out.println("Chi-squared probability: " +
		       ContingencyTables.chiSquared(matrix, false));
    System.out.println("Chi-squared value: " +
		       ContingencyTables.chiVal(matrix, false));
    System.out.println("Cochran's criterion fullfilled: " +
		       ContingencyTables.cochransCriterion(matrix));
    System.out.println("Cramer's V: " +
		       ContingencyTables.CramersV(matrix));
    System.out.println("Entropy of first row: " +
		       ContingencyTables.entropy(firstRow));
    System.out.println("Entropy conditioned on columns: " +
		       ContingencyTables.entropyConditionedOnColumns(matrix));
    System.out.println("Entropy conditioned on rows: " +
		       ContingencyTables.entropyConditionedOnRows(matrix));
    System.out.println("Entropy conditioned on rows (with Laplace): " +
		       ContingencyTables.entropyConditionedOnRows(matrix, matrix, 3));
    System.out.println("Entropy of rows: " +
		       ContingencyTables.entropyOverRows(matrix));
    System.out.println("Entropy of columns: " +
		       ContingencyTables.entropyOverColumns(matrix));
    System.out.println("Gain ratio: " +
		       ContingencyTables.gainRatio(matrix));
    System.out.println("Negative log2 of multiple hypergeometric probability: " +
		       ContingencyTables.log2MultipleHypergeometric(matrix));
    System.out.println("Symmetrical uncertainty: " +
		       ContingencyTables.symmetricalUncertainty(matrix));
    System.out.println("Tau value: " +
		       ContingencyTables.tauVal(matrix));
    double[][] newMatrix = new double[3][3];
    newMatrix[0][0] = 1; newMatrix[0][1] = 0; newMatrix[0][2] = 1;
    newMatrix[1][0] = 0; newMatrix[1][1] = 0; newMatrix[1][2] = 0;
    newMatrix[2][0] = 1; newMatrix[2][1] = 0; newMatrix[2][2] = 1;
    System.out.println("Matrix with empty row and column: ");
    for (int i = 0; i < newMatrix.length; i++) {
      for (int j = 0; j < newMatrix[i].length; j++) {
	System.out.print(newMatrix[i][j] + " ");
      }
      System.out.println();
    }
    System.out.println("Reduced matrix: ");
    newMatrix = ContingencyTables.reduceMatrix(newMatrix);
    for (int i = 0; i < newMatrix.length; i++) {
      for (int j = 0; j < newMatrix[i].length; j++) {
	System.out.print(newMatrix[i][j] + " ");
      }
      System.out.println();
    }
  }
}









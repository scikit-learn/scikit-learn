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
 *    CostMatrix.java
 *    Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.core.AttributeExpression;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Matrix;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.LineNumberReader;
import java.io.Reader;
import java.io.Serializable;
import java.io.StreamTokenizer;
import java.io.Writer;
import java.util.Random;
import java.util.StringTokenizer;

/**
 * Class for storing and manipulating a misclassification cost matrix.
 * The element at position i,j in the matrix is the penalty for classifying
 * an instance of class j as class i. Cost values can be fixed or
 * computed on a per-instance basis (cost sensitive evaluation only)
 * from the value of an attribute or an expression involving
 * attribute(s).
 *
 * @author Mark Hall
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 6041 $
 */
public class CostMatrix implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -1973792250544554965L;

  private int m_size;

  /** [rows][columns] */
  protected Object [][] m_matrix;

  /** The deafult file extension for cost matrix files */
  public static String FILE_EXTENSION = ".cost";

  /**
   * Creates a default cost matrix of a particular size.
   * All diagonal values will be 0 and all non-diagonal values 1.
   *
   * @param numOfClasses the number of classes that the cost matrix holds.
   */
  public CostMatrix(int numOfClasses) {
    m_size = numOfClasses;
    initialize();
  }

  /**
   * Creates a cost matrix that is a copy of another.
   *
   * @param toCopy the matrix to copy.
   */
  public CostMatrix(CostMatrix toCopy) {
    this(toCopy.size());

    for (int i = 0; i < m_size; i++) {
      for (int j = 0; j < m_size; j++) {
        setCell(i, j, toCopy.getCell(i, j));
      }
    }
  }

  /**
   * Initializes the matrix
   */
  public void initialize() {
    m_matrix = new Object[m_size][m_size];
    for (int i = 0; i < m_size; i++) {
      for (int j = 0; j < m_size; j++) {
        setCell(i, j, i == j ? new Double(0.0) : new Double(1.0));
      }
    }
  }

  /**
   * The number of rows (and columns)
   * @return the size of the matrix
   */
  public int size() {
    return m_size;
  }

  /**
   * Same as size
   * @return the number of columns
   */
  public int numColumns() {
    return size();
  }

  /**
   * Same as size
   * @return the number of rows
   */
  public int numRows() {
    return size();
  }

  private boolean replaceStrings() throws Exception {
    boolean nonDouble = false;

    for (int i = 0; i < m_size; i++) {
      for (int j = 0; j < m_size; j++) {
        if (getCell(i, j) instanceof String) {
          AttributeExpression temp = new AttributeExpression();
          temp.convertInfixToPostfix((String)getCell(i, j));
          setCell(i, j, temp);
          nonDouble = true;
        } else if (getCell(i, j) instanceof AttributeExpression) {
          nonDouble = true;
        }
      }
    }

    return nonDouble;
  }

  /**
   * Applies the cost matrix to a set of instances. If a random number generator is
   * supplied the instances will be resampled, otherwise they will be rewighted.
   * Adapted from code once sitting in Instances.java
   *
   * @param data the instances to reweight.
   * @param random a random number generator for resampling, if null then instances are
   * rewighted.
   * @return a new dataset reflecting the cost of misclassification.
   * @exception Exception if the data has no class or the matrix in inappropriate.
   */
  public Instances applyCostMatrix(Instances data, Random random)
    throws Exception {

    double sumOfWeightFactors = 0, sumOfMissClassWeights,
      sumOfWeights;
    double [] weightOfInstancesInClass, weightFactor, weightOfInstances;
    Instances newData;

    if (data.classIndex() < 0) {
      throw new Exception("Class index is not set!");
    }

    if (size() != data.numClasses()) {
      throw new Exception("Misclassification cost matrix has wrong format!");
    }

    // are there any non-fixed, per-instance costs defined in the matrix?
    if (replaceStrings()) {
      // could reweight in the two class case
      if (data.classAttribute().numValues() > 2) {
        throw new Exception("Can't resample/reweight instances using "
                            +"non-fixed cost values when there are more "
                            +"than two classes!");
      } else {
        // Store new weights
        weightOfInstances = new double[data.numInstances()];
        for (int i = 0; i < data.numInstances(); i++) {
          Instance inst = data.instance(i);
          int classValIndex = (int)inst.classValue();
          double factor = 1.0;
          Object element = (classValIndex == 0)
            ? getCell(classValIndex, 1)
            : getCell(classValIndex, 0);
          if (element instanceof Double) {
            factor = ((Double)element).doubleValue();
          } else {
            factor = ((AttributeExpression)element).evaluateExpression(inst);
          }
          weightOfInstances[i] = inst.weight() * factor;
          /*          System.err.println("Multiplying " + inst.classAttribute().value((int)inst.classValue())
                      +" by factor " + factor); */
        }

        // Change instances weight or do resampling
        if (random != null) {
          return data.resampleWithWeights(random, weightOfInstances);
        } else {
          Instances instances = new Instances(data);
          for (int i = 0; i < data.numInstances(); i++) {
            instances.instance(i).setWeight(weightOfInstances[i]);
          }
          return instances;
        }
      }
    }

    weightFactor = new double[data.numClasses()];
    weightOfInstancesInClass = new double[data.numClasses()];
    for (int j = 0; j < data.numInstances(); j++) {
      weightOfInstancesInClass[(int)data.instance(j).classValue()] +=
        data.instance(j).weight();
    }
    sumOfWeights = Utils.sum(weightOfInstancesInClass);

    // normalize the matrix if not already
    for (int i=0; i< m_size; i++) {
      if (!Utils.eq(((Double)getCell(i, i)).doubleValue(), 0)) {
        CostMatrix normMatrix = new CostMatrix(this);
        normMatrix.normalize();
        return normMatrix.applyCostMatrix(data, random);
      }
    }

    for (int i = 0; i < data.numClasses(); i++) {
      // Using Kai Ming Ting's formula for deriving weights for
      // the classes and Breiman's heuristic for multiclass
      // problems.

      sumOfMissClassWeights = 0;
      for (int j = 0; j < data.numClasses(); j++) {
        if (Utils.sm(((Double)getCell(i,j)).doubleValue(),0)) {
          throw new Exception("Neg. weights in misclassification "+
              "cost matrix!");
        }
        sumOfMissClassWeights
          += ((Double)getCell(i,j)).doubleValue();
      }
      weightFactor[i] = sumOfMissClassWeights * sumOfWeights;
      sumOfWeightFactors += sumOfMissClassWeights *
        weightOfInstancesInClass[i];
    }
    for (int i = 0; i < data.numClasses(); i++) {
      weightFactor[i] /= sumOfWeightFactors;
    }

    // Store new weights
    weightOfInstances = new double[data.numInstances()];
    for (int i = 0; i < data.numInstances(); i++) {
      weightOfInstances[i] = data.instance(i).weight()*
        weightFactor[(int)data.instance(i).classValue()];
    }

    // Change instances weight or do resampling
    if (random != null) {
      return data.resampleWithWeights(random, weightOfInstances);
    } else {
      Instances instances = new Instances(data);
      for (int i = 0; i < data.numInstances(); i++) {
        instances.instance(i).setWeight(weightOfInstances[i]);
      }
      return instances;
    }
  }

  /**
   * Calculates the expected misclassification cost for each possible class value,
   * given class probability estimates.
   *
   * @param classProbs the class probability estimates.
   * @return the expected costs.
   * @exception Exception if the wrong number of class probabilities is supplied.
   */
  public double[] expectedCosts(double[] classProbs) throws Exception {

    if (classProbs.length != m_size) {
      throw new Exception("Length of probability estimates don't "
                          +"match cost matrix");
    }

    double[] costs = new double[m_size];

    for (int x = 0; x < m_size; x++) {
      for (int y = 0; y < m_size; y++) {
        Object element = getCell(y, x);
        if (!(element instanceof Double)) {
          throw new Exception("Can't use non-fixed costs in "
                              +"computing expected costs.");
        }
        costs[x] += classProbs[y] * ((Double)element).doubleValue();
      }
    }

    return costs;
  }

  /**
   * Calculates the expected misclassification cost for each possible class value,
   * given class probability estimates.
   *
   * @param classProbs the class probability estimates.
   * @param inst the current instance for which the class probabilites
   * apply. Is used for computing any non-fixed cost values.
   * @return the expected costs.
   * @exception Exception if something goes wrong
   */
  public double[] expectedCosts(double [] classProbs,
                                Instance inst) throws Exception {

    if (classProbs.length != m_size) {
      throw new Exception("Length of probability estimates don't "
                          +"match cost matrix");
    }

    if (!replaceStrings()) {
      return expectedCosts(classProbs);
    }

    double[] costs = new double[m_size];

    for (int x = 0; x < m_size; x++) {
      for (int y = 0; y < m_size; y++) {
        Object element = getCell(y, x);
        double costVal;
        if (!(element instanceof Double)) {
          costVal =
            ((AttributeExpression)element).evaluateExpression(inst);
        } else {
          costVal = ((Double)element).doubleValue();
        }
        costs[x] += classProbs[y] * costVal;
      }
    }

    return costs;
  }

  /**
   * Gets the maximum cost for a particular class value.
   *
   * @param classVal the class value.
   * @return the maximum cost.
   * @exception Exception if cost matrix contains non-fixed
   * costs
   */
  public double getMaxCost(int classVal) throws Exception {

    double maxCost = Double.NEGATIVE_INFINITY;

    for (int i = 0; i < m_size; i++) {
      Object element = getCell(classVal, i);
      if (!(element instanceof Double)) {
          throw new Exception("Can't use non-fixed costs when "
                              +"getting max cost.");
      }
      double cost = ((Double)element).doubleValue();
      if (cost > maxCost) maxCost = cost;
    }

    return maxCost;
  }

  /**
   * Gets the maximum cost for a particular class value.
   *
   * @param classVal the class value.
   * @return the maximum cost.
   * @exception Exception if cost matrix contains non-fixed
   * costs
   */
  public double getMaxCost(int classVal, Instance inst)
    throws Exception {

    if (!replaceStrings()) {
      return getMaxCost(classVal);
    }

    double maxCost = Double.NEGATIVE_INFINITY;
    double cost;
    for (int i = 0; i < m_size; i++) {
      Object element = getCell(classVal, i);
      if (!(element instanceof Double)) {
        cost =
          ((AttributeExpression)element).evaluateExpression(inst);
      } else {
        cost = ((Double)element).doubleValue();
      }
      if (cost > maxCost) maxCost = cost;
    }

    return maxCost;
  }


  /**
   * Normalizes the matrix so that the diagonal contains zeros.
   *
   */
  public void normalize() {

    for (int y=0; y<m_size; y++) {
      double diag = ((Double)getCell(y, y)).doubleValue();
      for (int x=0; x<m_size; x++) {
        setCell(x, y, new Double(((Double)getCell(x, y)).
                                    doubleValue() - diag));
      }
    }
  }

  /**
   * Loads a cost matrix in the old format from a reader. Adapted from code once sitting
   * in Instances.java
   *
   * @param reader the reader to get the values from.
   * @exception Exception if the matrix cannot be read correctly.
   */
  public void readOldFormat(Reader reader) throws Exception {

    StreamTokenizer tokenizer;
    int currentToken;
    double firstIndex, secondIndex, weight;

    tokenizer = new StreamTokenizer(reader);

    initialize();

    tokenizer.commentChar('%');
    tokenizer.eolIsSignificant(true);
    while (StreamTokenizer.TT_EOF != (currentToken = tokenizer.nextToken())) {

      // Skip empty lines
      if (currentToken == StreamTokenizer.TT_EOL) {
        continue;
      }

      // Get index of first class.
      if (currentToken != StreamTokenizer.TT_NUMBER) {
        throw new Exception("Only numbers and comments allowed "+
            "in cost file!");
      }
      firstIndex = tokenizer.nval;
      if (!Utils.eq((double)(int)firstIndex,firstIndex)) {
        throw new Exception("First number in line has to be "+
            "index of a class!");
      }
      if ((int)firstIndex >= size()) {
        throw new Exception("Class index out of range!");
      }

      // Get index of second class.
      if (StreamTokenizer.TT_EOF == (currentToken = tokenizer.nextToken())) {
        throw new Exception("Premature end of file!");
      }
      if (currentToken == StreamTokenizer.TT_EOL) {
        throw new Exception("Premature end of line!");
      }
      if (currentToken != StreamTokenizer.TT_NUMBER) {
        throw new Exception("Only numbers and comments allowed "+
            "in cost file!");
      }
      secondIndex = tokenizer.nval;
      if (!Utils.eq((double)(int)secondIndex,secondIndex)) {
        throw new Exception("Second number in line has to be "+
            "index of a class!");
      }
      if ((int)secondIndex >= size()) {
        throw new Exception("Class index out of range!");
      }
      if ((int)secondIndex == (int)firstIndex) {
        throw new Exception("Diagonal of cost matrix non-zero!");
      }

      // Get cost factor.
      if (StreamTokenizer.TT_EOF == (currentToken = tokenizer.nextToken())) {
        throw new Exception("Premature end of file!");
      }
      if (currentToken == StreamTokenizer.TT_EOL) {
        throw new Exception("Premature end of line!");
      }
      if (currentToken != StreamTokenizer.TT_NUMBER) {
        throw new Exception("Only numbers and comments allowed "+
            "in cost file!");
      }
      weight = tokenizer.nval;
      if (!Utils.gr(weight,0)) {
        throw new Exception("Only positive weights allowed!");
      }
      setCell((int)firstIndex, (int)secondIndex,
          new Double(weight));
    }
  }

  /**
   * Reads a matrix from a reader. The first line in the file should
   * contain the number of rows and columns. Subsequent lines
   * contain elements of the matrix.
   * (FracPete: taken from old weka.core.Matrix class)
   *
   * @param     reader the reader containing the matrix
   * @throws    Exception if an error occurs
   * @see       #write(Writer)
   */
  public CostMatrix(Reader reader) throws Exception {
    LineNumberReader lnr = new LineNumberReader(reader);
    String line;
    int currentRow = -1;

    while ((line = lnr.readLine()) != null) {

      // Comments
      if (line.startsWith("%")) {
        continue;
      }

      StringTokenizer st = new StringTokenizer(line);
      // Ignore blank lines
      if (!st.hasMoreTokens()) {
        continue;
      }

      if (currentRow < 0) {
        int rows = Integer.parseInt(st.nextToken());
        if (!st.hasMoreTokens()) {
          throw new Exception("Line " + lnr.getLineNumber()
              + ": expected number of columns");
        }

        int cols = Integer.parseInt(st.nextToken());
        if (rows != cols) {
          throw new Exception("Trying to create a non-square cost "
                              +"matrix");
        }
        //        m_matrix = new Object[rows][cols];
        m_size = rows;
        initialize();
        currentRow++;
        continue;

      } else {
        if (currentRow == m_size) {
          throw new Exception("Line " + lnr.getLineNumber()
              + ": too many rows provided");
        }

        for (int i = 0; i < m_size; i++) {
          if (!st.hasMoreTokens()) {
            throw new Exception("Line " + lnr.getLineNumber()
                + ": too few matrix elements provided");
          }

          String nextTok = st.nextToken();
          // try to parse as a double first
          Double val = null;
          try {
            val = new Double(nextTok);
            double value = val.doubleValue();
          } catch (Exception ex) {
            val = null;
          }
          if (val == null) {
            setCell(currentRow, i, nextTok);
          } else {
            setCell(currentRow, i, val);
          }
        }
        currentRow++;
      }
    }

    if (currentRow == -1) {
      throw new Exception("Line " + lnr.getLineNumber()
                          + ": expected number of rows");
    } else if (currentRow != m_size) {
      throw new Exception("Line " + lnr.getLineNumber()
                          + ": too few rows provided");
    }
  }

  /**
   * Writes out a matrix. The format can be read via the
   * CostMatrix(Reader) constructor.
   * (FracPete: taken from old weka.core.Matrix class)
   *
   * @param     w the output Writer
   * @throws    Exception if an error occurs
   */
  public void write(Writer w) throws Exception {
    w.write("% Rows\tColumns\n");
    w.write("" + m_size + "\t" + m_size + "\n");
    w.write("% Matrix elements\n");
    for(int i = 0; i < m_size; i++) {
      for(int j = 0; j < m_size; j++) {
        w.write("" + getCell(i, j) + "\t");
      }
      w.write("\n");
    }
    w.flush();
  }

  /**
   * converts the Matrix into a single line Matlab string: matrix is enclosed
   * by parentheses, rows are separated by semicolon and single cells by
   * blanks, e.g., [1 2; 3 4].
   * @return      the matrix in Matlab single line format
   */
  public String toMatlab() {
    StringBuffer      result;
    int               i;
    int               n;

    result = new StringBuffer();

    result.append("[");

    for (i = 0; i < m_size; i++) {
      if (i > 0) {
        result.append("; ");
      }

      for (n = 0; n < m_size; n++) {
        if (n > 0) {
          result.append(" ");
        }
        result.append(getCell(i, n));
      }
    }

    result.append("]");

    return result.toString();
  }

  /**
   * Set the value of a particular cell in the matrix
   *
   * @param rowIndex the row
   * @param columnIndex the column
   * @param value the value to set
   */
  public final void setCell(int rowIndex, int columnIndex,
                               Object value) {
    m_matrix[rowIndex][columnIndex] = value;
  }

  /**
   * Return the contents of a particular cell. Note: this
   * method returns the Object stored at a particular cell.
   *
   * @param rowIndex the row
   * @param columnIndex the column
   * @return the value at the cell
   */
  public final Object getCell(int rowIndex, int columnIndex) {
    return m_matrix[rowIndex][columnIndex];
  }

  /**
   * Return the value of a cell as a double (for legacy code)
   *
   * @param rowIndex the row
   * @param columnIndex the column
   * @return the value at a particular cell as a double
   * @exception Exception if the value is not a double
   */
  public final double getElement(int rowIndex, int columnIndex)
    throws Exception {
    if (!(m_matrix[rowIndex][columnIndex] instanceof Double)) {
      throw new Exception("Cost matrix contains non-fixed costs!");
    }
    return ((Double)m_matrix[rowIndex][columnIndex]).doubleValue();
  }

  /**
   * Return the value of a cell as a double. Computes the
   * value for non-fixed costs using the supplied Instance
   *
   * @param rowIndex the row
   * @param columnIndex the column
   * @return the value from a particular cell
   * @exception Exception if something goes wrong
   */
  public final double getElement(int rowIndex, int columnIndex,
                                 Instance inst) throws Exception {

    if (m_matrix[rowIndex][columnIndex] instanceof Double) {
      return ((Double)m_matrix[rowIndex][columnIndex]).doubleValue();
    } else if (m_matrix[rowIndex][columnIndex] instanceof String) {
      replaceStrings();
    }

    return ((AttributeExpression)m_matrix[rowIndex][columnIndex]).
      evaluateExpression(inst);
  }

  /**
   * Set the value of a cell as a double
   *
   * @param rowIndex the row
   * @param columnIndex the column
   * @param value the value (double) to set
   */
  public final void setElement(int rowIndex, int columnIndex,
                               double value) {
    m_matrix[rowIndex][columnIndex] = new Double(value);
  }

  /**
   * creates a matrix from the given Matlab string.
   * @param matlab  the matrix in matlab format
   * @return        the matrix represented by the given string
   */
  public static Matrix parseMatlab(String matlab) throws Exception {
    return Matrix.parseMatlab(matlab);
  }

  /**
   * Converts a matrix to a string.
   * (FracPete: taken from old weka.core.Matrix class)
   *
   * @return    the converted string
   */
  public String toString() {
    // Determine the width required for the maximum element,
    // and check for fractional display requirement.
    double maxval = 0;
    boolean fractional = false;
    Object element = null;
    int widthNumber = 0;
    int widthExpression = 0;
    for (int i = 0; i < size(); i++) {
      for (int j = 0; j < size(); j++) {
        element = getCell(i, j);
        if (element instanceof Double) {
          double current = ((Double)element).doubleValue();

          if (current < 0)
            current *= -11;
          if (current > maxval)
            maxval = current;
          double fract = Math.abs(current - Math.rint(current));
          if (!fractional
              && ((Math.log(fract) / Math.log(10)) >= -2)) {
            fractional = true;
          }
        } else {
          if (element.toString().length() > widthExpression) {
            widthExpression = element.toString().length();
          }
        }
      }
    }
    if (maxval > 0) {
      widthNumber = (int)(Math.log(maxval) / Math.log(10)
                          + (fractional ? 4 : 1));
    }

    int width = (widthNumber > widthExpression)
      ? widthNumber
      : widthExpression;

    StringBuffer text = new StringBuffer();
    for (int i = 0; i < size(); i++) {
      for (int j = 0; j < size(); j++) {
        element = getCell(i, j);
        if (element instanceof Double) {
          text.append(" ").
            append(Utils.doubleToString(((Double)element).
                                        doubleValue(),
                                        width, (fractional ? 2 : 0)));
        } else {
          int diff = width - element.toString().length();
          if (diff > 0) {
            int left = diff % 2;
            left += diff / 2;
            String temp = Utils.padLeft(element.toString(),
                            element.toString().length()+left);
            temp = Utils.padRight(temp, width);
            text.append(" ").append(temp);
          } else {
            text.append(" ").
              append(element.toString());
          }
        }
      }
      text.append("\n");
    }

    return text.toString();
  }

  /**
   * Returns the revision string.
   *
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6041 $");
  }
}

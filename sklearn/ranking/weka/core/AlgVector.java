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
 *    AlgVector.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.Serializable;
import java.util.Random;

/**
 * Class for performing operations on an algebraic vector
 * of floating-point values.
 *
 * @author  Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class AlgVector 
  implements Cloneable, Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -4023736016850256591L;

  /** The values of the matrix */
  protected double[] m_Elements;

  /**
   * Constructs a vector and initializes it with default values.
   *
   * @param n 		the number of elements
   */
  public AlgVector(int n) {

    m_Elements = new double[n];
    initialize();
  }

 /**
   * Constructs a vector using a given array.
   *
   * @param array 	the values of the matrix
   */
  public AlgVector(double[] array) {
    
    m_Elements = new double[array.length];
    for (int i = 0; i < array.length; i++) {
      m_Elements[i] = array[i];
    }
  }

  /**
   * Constructs a vector using a given data format.
   * The vector has an element for each numerical attribute.
   * The other attributes (nominal, string) are ignored.
   * Random is used to initialize the attributes.
   *
   * @param format 	the data format to use
   * @param random 	for initializing the attributes
   * @throws Exception	if something goes wrong
   */
  public AlgVector(Instances format, Random random) throws Exception {
    
    int len = format.numAttributes();
    for (int i = 0; i < format.numAttributes(); i++) {
      if (!format.attribute(i).isNumeric()) len--;
    }
    if (len > 0) {
      m_Elements = new double[len];
      initialize(random);
    }
  }

  /**
   * Constructs a vector using an instance.
   * The vector has an element for each numerical attribute.
   * The other attributes (nominal, string) are ignored.
   *
   * @param instance 	with numeric attributes, that AlgVector gets build from
   * @throws Exception 	if instance doesn't have access to the data format or
   * 			no numeric attributes in the data
   */
  public AlgVector(Instance instance) throws Exception {
    
    int len = instance.numAttributes();
    for (int i = 0; i < instance.numAttributes(); i++) {
      if (!instance.attribute(i).isNumeric())
	len--;
    }
    if (len > 0) {
      m_Elements = new double[len];
      int n = 0;
      for (int i = 0; i < instance.numAttributes(); i++) {
	if (!instance.attribute(i).isNumeric())
	  continue;
	m_Elements[n] = instance.value(i);
	n++;
      }
    }
    else {
      throw new IllegalArgumentException("No numeric attributes in data!");
    }
  }

  /**
   * Creates and returns a clone of this object.
   *
   * @return 		a clone of this instance.
   * @throws CloneNotSupportedException if an error occurs
   */
  public Object clone() throws CloneNotSupportedException {

    AlgVector v = (AlgVector)super.clone();
    v.m_Elements = new double[numElements()];
    for (int i = 0; i < numElements(); i++) {
        v.m_Elements[i] = m_Elements[i];
      }
    
    return v;
  }

  /**
   * Resets the elements to the default value which is 0.0.
   */
  protected void initialize() {

    for (int i = 0; i < m_Elements.length; i++) {
      m_Elements[i] = 0.0;
    }
  }

  /**
   * Initializes the values with random numbers between 0 and 1.
   * 
   * @param random	the random number generator to use for initializing
   */
  protected void initialize(Random random) {

    for (int i = 0; i < m_Elements.length; i++) {
      m_Elements[i] = random.nextDouble();
    }
  }

  /**
   * Returns the value of a cell in the matrix.
   *
   * @param index 	the row's index
   * @return 		the value of the cell of the vector
   */
  public final double getElement(int index) {
    return m_Elements[index];
  }


  /**
   * Returns the number of elements in the vector.
   *
   * @return 		the number of rows
   */
  public final int numElements() {
  
    return m_Elements.length;
  }


  /**
   * Sets an element of the matrix to the given value.
   *
   * @param index 	the elements index
   * @param value 	the new value
   */
  public final void setElement(int index, double value) {
    
    m_Elements[index] = value;
  }

  /**
   * Sets the elements of the vector to values of the given array.
   * Performs a deep copy.
   *
   * @param elements 	an array of doubles
   */
  public final void setElements(double[] elements) {

    for (int i = 0; i < elements.length; i++) {
      m_Elements[i] = elements[i];
    }
  }
  
  /**
   * Gets the elements of the vector and returns them as double array.
   *
   * @return 		an array of doubles
   */
  public double[] getElements() {

    double [] elements = new double[this.numElements()];
    for (int i = 0; i < elements.length; i++) {
      elements[i] = m_Elements[i];
    }
    return elements;
  }

  /**
   * Gets the elements of the vector as an instance.
   * !! NON-numeric data is ignored sofar
   * 
   * @param model 	the dataset structure to fit the data to
   * @param random 	in case of nominal values a random label is taken
   * @return 		an array of doubles
   * @throws Exception	if length of vector is not number of numerical attributes
   */
  public Instance getAsInstance(Instances model, Random random) 
    throws Exception {

    Instance newInst = null;

    if (m_Elements != null) {
      newInst = new DenseInstance(model.numAttributes());
      newInst.setDataset(model);
      
      for (int i = 0, j = 0; i < model.numAttributes(); i++) {
	if (model.attribute(i).isNumeric()) {
	  if (j >= m_Elements.length)
	    throw new Exception("Datatypes are not compatible."); 
	  newInst.setValue(i, m_Elements[j++]);
	}
	if (model.attribute(i).isNominal() || model.attribute(i).isRanking()) {
	  int newVal = (int) 
	    (random.nextDouble() * (double) (model.attribute(i).numValues()));
	  if (newVal == (int) model.attribute(i).numValues())
	    newVal -= 1;
	  newInst.setValue(i, newVal);
	}
      }
    }
    return newInst;
  }
    
  /**
   * Returns the sum of this vector with another.
   *
   * @param other 	the vector to add
   * @return 		a vector containing the sum.
   */
  public final AlgVector add(AlgVector other) {
  
    AlgVector b = null;

    if (m_Elements != null) {
      int n = m_Elements.length;
       try {
	b = (AlgVector)clone();
      } catch (CloneNotSupportedException ex) {
	b = new AlgVector(n);
      }
    
      for(int i = 0; i < n; i++) {
	b.m_Elements[i] = m_Elements[i] + other.m_Elements[i];
      }
    }
    
    return b;
  }

  /**
   * Returns the difference of this vector minus another.
   *
   * @param other 	the vector to subtract
   * @return 		a vector containing the difference vector.
   */
  public final AlgVector substract(AlgVector other) {
  
    int n = m_Elements.length;
    AlgVector b;
    try {
      b = (AlgVector)clone();
    } catch (CloneNotSupportedException ex) {
      b = new AlgVector(n);
    }
    
    for(int i = 0; i < n; i++) {
      b.m_Elements[i] = m_Elements[i] - other.m_Elements[i];
    }
    
    return b;
  }
 
  /**
   * Returns the inner (or dot) product of two vectors
   *
   * @param b 		the multiplication matrix
   * @return 		the double representing the dot product
   */
  public final double dotMultiply(AlgVector b) {
   
    double sum = 0.0;
 
    if (m_Elements != null) {
      int n = m_Elements.length;
      
      for(int i = 0; i < n; i++) {
	sum += m_Elements[i] * b.m_Elements[i];
      }
    }
    
    return sum;
  }

  /**
   * Computes the scalar product of this vector with a scalar
   *
   * @param s 		the scalar
   */
  public final void scalarMultiply(double s) {
   
    if (m_Elements != null) {
      int n = m_Elements.length;
      
      for(int i = 0; i < n; i++) {
	m_Elements[i] = s * m_Elements[i];
      }
    }
  }

  /**
   * Changes the length of a vector.
   *
   * @param len 	the new length of the vector
   */
  public void changeLength(double len) {
   
    double factor = this.norm();
    factor = len / factor;
    scalarMultiply(factor);
  }

  /**
   * Returns the norm of the vector
   *
   * @return 		the norm of the vector
   */
  public double norm() {
   
    if (m_Elements != null) {
      int n = m_Elements.length;
      double sum = 0.0;
      
      for(int i = 0; i < n; i++) {
	sum += m_Elements[i] * m_Elements[i];
      }
    return Math.pow(sum, 0.5);
    }
    else return 0.0;
  }

  /**
   * Norms this vector to length 1.0
   */
  public final void normVector() {
   
    double len = this.norm();
    this.scalarMultiply(1 / len);
  }

  /** 
   * Converts a vector to a string
   *
   * @return 		the converted string
   */
  public String toString() {

    StringBuffer text = new StringBuffer();
    for (int i = 0; i < m_Elements.length; i++) {
      if (i > 0) text.append(",");
      text.append(Utils.doubleToString(m_Elements[i],6));
    }

    text.append("\n");
    return text.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  /**
   * Main method for testing this class, can take an ARFF file as first argument.
   * 
   * @param args 	commandline options
   * @throws Exception	if something goes wrong in testing
   */
  public static void main(String[] args) throws Exception {
    
    double[] first = {2.3, 1.2, 5.0};
    
    try {
      AlgVector test = new AlgVector(first);
      System.out.println("test:\n " + test);
     
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}

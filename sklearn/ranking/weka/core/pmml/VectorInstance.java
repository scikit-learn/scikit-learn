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
 *    VectorInstance.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * Class encapsulating a PMML VectorInstance construct
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6466 $
 */
public class VectorInstance implements Serializable {
  
  /** For serialization */
  private static final long serialVersionUID = -7543200367512646290L;

  /** The ID of this instance */
  protected String m_ID;
  
  /** The usually sparse elements of this vector */
  protected Array m_values;
  
  /** The fields indexed by this VectorInstance */
  protected List<FieldRef> m_vectorFields;
  
  /**
   * Constructor
   * 
   * @param values the Array of values for this vector instance
   * @param vectorFields the mining fields indexed by this vector instance
   */
  public VectorInstance(Array values, List<FieldRef> vectorFields) {
    m_values = values;
    m_vectorFields = vectorFields;
  }
  
  /**
   * Constructor
   * 
   * @param vecElement PMML element containing this vector instance
   * @param vectorFields the mining fields indexed by this vector instance
   * @throws Exception if something goes wrong
   */
  public VectorInstance(Element vecElement, List<FieldRef> vectorFields) 
    throws Exception {
    
    m_vectorFields = vectorFields;
    
    // get the ID
    
    String id = vecElement.getAttribute("id");
    if (id == null || id.length() == 0) {
      throw new Exception("[VectorInstance] no ID attribute defined!");
    }
    
    m_ID = id;
    
    // check for both types of array
    NodeList s_arrL = vecElement.getElementsByTagName("REAL-SparseArray");
    NodeList d_arrL = vecElement.getElementsByTagName("REAL-ARRAY");
    
    if (s_arrL.getLength() == 0 && d_arrL.getLength() == 0) {
      throw new Exception("[VectorInstance] no arrays defined!");
    }
    
    NodeList arrL = (s_arrL.getLength() > 0)
      ? s_arrL
      : d_arrL;
    
    // should be just one array per vector instance
    Element theArray = (Element)arrL.item(0);
    
    m_values = Array.create(theArray);
  }
  
  /**
   * Get the ID of this vector instance
   * 
   * @return the ID of this vector instance
   */
  public String getID() {
    return m_ID;
  }
  
  /**
   * Get the Array of values encapsulated in this vector instance
   * 
   * @return the Array of values encapsulated in this vector instance
   */
  public Array getValues() {
    return m_values;
  }
  
  /**
   * Get the mining fields that are indexed by this vector instance
   * 
   * @return the mining fields that are indexed by this vector instance
   */
  public List<FieldRef> getVectorFields() {
    return m_vectorFields;
  }
  
  /**
   * Subtract the values in the supplied array from this vector instance
   * 
   * @param other an array of values
   * @return a new VectorInstance containing the result of the subtraction
   * @throws Exception if something goes wrong
   */
  public VectorInstance subtract(double[] other) throws Exception {
    // other is a non-sparse array of values
    
    ArrayList<Object> diffVals = new ArrayList<Object>();
    for (int i = 0; i < other.length; i++) {
      double x = m_values.valueDouble(i);
      double y = other[i];
//      System.err.println("x: " + x + " y: " + y);
      double result = x - y;
      diffVals.add(new Double(result));
    }
    
    Array newArray = Array.create(diffVals, null);
    
    return new VectorInstance(newArray, m_vectorFields);
  }
  
  /**
   * Subtract the supplied VectorInstance from this one and return the
   * result as a new VectorInstance
   * 
   * @param other the vector instance to subtract
   * @return a new VectorInstance containing the result of the subtraction
   * @throws Exception if something goes wrong
   */
  public VectorInstance subtract(VectorInstance other) throws Exception {
    // IMPORTANT: assumes that other has the same vector fields
    // as this vector instance. Otherwise results will be meaningless
    if (m_vectorFields.size() != other.getVectorFields().size()) {
      throw new Exception("[VectorInstance.dotProduct] supplied vector instance does " +
                "not have the same number of vector fields as this vector instance!");
    }
    
    ArrayList<Object> diffVals = new ArrayList<Object>();
    for (int i = 0; i < m_vectorFields.size(); i++) {
      double x = m_values.valueDouble(i);
      double y = other.getValues().valueDouble(i);
      double result = x - y;
      diffVals.add(new Double(result));
    }
    
    Array newArray = Array.create(diffVals, null);
    
    return new VectorInstance(newArray, m_vectorFields);
  }
  
  /**
   * Computes the dot product between this vector instance and the
   * argument
   * 
   * @param other the vector instance to perform the dot product with
   * @return the dot product of this and the supplied vector instance
   * @throws Exception if something goes wrong
   */
  public double dotProduct(VectorInstance other) throws Exception {
    // IMPORTANT: assumes that other has the same vector fields
    // as this vector instance. Otherwise results will be meaningless
    if (m_vectorFields.size() != other.getVectorFields().size()) {
      throw new Exception("[VectorInstance.dotProduct] supplied vector instance does " +
      		"not have the same number of vector fields as this vector instance!");
    }
    double result = 0;
    
    Array otherValues = other.getValues();
    
    // do a fast dot product
    int n1 = m_values.numValues();
    int n2 = otherValues.numValues();
    
    for (int p1 = 0, p2 = 0; p1 < n1 && p2 < n2;) {
      int ind1 = m_values.index(p1);
      int ind2 = otherValues.index(p2);
      
      if (ind1 == ind2) {
//        System.err.println("Here..." + m_values.valueSparseDouble(p1) + " " + otherValues.valueSparseDouble(p2));
        result += m_values.valueSparseDouble(p1) * otherValues.valueSparseDouble(p2);
        p1++;
        p2++;
      } else if (ind1 > ind2) {
        p2++;
      } else {
        p1 ++;
      }      
    }
    
    return result;
  }
  
  /**
   * Computes the dot product between this vector instance and the
   * supplied array of values
   * 
   * @param other an array of values to perform the dot product with
   * @return the dot product of this vector instance with the argument
   * @throws Exception if something goes wrong
   */
  public double dotProduct(double[] other) throws Exception {
    // other is a non-sparse array of values
    
    double result = 0;
    // do a fast dot product
    int n1 = m_values.numValues();
    for (int i = 0; i < n1; i++) {
      int ind1 = m_values.index(i);
      
      result += m_values.valueSparseDouble(i) * other[ind1];
    }
    
    return result;
  }
}

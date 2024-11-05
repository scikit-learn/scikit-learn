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
 *    Array.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Element;

/**
 * Class for encapsulating a PMML Array element.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6467 $
 */
public class Array implements Serializable {
  
  /**
   * Utility method to check if an XML element is an array.
   * 
   * @param arrayE the XML element to check
   * @return returns true if the XML element is an array
   */
  public static boolean isArray(Element arrayE) {
    String name = arrayE.getTagName();

    if (name.equals("Array") || name.equals("NUM-ARRAY") || name.equals("INT-ARRAY")
        || name.equals("REAL-ARRAY") || name.equals("STRING-ARRAY") 
        || isSparseArray(arrayE)) {
      return true;
    }
    return false;
  }
  
  /**
   * Utility method to check if an XML element is a sparse array.
   * 
   * @param arrayE the XML element to check.
   * @return true if the XML element is a sparse array.
   */
  private static boolean isSparseArray(Element arrayE) {
    String name = arrayE.getTagName();
    
    if (name.equals("INT-SparseArray") || name.equals("REAL-SparseArray")) {
      return true;
    }
    
    return false;
  }
  
  public static Array create(List<Object> values, 
      List<Integer> indices) throws Exception {
    
    ArrayType type = null;
    
    Object first = values.get(0);
    if ((first instanceof Double) || (first instanceof Float)) {
      type = ArrayType.REAL;
    } else if ((first instanceof Integer) || (first instanceof Long)) {
      type = ArrayType.INT;
    } else if ((first instanceof String)) {
      type = ArrayType.STRING;
    } else {
      throw new Exception("[Array] unsupport type!");
    }
    
    if (indices != null) {
      // array is sparse
      
      if (indices.size() != values.size()) {
        throw new Exception("[Array] num values is not equal to num indices!!");
      }
      
      if (type == ArrayType.REAL) {
        type = ArrayType.REAL_SPARSE;
      } else if (type == ArrayType.INT) {
        type = ArrayType.INT_SPARSE;
      } else {
        throw new Exception("[Array] sparse arrays can only be integer, long, float or double!");
      }
      
      return new SparseArray(type, values, indices);
    }
    
    return new Array(type, values);
  }
  
  /**
   * Static factory method for creating non-sparse or sparse
   * array types as needed.
   * 
   * @param arrayE the XML element encapsulating the array
   * @return an appropriate Array type
   * @throws Exception if there is a problem when constructing the array
   */
  public static Array create(Element arrayE) throws Exception {
    if (!isArray(arrayE)) {
      throw new Exception("[Array] the supplied element does not contain an array!");
    }
    
    if (isSparseArray(arrayE)) {
      return new SparseArray(arrayE);
    }
     
    return new Array(arrayE);
  }
  
  public static enum ArrayType {
    NUM("NUM-ARRAY"),
    INT("INT-ARRAY"),
    REAL("REAL-ARRAY"),
    STRING("STRING-ARRAY"),
    REAL_SPARSE("REAL-SparseArray"),
    INT_SPARSE("INT-SparseArray");
    
    private final String m_stringVal;
    
    ArrayType(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  /** The values of the array */
  protected ArrayList<String> m_values = new ArrayList<String>();
  
  /** The type of the array */
  protected ArrayType m_type = ArrayType.NUM;
  
  protected void initialize(Element arrayE) throws Exception {
    String arrayS = arrayE.getTagName();
    
    // get the type of the array
    if (arrayS.equals("Array")) {
      String type = arrayE.getAttribute("type"); 
      if (type.equals("int")) {
        m_type = ArrayType.INT;
      } else if (type.equals("real")) {
        m_type = ArrayType.REAL;
      } else if (type.equals("string")) {
        m_type = ArrayType.STRING;
      }
    } else {
      for (ArrayType a : ArrayType.values()) {
        if (a.toString().equals(arrayS)) {
          m_type = a;
          break;
        }
      }
    }
    // now read the values
    String contents = arrayE.getChildNodes().item(0).getNodeValue();
    StringReader sr = new StringReader(contents);
    StreamTokenizer st = new StreamTokenizer(sr);
    st.resetSyntax();
    st.whitespaceChars(0, ' ');
    st.wordChars(' '+1,'\u00FF');
    st.whitespaceChars(' ',' ');
    st.quoteChar('"');
    st.quoteChar('\'');
    //m_Tokenizer.eolIsSignificant(true);
    
    st.nextToken();
    while (st.ttype != StreamTokenizer.TT_EOF && 
        st.ttype != StreamTokenizer.TT_EOL) {
      m_values.add(st.sval);
      st.nextToken();
    }
  }
  
  /**
   * Construct an array from an XML node
   * 
   * @param arrayE the Element containing the XML
   * @throws Exception if something goes wrong
   */
  protected Array(Element arrayE) throws Exception {
    initialize(arrayE);
  }
  
  /**
   * Construct an array from the given values.
   * 
   * @param type the type of the elements.
   * @param values the values of the array.
   */
  protected Array(ArrayType type, List<Object> values) {
    m_values = new ArrayList<String>();
    m_type = type;
    
    for (Object o : values) {
      m_values.add(o.toString());
    }
  }
  
  /**
   * Get the type of this array.
   * 
   * @return the type of the array.
   */
  public ArrayType getType() {
    return m_type;
  }
  
  /**
   * Is this array a SparseArray?
   * 
   * @return true if this array is sparse.
   */
  public boolean isSparse() {
    return false;
  }
  
  /**
   * Get the number of values in this array.
   * 
   * @return the number of values in this array.
   */
  public int numValues() {
    return m_values.size();
  }
  
  /**
   * Returns true if the array contains this string value.
   * 
   * @param value the value to check for.
   * @return true if the array contains this string value
   */
  public boolean contains(String value) {
    return m_values.contains(value);
  }
  
  /**
   * Returns true if the array contains this integer value.
   * 
   * @param value the value to check for
   * @return true if the array contains this integer value
   */
  public boolean contains(int value) {
    return contains(new Integer(value).toString());
  }
  
  /**
   * Returns true if the array contains this real value.
   * 
   * @param value the value to check for
   * @return true if the array contains this real value
   */
  public boolean contains(double value) {
    return contains(new Double(value).toString());
  }
  
  /**
   * Returns true if the array contains this real value.
   * 
   * @param value the value to check for
   * @return true if the array contains this real value
   */
  public boolean contains(float value) {
    return contains(new Float(value).toString());
  }
  
  private void checkInRange(int index) throws Exception {
    if (index >= m_values.size() || index < 0) {
      throw new IllegalArgumentException("[Array] index out of range " + index);
    }
  }
  
  /**
   * Returns the index of the value stored at the given position
   * 
   * @param position the position
   * @return the index of the value stored at the given position
   */
  public int index(int position) {
    return position; // position is the index for dense arrays
  }
  
  /**
   * Gets the value at index from the array.
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as as String.
   * @throws Exception if index is out of bounds.
   */
  public String value(int index) throws Exception {
    return actualValue(index);
  }
  
  /**
   * Gets the value at index from the array
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as as String.
   * @throws Exception if index is out of bounds.
   */
  protected String actualValue(int index) throws Exception {
    checkInRange(index);
    return m_values.get(index);
  }
  
  /**
   * Gets the value at index from the array as a String. Calls
   * value().
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as a String.
   * @throws Exception if index is out of bounds.
   */
  public String valueString(int index) throws Exception {
    return value(index);
  }
  
  /**
   * Gets the value at index from the array as a double.
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as a double.
   * @throws Exception if index is out of bounds.
   */
  public double valueDouble(int index) throws Exception {
    if (m_type == ArrayType.STRING) {
      throw new Exception("[Array] Array does not contain numbers!");
    }
    return Double.parseDouble(value(index));
  }
  
  /**
   * Gets the value at index from the array as a float.
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as a float.
   * @throws Exception if index is out of bounds.
   */
  public float valueFloat(int index) throws Exception {
    if (m_type == ArrayType.STRING) {
      throw new Exception("[Array] Array does not contain numbers!");
    }
    return Float.parseFloat(value(index));
  }
  
  /**
   * Gets the value at index from the array as an int.
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as an int.
   * @throws Exception if index is out of bounds.
   */
  public int valueInt(int index) throws Exception {
    if (m_type != ArrayType.INT && m_type != ArrayType.INT_SPARSE) {
      throw new Exception("[Array] Array does not contain integers!");
    }
    return Integer.parseInt(value(index));
  }
  
  /**
   * Gets the value at indexOfIndex from the array. Does the
   * same as value() if this array is not sparse.
   * 
   * @param indexOfIndex the index of the index of the value to get.
   * @return a value from the array as a String.
   * @throws Exception if indexOfIndex is out of bounds.
   */
  public String valueSparse(int indexOfIndex) throws Exception {
    return actualValue(indexOfIndex);
  }
  
  /**
   * Gets the value at indexOfIndex from the array. Does the
   * same as value() if this array is not sparse.
   * 
   * @param indexOfIndex the index of the index of the value to get.
   * @return a value from the array as a String.
   * @throws Exception if indexOfIndex is out of bounds.
   */
  public String valueSparseString(int indexOfIndex) throws Exception {
    return valueSparse(indexOfIndex);
  }
  
  /**
   * Gets the value at indexOfIndex from the array. Does the
   * same as value() if this array is not sparse.
   * 
   * @param indexOfIndex the index of the index of the value to get.
   * @return a value from the array as a double.
   * @throws Exception if indexOfIndex is out of bounds.
   */
  public double valueSparseDouble(int indexOfIndex) throws Exception {
    return Double.parseDouble(actualValue(indexOfIndex));
  }
  
  /**
   * Gets the value at indexOfIndex from the array. Does the
   * same as value() if this array is not sparse.
   * 
   * @param indexOfIndex the index of the index of the value to get.
   * @return a value from the array as a float.
   * @throws Exception if indexOfIndex is out of bounds.
   */
  public float valueSparseFloat(int indexOfIndex) throws Exception {
    return Float.parseFloat(actualValue(indexOfIndex));
  }
  
  /**
   * Gets the value at indexOfIndex from the array. Does the
   * same as value() if this array is not sparse.
   * 
   * @param indexOfIndex the index of the index of the value to get.
   * @return a value from the array as an int.
   * @throws Exception if indexOfIndex is out of bounds.
   */
  public int valueSparseInt(int indexOfIndex) throws Exception {
    return Integer.parseInt(actualValue(indexOfIndex));
  }
  
  public String toString() {
    StringBuffer text = new StringBuffer();
    
    text.append("<");
    for (int i = 0; i < m_values.size(); i++) {
      text.append(m_values.get(i));
      if (i < m_values.size() - 1) {
        text.append(",");
      }
    }
    
    text.append(">");
    return text.toString();
  }
}

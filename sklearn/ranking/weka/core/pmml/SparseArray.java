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
 *    SparseArray.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * Implementation of a sparse array. Extends Array.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6466 $
 */
public class SparseArray extends Array implements Serializable {
  
  /** For serialization */
  private static final long serialVersionUID = 8129550573612673674L;

  /** The size of the array if known */
  protected int m_numValues;
  
  /** The number of non-zero elements */
  protected int m_numNonZero;
  
  /** The indices of the sparse array */
  protected List<Integer> m_indices;// = new ArrayList<Integer>();
  
  /**
   * Initializes the array from the supplied XML element.
   * 
   * @param arrayE the XML containing the array
   * @throws Exception if something goes wrong
   */
  protected void initialize(Element arrayE) throws Exception {
    m_indices = new ArrayList<Integer>();
    
    String arrayS = arrayE.getTagName();
    String entriesName = null;
    
    if (arrayS.equals(ArrayType.REAL_SPARSE.toString())) {
      m_type = ArrayType.REAL_SPARSE;
      entriesName = "REAL-Entries";
    } else {
      m_type = ArrayType.INT_SPARSE;
      entriesName = "INT-Entries";
    }
    
    // see if we can get the "n" attribute to determine the
    // size
    String N = arrayE.getAttribute("n");
    if (N != null && N.length() > 0) {
      m_numValues = Integer.parseInt(N);
    }
    
    // get the values
    NodeList v = arrayE.getElementsByTagName(entriesName);
    if (v == null || v.getLength() == 0) {
      // there are no entries (and indices), so this
      // array must contain all zeros  
      m_numNonZero = 0;
    } else {
      Element entries = (Element)v.item(0);
      String contents = entries.getChildNodes().item(0).getNodeValue();
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
      
      // get the indices
      NodeList i = arrayE.getElementsByTagName("Indices");
      Element indices = (Element)i.item(0);
      contents = indices.getChildNodes().item(0).getNodeValue();
      sr = new StringReader(contents);
      st = new StreamTokenizer(sr);
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
        Integer newInt = new Integer(Integer.parseInt(st.sval) - 1);
        m_indices.add(newInt);
        st.nextToken();
      }
      
      m_numNonZero = m_indices.size();
    }
    
  }

  /**
   * Construct a sparse array from an XML node
   * 
   * @param arrayE the Element containing the XML
   * @throws Exception if something goes wrong
   */
  protected SparseArray(Element arrayE) throws Exception {
    super(arrayE);
  }
  
  /**
   * Construct a sparse array from the given values and indices
   * 
   * @param type the type of the array elements
   * @param values the values of the array
   * @param indices the indices of the array
   */
  protected SparseArray(ArrayType type, List<Object> values, 
      List<Integer> indices) {
    super(type, values);
    
    m_indices = indices;
  }
  
  /**
   * Overrides isSparse() in Array and always returns true.
   * 
   * @return true always
   */
  public boolean isSparse() {
    return true;
  }
  
  /**
   * Get the number of values in this array.
   * 
   * @return the number of values in this array.
   */
  public int numValues() {
    return m_numValues;
  }
  
  /**
   * Get the number of non-zero values in this sparse array
   * 
   * @return the number of values that are non-zero
   */
  public int numNonZero() {
    return m_numNonZero;
  }
  
  /**
   * Returns the index of the value stored at the given position
   * 
   * @param position the position
   * @return the index of the value stored at the given position
   */
  public int index(int position) {
    return m_indices.get(position);
  }
  
  /**
   * Locates the greatest index that is not greater than the
   * given index.
   *
   * @return the internal index of the attribute index. Returns
   * -1 if no index with this property could be found
   */
  public int locateIndex(int index) {

    int min = 0, max = m_indices.size() - 1;

    if (max == -1) {
      return -1;
    }

    // Binary search
    while ((m_indices.get(min) <= index) 
        && (m_indices.get(max) >= index)) {
      int current = (max + min) / 2;
      if (m_indices.get(current) > index) {
        max = current - 1;
      } else if (m_indices.get(current) < index) {
        min = current + 1;
      } else {
        return current;
      }
    }
    if (m_indices.get(max) < index) {
      return max;
    } else {
      return min - 1;
    }
  }
  
  /**
   * Gets the value at index from the array.
   * 
   * @param index the index of the value to get.
   * @return the value at index in the array as as String.
   * @throws Exception if index is out of bounds.
   */
  public String value(int arrIndex) throws Exception {
    int index = locateIndex(arrIndex);
    if (index >= 0 && (m_indices.get(index) == arrIndex)) {

      return m_values.get(index);
    } else {
      return "0";
    }
  }
  
  public String toString() {
    StringBuffer text = new StringBuffer();
    
    text.append("<");
    for (int i = 0; i < m_indices.size(); i++) {
      text.append(m_indices.get(i).toString() + " ");
      text.append(m_values.get(i));
      if (i < m_indices.size() - 1) {
        text.append(",");
      }
    }
    
    text.append(">");
    return text.toString();
  }
}

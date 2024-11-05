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
 *    NominalItem.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.io.Serializable;

import weka.core.Attribute;


/**
 * Class that encapsulates a nominal item.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6543 $
 */
public class NominalItem extends Item implements Serializable {

  /** For serialization */
  private static final long serialVersionUID = 2182122099990462066L;
  
  /** The index of the value considered to be positive */
  protected int m_valueIndex;
  
  /**
   * Constructs a new NominalItem.
   * 
   * @param att the attribute that backs the item.
   * @param valueIndex the index of the value for this item.
   * @throws Exception if the NominalItem can't be constructed.
   */
  public NominalItem(Attribute att, int valueIndex) throws Exception {
    
    super(att);
    
    if (att.isNumeric()) {
      throw new Exception("NominalItem must be constructed using a nominal attribute");
    }
    m_attribute = att;
    if (m_attribute.numValues() == 1) {
      m_valueIndex = 0; // unary attribute (? used to indicate absence from a basket)
    } else {
      m_valueIndex = valueIndex;
    }
  }
  
  /**
   * Get the value index for this item.
   * 
   * @return the value index.
   */
  public int getValueIndex() {
    return m_valueIndex;
  }
  
  /**
   * Get this item's value as a String.
   * 
   * @return this item's value as a String.
   */
  public String getItemValueAsString() {
    return m_attribute.value(m_valueIndex);
  }
  
  /**
   * Get this item's comparison operator as a String.
   * 
   * @return this item's comparison operator as a String.
   */
  public String getComparisonAsString() {
    return "=";
  }
    
  /**
   * A string representation of this item, (i.e.
   * <attribute name> <comparison operator> <item value>).
   * This default implementation just prints the attribute
   * name and (optionally) frequency information.
   * 
   * @param freq true if the frequency should be included.
   * @return a string representation of this item. 
   */
  public String toString(boolean freq) {
    String result = m_attribute.name() + "=" + m_attribute.value(m_valueIndex);
    if (freq) {
      result += ":" + m_frequency;
    }
    return result;
  }
  
  /**
   * Equals. Just compares attribute and valueIndex.
   * @return true if this NominalItem is equal to the argument.
   */
  public boolean equals(Object compareTo) {
    if (!(compareTo instanceof NominalItem)) {
      return false;
    }
    
    NominalItem b = (NominalItem)compareTo;
    if (m_attribute.equals(b.getAttribute()) && 
//        m_frequency == b.getFrequency() && 
        m_valueIndex == b.getValueIndex()) {
      return true;
    }
    
    return false;
  }
}

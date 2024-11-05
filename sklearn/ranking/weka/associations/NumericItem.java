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
 *    NumericItem.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.io.Serializable;

import weka.core.Attribute;
import weka.core.Utils;

/**
 * Class that encapsulates a numeric item.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6543 $
 */
public class NumericItem extends Item implements Serializable {
  
  /** For serialization */
  private static final long serialVersionUID = -7869433770765864800L;

  public static enum Comparison {
    NONE,
    EQUAL,
    LESS_THAN_OR_EQUAL_TO,
    GREATER_THAN;
  }
  
  /** The numeric test */
  protected double m_splitPoint;
  
  /** The comparison operator */
  protected Comparison m_comparison;

  /**
   * Constructs a new <code>NumericItem</code>
   * 
   * @param att the attribute that backs the item.
   * @param splitPoint the numeric test value.
   * @param comp the comparison operator.
   * @throws Exception if the NumericItem can't be constructed.
   */
  public NumericItem(Attribute att, double splitPoint, Comparison comp) 
    throws Exception {
    super(att);
    
    if (!att.isNumeric()) {
      throw new Exception("NumericItem must be constructed using a numeric attribute");
    }
    
    m_comparison = comp;
    m_splitPoint = splitPoint;
  }
  
  /**
   * Gets the numeric test.
   * 
   * @return the numeric test value for this item.
   */
  public double getSplitPoint() {
    return m_splitPoint;
  }
  
  /**
   * Gets the comparison operator for this item.
   * 
   * @return the comparison operator for this item.
   */
  public Comparison getComparison() {
    return m_comparison;
  }
  
  /**
   * Get this item's value as a String.
   * 
   * @return this item's value as a String.
   */
  public String getItemValueAsString() {
    return Utils.doubleToString(m_splitPoint, 3);
  }
  
  /**
   * Get this item's comparison operator as a String.
   * 
   * @return this item's comparison operator as a String.
   */
  public String getComparisonAsString() {
    String result = null;
    
    switch (m_comparison) {
    case EQUAL:
      result = "=";
      break;
    case LESS_THAN_OR_EQUAL_TO:
      result = "<=";
      break;
    case GREATER_THAN:
      result = ">";
      break;
    }
    
    return result;
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
    StringBuffer result = new StringBuffer();
    result.append(m_attribute.name() + " ");
    switch (m_comparison) {
    case EQUAL:
      result.append("=");
      break;
    case LESS_THAN_OR_EQUAL_TO:
      result.append("<=");
      break;
    case GREATER_THAN:
      result.append(">");
      break;
    }
    
    result.append(" " + Utils.doubleToString(m_splitPoint, 4));
    if (freq) {
      result.append(":" + m_frequency);
    }
    
    return result.toString();
  }
  
  /**
   * Equals. Compares the attribute, numeric test and comparison
   * operator
   * 
   * @return true if this NumericItem is equal to the argument.
   */
  public boolean equals(Object compareTo) {
    if (!(compareTo instanceof NumericItem)) {
      return false;
    }
    
    NumericItem b = (NumericItem)compareTo;
    if (m_attribute.equals(b.getAttribute()) && 
        m_comparison == b.getComparison() && 
        (new Double(m_splitPoint).equals(new Double(b.getSplitPoint())))) {
      return true;
    }
    
    return false;
  }
}

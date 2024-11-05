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
 * StringLocator.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.io.Serializable;
import java.util.Vector;

/**
 * This class locates and records the indices of a certain type of attributes, 
 * recursively in case of Relational attributes.
 * 
 * @author fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see Attribute#RELATIONAL
 */
public class AttributeLocator 
  implements Serializable, Comparable<AttributeLocator>, RevisionHandler {
  
  /** for serialization */
  private static final long serialVersionUID = -2932848827681070345L;

  /** the attribute indices that may be inspected */
  protected int[] m_AllowedIndices = null;
  
  /** contains the attribute locations, either true or false Boolean objects */
  protected Vector<Boolean> m_Attributes = null;
  
  /** contains the locator locations, either null or a AttributeLocator reference */
  protected Vector<AttributeLocator> m_Locators = null;

  /** the type of the attribute */
  protected int m_Type = -1;
  
  /** the referenced data */
  protected Instances m_Data = null;

  /** the indices */
  protected int[] m_Indices = null;

  /** the indices of locator objects */
  protected int[] m_LocatorIndices = null;
  
  /**
   * Initializes the AttributeLocator with the given data for the specified
   * type of attribute. Checks all attributes.
   * 
   * @param data	the data to work on
   * @param type	the type of attribute to locate
   */
  public AttributeLocator(Instances data, int type) {
    this(data, type, 0, data.numAttributes() - 1);
  }
  
  /**
   * Initializes the AttributeLocator with the given data for the specified
   * type of attribute. Checks only the given range.
   * 
   * @param data	the data to work on
   * @param type	the type of attribute to locate
   * @param fromIndex	the first index to inspect (including)
   * @param toIndex	the last index to check (including)
   */
  public AttributeLocator(Instances data, int type, int fromIndex, int toIndex) {
    super();

    int[] indices = new int[toIndex - fromIndex + 1];
    for (int i = 0; i < indices.length; i++)
      indices[i] = fromIndex + i;
    
    initialize(data, type, indices);
  }
  
  /**
   * initializes the AttributeLocator with the given data for the specified
   * type of attribute. Checks only the given attribute indices.
   * 
   * @param data	the data to work on
   * @param type	the type of attribute to locate
   * @param indices	the attribute indices to check
   */
  public AttributeLocator(Instances data, int type, int[] indices) {
    super();

    initialize(data, type, indices);
  }
  
  /**
   * initializes the AttributeLocator
   * 
   * @param data	the data to base the search for attributes on
   * @param type	the type of attribute to look for
   * @param indices	the indices that are allowed to check
   */
  protected void initialize(Instances data, int type, int[] indices) {
    m_Data = new Instances(data, 0);
    m_Type = type;
    
    m_AllowedIndices = new int[indices.length];
    System.arraycopy(indices, 0, m_AllowedIndices, 0, indices.length);
    
    locate();

    m_Indices        = find(true);
    m_LocatorIndices = find(false);
  }
  
  /**
   * returns the type of attribute that is located
   * 
   * @return		the type of attribute
   */
  public int getType() {
    return m_Type;
  }
  
  /**
   * returns the indices that are allowed to check for the attribute type
   * 
   * @return 		the indices that are checked for the attribute type
   */
  public int[] getAllowedIndices() {
    return m_AllowedIndices;
  }
  
  /**
   * sets up the structure
   */
  protected void locate() {
    int         i;
    
    m_Attributes = new Vector<Boolean>();
    m_Locators   = new Vector<AttributeLocator>();
    
    for (i = 0; i < m_AllowedIndices.length; i++) {
      if (m_Data.attribute(m_AllowedIndices[i]).type() == Attribute.RELATIONAL)
	m_Locators.add(new AttributeLocator(m_Data.attribute(m_AllowedIndices[i]).relation(), getType()));
      else
	m_Locators.add(null);
      
      if (m_Data.attribute(m_AllowedIndices[i]).type() == getType())
        m_Attributes.add(new Boolean(true));
      else
        m_Attributes.add(new Boolean(false));
    }
  }
  
  /**
   * returns the underlying data
   * 
   * @return      the underlying Instances object
   */
  public Instances getData() {
    return m_Data;
  }
  
  /**
   * returns the indices of the searched-for attributes (if TRUE) or the indices
   * of AttributeLocator objects (if FALSE)
   * 
   * @param findAtts      if true the indices of attributes are located,
   *                      otherwise the ones of AttributeLocator objects
   * @return              the indices of the attributes or the AttributeLocator objects
   */
  protected int[] find(boolean findAtts) {
    int		i;
    int[]	result;
    Vector<Integer>	indices;

    // determine locations
    indices = new Vector<Integer>();
    if (findAtts) {
      for (i = 0; i < m_Attributes.size(); i++) {
	if (((Boolean) m_Attributes.get(i)).booleanValue())
	  indices.add(new Integer(i));
      }
    }
    else {
      for (i = 0; i < m_Locators.size(); i++) {
	if (m_Locators.get(i) != null)
	  indices.add(new Integer(i));
      }
    }
    
    // fill array
    result = new int[indices.size()];
    for (i = 0; i < indices.size(); i++)
      result[i] = ((Integer) indices.get(i)).intValue();
    
    return result;
  }

  /**
   * returns actual index in the Instances object.
   * 
   * @param index	the index in the m_AllowedIndices array
   * @return		the actual index in the instances object
   */
  public int getActualIndex(int index) {
    return m_AllowedIndices[index];
  }
  
  /**
   * Returns the indices of the attributes. These indices are referring
   * to the m_AllowedIndices array, not the actual indices in the Instances
   * object.
   * 
   * @return	the indices of the attributes
   * @see	#getActualIndex(int)
   */
  public int[] getAttributeIndices() {
    return m_Indices;
  }
  
  /**
   * Returns the indices of the AttributeLocator objects.  These indices are 
   * referring to the m_AllowedIndices array, not the actual indices in the 
   * Instances object.
   * 
   * @return	the indices of the AttributeLocator objects
   * @see	#getActualIndex(int)
   */
  public int[] getLocatorIndices() {
    return m_LocatorIndices;
  }
  
  /**
   * Returns the AttributeLocator at the given index. This index refers to
   * the index of the m_AllowedIndices array, not the actual Instances object.
   * 
   * @param index   the index of the locator to retrieve
   * @return        the AttributeLocator at the given index
   */
  public AttributeLocator getLocator(int index) {
    return (AttributeLocator) m_Locators.get(index);
  }
  
  /**
   * Compares this object with the specified object for order. Returns a 
   * negative integer, zero, or a positive integer as this object is less 
   * than, equal to, or greater than the specified object. Only type and
   * indices are checked.
   * 
   * @param o		the object to compare with
   * @return		-1 if less than, 0 if equal, +1 if greater than the 
   * 			given object
   */
  public int compareTo(AttributeLocator o) {
    int		result;
    int		i;
    
    result = 0;
    
    // 1. check type
    if (this.getType() < o.getType()) {
      result = -1;
    }
    else if (this.getType() > o.getType()) {
      result = 1;
    }
    else {
      // 2. check indices
      if (this.getAllowedIndices().length < o.getAllowedIndices().length) {
	result = -1;
      }
      else if (this.getAllowedIndices().length > o.getAllowedIndices().length) {
	result = 1;
      }
      else {
	for (i = 0; i < this.getAllowedIndices().length; i++) {
	  if (this.getAllowedIndices()[i] < o.getAllowedIndices()[i]) {
	    result = -1;
	    break;
	  }
	  else if (this.getAllowedIndices()[i] > o.getAllowedIndices()[i]) {
	    result = 1;
	    break;
	  }
	  else {
	    result = 0;
	  }
	}
      }
    }
    
    return result;
  }
  
  /**
   * Indicates whether some other object is "equal to" this one. Only type
   * and indices are checked.
   * 
   * @param o		the AttributeLocator to check for equality
   * @return		true if the AttributeLocators have the same type and 
   * 			indices
   */
  public boolean equals(Object o) {
    return (compareTo((AttributeLocator) o) == 0);
  }
  
  /**
   * returns a string representation of this object
   * 
   * @return 		a string representation
   */
  public String toString() {
    return m_Attributes.toString();
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

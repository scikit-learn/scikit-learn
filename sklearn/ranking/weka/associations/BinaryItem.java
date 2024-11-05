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
 *    BinaryItem.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import java.io.Serializable;

import weka.core.Attribute;

  
  /**
   * Class that encapsulates an item whose backing Attribute is
   * binary or unary.
   * 
   * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
   * @version $Revision; $
   */
  public class BinaryItem extends NominalItem 
    implements Serializable {
    
    /** For serialization */
    private static final long serialVersionUID = -3372941834914147669L;    
    
    /**
     * Constructor.
     * 
     * @param att the attribute that backs this item.
     * @param valueIndex the index of the value for this item.
     * @throws Exception if the backing attribute is not binary or unary.
     */
    public BinaryItem(Attribute att, int valueIndex) throws Exception {
      super(att, valueIndex);
      
      if (att.isNumeric() || ((att.isNominal() || att.isRanking()) && att.numValues() > 2)) {
        throw new Exception("BinaryItem must be constructed using a nominal attribute" +
        		" with at most 2 values!");
      }
    }
       
    /**
     * Equals. Just compares attribute and valueIndex.
     * @return true if this BinaryItem is equal to the argument.
     */
    public boolean equals(Object compareTo) {
      if (!(compareTo instanceof BinaryItem)) {
        return false;
      }
      
      BinaryItem b = (BinaryItem)compareTo;
      if (m_attribute.equals(b.getAttribute()) && 
//          m_frequency == b.getFrequency() &&
          m_valueIndex == b.getValueIndex()) {
        return true;
      }
      
      return false;
    }
    
    public int hashCode() {
      return (m_attribute.name().hashCode() ^ 
          m_attribute.numValues()) * m_frequency;
    }
  }
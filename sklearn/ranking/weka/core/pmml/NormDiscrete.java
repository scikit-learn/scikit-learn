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
 *    NormDiscrete.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.util.ArrayList;

import org.w3c.dom.Element;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Utils;

/**
 * Class encapsulating a NormDiscrete Expression. Creates an
 * indicator for a particular discrete value.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public class NormDiscrete extends Expression {
  
  /**
   * For serialization 
   */
  private static final long serialVersionUID = -8854409417983908220L;

  /** The name of the field to lookup our value in */
  protected String m_fieldName;
  
  /** The actual attribute itself */
  protected Attribute m_field;
  
  /** The index of the attribute */
  protected int m_fieldIndex = -1;
  
  /** The actual value (as a String) that will correspond to an output of 1 */
  protected String m_fieldValue;
  
  /** True if a replacement for missing values has been specified */
  protected boolean m_mapMissingDefined = false;
  
  /** The value of the missing value replacement (if defined) */
  protected double m_mapMissingTo;
  
  /**
   *  If we are referring to a nominal (rather than String) attribute
   * then this holds the index of the value in question. Will be faster
   * than searching for the value each time.
   */
  protected int m_fieldValueIndex = -1;
  
  /**
   * Constructor. Reads the field name and field value for this NormDiscrete
   * Expression.
   * 
   * @param normDisc the Element encapsulating this NormDiscrete
   * @param opType the optype for this expression (taken from either the
   * enclosing DefineFunction or DerivedField)
   * @param fieldDefs an ArrayList of Attributes for the fields that this
   * Expression might need to access
   * enclosing DefineFunction or DerivedField)
   * @throws Exception if there is a problem parsing this Apply Expression
   */
  public NormDiscrete(Element normDisc, FieldMetaInfo.Optype opType, ArrayList<Attribute> fieldDefs)
    throws Exception {
    super(opType, fieldDefs);
    
    if (opType != FieldMetaInfo.Optype.CONTINUOUS) {
      throw new Exception("[NormDiscrete] can only have a continuous optype");
    }
    
    m_fieldName = normDisc.getAttribute("field");
    m_fieldValue = normDisc.getAttribute("value");
    
    String mapMissing = normDisc.getAttribute("mapMissingTo");
    if (mapMissing != null && mapMissing.length() > 0) {
      m_mapMissingTo = Double.parseDouble(mapMissing);
      m_mapMissingDefined = true;
    }
    
    if (fieldDefs != null) {
      setUpField();
    }
  }
  
  /**
   * Set the field definitions for this Expression to use
   * 
   * @param fieldDefs the field definitions to use
   * @throws Exception if there is a problem setting the field definitions
   */
  public void setFieldDefs(ArrayList<Attribute> fieldDefs) throws Exception {
    super.setFieldDefs(fieldDefs);
    setUpField();
  }
  
  /**
   * Find the named field, set up the index(es) etc.
   * 
   * @throws Exception if a problem occurs.
   */
  private void setUpField() throws Exception {
    m_fieldIndex = -1;
    m_fieldValueIndex = -1;
    m_field = null;
    
    if (m_fieldDefs != null) {
      m_fieldIndex = getFieldDefIndex(m_fieldName);

      if (m_fieldIndex < 0) {
        throw new Exception("[NormDiscrete] Can't find field " + m_fieldName
            + " in the supplied field definitions.");
      }
      m_field = m_fieldDefs.get(m_fieldIndex);
      
      if (!(m_field.isString() || m_field.isNominal() || m_field.isRanking())) {
        throw new Exception("[NormDiscrete] reference field " + m_fieldName
            +" must be categorical");
      }
      
      if (m_field.isNominal() || m_field.isRanking()) {
        // set up the value index
        m_fieldValueIndex = m_field.indexOfValue(m_fieldValue);
        if (m_fieldValueIndex < 0) {
          throw new Exception("[NormDiscrete] Unable to find value " + m_fieldValue
              + " in nominal attribute " + m_field.name());
        }
      } else if (m_field.isString()) {
        // add our value to this attribute (if it is already there
        // then this will have no effect).
        m_fieldValueIndex = m_field.addStringValue(m_fieldValue);
      }
    }
  }

  /**
   * Return the structure of the result of applying this Expression
   * as an Attribute.
   * 
   * @return the structure of the result of applying this Expression as an
   * Attribute.
   */
  protected Attribute getOutputDef() {    
    return new Attribute(m_fieldName + "=" + m_fieldValue);
  }

  /**
   * Get the result of evaluating the expression. In the case
   * of a continuous optype, a real number is returned; in
   * the case of a categorical/ordinal optype, the index of the nominal
   * value is returned as a double.
   * 
   * @param incoming the incoming parameter values
   * @return the result of evaluating the expression
   * @throws Exception if there is a problem computing the result
   */
  public double getResult(double[] incoming) throws Exception {
    
    double result = 0.0;
    if (Utils.isMissingValue(incoming[m_fieldIndex])) {
      if (m_mapMissingDefined) {
        result = m_mapMissingTo; // return the replacement
      } else {
        result = incoming[m_fieldIndex]; // just return the missing value
      }
    } else {
      if (m_fieldValueIndex == (int)incoming[m_fieldIndex]) {
        result = 1.0;
      }
    }
    
    return result;
  }

  /**
   * Always throws an Exception since the result of NormDiscrete must
   * be continuous.
   * 
   * @param incoming the incoming parameter values
   * @throws Exception always
   */
  public String getResultCategorical(double[] incoming) throws Exception {
    throw new Exception("[NormDiscrete] Can't return the result as a categorical value!");
  }
  
  public String toString(String pad) {
    StringBuffer buff = new StringBuffer();
    buff.append("NormDiscrete: " + m_fieldName + "=" + m_fieldValue);
    if (m_mapMissingDefined) {
      buff.append("\n" + pad + "map missing values to: " + m_mapMissingTo);
    }
    
    return buff.toString();
  }
}

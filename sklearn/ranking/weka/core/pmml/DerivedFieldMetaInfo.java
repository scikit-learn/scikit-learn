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
 *    DerivedFieldMetaInfo.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;
import weka.core.Instances;

public class DerivedFieldMetaInfo extends FieldMetaInfo implements Serializable {
  
  /** for serialization */
  
  /** display name */
  protected String m_displayName = null;
  
  /** 
   * the list of values (if the field is ordinal) - may be of size zero if none are specified.
   * If none are specified, we may be able to construct this by querying the Expression in
   * this derived field 
   */
  protected ArrayList<Value> m_values = new ArrayList<Value>();
  
  /** the single expression that defines the value of this field */
  protected Expression m_expression;
  
  public DerivedFieldMetaInfo(Element derivedField, ArrayList<Attribute> fieldDefs,
                              TransformationDictionary transDict) throws Exception {
    super(derivedField);
    // m_fieldName = derivedField.getAttribute("name");
    String displayName = derivedField.getAttribute("displayName");
    if (displayName != null && displayName.length() > 0) {
      m_displayName = displayName;
    }
    
    // get any values
    NodeList valL = derivedField.getElementsByTagName("Value");
    if (valL.getLength() > 0) {
      for (int i = 0; i < valL.getLength(); i++) {
        Node valueN = valL.item(i);
        if (valueN.getNodeType() == Node.ELEMENT_NODE) {
          Value v = new Value((Element)valueN);
          m_values.add(v);
        }
      }
    }
    
    // now get the expression
    m_expression = Expression.getExpression(derivedField, m_optype, fieldDefs, transDict);
  }
  
  /**
   * Upadate the field definitions for this derived field
   * 
   * @param fieldDefs the fields as an ArrayList of Attributes
   * @throws Exception if there is a problem setting the field definitions
   */
  public void setFieldDefs(ArrayList<Attribute> fieldDefs) throws Exception {
    m_expression.setFieldDefs(fieldDefs);
  }
  
  /**
   * Upadate the field definitions for this derived field
   * 
   * @param fields the fields as an Instances object
   * @throws Exception if there is a problem setting the field definitions
   */
  public void setFieldDefs(Instances fields) throws Exception {
    ArrayList<Attribute> tempDefs = new ArrayList<Attribute>();
    for (int i = 0; i < fields.numAttributes(); i++) {
      tempDefs.add(fields.attribute(i));
    }
    setFieldDefs(tempDefs);
  }
  
  /**
   * Get this derived field as an Attribute.
   * 
   * @return an Attribute for this derived field.
   */
  public Attribute getFieldAsAttribute() {
    return m_expression.getOutputDef().copy(m_fieldName);
  }
  
  /**
   * Get the derived field value for the given incoming vector of
   * values. Incoming values are assumed to be in the same order
   * as the attributes supplied in the field definitions ArrayList
   * used to construct this DerivedField.
   * 
   * If the optype of this derived field is continuous, then a real
   * number is returned. Otherwise, the number returned is the index
   * of the categorical/ordinal value corresponding to result of computing
   * the derived field value.
   * 
   * @param incoming the incoming parameter values
   * @return the result of computing the derived value
   * @throws Exception if there is a problem computing the value
   */
  public double getDerivedValue(double[] incoming) throws Exception {
    return m_expression.getResult(incoming);
  }
  
  public String toString() {
    StringBuffer buff = new StringBuffer();
    buff.append(getFieldAsAttribute() + "\nexpression:\n");
    buff.append(m_expression + "\n");
    
    return buff.toString();
  }
}

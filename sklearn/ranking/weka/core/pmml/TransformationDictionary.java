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
 *    TransformationDictionary.java
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
import weka.core.SerializedObject;

/**
 * Class encapsulating the TransformationDictionary element. Contains
 * a list of DefineFunctions and DerivedFields (if any).
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com
 * @version $Revision 1.0 $
 */
class TransformationDictionary implements Serializable {
  
  /** The defined functions (if any) */
  protected ArrayList<DefineFunction> m_defineFunctions = new ArrayList<DefineFunction>();
  
  /** The derived fields (if any) */
  protected ArrayList<DerivedFieldMetaInfo> m_derivedFields = 
    new ArrayList<DerivedFieldMetaInfo>();
  
  /**
   * Construct a new TransformationDictionary
   * 
   * @param dictionary the Element containing the dictionary
   * @param dataDictionary the data dictionary as an Instances object
   * @throws Exception if there is a problem constructing the transformation
   * dictionary
   */
  protected TransformationDictionary(Element dictionary, 
                                  Instances dataDictionary) throws Exception {
    
    // set up incoming field definitions
/*    ArrayList<Attribute> incomingFieldDefs = new ArrayList<Attribute>();
    for (int i = 0; i < dataDictionary.numAttributes(); i++) {
      incomingFieldDefs.add(dataDictionary.attribute(i));
    } */
    
    // get any derived fields and DefineFunctions
    NodeList derivedL = dictionary.getChildNodes();
    for (int i = 0; i < derivedL.getLength(); i++) {
      Node child = derivedL.item(i);
      if (child.getNodeType() == Node.ELEMENT_NODE) {
        String tagName = ((Element)child).getTagName();
        if (tagName.equals("DerivedField")) {
          DerivedFieldMetaInfo df = new DerivedFieldMetaInfo((Element)child, null /*incomingFieldDefs*/, null);
          m_derivedFields.add(df);
        } else if (tagName.equals("DefineFunction")) {
          DefineFunction defF = new DefineFunction((Element)child, null);
          m_defineFunctions.add(defF);
        }
      }
    }
  }
  
  /**
   * Set the field definitions for the derived fields. Usually called once the
   * structure of the mining schema + derived fields has been determined. 
   * Calling this method with an array list of field definitions in the order of 
   * attributes in the mining schema + derived fields will allow the expressions 
   * used in the derived fields to access the correct attribute values from the 
   * incoming instance (also allows for derived fields to reference other
   * derived fields). This is necessary because construction of the TransformationDictionary uses the
   * data dictionary to reference fields (the order of fields in the data dictionary
   * is not guaranteed to be the same as the order in the mining schema).
   * 
   * IMPORTANT: for derived field x to be able to reference derived field y, y must
   * have been declared before x in the PMML file. This is because the process
   * of constructing an input vector of values to the model proceeds in a linear
   * left-to-right fashion - so any referenced derived field (e.g. field y),
   * must have already computed its value when x is evaluated.
   * 
   * @param fieldDefs the definition of the incoming fields as an array list of attributes
   * @throws Exception if a problem occurs
   */
  protected void setFieldDefsForDerivedFields(ArrayList<Attribute> fieldDefs) throws Exception {
    for (int i = 0; i < m_derivedFields.size(); i++) {
      m_derivedFields.get(i).setFieldDefs(fieldDefs);
    }
    
    // refresh the define functions - force them to pass on their parameter
    // definitions as field defs to their encapsulated expression. Parameter
    // defs were not passed on by expressions encapsulated in DefineFunctions
    // at construction time because the encapsulated expression does not know
    // whether it is contained in a DefineFunction or a DerivedField. Since
    // we delay passing on field definitions until all derived fields are 
    // loaded (in order to allow derived fields to reference other derived fields),
    // we must tell DefineFunctions to pass on their parameter definitions
    for (int i = 0; i < m_defineFunctions.size(); i++) {
      m_defineFunctions.get(i).pushParameterDefs();
    }
  }
  
  /**
   * Set the field definitions for the derived fields. Usually called once the
   * structure of the mining schema has been determined. Calling this method
   * with an array list of field definitions in the order of attributes in the
   * mining schema will allow the expressions used in the derived fields to
   * access the correct attribute values from the incoming instances. This is
   * necessary because construction of the TransformationDictionary uses the
   * data dictionary to reference fields (the order of fields in the data dictionary
   * is not guaranteed to be the same as the order in the mining schema).
   * 
   * @param fieldDefs the definition of the incoming fields as an Instances object
   * @throws Exception if a problem occurs
   */
  protected void setFieldDefsForDerivedFields(Instances fieldDefs) throws Exception {
    ArrayList<Attribute> tempDefs = new ArrayList<Attribute>();
    for (int i = 0; i < fieldDefs.numAttributes(); i++) {
      tempDefs.add(fieldDefs.attribute(i));
    }
    setFieldDefsForDerivedFields(tempDefs);
  }
  
  protected ArrayList<DerivedFieldMetaInfo> getDerivedFields() {
    return new ArrayList<DerivedFieldMetaInfo>(m_derivedFields);
  }
  
  /**
   * Get a named DefineFunction. Returns a deep copy of the
   * function.
   * 
   * @param functionName the name of the function to get
   * @return the named function or null if it cannot be found
   * @throws Exception if there is a problem deep copying the function
   */
  protected DefineFunction getFunction(String functionName) throws Exception {

    DefineFunction copy = null;
    DefineFunction match = null;
    for (DefineFunction f : m_defineFunctions) {
      if (f.getName().equals(functionName)) {
        match = f;
        //System.err.println("Found a match!!!");
        break;
      }
    }
    
    if (match != null) {
      SerializedObject so = new SerializedObject(match, false);
      copy = (DefineFunction)so.getObject();
      //System.err.println(copy);
    }
    
    return copy;
  }
  
  public String toString() {
    StringBuffer buff = new StringBuffer();
    
    buff.append("Transformation dictionary:\n");
    
    if (m_derivedFields.size() > 0) {
      buff.append("derived fields:\n");
      for (DerivedFieldMetaInfo d : m_derivedFields) {
        buff.append("" + d.getFieldAsAttribute() + "\n");
      }
    }
    
    if (m_defineFunctions.size() > 0) {
      buff.append("\nfunctions:\n");
      for (DefineFunction f : m_defineFunctions) {
        buff.append(f.toString("  ") + "\n");
      }
    }
    
    buff.append("\n");
    
    return buff.toString();
  }
}

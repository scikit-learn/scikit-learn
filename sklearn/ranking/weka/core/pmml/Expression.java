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
 *    Expression.java
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

public abstract class Expression  implements Serializable {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = 4448840549804800321L;
  
  /** The optype of this Expression */
  protected FieldMetaInfo.Optype m_opType;

  /** The field defs */
  protected ArrayList<Attribute> m_fieldDefs = null;
  
  // NOTE - might need to pass in mining schema in order
  // to determine values for nominal optypes
  public Expression(FieldMetaInfo.Optype opType, ArrayList<Attribute> fieldDefs) {
    m_opType = opType;
    m_fieldDefs = fieldDefs;
  }
  
  /**
   * Set the field definitions for this Expression to use
   * 
   * @param fieldDefs the field definitions to use
   * @throws Exception if there is a problem setting the field definitions
   */
  public void setFieldDefs(ArrayList<Attribute> fieldDefs) throws Exception {
    m_fieldDefs = fieldDefs;
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
  public abstract double getResult(double[] incoming) throws Exception;
  
  /**
   * Get the result of evaluating the expression for continuous
   * optype. Is the same as calling getResult() when the optype
   * is continuous.
   * 
   * @param incoming the incoming parameter values
   * mining schema
   * @return the result of evaluating the expression.
   * @throws Exception if the optype is not continuous.
   */
  public double getResultContinuous(double[] incoming) throws Exception {
    if (!(m_opType == FieldMetaInfo.Optype.CONTINUOUS)) {
      throw new Exception("[Expression] Can't return continuous result "
          + "as optype is not continuous");
    }
    return getResult(incoming);
  }
  
  /**
   * Gets the result of evaluating the expression when the
   * optype is categorical or ordinal as the actual String
   * value.
   * 
   * @param incoming the incoming parameter values 
   * @return the result of evaluating the expression
   * @throws Exception if the optype is continuous
   */
  public abstract String getResultCategorical(double[] incoming) 
    throws Exception;
  
  /**
   * Return the structure of the result of applying this Expression
   * as an Attribute.
   * 
   * @return the structure of the result of applying this Expression as an
   * Attribute.
   */
  protected abstract Attribute getOutputDef(); 

  /**
   * Static factory method that returns a subclass of Expression that
   * encapsulates the type of expression contained in the Element
   * supplied. Assumes that there is just one expression contained
   * in the supplied Element.
   * 
   * @param container the Node containing the expression
   * @param opType the optype of the value returned by this Expression.
   * @param fieldDefs an ArrayList of Attributes for the fields that this
   * Expression may need to access
   * Since Expressions are children of either DerivedFields or
   * DefineFuntions, they will have the same optype as their parent.
   * @param transDict the TransformationDictionary (may be null if there
   * is no dictionary)
   * @return an Expression object or null if there is no known expression in
   * the container 
   * @throws Exception for unsupported Expression types 
   */
  public static Expression getExpression(Node container, 
      FieldMetaInfo.Optype opType,
      ArrayList<Attribute> fieldDefs,
      TransformationDictionary transDict) throws Exception {
    
    // we need to examine children of this Node to find an expression,
    // not the entire subtree (as would be returned by Element.getElementsByTagName()
    
    Expression result = null;
    String tagName = "";
    
    NodeList children = container.getChildNodes();
    if (children.getLength() == 0) {
      throw new Exception("[Expression] container has no children!");
    }
    
    // at this level in the tree there should be only one expression type
    // specified - look for it here.
    for (int i = 0; i < children.getLength(); i++) {
      Node child = children.item(i);
      if (child.getNodeType() == Node.ELEMENT_NODE) {
        tagName = ((Element)child).getTagName();
        result = getExpression(tagName, child, opType, fieldDefs, transDict);
        if (result != null) {
          break;
        }
      }
    }
    
    return result;
  }
  
  /**
   * Static factory method that returns a subclass of Expression that
   * encapsulates the type of expression supplied as an argument.
   * 
   * @param name the name of the Expression to get
   * @param expression the Node containing the expression
   * @param opType the optype of the value returned by this Expression.
   * @param fieldDefs an ArrayList of Attributes for the fields that this
   * Expression may need to access
   * Since Expressions are children of either DerivedFields or
   * DefineFuntions, they will have the same optype as their parent.
   * @param transDict the TransformationDictionary (may be null if there
   * is no dictionary)
   * @return an Expression object or null if there is no known expression in
   * the container 
   * @throws Exception for unsupported Expression types 
   */
  public static Expression getExpression(String name, 
      Node expression,
      FieldMetaInfo.Optype opType,
      ArrayList<Attribute> fieldDefs,
      TransformationDictionary transDict) throws Exception {
   
    Expression result = null;
    
    if (name.equals("Constant")) {
      // construct a Constant expression
      result = new Constant((Element)expression, opType, fieldDefs);
    } else if (name.equals("FieldRef")) {
      // construct a FieldRef expression
      result = new FieldRef((Element)expression, opType, fieldDefs);
    } else if (name.equals("Apply")) {
      // construct an Apply expression
      result = new Apply((Element)expression, opType, fieldDefs, transDict);
    } else if (name.equals("NormDiscrete")) {
      result = new NormDiscrete((Element)expression, opType, fieldDefs);
    } else if (name.equals("NormContinuous")) {
      result = new NormContinuous((Element)expression, opType, fieldDefs);
    } else if (name.equals("Discretize")) {
      result = new Discretize((Element)expression, opType, fieldDefs);
    } else if (name.equals("MapValues") ||
        name.equals("Aggregate")) {
      throw new Exception("[Expression] Unhandled Expression type " + name);
    }
    return result;
  }
  
  /**
   * Return the named attribute from the list of reference fields.
   * 
   * @param attName the name of the attribute to retrieve
   * @return the named attribute (or null if it can't be found).
   */
  public Attribute getFieldDef(String attName) {
    Attribute returnV = null;
    for (int i = 0; i < m_fieldDefs.size(); i++) {
      if (m_fieldDefs.get(i).name().equals(attName)) {
        returnV = m_fieldDefs.get(i);
        break;
      }
    }
    return returnV;
  }
  
  public int getFieldDefIndex(String attName) {
    int returnV = -1;
    for (int i = 0; i < m_fieldDefs.size(); i++) {
      if (m_fieldDefs.get(i).name().equals(attName)) {
        returnV = i;
        break;
      }
    }
    return returnV;
  }
  
  /**
   * Get the optype of the result of applying this Expression.
   * 
   * @return the optype of the result of applying this Expression
   */
  public FieldMetaInfo.Optype getOptype() {
    return m_opType;
  }
  
  public String toString() {
    return toString("");
  }
  
  public String toString(String pad) {
    return pad + this.getClass().getName();
  }
}

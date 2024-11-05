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
 *    DefineFunction.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.util.ArrayList;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import weka.core.Attribute;

/**
 * Class encapsulating DefineFunction (used in TransformationDictionary).
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com
 * @version $Revision 1.0 $
 */
public class DefineFunction extends Function {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = -1976646917527243888L;

  /**
   * Inner class for handling Parameters
   */
  protected class ParameterField extends FieldMetaInfo {
    /**
     * For serialization
     */
    private static final long serialVersionUID = 3918895902507585558L;

    protected ParameterField(Element field) {
      super(field);
    }
    
    public Attribute getFieldAsAttribute() {
      if (m_optype == Optype.CONTINUOUS) {
        return new Attribute(m_fieldName);
      }
      // return a string attribute for categorical/ordinal optypes
      return new Attribute(m_fieldName, (ArrayList<String>)null);
    }
  }
  
  /** 
   * The list of parameters expected by this function. We can use this to do 
   * some error/type checking when users call setParameterDefs() on us 
   */
  protected ArrayList<ParameterField> m_parameters = new ArrayList<ParameterField>();
  
  /** The optype for this function */
  FieldMetaInfo.Optype m_optype = FieldMetaInfo.Optype.NONE;
  
  /** The Expression for this function to use */
  protected Expression m_expression = null;
  
  public DefineFunction(Element container, TransformationDictionary transDict) throws Exception {
    
    m_functionName = container.getAttribute("name");
    
    // get the optype for this function
    String opType = container.getAttribute("optype");
    if (opType != null && opType.length() > 0) {
      for (FieldMetaInfo.Optype o : FieldMetaInfo.Optype.values()) {
        if (o.toString().equals(opType)) {
          m_optype = o;
          break;
        }
      }
    } else {
      throw new Exception("[DefineFunction] no optype specified!!");
    }
    
    m_parameterDefs = new ArrayList<Attribute>();
    
    // get all the parameters
    NodeList paramL = container.getElementsByTagName("ParameterField");
    for (int i = 0; i < paramL.getLength(); i++) {
      Node paramN = paramL.item(i);
      if (paramN.getNodeType() == Node.ELEMENT_NODE) {
        ParameterField newP = new ParameterField((Element)paramN);
        m_parameters.add(newP);
        
        // set up default parameter definitions - these will probably get replaced
        // by more informative ones (i.e. possibly nominal attributes instead of
        // string attributes) later
        m_parameterDefs.add(newP.getFieldAsAttribute());
      }
    }
    
    m_expression = Expression.getExpression(container, m_optype, m_parameterDefs, transDict);
    
    // check that the optype of the Expression is compatible with ours
    if (m_optype == FieldMetaInfo.Optype.CONTINUOUS && 
        m_expression.getOptype() != m_optype) {
      throw new Exception("[DefineFunction] optype is continuous but our Expression's optype "
          + "is not.");
    }
    
    if ((m_optype == FieldMetaInfo.Optype.CATEGORICAL || m_optype == FieldMetaInfo.Optype.ORDINAL) !=
        (m_expression.getOptype() == FieldMetaInfo.Optype.CATEGORICAL || 
         m_expression.getOptype() == FieldMetaInfo.Optype.ORDINAL)) {
      throw new Exception("[DefineFunction] optype is categorical/ordinal but our Expression's optype "
          + "is not."); 
    }
  }
  
  public void pushParameterDefs() throws Exception {
    if (m_parameterDefs == null) {
      throw new Exception("[DefineFunction] parameter definitions are null! Can't "
          + "push them to encapsulated expression.");
    }
    
    m_expression.setFieldDefs(m_parameterDefs);
  }

  /**
   * Get the structure of the result produced by this function.
   * 
   * @return the structure of the result produced by this function.
   */
  public Attribute getOutputDef() {
    return m_expression.getOutputDef();
  }

  /**
   * Returns an array of the names of the parameters expected
   * as input by this function. May return null if this function
   * can take an unbounded number of parameters (i.e. min, max, etc.).
   * 
   * @return an array of the parameter names or null if there are an
   * unbounded number of parameters.
   */
  public String[] getParameterNames() {
    String[] result = new String[m_parameters.size()];
    for (int i = 0; i < m_parameters.size(); i++) {
      result[i] = m_parameters.get(i).getFieldName();
    }
    
    return result;
  }

  /**
   * Get the result of applying this function.
   * 
   * @param incoming the arguments to this function (supplied in order to match that
   * of the parameter definitions
   * @return the result of applying this function. When the optype is
   * categorical or ordinal, an index into the values of the output definition
   * is returned.
   * @throws Exception if there is a problem computing the result of this function
   */
  public double getResult(double[] incoming) throws Exception {
    
    if (incoming.length != m_parameters.size()) {
      throw new IllegalArgumentException("[DefineFunction] wrong number of arguments: expected "
          + m_parameters.size() + ", recieved " + incoming.length);
    }

    return m_expression.getResult(incoming);
  }

  /**
   * Set the structure of the parameters that are expected as input by
   * this function. This must be called before getOutputDef() is called.
   * 
   * @param paramDefs the structure of the input parameters
   * @throws Exception if the number or types of parameters are not acceptable by
   * this function
   */
  public void setParameterDefs(ArrayList<Attribute> paramDefs) throws Exception {
    if (paramDefs.size() != m_parameters.size()) {
      throw new Exception("[DefineFunction] number of parameter definitions does not match "
          + "number of parameters!");
    }
    
    // check these defs against the optypes of the parameters
    for (int i = 0; i < m_parameters.size(); i++) {
      if (m_parameters.get(i).getOptype() == FieldMetaInfo.Optype.CONTINUOUS) {
        if (!paramDefs.get(i).isNumeric()) {
          throw new Exception("[DefineFunction] parameter "
              + m_parameters.get(i).getFieldName() + " is continuous, but corresponding "
              + "supplied parameter def " + paramDefs.get(i).name() + " is not!");
        }
      } else {
        if (!paramDefs.get(i).isNominal() && !paramDefs.get(i).isRanking() && !paramDefs.get(i).isString()) {
          throw new Exception("[DefineFunction] parameter "
              + m_parameters.get(i).getFieldName() + " is categorical/ordinal, but corresponding "
              + "supplied parameter def " + paramDefs.get(i).name() + " is not!");
        }
      }
    }
    
    // now we need to rename these argument definitions to match the names of
    // the actual parameters
    ArrayList<Attribute> newParamDefs = new ArrayList<Attribute>();
    for (int i = 0; i < paramDefs.size(); i++) {
      Attribute a = paramDefs.get(i);
      newParamDefs.add(a.copy(m_parameters.get(i).getFieldName()));
    }
    
    m_parameterDefs = newParamDefs;
    
    // update the Expression
    m_expression.setFieldDefs(m_parameterDefs);
  }
  
  public String toString() {
    return toString("");
  }
  
  public String toString(String pad) {
    StringBuffer buff = new StringBuffer();
    
    buff.append(pad + "DefineFunction (" + m_functionName + "):\n"
        + pad + "nparameters:\n");
    
    for (ParameterField p : m_parameters) {
      buff.append(pad + p.getFieldAsAttribute() + "\n");
    }
    
    buff.append(pad + "expression:\n" + m_expression.toString(pad + "  "));
    return buff.toString();
  }
}

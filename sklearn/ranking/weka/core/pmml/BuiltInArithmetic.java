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
 *    Arithmetic.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.util.ArrayList;

import weka.core.Attribute;

/**
 * Built-in function for +, -, *, /.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public class BuiltInArithmetic extends Function {
  
  /**
   * For serialization.
   */
  private static final long serialVersionUID = 2275009453597279459L;

  /**
   * Enumerated type for the operator.
   */
  enum Operator {
    ADDITION (" + ") {
      double eval(double a, double b) {
        return a + b;
      }
    },
    SUBTRACTION (" - ") {
      double eval(double a, double b) {
        return a - b;
      }
    },
    MULTIPLICATION (" * ") {
      double eval(double a, double b) {
        return a * b;
      }
    },
    DIVISION (" / ") {
      double eval(double a, double b) {
        return a / b;
      }
    };
    
    abstract double eval(double a, double b);
    
    private final String m_stringVal;
    
    Operator(String opName) {
      m_stringVal = opName;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  /** The operator for this function */
  protected Operator m_operator = Operator.ADDITION;
  
  /**
   * Construct a new Arithmetic built-in pmml function.
   * @param op the operator to use.
   */
  public BuiltInArithmetic(Operator op) {
    m_operator = op;
    m_functionName = m_operator.toString();
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
    m_parameterDefs = paramDefs;
    
    if (m_parameterDefs.size() != 2) {
      throw new Exception("[Arithmetic] wrong number of parameters. Recieved " 
          + m_parameterDefs.size() + ", expected 2.");
    }
  }
  
  /**
   * Returns an array of the names of the parameters expected
   * as input by this function
   * 
   * @return an array of the parameter names
   */
  public String[] getParameterNames() {
    String[] result = {"A", "B"};
    return result;
  }
  
  /**
   * Get the structure of the result produced by this function.
   * Subclasses must implement.
   * 
   * @return the structure of the result produced by this function.
   */
  public Attribute getOutputDef() {
    return new Attribute("BuiltInArithmeticResult:" + m_operator.toString());
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
    if (m_parameterDefs == null) {
      throw new Exception("[BuiltInArithmetic] incoming parameter structure has not been set!");
    }
    
    if (m_parameterDefs.size() != 2 || incoming.length != 2) {
      throw new Exception("[BuiltInArithmetic] wrong number of parameters!");
    }
    
    double result = m_operator.eval(incoming[0], incoming[1]);
    
    return result;
  }
  
  public String toString() {
    return toString("");
  }
  
  public String toString(String pad) {
    return pad + m_parameterDefs.get(0).name() + m_functionName
      + m_parameterDefs.get(1).name();
  }
}

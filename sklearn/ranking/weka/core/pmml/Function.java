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
 *    Function.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import java.io.Serializable;
import java.util.ArrayList;

import weka.core.Attribute;

/**
 * Abstract superclass for PMML built-in and DefineFunctions.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public abstract class Function implements Serializable {
  
  /**
   * For serialization 
   */
  private static final long serialVersionUID = -6997738288201933171L;
  
  /** The name of this function */
  protected String m_functionName;
  
  /** The structure of the parameters to this function */
  protected ArrayList<Attribute> m_parameterDefs = null;
  
    
  public String getName() {
    return m_functionName;
  }
  
  /**
   * Returns an array of the names of the parameters expected
   * as input by this function. May return null if this function
   * can take an unbounded number of parameters (i.e. min, max, etc.).
   * 
   * @return an array of the parameter names or null if there are an
   * unbounded number of parameters.
   */
  public abstract String[] getParameterNames();
  
  /**
   * Set the structure of the parameters that are expected as input by
   * this function. This must be called before getOutputDef() is called.
   * 
   * @param paramDefs the structure of the input parameters
   * @throws Exception if the number or types of parameters are not acceptable by
   * this function
   */
  public abstract void setParameterDefs(ArrayList<Attribute> paramDefs) throws Exception;
  
  /**
   * Get the structure of the result produced by this function.
   * 
   * @return the structure of the result produced by this function.
   */
  public abstract Attribute getOutputDef();
  
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
  public abstract double getResult(double[] incoming) throws Exception;
  
  /**
   * Get the result of applying this function. Subclasses should overide this
   * method when they might produce categorical values where the legal set of
   * values can't be determined apriori (i.e. by just using the input parameter
   * definitions). An example is the substring function - in this case there
   * is no way of knowing apriori what all the legal values will be because the
   * start position and length parameters are not known until the function is
   * invoked. In this scenario, a client could call getResultCategorical()
   * repeatedly (in an initialization routine) in order to manually build the
   * list of legal values and then call this method at processing time, passing
   * in the pre-computed output structure.
   * 
   *  This default implementation ignores the supplied output definition argument
   *  and simply invokes getResult(incoming).
   * 
   * @param incoming the arguments to this function (supplied in order to match that
   * of the parameter definitions
   * @param outputDef the output definition to use for looking up the index of
   * result values (in the case of categorical results)
   * @return the result of applying this function. When the optype is
   * categorical or ordinal, an index into the values of the output definition
   * is returned.
   * @throws Exception if there is a problem computing the result of this function
   *
  public double getResult(double[] incoming, Attribute outputDef) throws Exception {
    if (outputDef.isString()) {
      throw new Exception("[Function] outputDef argument must not be a String attribute!");
    }
    return getResult(incoming);
  }*/
  
  /**
   * Get the result of applying this function when the output type categorical.
   * Will throw an exception for numeric output. If subclasses output definition
   * is a string attribute (i.e. because all legal values can't be computed apriori),
   * then the subclass will need to overide this method and return something sensible
   * in this case.
   * 
   * @param incoming the incoming arguments to this function (supplied in order to match
   * that of the parameter definitions
   * @return the result of applying this function as a String.
   * @throws Exception if this method is not applicable because the optype is not
   * categorical/ordinal
   *
  public String getResultCategorical(double[] incoming) throws Exception {
    if (getOutputDef().isNumeric()) {
      throw new Exception("[Function] can't return nominal value, output is numeric!!");
    }
    
    if (getOutputDef().isString()) {
      throw new Exception("[Function] subclass neeeds to overide this method and do "
          + "something sensible when the output def is a string attribute.");
    }
    
    return getOutputDef().value((int)getResult(incoming));
  } */
  
  
  //public static FieldMetaInfo.Optype
  
  /**
   * Get a built-in PMML Function.
   * 
   * @param name the name of the function to get.
   * @return a built-in Function or null if the named function is not 
   * known/supported.
   */
  public static Function getFunction(String name) {
    Function result = null;
    
    name = name.trim();
    if (name.equals("+")) {
      result = new BuiltInArithmetic(BuiltInArithmetic.Operator.ADDITION);
    } else if (name.equals("-")) {
      result = new BuiltInArithmetic(BuiltInArithmetic.Operator.SUBTRACTION);
    } else if (name.equals("*")) {
      result = new BuiltInArithmetic(BuiltInArithmetic.Operator.MULTIPLICATION);
    } else if (name.equals("/")) {
      result = new BuiltInArithmetic(BuiltInArithmetic.Operator.DIVISION);
    } else if (name.equals(BuiltInMath.MathFunc.MIN.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.MIN);
    } else if (name.equals(BuiltInMath.MathFunc.MAX.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.MAX);
    } else if (name.equals(BuiltInMath.MathFunc.SUM.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.SUM);
    } else if (name.equals(BuiltInMath.MathFunc.AVG.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.AVG);
    } else if (name.equals(BuiltInMath.MathFunc.LOG10.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.LOG10);
    } else if (name.equals(BuiltInMath.MathFunc.LN.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.LN);
    } else if (name.equals(BuiltInMath.MathFunc.SQRT.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.SQRT);
    } else if (name.equals(BuiltInMath.MathFunc.ABS.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.ABS);
    } else if (name.equals(BuiltInMath.MathFunc.EXP.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.EXP);
    } else if (name.equals(BuiltInMath.MathFunc.POW.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.POW);
    } else if (name.equals(BuiltInMath.MathFunc.THRESHOLD.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.THRESHOLD);
    } else if (name.equals(BuiltInMath.MathFunc.FLOOR.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.FLOOR);
    } else if (name.equals(BuiltInMath.MathFunc.CEIL.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.CEIL);
    } else if (name.equals(BuiltInMath.MathFunc.ROUND.toString())) {
      result = new BuiltInMath(BuiltInMath.MathFunc.ROUND);
    }
    
    return result;
  }
  
  
  /**
   * Get either a function. Built-in functions are queried first, and then
   * DefineFunctions in the TransformationDictionary (if any).
   * 
   * @param name the name of the function to get.
   * @param transDict the TransformationDictionary (may be null if there is
   * no dictionary).
   * @return the function
   * @throws Exception if the named function is not known/supported.
   */
  public static Function getFunction(String name, TransformationDictionary transDict)
    throws Exception {
    
    Function result = getFunction(name);
    
    // try the defined functions in the TransformationDictionary (if any)
    if (result == null && transDict != null) {
      result = transDict.getFunction(name);
    }
    
    if (result == null) {
      throw new Exception("[Function] unknown/unsupported function " + name);
    }
    
    return result;
  }
  
  public String toString() {
    return toString("");
  }
  
  public String toString(String pad) {
    return pad + this.getClass().getName();
  }
}

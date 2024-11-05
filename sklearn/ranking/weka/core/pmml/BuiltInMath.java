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
import weka.core.Utils;

/**
 * Built-in function for min, max, sum, avg, log10,
 * ln, sqrt, abs, exp, pow, threshold, floor, ceil and round.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision 1.0 $
 */
public class BuiltInMath extends Function {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = -8092338695602573652L;

  /**
   * Enum for the math functions.
   */
  enum MathFunc {
    MIN ("min") {
      double eval(double[] args) {
        return args[Utils.minIndex(args)];
      }
      
      boolean legalNumParams(int num) {
        return (num > 0);
      }
      
      String[] getParameterNames() {
        return null; // unbounded number of parameters
      }
    },
    MAX ("max") {
      double eval(double[] args) {
        return args[Utils.maxIndex(args)];
      }
      
      boolean legalNumParams(int num) {
        return (num > 0);
      }
      
      String[] getParameterNames() {
        return null; // unbounded number of parameters
      }
    },
    SUM ("sum") {
      double eval(double[] args) {
        return Utils.sum(args);
      }
      
      boolean legalNumParams(int num) {
        return (num > 0);
      }
      
      String[] getParameterNames() {
        return null; // unbounded number of parameters
      }
    },
    AVG ("avg") {
      double eval(double[] args) {
        return Utils.mean(args);
      }
      
      boolean legalNumParams(int num) {
        return (num > 0);
      }
      
      String[] getParameterNames() {
        return null; // unbounded number of parameters
      }
    },
    LOG10 ("log10") {
      double eval(double[] args) {
        return Math.log10(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    LN ("ln") {
      double eval(double[] args) {
        return Math.log(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    SQRT ("sqrt") {
      double eval(double[] args) {
        return Math.sqrt(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    ABS ("abs") {
      double eval(double[] args) {
        return Math.abs(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    EXP ("exp") {
      double eval(double[] args) {
        return Math.exp(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    POW ("pow") {
      double eval(double[] args) {
        return Math.pow(args[0], args[1]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 2);
      }
      
      String[] getParameterNames() {
        return new String[] {"A", "B"};
      }
    },
    THRESHOLD ("threshold") {
      double eval(double[] args) {
        if (args[0] > args[1]) {
          return 1.0;
        } else {
          return 0.0;
        }
      }
      
      boolean legalNumParams(int num) {
        return (num == 2);
      }
      
      String[] getParameterNames() {
        return new String[] {"A", "B"};
      }
    },
    FLOOR ("floor") {
      double eval(double[] args) {
        return Math.floor(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    CEIL ("ceil") {
      double eval(double[] args) {
        return Math.ceil(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    },
    ROUND ("round") {
      double eval(double[] args) {
        return Math.round(args[0]);
      }
      
      boolean legalNumParams(int num) {
        return (num == 1);
      }
      
      String[] getParameterNames() {
        return new String[] {"A"};
      }
    };
    
    abstract double eval(double[] args);
    abstract boolean legalNumParams(int num);
    abstract String[] getParameterNames();
    
    private final String m_stringVal;
    
    MathFunc(String funcName) {
      m_stringVal = funcName;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  /** The function to apply */
  protected MathFunc m_func = MathFunc.ABS;
  
  /**
   * Construct a new built-in pmml Math function.
   * @param func the math function to use
   */
  public BuiltInMath(MathFunc func) {
    m_func = func;
    m_functionName = m_func.toString();
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
    
    if (!m_func.legalNumParams(m_parameterDefs.size())) {
      throw new Exception("[BuiltInMath] illegal number of parameters for function: " 
          + m_functionName);
    }
  }

  /**
   * Get the structure of the result produced by this function.
   * Subclasses must implement.
   * 
   * @return the structure of the result produced by this function.
   */
  public Attribute getOutputDef() {
    return new Attribute("BuiltInMathResult:" + m_func.toString());
  }

  /**
   * Returns an array of the names of the parameters expected
   * as input by this function. May return null if the function
   * can accept an unbounded number of arguments.
   * 
   * @return an array of the parameter names (or null if the function
   * can accept any number of arguments).
   */
  public String[] getParameterNames() {
    return m_func.getParameterNames();
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
      throw new Exception("[BuiltInMath] incoming parameter structure has not been set");
    }
    
    if (!m_func.legalNumParams(incoming.length)) {
      throw new Exception("[BuiltInMath] wrong number of parameters!");
    }
    
    double result = m_func.eval(incoming);
    
    return result;
  }
  
  public String toString() {
    String result = m_func.toString() + "(";
    for (int i = 0; i < m_parameterDefs.size(); i++) {
      result += m_parameterDefs.get(i).name();
      if (i != m_parameterDefs.size() - 1) {
        result += ", ";
      } else {
        result += ")";
      }
    }
    return result;
  }
}

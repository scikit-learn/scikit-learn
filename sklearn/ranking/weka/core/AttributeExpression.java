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
 *    AttributeExpression.java
 *    Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.io.Serializable;
import java.util.Stack;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * A general purpose class for parsing mathematical expressions
 * involving attribute values. Values can be provided in an array
 * or in an Instance. Values are accessed in the expression by
 * prefixing their index (starting at 1) with the character 'a'.
 *
 * <pre> Example expression: a1^2*a5/log(a7*4.0) </pre>
 *
 * Supported opperators: +, -, *, /, ^, log, abs, cos, exp, sqrt,
 * floor, ceil, rint, tan, sin, (, ).
 *
 * @author Mark Hall
 * @version $Revision: 5988 $
 */
public class AttributeExpression
  implements Serializable, RevisionHandler {

  /** for serialization */
  static final long serialVersionUID = 402130123261736245L;
  
  /**
   * Interface implemented by operators and operants.
   */
  private interface ExpressionComponent {};

  /**
   * Inner class handling an attribute index as an operand
   */
  private class AttributeOperand  
    implements ExpressionComponent, Serializable, RevisionHandler {
    
    /** for serialization */
    static final long serialVersionUID = -7674280127286031105L;

    /** the index of the attribute */
    protected int m_attributeIndex;

    /** true if the value of the attribute are to be multiplied by -1 */
    protected boolean m_negative;

    /**
     * Constructor
     * 
     * @param operand
     * @param sign
     * @throws Exception
     */
    public AttributeOperand(String operand, boolean sign) throws Exception {
      // strip the leading 'a' and set the index
      m_attributeIndex = (Integer.parseInt(operand.substring(1)))-1;
      m_negative = sign;
    }

    /**
     * Return a string describing this object
     * @return a string descibing the attribute operand
     */
    public String toString() {
      String result = "";
      if (m_negative) {
	result += '-';
      }
      return result+"a"+(m_attributeIndex+1);
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5988 $");
    }
  }

  /**
   * Inner class for storing numeric constant opperands
   */
  private class NumericOperand 
    implements ExpressionComponent, Serializable, RevisionHandler {
    
    /** for serialization */
    static final long serialVersionUID = 9037007836243662859L;

    /** numeric constant */
    protected double m_numericConst;

    /**
     * Constructor
     * 
     * @param operand
     * @param sign
     * @throws Exception
     */
    public NumericOperand(String operand, boolean sign) throws Exception {
      m_numericConst = Double.valueOf(operand).doubleValue();
      if (sign) {
	m_numericConst *= -1.0;
      }
    }
    
    /**
     * Return a string describing this object
     * @return a string descibing the numeric operand
     */
    public String toString() {
      return ""+m_numericConst;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5988 $");
    }
  }

  /**
   * Inner class for storing operators
   */
  private class Operator 
    implements ExpressionComponent, Serializable, RevisionHandler {
    
    /** for serialization */
    static final long serialVersionUID = -2760353522666004638L;

    /** the operator */
    protected char m_operator;

    /**
     * Constructor
     * 
     * @param opp the operator
     */
    public Operator(char opp) {
      if (!isOperator(opp)) {
	throw new IllegalArgumentException("Unrecognized operator:" + opp);
      }
      m_operator = opp;
    }

    /**
     * Apply this operator to the supplied arguments
     * @param first the first argument
     * @param second the second argument
     * @return the result
     */
    protected double applyOperator(double first, double second) {
      switch (m_operator) {
	case '+' :
	  return (first+second);
	case '-' :
	  return (first-second);
	case '*' :
	  return (first*second);
	case '/' :
	  return (first/second);
	case '^' :
	  return Math.pow(first,second);
      }
      return Double.NaN;
    }

    /**
     * Apply this operator (function) to the supplied argument
     * @param value the argument
     * @return the result
     */
    protected double applyFunction(double value) {
      switch (m_operator) {
	case 'l' :
	  return Math.log(value);
	case 'b' :
	  return Math.abs(value);
	case 'c' :
	  return Math.cos(value);
	case 'e' :
	  return Math.exp(value);
	case 's' :
	  return Math.sqrt(value);
	case 'f' :
	  return Math.floor(value);
	case 'h' :
	  return Math.ceil(value);
	case 'r' :
	  return Math.rint(value);
	case 't' :
	  return Math.tan(value);
	case 'n' :
	  return Math.sin(value);
      }
      return Double.NaN;
    }

    /**
     * Return a string describing this object
     * @return a string descibing the operator
     */
    public String toString() {
      return ""+m_operator;
    }
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5988 $");
    }
  }

  /** Operator stack */
  private Stack<String> m_operatorStack = new Stack<String>();

  /** Supported operators. l = log, b = abs, c = cos, e = exp, s = sqrt, 
      f = floor, h = ceil, r = rint, t = tan, n = sin */
  private static final String OPERATORS = "+-*/()^lbcesfhrtn";
  /** Unary functions. l = log, b = abs, c = cos, e = exp, s = sqrt, 
      f = floor, h = ceil, r = rint, t = tan, n = sin */
  private static final String UNARY_FUNCTIONS = "lbcesfhrtn";

  /** Holds the original infix expression */
  private String m_originalInfix;
  
  /** Holds the expression in postfix form */
  private Vector<ExpressionComponent> m_postFixExpVector;

  /** True if the next numeric constant or attribute index is negative */
  private boolean m_signMod = false;

  /** Holds the previous token */
  private String m_previousTok = "";

    /**
   * Handles the processing of an infix operand to postfix
   * @param tok the infix operand
   * @throws Exception if there is difficulty parsing the operand
   */
  private void handleOperand(String tok) throws Exception {
    // if it contains an 'a' then its an attribute index
    if (tok.indexOf('a') != -1) {
      m_postFixExpVector.addElement(new AttributeOperand(tok,m_signMod));
    } else {
      try {
	// should be a numeric constant
	m_postFixExpVector.addElement(new NumericOperand(tok, m_signMod));
      } catch (NumberFormatException ne) {
	throw new Exception("Trouble parsing numeric constant");
      }
    }
    m_signMod = false;
  }

  /**
   * Handles the processing of an infix operator to postfix
   * @param tok the infix operator
   * @throws Exception if there is difficulty parsing the operator
   */
  private void handleOperator(String tok) throws Exception {
    boolean push = true;

    char tokchar = tok.charAt(0);
    if (tokchar == ')') {
      String popop = " ";
      do {
	popop = (String)(m_operatorStack.pop());
	if (popop.charAt(0) != '(') {
	  m_postFixExpVector.addElement(new Operator(popop.charAt(0)));
	}
      } while (popop.charAt(0) != '(');
    } else {
      int infixToc = infixPriority(tok.charAt(0));
      while (!m_operatorStack.empty() && 
	     stackPriority(((String)(m_operatorStack.peek())).charAt(0)) 
	     >= infixToc) {
	
	// try an catch double operators and see if the current one can
	// be interpreted as the sign of an upcoming number
	if (m_previousTok.length() == 1 && 
	    isOperator(m_previousTok.charAt(0)) &&
	    m_previousTok.charAt(0) != ')') {
	  if (tok.charAt(0) == '-') {
	    m_signMod = true;
	  } else {
	    m_signMod = false;
	  }
	  push = false;
	  break;
	} else {
	  String popop = (String)(m_operatorStack.pop());
	  m_postFixExpVector.addElement(new Operator(popop.charAt(0)));
	}
      }
      if (m_postFixExpVector.size() == 0) {
	if (tok.charAt(0) == '-') {
	  m_signMod = true;
	  push = false;
	}
      }

      if (push) {
	m_operatorStack.push(tok);
      }
    }
  }

  /**
   * Converts a string containing a mathematical expression in infix form
   * to postfix form. The result is stored in the vector m_postfixExpVector
   *
   * @param infixExp the infix expression to convert
   * @throws Exception if something goes wrong during the conversion
   */
  public void convertInfixToPostfix(String infixExp) throws Exception {
    m_originalInfix = infixExp;

    infixExp = Utils.removeSubstring(infixExp, " ");
    infixExp = Utils.replaceSubstring(infixExp,"log","l");
    infixExp = Utils.replaceSubstring(infixExp,"abs","b");
    infixExp = Utils.replaceSubstring(infixExp,"cos","c");
    infixExp = Utils.replaceSubstring(infixExp,"exp","e");
    infixExp = Utils.replaceSubstring(infixExp,"sqrt","s");
    infixExp = Utils.replaceSubstring(infixExp,"floor","f");
    infixExp = Utils.replaceSubstring(infixExp,"ceil","h");
    infixExp = Utils.replaceSubstring(infixExp,"rint","r");
    infixExp = Utils.replaceSubstring(infixExp,"tan","t");
    infixExp = Utils.replaceSubstring(infixExp,"sin","n");

    StringTokenizer tokenizer = new StringTokenizer(infixExp, OPERATORS, true);
    m_postFixExpVector = new Vector<ExpressionComponent>();

    while (tokenizer.hasMoreTokens()) {
      String tok = tokenizer.nextToken();
      
      if (tok.length() > 1) {
	handleOperand(tok);
      } else {
	// probably an operator, but could be a single char operand
	if (isOperator(tok.charAt(0))) {
	  handleOperator(tok);
	} else {
	  // should be a numeric constant
	  handleOperand(tok);
	}
      }
      m_previousTok = tok;
    }
    while (!m_operatorStack.empty()) {
      String popop = (String)(m_operatorStack.pop());
      if (popop.charAt(0) == '(' || popop.charAt(0) == ')') {
	throw new Exception("Mis-matched parenthesis!");
      }
      m_postFixExpVector.addElement(new Operator(popop.charAt(0)));
    }
  }

  /**
   * Evaluate the expression using the supplied Instance.
   * Assumes that the infix expression has been converted to 
   * postfix and stored in m_postFixExpVector
   *
   * @param instance the Instance containing values to apply
   * the expression to
   * @throws Exception if something goes wrong
   */
  public double evaluateExpression(Instance instance)
    throws Exception {
    double [] vals = new double [instance.numAttributes()+1];
    for(int i = 0; i < instance.numAttributes(); i++) {
      if (instance.isMissing(i)) {
	vals[i] = Utils.missingValue();
      } else {
	vals[i] = instance.value(i);
      }
    }
    
    evaluateExpression(vals);
    return vals[vals.length - 1];
  }

  /**
   * Evaluate the expression using the supplied array of attribute values.
   * The result is stored in the last element of the array. Assumes that
   * the infix expression has been converted to postfix and stored in
   * m_postFixExpVector
   * @param vals the values to apply the expression to
   * @throws Exception if something goes wrong
   */
  public void evaluateExpression(double [] vals) throws Exception {
    Stack<Double> operands = new Stack<Double>();

    for (int i=0;i<m_postFixExpVector.size();i++) {
      Object nextob = m_postFixExpVector.elementAt(i);
      if (nextob instanceof NumericOperand) {
	operands.push(new Double(((NumericOperand)nextob).m_numericConst));
      } else if (nextob instanceof AttributeOperand) {
	double value = vals[((AttributeOperand)nextob).m_attributeIndex];
	/*if (Utils.isMissingValue(value)) {
	  vals[vals.length-1] = Utils.missingValue();
	  break;
	}*/
	if (((AttributeOperand)nextob).m_negative) {
	  value = -value;
	}
	operands.push(new Double(value));
      } else if (nextob instanceof Operator) {
	char op = ((Operator)nextob).m_operator;
	if (isUnaryFunction(op)) {
	  double operand = ((Double)operands.pop()).doubleValue();
	  double result = ((Operator)nextob).applyFunction(operand);
	  operands.push(new Double(result));
	} else {
	  double second = ((Double)operands.pop()).doubleValue();
	  double first = ((Double)operands.pop()).doubleValue();
	  double result = ((Operator)nextob).applyOperator(first,second);
	  operands.push(new Double(result));
	}
      } else {
	throw new Exception("Unknown object in postfix vector!");
      }
    }

    if (operands.size() != 1) {
      throw new Exception("Problem applying function");
    }

    Double result = ((Double)operands.pop());
    if (result.isNaN() || result.isInfinite()) {
      vals[vals.length-1] = Utils.missingValue();
    } else {
      vals[vals.length-1] = result.doubleValue();
    }
  }

  /**
   * Returns true if a token is an operator
   * @param tok the token to check
   * @return true if the supplied token is an operator
   */
  private boolean isOperator(char tok) {
    if (OPERATORS.indexOf(tok) == -1) {
      return false;
    }

    return true;
  }

  /**
   * Returns true if a token is a unary function
   * @param tok the token to check
   * @return true if the supplied token is a unary function
   */
  private boolean isUnaryFunction(char tok) {
    if (UNARY_FUNCTIONS.indexOf(tok) == -1) {
      return false;
    }

    return true;
  }

  /**
   * Return the infix priority of an operator
   * @param opp the operator
   * @return the infix priority
   */
  private int infixPriority(char opp) {
    switch (opp) {
      case 'l' : 
      case 'b' :
      case 'c' :
      case 'e' :
      case 's' :
      case 'f' :
      case 'h' :
      case 'r' :
      case 't' :
      case 'n' :
	return 3;
      case '^' :
	return 2;
      case '*' : 
	return 2;
      case '/' : 
	return 2;
      case '+' :
	return 1;
      case '-' :
	return 1;
      case '(' :
	return 4;
      case ')' :
	return 0;
      default :
	throw new IllegalArgumentException("Unrecognized operator:" + opp);
    }
  }

  /**
   * Return the stack priority of an operator
   * @param opp the operator
   * @return the stack priority
   */
  private int stackPriority(char opp) {
     switch (opp) {
       case 'l' :
       case 'b' :
       case 'c' :
       case 'e' :
       case 's' :
       case 'f' :
       case 'h' :
       case 'r' :
       case 't' :
       case 'n' :
	 return 3;
       case '^' :
	 return 2;
       case '*' : 
	 return 2;
       case '/' : 
	 return 2;
       case '+' :
	 return 1;
       case '-' :
	 return 1;
       case '(' :
	 return 0;
       case ')' :
	 return -1;
       default :
	 throw new IllegalArgumentException("Unrecognized operator:" + opp);
    }
  }

  /**
   * Return the postfix expression
   *
   * @return the postfix expression as a String
   */
  public String getPostFixExpression() {
    return m_postFixExpVector.toString();
  }

  public String toString() {
    return m_originalInfix;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5988 $");
  }
}

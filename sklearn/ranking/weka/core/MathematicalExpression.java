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
 *    MathematicalExpression.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import weka.core.mathematicalexpression.Parser;
import weka.core.mathematicalexpression.Scanner;
import java_cup.runtime.DefaultSymbolFactory;
import java_cup.runtime.SymbolFactory;

import java.io.ByteArrayInputStream;
import java.util.HashMap;

/** 
 * Class for evaluating a string adhering the following grammar:<br/>
 * 
 * <pre>
 * expr_list ::= expr_list expr_part | expr_part ;
 * expr_part ::= expr ;
 * expr      ::=   NUMBER
 *               | ( expr )
 *               | opexpr
 *               | varexpr
 *               | funcexpr
 *               ;
 * 
 * opexpr    ::=   expr + expr
 *               | expr - expr
 *               | expr * expr
 *               | expr / expr
 *               ;
 * 
 * varexpr  ::=  VARIABLE ;
 * 
 * funcexpr ::=    abs ( expr )
 *               | sqrt ( expr )
 *               | log ( expr )
 *               | exp ( expr )
 *               | sin ( expr )
 *               | cos ( expr )
 *               | tan ( expr )
 *               | rint ( expr )
 *               | floor ( expr )
 *               | pow ( expr , expr )
 *               | ceil ( expr )
 *               | ifelse ( boolexpr , expr (if true) , expr (if false) )
 *               ;
 * 
 * boolexpr ::=    BOOLEAN
 *               | true
 *               | false
 *               | expr &lt; expr
 *               | expr &lt;= expr
 *               | expr &gt; expr
 *               | expr &gt;= expr
 *               | expr = expr
 *               | ( boolexpr )
 *               | ! boolexpr
 *               | boolexpr & boolexpr
 *               | boolexpr | boolexpr
 *               ;
 * </pre>
 *
 * Code example 1:
 * <pre>
 * String expr = "pow(BASE,EXPONENT)*MULT";
 * HashMap symbols = new HashMap();
 * symbols.put("BASE", new Double(2));
 * symbols.put("EXPONENT", new Double(9));
 * symbols.put("MULT", new Double(0.1));
 * double result = MathematicalExpression.evaluate(expr, symbols);
 * System.out.println(expr + " and " + symbols + " = " + result);
 * </pre>
 * 
 * Code Example 2 (uses the "ifelse" construct):
 * <pre>
 * String expr = "ifelse(I<0,pow(BASE,I*0.5),pow(BASE,I))";
 * MathematicalExpression.TreeNode tree = MathematicalExpression.parse(expr);
 * HashMap symbols = new HashMap();
 * symbols.put("BASE", new Double(2));
 * for (int i = -10; i <= 10; i++) {
 *   symbols.put("I", new Double(i));
 *   double result = MathematicalExpression.evaluate(expr, symbols);
 *   System.out.println(expr + " and " + symbols + " = " + result);
 * }
 * </pre>
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4939 $
 */
public class MathematicalExpression
  implements RevisionHandler {
  
  /**
   * Parses and evaluates the given expression.
   * Returns the result of the mathematical expression, based on the given 
   * values of the symbols.
   * 
   * @param expr	the expression to evaluate
   * @param symbols	the symbol/value mapping
   * @return		the evaluated result
   * @throws Exception	if something goes wrong
   */
  public static double evaluate(String expr, HashMap symbols) throws Exception {
    SymbolFactory 		sf;
    ByteArrayInputStream 	parserInput;
    Parser 			parser;
    
    sf          = new DefaultSymbolFactory();
    parserInput = new ByteArrayInputStream(expr.getBytes());
    parser      = new Parser(new Scanner(parserInput, sf), sf);
    parser.setSymbols(symbols);
    parser.parse();
    
    return parser.getResult();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 4939 $");
  }
}

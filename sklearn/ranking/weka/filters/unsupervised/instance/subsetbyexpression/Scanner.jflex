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
 * Scanner.java
 * Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.instance.subsetbyexpression;

import java_cup.runtime.SymbolFactory;
import java.io.*;

/**
 * A scanner for evaluating whether an Instance is to be included in a subset
 * or not.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4939 $
 */
%%
%unicode
%char
%cup
%public
%class Scanner
%{
  // Author: FracPete (fracpete at waikato dot ac dot nz)
  // Version: $Revision: 4939 $
  protected SymbolFactory m_SymFactory;

  protected StringBuffer m_String = new StringBuffer();

  public Scanner(InputStream r, SymbolFactory sf) {
    this(r);
    m_SymFactory = sf;
  }
%}
%eofval{
    return m_SymFactory.newSymbol("EOF",sym.EOF);
%eofval}

%state STRING

%%
<YYINITIAL> {
  // operands
  "-" { return m_SymFactory.newSymbol("Minus", sym.MINUS); }
  "+" { return m_SymFactory.newSymbol("Plus", sym.PLUS); }
  "*" { return m_SymFactory.newSymbol("Times", sym.TIMES); }
  "/" { return m_SymFactory.newSymbol("Division", sym.DIVISION); }

  // boolean stuff
  "<" { return m_SymFactory.newSymbol("Less than", sym.LT); }
  "<=" { return m_SymFactory.newSymbol("Less or equal than", sym.LE); }
  ">" { return m_SymFactory.newSymbol("Greater than", sym.GT); }
  ">=" { return m_SymFactory.newSymbol("Greater or equal than", sym.GE); }
  "=" { return m_SymFactory.newSymbol("Equals", sym.EQ); }
  "is" { return m_SymFactory.newSymbol("Is", sym.IS); }
  "not" { return m_SymFactory.newSymbol("Not", sym.NOT); }
  "and" { return m_SymFactory.newSymbol("And", sym.AND); }
  "or" { return m_SymFactory.newSymbol("Or", sym.OR); }
  "true" { return m_SymFactory.newSymbol("True", sym.TRUE); }
  "false" { return m_SymFactory.newSymbol("False", sym.FALSE); }

  // functions
  "abs" { return m_SymFactory.newSymbol("Abs", sym.ABS); }
  "sqrt" { return m_SymFactory.newSymbol("Sqrt", sym.SQRT); }
  "log" { return m_SymFactory.newSymbol("Log", sym.LOG); }
  "exp" { return m_SymFactory.newSymbol("Exp", sym.EXP); }
  "sin" { return m_SymFactory.newSymbol("Sin", sym.SIN); }
  "cos" { return m_SymFactory.newSymbol("Cos", sym.COS); }
  "tan" { return m_SymFactory.newSymbol("Tan", sym.TAN); }
  "rint" { return m_SymFactory.newSymbol("Rint", sym.RINT); }
  "floor" { return m_SymFactory.newSymbol("Floor", sym.FLOOR); }
  "pow" { return m_SymFactory.newSymbol("Pow", sym.POW); }
  "ceil" { return m_SymFactory.newSymbol("Ceil", sym.CEIL); }

  // numbers and variables
  "'" { yybegin(STRING); m_String.setLength(0); }
  [0-9][0-9]*\.?[0-9]* { return m_SymFactory.newSymbol("Number", sym.NUMBER, new Double(yytext())); }
  -[0-9][0-9]*\.?[0-9]* { return m_SymFactory.newSymbol("Number", sym.NUMBER, new Double(yytext())); }
  [A][T][T][0-9][0-9]* { return m_SymFactory.newSymbol("Attribute", sym.ATTRIBUTE, new String(yytext())); }
  "CLASS" { return m_SymFactory.newSymbol("Class", sym.ATTRIBUTE, new String(yytext())); }

  // whitespaces
  [ \r\n\t\f] { /* ignore white space. */ }

  // various
  "," { return m_SymFactory.newSymbol("Comma", sym.COMMA); }
  "(" { return m_SymFactory.newSymbol("Left Bracket", sym.LPAREN); }
  ")" { return m_SymFactory.newSymbol("Right Bracket", sym.RPAREN); }
  "ismissing" { return m_SymFactory.newSymbol("Missing", sym.ISMISSING); }
}

<STRING> {
  "'" { yybegin(YYINITIAL); return m_SymFactory.newSymbol("String", sym.STRING, new String(m_String.toString())); }
  . { m_String.append(yytext()); }
}

// catch all
. { System.err.println("Illegal character: " + yytext()); }

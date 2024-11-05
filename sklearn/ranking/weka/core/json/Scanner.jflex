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
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.core.json;

import java_cup.runtime.SymbolFactory;
import java.io.*;

/**
 * A scanner for JSON data files.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5786 $
 */
%%
%cup
%public
%class Scanner
%{
  // Author: FracPete (fracpete at waikato dot ac dot nz)
  // Version: $Revision: 5786 $
  protected SymbolFactory m_SF;

  protected StringBuffer m_String = new StringBuffer();

  public Scanner(InputStream r, SymbolFactory sf) {
    this(r);
    m_SF = sf;
  }

  public Scanner(Reader r, SymbolFactory sf) {
    this(r);
    m_SF = sf;
  }
%}
%eofval{
    return m_SF.newSymbol("EOF", sym.EOF);
%eofval}

%state STRING

%%
<YYINITIAL> "{"                  { return m_SF.newSymbol("Left curly bracket", sym.LCURLY); }
<YYINITIAL> "}"                  { return m_SF.newSymbol("Right curly bracket", sym.RCURLY); }

<YYINITIAL> {
  "["                            { return m_SF.newSymbol("Left square bracket", sym.LSQUARE); }
  "]"                            { return m_SF.newSymbol("Right square bracket", sym.RSQUARE); }
  ","                            { return m_SF.newSymbol("Comma", sym.COMMA); }
  ":"                            { return m_SF.newSymbol("Colon", sym.COLON); }
  "null"                         { return m_SF.newSymbol("Null", sym.NULL); }
  "true"                         { return m_SF.newSymbol("Boolean", sym.BOOLEAN, new Boolean(yytext())); }
  "false"                        { return m_SF.newSymbol("Boolean", sym.BOOLEAN, new Boolean(yytext())); }
  [0-9][0-9]*                    { return m_SF.newSymbol("Integer", sym.INTEGER, new Integer(yytext())); }
  [0-9][0-9]*\.?[0-9]*           { return m_SF.newSymbol("Double", sym.DOUBLE, new Double(yytext())); }
  -[0-9][0-9]*\.?[0-9]*          { return m_SF.newSymbol("Double", sym.DOUBLE, new Double(yytext())); }
  \"                             { m_String.setLength(0); yybegin(STRING); }
  [ \r\n\t\f]                    { /* ignore white space. */ }
}

<STRING> {
  \"                             { yybegin(YYINITIAL); return m_SF.newSymbol("String", sym.STRING, m_String.toString()); }
  [^\n\r\"\\]+                   { m_String.append(yytext()); }
  \\\"                           { m_String.append('\"'); }
  \\b                            { m_String.append('\b'); }
  \\f                            { m_String.append('\f'); }
  \\n                            { m_String.append('\n'); }
  \\r                            { m_String.append('\r'); }
  \\t                            { m_String.append('\t'); }
  \\                             { m_String.append('\\'); }
}

// catch all
.                                { System.err.println("Illegal character: " + yytext()); }

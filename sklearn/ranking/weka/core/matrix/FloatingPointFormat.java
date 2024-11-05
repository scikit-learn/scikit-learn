/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or (at
 *    your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/*
 *    FloatingPoint.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.matrix;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.text.DecimalFormat;
import java.text.FieldPosition;

/**
 * Class for the format of floating point numbers
 *
 * @author Yong Wang
 * @version $Revision: 5953 $
 */
public class FloatingPointFormat
  extends DecimalFormat
  implements RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 4500373755333429499L;
    
  protected DecimalFormat nf ;
  protected int width;
  protected int decimal;
  protected boolean trailing = true;

  /**
   * Default constructor
   */
  public FloatingPointFormat () {
    this( 8, 5 );
  }

  public FloatingPointFormat ( int digits ) {
    this( 8, 2 );
  }

  public FloatingPointFormat( int w, int d ) {
    width = w;
    decimal = d;
    nf = new DecimalFormat( pattern(w, d) );
    nf.setPositivePrefix(" ");
    nf.setNegativePrefix("-");
  }

  public FloatingPointFormat( int w, int d, boolean trailingZeros ) {
    this( w, d );
    this.trailing = trailingZeros;
  }

  public StringBuffer format(double number, StringBuffer toAppendTo, 
			     FieldPosition pos) {
    StringBuffer s = new StringBuffer( nf.format(number) );
    if( s.length() > width ) {
      if( s.charAt(0) == ' ' && s.length() == width + 1 ) {
	s.deleteCharAt(0);
      }
      else {
	s.setLength( width );
	for( int i = 0; i < width; i ++ )
	  s.setCharAt(i, '*');
      }
    }
    else {
      for (int i = 0; i < width - s.length(); i++)  // padding
	s.insert(0,' ');
    }
    if( !trailing && decimal > 0 ) { // delete trailing zeros
      while( s.charAt( s.length()-1 ) == '0' )
	s.deleteCharAt( s.length()-1 );
      if( s.charAt( s.length()-1 ) == '.' )
	s.deleteCharAt( s.length()-1 );
    }
	
    return toAppendTo.append( s );
  }

  public static String  pattern( int w, int d ) {
    StringBuffer s = new StringBuffer();      // "-##0.00"   // fw.d
    s.append( padding(w - d - 3, '#') );
    if( d == 0) s.append('0');
    else {
      s.append("0.");
      s.append( padding( d, '0') );
    }
    return s.toString();
  }

  private static StringBuffer  padding( int n, char c ) {
    StringBuffer text = new StringBuffer();
	
    for(int i = 0; i < n; i++ ){
      text.append( c );
    }

    return text;
  }

  public int width () {
    if( !trailing ) throw new RuntimeException( "flexible width" );
    return width;
  }

  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
}

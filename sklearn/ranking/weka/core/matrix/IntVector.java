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
 *    IntVector.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.matrix;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.util.Arrays;

/**
 * A vector specialized on integers.
 * 
 * @author Yong Wang
 * @version $Revision: 5953 $
 */
public class  IntVector
  implements Cloneable, RevisionHandler {

  /** Array for internal storage of elements. */
  int[]  V;

  /** size of the vector */
  private int  sizeOfVector;


  /* ------------------------
     Constructors
     * ------------------------ */

  /** Constructs a null vector.
   */
  public IntVector(){
    V = new int[ 0 ];
    setSize( 0 );
  }
    
  /** Constructs an n-vector of zeros. 
   *  @param n    Length.
  */
  public IntVector( int n ){
    V = new int[ n ];
    setSize( n );
  }
    
  /** Constructs an n-vector of a constant
   *  @param n    Length.
  */
  public IntVector( int n, int s ){
    this(n);
    set( s );
  }
    
  /** Constructs a vector given an int array
   *  @param v the int array
  */
  public IntVector( int v[] ){
    if( v == null ) {
      V = new int[ 0 ];
      setSize( 0 );
    }
    else {
      V = new int[ v.length ];
      setSize( v.length );
      set(0, size() - 1, v, 0);
    }
  }
    
  /* ------------------------
     Public Methods
     * ------------------------ */
    
  /** Gets the size of the vector.
   *  @return Size.  */
  public int  size(){
    return sizeOfVector;
  }
    
  /** 
   * Sets the size of the vector. The provided size can't be greater than
   * the capacity of the vector.
   * @param size the new Size.
   */
  public void  setSize( int size ){
    if( size > capacity() ) 
      throw new IllegalArgumentException("insufficient capacity");
    sizeOfVector = size;
  }
  
  /** Sets the value of an element.
   *  @param s the value for the element */
  public void  set( int s ) {
    for( int i = 0; i < size(); i++ )
      set(i, s);
  }

  /** Sets the values of elements from an int array.
   *  @param i0 the index of the first element
   *  @param i1 the index of the last element
   *  @param v the int array that stores the values
   *  @param j0 the index of the first element in the int array */
  public void  set( int i0, int i1, int [] v, int j0){
    for(int i = i0; i<= i1; i++)
      set( i, v[j0 + i - i0] );
  }

  /** Sets the values of elements from another IntVector.
   *  @param i0 the index of the first element
   *  @param i1 the index of the last element
   *  @param v the IntVector that stores the values
   *  @param j0 the index of the first element in the IntVector */
  public void  set( int i0, int i1, IntVector v, int j0){
    for(int i = i0; i<= i1; i++)
      set( i, v.get(j0 + i - i0) );
  }

  /** Sets the values of elements from another IntVector.
   *  @param v the IntVector that stores the values 
   */
  public void  set( IntVector v ){
    set( 0, v.size() - 1, v, 0);
  }

  /** Generates an IntVector that stores all integers inclusively between
   *  two integers.
   *  @param i0 the first integer
   *  @param i1 the second integer 
   */
  public static IntVector  seq( int i0, int i1 ) {
    if( i1 < i0 ) throw new IllegalArgumentException("i1 < i0 ");
    IntVector v = new IntVector( i1 - i0 + 1 );
    for( int i = 0; i < i1 - i0 + 1; i++ ) {
      v.set(i, i + i0);
    }
    return v; 
  } 
  
  /** Access the internal one-dimensional array.
      @return Pointer to the one-dimensional array of vector elements. */
  public int []  getArray() {
    return V;
  }
    
  /** Sets the internal one-dimensional array.
      @param a Pointer to the one-dimensional array of vector elements. */
  protected void  setArray( int [] a ) {
    V = a;
  }
    
  /** Sorts the elements in place 
   */
  public void  sort() {
    Arrays.sort( V, 0, size() );
  }

  /** Returns a copy of the internal one-dimensional array.
      @return One-dimensional array copy of vector elements.  */
  public int[]  getArrayCopy() {
    int [] b = new int[ size() ];
    for( int i = 0; i <= size() - 1; i++ ) {
      b[i] = V[i];
    }
    return b;
  }

  /** Returns the capacity of the vector 
   */
  public int capacity() {
    return V.length;
  }

  /** Sets the capacity of the vector 
   *  @param capacity the new capacity of the vector
   */
  public void  setCapacity( int capacity ) {
    if( capacity == capacity() ) return;
    int [] old_V = V;
    int m = Math.min( capacity, size() );
    V = new int[ capacity ];
    setSize( capacity );
    set(0, m-1, old_V, 0);
  }

  /** Sets a single element.
   *  @param i    the index of the element
   *  @param s    the new value
  */
  public void  set( int i, int s ) {
    V[i] = s;
  }
    
  /** Gets the value of an element.
   *  @param i    the index of the element
   *  @return     the value of the element
  */
  public int  get( int i ) {
    return V[i];
  }
  
  /** Makes a deep copy of the vector
   */
  public IntVector  copy() { 
    return (IntVector) clone();
  }
    
  /** Clones the IntVector object.
   */
  public Object  clone() { 
    IntVector u = new IntVector( size() );
    for( int i = 0; i < size(); i++) 
      u.V[i] = V[i];
    return u;
  }
  
  /** Returns a subvector.
   *  @param i0   the index of the first element
   *  @param i1   the index of the last element
   *  @return the subvector
  */
  public IntVector  subvector( int i0, int i1 ) 
  {
    IntVector v = new IntVector( i1-i0+1 );
    v.set(0, i1 - i0, this, i0);
    return v;
  }

  /** Returns a subvector as indexed by an IntVector.
   *  @param index   the index
   *  @return the subvector
  */
  public IntVector  subvector( IntVector index ) {
    IntVector v = new IntVector( index.size() );
    for( int i = 0; i < index.size(); i++ )
      v.V[i] = V[index.V[i]];
    return v;
  }

  /**
   *  Swaps the values stored at i and j
   *  @param i the index i
   *  @param j the index j
   */
  public void  swap( int i, int j ){
    if( i == j ) return;
    int t = get( i );
    set( i, get(j) );
    set( j, t );
  }
  
  /** 
   *  Shifts an element to another position. Elements between them are
   *  shifted one position left.
   *  @param i the index of the element
   *  @param j the index of the new position */
  public void  shift( int i, int j ){
    if( i == j ) return;
    if( i < j ) {
      int t = V[i];
      for( int k = i; k <= j-1; k++ )
  	V[k] = V[k+1];
      V[j] = t;
    }
    else shift( j, i );
  }
  
  /** 
   *  Shifts an element to the end of the vector. Elements between them are
   *  shifted one position left.
   *  @param j the index of the element
   */
  public void  shiftToEnd( int j ){
    shift( j, size()-1 );
  }
  
  /** 
   * Returns true if the vector is empty
   */
  public boolean  isEmpty() {
    if( size() == 0 ) return true;
    return false;
  }

  /** Converts the IntVecor to a string
   */ 
  public String  toString() {
    return toString( 5, false );
  }
    
  /** Convert the IntVecor to a string
   *  @param digits number of digits to be shown
   *  @param trailing true if trailing zeros are to be shown
   */ 
  public String  toString( int digits, boolean trailing ) {
    if( isEmpty() ) return "null vector";

    StringBuffer text = new StringBuffer();
    FlexibleDecimalFormat nf = new FlexibleDecimalFormat( digits, 
							  trailing );
    nf.grouping( true );
    for( int i = 0; i < size(); i ++ ) nf.update( get(i) );
    int count = 0;
    int width = 80;
    String number;
    for( int i = 0; i < size(); i++ ) {
      number = nf.format(get(i));
      count += 1 + number.length();
      if( count > width-1 ) { 
	text.append('\n'); 
	count = 1 + number.length();
      }
      text.append( " " + number );
    }
	
    return text.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }

  /** 
   *  Tests the IntVector class
   */
  public static void  main( String args[] ) {
    
    IntVector u = new IntVector();
    System.out.println( u );

    IntVector v = IntVector.seq(10, 25);
    System.out.println( v );

    IntVector w = IntVector.seq(25, 10);
    System.out.println( w );

  }
}

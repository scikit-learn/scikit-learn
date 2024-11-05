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
 *    DoubleVector.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.matrix;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.lang.reflect.Method;
import java.util.Arrays;

/**
 * A vector specialized on doubles.
 * 
 * @author Yong Wang
 * @version $Revision: 5953 $
 */
public class  DoubleVector
  implements Cloneable, RevisionHandler {

  double[] V; // array for internal storage of elements.

  private int  sizeOfVector;      // size of the vector
  
  /* ------------------------
     Constructors
     * ------------------------ */

  /** Constructs a null vector.
   */
  public DoubleVector() {
    this( 0 );
  }
    
  /** Constructs an n-vector of zeros. 
      @param n    length.
  */
  public DoubleVector( int n ){
    V = new double[ n ];
    setSize( n );
  }
    
  /** Constructs a constant n-vector.
      @param n    length.
      @param s    the scalar value used to fill the vector
  */
  public DoubleVector( int n, double s ){
    this( n );
    set( s );
  }
    
  /** Constructs a vector directly from a double array
   *  @param v   the array
   */
  public DoubleVector( double v[] ){
    if( v == null ) {
      V = new double[0];
      setSize( 0 );
    }
    else {
      V = v;
      setSize( v.length );
    }
  }
    
  /* ------------------------
   *  Public Methods
   * ------------------------ */
    
  /** Set a single element.
   *  @param i    Index.
   *  @param s    a[i].
   */
  public void  set( int i, double s ) {
    
    V[i] = s;
  }
    
  /** Set all elements to a value
   *  @param s    the value
   */
  public void  set( double s ) {
    set(0, size()-1, s);
  }

  /** Set some elements to a value
   *  @param i0 the index of the first element
   *  @param i1 the index of the second element
   *  @param s the value 
   */
  public void set( int i0, int i1, double s ) {

    for(int i = i0; i <= i1; i++ )
      V[i] = s;
  }

  /** Set some elements using a 2-D array
   *  @param i0 the index of the first element
   *  @param i1 the index of the second element
   *  @param j0 the index of the starting element in the 2-D array
   *  @param v the values
   */
  public void  set( int i0, int i1, double [] v, int j0){
    for(int i = i0; i<= i1; i++)
      V[i] = v[j0 + i - i0];
  }
    
  /** Set the elements using a DoubleVector
   *  @param v the DoubleVector
   */
  public void  set( DoubleVector v ){
    set( 0, v.size() - 1, v, 0);
  }
  
  /** Set some elements using a DoubleVector.
   *  @param i0 the index of the first element
   *  @param i1 the index of the second element
   *  @param v the DoubleVector
   *  @param j0 the index of the starting element in the DoubleVector
   */
  public void  set( int i0, int i1, DoubleVector v, int j0){
    for(int i = i0; i<= i1; i++)
      V[i] = v.V[j0 + i - i0];
  }
  
  /** Access the internal one-dimensional array.
      @return     Pointer to the one-dimensional array of vector elements.
  */
  public double []  getArray() {
    return V;
  }
    
  void  setArray( double [] a ) {
    V = a;
  }

  /** Returns a copy of the DoubleVector usng a double array.
      @return the one-dimensional array.  */
  public double[] getArrayCopy() {
    double v[] = new double[size()];
    
    for(int i= 0; i < size(); i++ ) 
      v[i] = V[i];
    
    return v;
  }
    
  /** Sorts the array in place */
  public void  sort() {
    Arrays.sort( V, 0, size() );
  }

  /** Sorts the array in place with index returned */
  public IntVector  sortWithIndex() {
    IntVector index = IntVector.seq( 0, size()-1 );
    sortWithIndex( 0, size()-1, index );
    return index;
  }
  
  /** Sorts the array in place with index changed 
   *  @param xi   first index
   *  @param xj   last index
   *  @param index array that stores all indices
   */
  public void  sortWithIndex( int xi, int xj, IntVector index ) {
    if( xi < xj ) { 
      double x, f, k;
      int xm = (int) (xi + xj) / 2; // median index
      x = Math.min( V[xi],             // median of three
		    Math.max( V[xm], V[xj])); 
      int i = xi;
      int j = xj;
      while( i < j ) {
	while( V[i] < x && i < xj ) i++;
	while( V[j] > x && j > xi ) j--;
	if( i <= j ){
	  swap(i, j);
	  index.swap(i, j);
	  i++;
	  j--;
	}
      }
      sortWithIndex(xi, j, index);
      sortWithIndex(i, xj, index);
    }
  }
  
  /** Gets the size of the vector.
      @return     the size
  */
  public int  size(){
    return sizeOfVector;
  }
    
  /** 
   *  Sets the size of the vector
   *  @param m the size
   */ 
  public void  setSize( int m ){
    if( m > capacity() ) 
      throw new IllegalArgumentException("insufficient capacity");
    sizeOfVector = m;
  }
    
  /** Gets the capacity of the vector.
   *  @return     the capacity.
   */
  public int  capacity() {
    if( V == null ) return 0;
    return V.length;
  } 

  /** Sets the capacity of the vector
   *  @param n the capacity.  
   */
  public void  setCapacity ( int n ) {
    if( n == capacity() ) return;
    double [] oldV = V;
    int m = Math.min( n, size() );
    V = new double[ n ];
    setSize( m );
    set(0, m-1, oldV, 0);
  }

  /** Gets a single element.
   *  @param i    Index.
   *  @return     the value of the i-th element
   */
  public double  get( int i ) {
    return V[i];
  }
    
  /** 
   *  Adds a value to an element 
   *  @param i  the index of the element 
   *  @param s the value
   */
  public void  setPlus( int i, double s ) {
    V[i] += s;
  }
    
  /** 
   *  Multiplies a value to an element 
   *  @param i  the index of the element 
   *  @param s the value
   */
  public void  setTimes( int i, double s ) {
    V[i] *= s;
  }
    
  /**
   *  Adds an element into the vector
   *  @param x  the value of the new element
   */
  public void addElement( double x ) {
    if( capacity() == 0 ) setCapacity( 10 );
    if( size() == capacity() ) setCapacity( 2 * capacity() );
    V[size()] = x;
    setSize( size() + 1 );
  }
  
  /**
   *  Returns the squared vector 
   */
  public DoubleVector square() {
    DoubleVector v = new DoubleVector( size() ); 
    for(int i = 0; i < size(); i++ ) v.V[i] = V[i] * V[i];
    return v;
  }

  /**
   *  Returns the square-root of all the elements in the vector 
   */
  public DoubleVector sqrt() {
    DoubleVector v = new DoubleVector( size() ); 
    for(int i = 0; i < size(); i++ ) v.V[i] = Math.sqrt(V[i]);
    return v;
  }

  /** Makes a deep copy of the vector
   */
  public DoubleVector  copy() { 
    return (DoubleVector) clone();
  }
    
  /** Clones the DoubleVector object.
   */
  public Object  clone() { 
    int n = size();
    DoubleVector u = new DoubleVector( n );
    for( int i = 0; i < n; i++) 
      u.V[i] = V[i];
    return u;
  }
    
  /** 
   * Returns the inner product of two DoubleVectors
   * @param v the second DoubleVector
   * @return the product
   */
  public double  innerProduct(DoubleVector v) {
    if(size() != v.size()) 
      throw new IllegalArgumentException("sizes unmatch");
    double p = 0;
    for (int i = 0; i < size(); i++) {
      p += V[i] * v.V[i];
    }
    return p;
  }
    
  /** 
   * Returns the signs of all elements in terms of -1, 0 and +1.
   */
  public DoubleVector sign() 
  {
    DoubleVector s = new DoubleVector( size() );
    for( int i = 0; i < size(); i++ ) {
      if( V[i] > 0 ) s.V[i] = 1;
      else if( V[i] < 0 ) s.V[i] = -1;
      else s.V[i] = 0;
    }
    return s;
  } 

  /** Returns the sum of all elements in the vector.
   */
  public double  sum() 
  {
    double s = 0;
    for( int i=0; i< size(); i++) s += V[i];
    return s;
  }
    
  /** Returns the squared sum of all elements in the vector.
   */
  public double  sum2()
  {
    double s2 = 0;
    for( int i=0; i< size(); i++) s2 += V[i] * V[i];
    return s2;
  }
  
  /** Returns the L1-norm of the vector
   */
  public double norm1()
  {
    double s = 0;
    for( int i=0; i< size(); i++) s += Math.abs(V[i]);
    return s;
  }
  
  /** Returns the L2-norm of the vector
   */
  public double norm2()
  {
    return Math.sqrt( sum2() );
  }

  /** Returns ||u-v||^2
   *  @param v the second vector
   */
  public double sum2( DoubleVector v ) 
  {
    return minus( v ).sum2();
  }
  
  /** Returns a subvector.
   *  @param i0   the index of the first element
   *  @param i1   the index of the last element
   *  @return     v[i0:i1]
   */
  public DoubleVector  subvector( int i0, int i1 ) 
  {
    DoubleVector v = new DoubleVector( i1-i0+1 );
    v.set(0, i1 - i0, this, i0);
    return v;
  }
  
  /** Returns a subvector.
   *  @param index stores the indices of the needed elements
   *  @return     v[index]
   */
  public DoubleVector  subvector( IntVector index ) {
    DoubleVector v = new DoubleVector( index.size() );
    for( int i = 0; i < index.size(); i++ )
      v.V[i] = V[index.V[i]];
    return v;
  }

  /** Returns a vector from the pivoting indices. Elements not indexed are
   *  set to zero.
   *  @param index stores the pivoting indices
   *  @param length the total number of the potential elements
   *  @return the subvector */
  public DoubleVector  unpivoting( IntVector index, int length ) {
    if( index.size() > length ) 
      throw new IllegalArgumentException("index.size() > length ");
    DoubleVector u = new DoubleVector( length );
    for( int i = 0; i < index.size(); i++ ) {
      u.V[index.V[i]] =  V[i];
    }
    return u;
  }
  
  /** Adds a value to all the elements 
   *  @param x the value
   */
  public DoubleVector  plus ( double x ) {
    return copy().plusEquals( x );	
  }
  
  /** Adds a value to all the elements in place
   *  @param x the value
   */
  public DoubleVector plusEquals ( double x ) {
    for( int i = 0; i < size(); i++ )
      V[i] += x;
    return this;
  }
  
  /** 
   *  Adds another vector element by element 
   *  @param v the second vector
   */
  public DoubleVector  plus( DoubleVector v ) {
    return copy().plusEquals( v );
  }
  
  /** 
   *  Adds another vector in place element by element 
   *  @param v the second vector
   */
  public DoubleVector  plusEquals( DoubleVector v ) {
    for(int i = 0; i < size(); i++ )
      V[i] += v.V[i];
    return this;
  }
  
  /** 
   *  Subtracts a value
   *  @param x the value
   */
  public DoubleVector  minus( double x ) {
    return plus( -x );
  }
  
  /** 
   *  Subtracts a value in place
   *  @param x the value
   */
  public DoubleVector  minusEquals( double x ) {
    plusEquals( -x );
    return this;
  }
  
  /** 
   *  Subtracts another DoubleVector element by element 
   *  @param v the second DoubleVector
   */
  public DoubleVector  minus( DoubleVector v ) {
    return copy().minusEquals( v );
  }
  
  /** 
   *  Subtracts another DoubleVector element by element in place
   *  @param v the second DoubleVector
   */
  public DoubleVector  minusEquals( DoubleVector v ) {
    for(int i = 0; i < size(); i++ )
      V[i] -=  v.V[i];
    return this;
  }
    
  /** Multiplies a scalar
      @param s    scalar
      @return     s * v
  */
  public DoubleVector  times( double s ) {
    return copy().timesEquals( s );
  }
    
  /** Multiply a vector by a scalar in place, u = s * u
      @param s    scalar
      @return     replace u by s * u
  */
  public DoubleVector  timesEquals( double s ) {
    for (int i = 0; i < size(); i++) {
      V[i] *= s;
    }
    return this;
  }
    
  /** 
   *  Multiplies another DoubleVector element by element
   *  @param v the second DoubleVector
   */
  public DoubleVector  times( DoubleVector v ) {
    return copy().timesEquals( v ); 
    
  }
    
  /** 
   *  Multiplies another DoubleVector element by element in place 
   *  @param v the second DoubleVector
   */
  public DoubleVector  timesEquals( DoubleVector v ) {
    for(int i = 0; i < size(); i++ )
      V[i] *= v.V[i];
    return this;
  }
  
  /** 
   *  Divided by another DoubleVector element by element
   *  @param v the second DoubleVector
   */
  public DoubleVector  dividedBy ( DoubleVector v ) {
    return copy().dividedByEquals( v );
  }
  
  /** 
   *  Divided by another DoubleVector element by element in place 
   *  @param v the second DoubleVector
   */
  public DoubleVector  dividedByEquals ( DoubleVector v ) {
    for( int i = 0; i < size(); i++ ) {
      V[i] /= v.V[i];
    }
    return this;
  }
  
  /** 
   *  Checks if it is an empty vector
   */
  public boolean  isEmpty() {
    if( size() == 0 ) return true;
    return false;
  }
  
  /** 
   * Returns a vector that stores the cumulated values of the original
   * vector */
  public DoubleVector cumulate() 
  {
    return copy().cumulateInPlace();
  }
	
  /** 
   * Cumulates the original vector in place 
   */
  public DoubleVector cumulateInPlace() 
  {
    for (int i = 1; i < size(); i++) {
      V[i] += V[i-1];
    }
    return this;
  }
	
  /** 
   * Returns the index of the maximum. <p>
   * If multiple maximums exist, the index of the first is returned.
   */
  public int  indexOfMax()
  {
    int index = 0;
    double ma = V[0];

    for( int i = 1; i < size(); i++ ){
      if( ma < V[i] ) {
	ma = V[i];
	index = i;
      }
    }
    return index;
  }
  

  /** 
   * Returns true if vector not sorted
   */
  public boolean unsorted () {
    if( size() < 2 ) return false;
    for( int i = 1; i < size(); i++ ) {
      if( V[i-1] > V[i] )
	return true;
    }
    return false;
  }
  
  /**
   *  Combine two vectors together
   *  @param v the second vector
   */
  public DoubleVector  cat( DoubleVector v ) {
    DoubleVector w = new DoubleVector( size() + v.size() );
    w.set(0, size() - 1, this, 0);
    w.set(size(), size() + v.size()-1, v, 0);
    return w;
  }
    
  /**
   *  Swaps the values stored at i and j
   *  @param i the index i
   *  @param j the index j
   */
  public void  swap( int i, int j ){
    if( i == j ) return;
    double t = V[i];
    V[i] = V[j];
    V[j] = t;
  }

  /**
   *  Returns the maximum value of all elements
   */
  public double max () {
    if( size() < 1 ) throw new IllegalArgumentException("zero size");
    double ma = V[0];
    if( size() < 2 ) return ma;
    for( int i = 1; i < size(); i++ ) {
      if( V[i] > ma ) ma = V[i];
    }
    return ma;
  }
  

  /**
   *  Applies a method to the vector
   *  @param className the class name
   *  @param method the method
   */
  public DoubleVector map( String className, String method ) {
    try {
      Class<?> c = Class.forName( className );
      Class [] cs = new Class[1]; 
      cs[ 0 ] = Double.TYPE;
      Method m = c.getMethod( method, cs );
      
      DoubleVector w = new DoubleVector( size() );
      Object [] obj = new Object[1];
      for( int i = 0; i < size(); i++ ) {
	obj[0] = new Double( V[i] );
	w.set( i, Double.parseDouble(m.invoke( null, obj ).toString()) ); 
      }
      return w;
    }
    catch ( Exception e ) {
      e.printStackTrace();
      System.exit(1);
    }
    return null;
  }

  /**
   * Returns the reverse vector
   */ 
  public DoubleVector  rev() {
    int n = size();
    DoubleVector w = new DoubleVector( n );
    for(int i = 0; i < n; i++ )
      w.V[i] = V[n-i-1];
    return w;
  }
  
  /**
   * Returns a random vector of uniform distribution
   * @param n the size of the vector
   */ 
  public static DoubleVector  random( int n ) {
    DoubleVector v = new DoubleVector( n );
    for (int i = 0; i < n; i++) {
      v.V[i] = Math.random();
    }
    return v;
  }

  /** Convert the DoubleVecor to a string
   */ 
  public String  toString() {
    return toString( 5, false );
  }
    
  /** Convert the DoubleVecor to a string
   *  @param digits the number of digits after decimal point
   *  @param trailing true if trailing zeros are to be shown
   */ 
  public String  toString( int digits, boolean trailing ) {
    if( isEmpty() ) return "null vector";

    StringBuffer text = new StringBuffer();
    FlexibleDecimalFormat nf = new FlexibleDecimalFormat( digits, 
							  trailing );
    nf.grouping( true );
    for( int i = 0; i < size(); i ++ ) nf.update( V[i] );
    int count = 0;
    int width = 80;
    String number;
    for( int i = 0; i < size(); i++ ) {
      number = nf.format(V[i]);
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

  public static void  main( String args[] ) {

    
    DoubleVector u = random(10);
    DoubleVector v = random(10);
    DoubleVector a = random(10);
    DoubleVector w = a; 

    System.out.println( random(10).plus(v).plus(w) );

  }
}


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
 *    FastVector.java
 *    Copyright (C) 1999, 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.util.Enumeration;
import java.util.Collection;
import java.util.ArrayList;

/**
 * Simple extension of ArrayList. Exists for legacy reasons.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
@Deprecated 
public class FastVector<E> extends ArrayList<E> implements Copyable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -2173635135622930169L;

  /**
   * Constructs an empty vector with initial
   * capacity zero.
   */
  public FastVector() {
    super();
  }

  /**
   * Constructs a vector with the given capacity.
   *
   * @param capacity the vector's initial capacity
   */
  public FastVector(int capacity) {
    super(capacity);
  }

  /**
   * Adds an element to this vector. Increases its
   * capacity if its not large enough.
   *
   * @param element the element to add
   */
  public final void addElement(E element) {
    add(element);
  }

  /**
   * Produces a shallow copy of this vector.
   *
   * @return the new vector
   */
  public final FastVector<E> copy() {
    return Utils.cast(clone());
  }

  /**
   * Clones the vector and shallow copies all its elements.
   * The elements have to implement the Copyable interface.
   * 
   * @return the new vector
   */
  public final FastVector<E> copyElements() {

    FastVector<E> copy = copy();
    for (int i = 0; i < size(); i++) {
      copy.set(i, Utils.<E>cast(((Copyable)get(i)).copy()));
    }
    return copy;
  }

  /**
   * Returns the element at the given position.
   *
   * @param index the element's index
   * @return the element with the given index
   */
  public final E elementAt(int index) {
    return get(index);
  }

  /**
   * Returns an enumeration of this vector.
   *
   * @return an enumeration of this vector
   */
  public final Enumeration elements() {
    return new WekaEnumeration(this);
  }

  /**
   * Returns an enumeration of this vector, skipping the
   * element with the given index.
   *
   * @param index the element to skip
   * @return an enumeration of this vector
   */
  public final Enumeration elements(int index) {
    return new WekaEnumeration(this, index);
  }

  /**
   * Returns the first element of the vector.
   *
   * @return the first element of the vector
   */
  public final  E firstElement() {
    return get(0);
  }

  /**
   * Inserts an element at the given position.
   *
   * @param element the element to be inserted
   * @param index the element's index
   */
  public final void insertElementAt(E element, int index) {
    add(index, element);
  }

  /**
   * Returns the last element of the vector.
   *
   * @return the last element of the vector
   */
  public final E lastElement() {
    return get(size() - 1);
  }

  /**
   * Deletes an element from this vector.
   *
   * @param index the index of the element to be deleted
   */
  public final void removeElementAt(int index) {
    remove(index);
  }

  /**
   * Removes all components from this vector and sets its 
   * size to zero. 
   */
  public final void removeAllElements() {
    clear();
  }

  /**
   * Appends all elements of the supplied vector to this vector.
   *
   * @param toAppend the FastVector containing elements to append.
   */
  public final void appendElements(Collection<? extends E> toAppend) {
    addAll(toAppend);
  }

  /**
   * Sets the vector's capacity to the given value.
   *
   * @param capacity the new capacity
   */
  public final void setCapacity(int capacity) {
    ensureCapacity(capacity);
  }

  /**
   * Sets the element at the given index.
   *
   * @param element the element to be put into the vector
   * @param index the index at which the element is to be placed
   */
  public final void setElementAt(E element, int index) {
    set(index, element);
  }

  /**
   * Swaps two elements in the vector.
   *
   * @param first index of the first element
   * @param second index of the second element
   */
  public final void swap(int first, int second) {
    
    E in = get(first);
    set(first, get(second));
    set(second, in);
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


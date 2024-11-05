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
 * Stack.java
 * Copyright (C) 2006 Alina Beygelzimer and Sham Kakade and John Langford
 */

package weka.core.neighboursearch.covertrees;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Class implementing a stack.
 * 
 * @param <T> The type of elements to be stored in 
 * the stack.
 * @author Alina Beygelzimer (original C++ code)
 * @author Sham Kakade (original C++ code)
 * @author John Langford (original C++ code)
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz) (Java port)
 * @version $Revision: 5953 $
 */
public class Stack<T>
  implements Serializable, RevisionHandler {

  /** for serialization. */
  private static final long serialVersionUID = 5604056321825539264L;
  
  /** The number of elements in the stack. */
  public int length;

  /** The elements inside the stack. */
  public ArrayList<T> elements;
  
  /** Constructor. */
  public Stack() { 
    length=0; 
    elements = new ArrayList<T>();
  }

  /**
   * Constructor. 
   * @param capacity The initial capacity of the stack.
   */
  public Stack(int capacity) { 
    length=0; 
    elements = new ArrayList<T>(capacity);
  }
  
  /**
   * Returns the last element in the stack.
   * @return The last element.
   */
  public T last() { 
    return elements.get(length-1);
  }
  
  /**
   * Returns the ith element in the stack.
   * @param i The index of the element to return.
   * @return The ith element. 
   */
  public T element(int i) { 
    return elements.get(i);
  }
  
  /**
   * Sets the ith element in the stack. 
   * @param i The index at which the element is
   * to be inserted. 
   * @param e The element to insert. 
   */
  public void set(int i, T e) {
    elements.set(i, e);
  }
  
  /**
   * Returns a sublist of the elements in the 
   * stack. 
   * @param beginIdx The start index of the 
   * sublist. 
   * @param uptoLength The length of the 
   * sublist.
   * @return The sublist starting from 
   * beginIdx and of length uptoLength.
   */
  public List subList(int beginIdx, int uptoLength) {
    return elements.subList(beginIdx, uptoLength);
  }
  
  /** Removes all the elements from the stack. */
  public void clear() {
    elements.clear();
    length=0;
  }
  
  /** 
   * Adds all the given elements in the stack.
   * @param c The collection of elements to add
   * in the stack.
   */ 
  public void addAll(Collection<? extends T> c) {
    elements.addAll(c);
    length = c.size();
  }
  
  /**
   * Replace all elements in the stack with 
   * the elements of another given stack.
   * It first removes all the elements 
   * currently in the stack, and then adds all
   * the elements of the provided stack. 
   * @param s The stack whose elements should 
   * be put in this stack.
   */
  public void replaceAllBy(Stack<T> s) {
    elements.clear();
    elements.addAll(s.elements);
    length = elements.size();
  }
  
  /** 
   * Pops (removes) the first (last added) 
   * element in the stack.
   * @return The poped element. 
   */
  public T pop() {
    length--;
    return elements.remove(length);    
  }
  
  /**
   * Pushes the given element to the stack.
   * @param new_ele The element to be pushed
   * to the stack. 
   */
  public void push(T new_ele) {
    length++;
    elements.add(new_ele);
  }
  
  /**
   * Pushes the given element onto the given 
   * stack. 
   * @param v The stack onto push the element.
   * @param new_ele The element to push.
   */
  public void push(Stack<T> v, T new_ele) {
    length++; 
    v.elements.add(new_ele);
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

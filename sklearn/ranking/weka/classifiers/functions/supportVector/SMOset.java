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
 *    SMOset.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions.supportVector;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/**
 * Stores a set of integer of a given size.
 * 
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class SMOset
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -8364829283188675777L;
  
  /** The current number of elements in the set */
  private int m_number;

  /** The first element in the set */
  private int m_first;

  /** Indicators */
  private boolean[] m_indicators;

  /** The next element for each element */
  private int[] m_next;

  /** The previous element for each element */
  private int[] m_previous;

  /**
   * Creates a new set of the given size.
   */
  public SMOset(int size) {
      
    m_indicators = new boolean[size];
    m_next = new int[size];
    m_previous = new int[size];
    m_number = 0;
    m_first = -1;
  }
 
  /**
   * Checks whether an element is in the set.
   */
  public boolean contains(int index) {

    return m_indicators[index];
  }

  /**
   * Deletes an element from the set.
   */
  public void delete(int index) {

    if (m_indicators[index]) {
      if (m_first == index) {
	m_first = m_next[index];
      } else {
	m_next[m_previous[index]] = m_next[index];
      }
      if (m_next[index] != -1) {
	m_previous[m_next[index]] = m_previous[index];
      }
      m_indicators[index] = false;
      m_number--;
    }
  }

  /**
   * Inserts an element into the set.
   */
  public void insert(int index) {

    if (!m_indicators[index]) {
      if (m_number == 0) {
	m_first = index;
	m_next[index] = -1;
	m_previous[index] = -1;
      } else {
	m_previous[m_first] = index;
	m_next[index] = m_first;
	m_previous[index] = -1;
	m_first = index;
      }
      m_indicators[index] = true;
      m_number++;
    }
  }

  /** 
   * Gets the next element in the set. -1 gets the first one.
   */
  public int getNext(int index) {

    if (index == -1) {
      return m_first;
    } else {
      return m_next[index];
    }
  }

  /**
   * Prints all the current elements in the set.
   */
  public void printElements() {

    for (int i = getNext(-1); i != -1; i = getNext(i)) {
      System.err.print(i + " ");
    }
    System.err.println();
    for (int i = 0; i < m_indicators.length; i++) {
      if (m_indicators[i]) {
	System.err.print(i + " ");
      }
    }
    System.err.println();
    System.err.println(m_number);
  }

  /** 
   * Returns the number of elements in the set.
   */
  public int numElements() {
      
    return m_number;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.5 $");
  }
}

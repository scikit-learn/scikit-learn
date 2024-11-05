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
 *    TestSetEvent.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;

import java.util.EventObject;

/**
 * Event encapsulating a test set
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 4761 $
 */
public class TestSetEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 8780718708498854231L;
  
  /**
   * The test set instances
   */
  protected Instances m_testSet;
  private boolean m_structureOnly;
  
  /**
   * What run number is this training set from. 
   */
  protected int m_runNumber = 1;
  
  
  /**
   * Maximum number of runs. 
   */
  protected int m_maxRunNumber = 1;

  /**
   * what number is this test set (ie fold 2 of 10 folds)
   */
  protected int m_setNumber;

  /**
   * Maximum number of sets (ie 10 in a 10 fold)
   */
  protected int m_maxSetNumber;

  /**
   * Creates a new <code>TestSetEvent</code>
   *
   * @param source the source of the event
   * @param testSet the test instances
   */
  public TestSetEvent(Object source, Instances testSet) {
    super(source);
    m_testSet = testSet;
    if (m_testSet != null && m_testSet.numInstances() == 0) {
      m_structureOnly = true;
    }
  }

  /**
   * Creates a new <code>TestSetEvent</code>
   *
   * @param source the source of the event
   * @param testSet the test instances
   * @param setNum the number of the test set
   * @param maxSetNum the maximum number of sets
   */
  public TestSetEvent(Object source, Instances testSet, 
                      int setNum, int maxSetNum) {
    this(source, testSet);
    m_setNumber = setNum;
    m_maxSetNumber = maxSetNum;
  }
  
  /**
   * Creates a new <code>TestSetEvent</code>
   * 
   * @param source the source of the event
   * @param testSet the test instances
   * @param runNum the run number that the test set belongs to
   * @param maxRunNum the maximum run number
   * @param setNum the number of the test set
   * @param maxSetNum the maximum number of sets
   */
  public TestSetEvent(Object source, Instances testSet,
      int runNum, int maxRunNum, int setNum, int maxSetNum) {
    this(source, testSet, setNum, maxSetNum);
    
    m_runNumber = runNum;
    m_maxRunNumber = maxRunNum;
  }

  /**
   * Get the test set instances
   *
   * @return an <code>Instances</code> value
   */
  public Instances getTestSet() {
    return m_testSet;
  }
  
  /**
   * Get the run number that this training set belongs to.
   * 
   * @return the run number for this training set.
   */
  public int getRunNumber() {
    return m_runNumber;
  }
  
  /**
   * Get the maximum number of runs.
   * 
   * @return return the maximum number of runs.
   */
  public int getMaxRunNumber() {
    return m_maxRunNumber;
  }

  /**
   * Get the test set number (eg. fold 2 of a 10 fold split)
   *
   * @return an <code>int</code> value
   */
  public int getSetNumber() {
    return m_setNumber;
  }

  /**
   * Get the maximum set number
   *
   * @return an <code>int</code> value
   */
  public int getMaxSetNumber() {
    return m_maxSetNumber;
  }

  /**
   * Returns true if the encapsulated instances
   * contain just header information
   *
   * @return true if only header information is
   * available in this DataSetEvent
   */
  public boolean isStructureOnly() {
    return m_structureOnly;
  }
}

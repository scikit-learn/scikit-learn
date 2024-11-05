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
 *    TrainingSetEvent.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;

import java.util.EventObject;

/**
 * Event encapsulating a training set
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 4761 $
 */
public class TrainingSetEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 5872343811810872662L;
  
  /**
   * The training instances
   */
  protected Instances m_trainingSet;
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
   * what number is this training set (ie fold 2 of 10 folds)
   */
  protected int m_setNumber;

  /**
   * Maximum number of sets (ie 10 in a 10 fold)
   */
  protected int m_maxSetNumber;

  /**
   * Creates a new <code>TrainingSetEvent</code>
   *
   * @param source the source of the event
   * @param trainSet the training instances
   */
  public TrainingSetEvent(Object source, Instances trainSet) {
    super(source);
    m_trainingSet = trainSet;
    if (m_trainingSet != null && m_trainingSet.numInstances() == 0) {
      m_structureOnly = true;
    }
  }

  /**
   * Creates a new <code>TrainingSetEvent</code>
   *
   * @param source the source of the event
   * @param trainSet the training instances
   * @param setNum the number of the training set
   * @param maxSetNum the maximum number of sets
   */
  public TrainingSetEvent(Object source, Instances trainSet, int setNum, int maxSetNum) {
    this(source, trainSet);
    m_setNumber = setNum;
    m_maxSetNumber = maxSetNum;
  }
  
  /**
   * Creates a new <code>TrainingSetEvent</code>
   * 
   * @param source the source of the event
   * @param trainSet the training instances
   * @param runNum the run number that the training set belongs to
   * @param maxRunNum the maximum run number
   * @param setNum the number of the training set
   * @param maxSetNum the maximum number of sets
   */
  public TrainingSetEvent(Object source, Instances trainSet, int runNum,
      int maxRunNum, int setNum, int maxSetNum) {
    this(source, trainSet, setNum, maxSetNum);
    
    m_runNumber = runNum;
    m_maxRunNumber = maxRunNum; 
  }

  /**
   * Get the training instances
   *
   * @return an <code>Instances</code> value
   */
  public Instances getTrainingSet() {
    return m_trainingSet;
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
   * Get the set number (eg. fold 2 of a 10 fold split)
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

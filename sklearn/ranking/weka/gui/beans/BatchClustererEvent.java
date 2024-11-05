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
 *    BatchClustererEvent.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.clusterers.Clusterer;

import java.util.EventObject;

/**
 * Class encapsulating a built clusterer and a batch of instances to
 * test on.
 *
 * @author <a href="mailto:mutter@cs.waikato.ac.nz">Stefan Mutter</a>
 * @version $Revision: 1.4 $
 * @since 1.0
 * @see EventObject
 */
public class BatchClustererEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 7268777944939129714L;

  /**
   * The clusterer
   */
  protected Clusterer m_clusterer;
  //  protected Instances m_trainingSet;

  /**
   * Training or Test Instances
   */
  //  protected Instances m_testSet;
  protected DataSetEvent m_testSet;

  /**
   * The set number for the test set
   */
  protected int m_setNumber;

  /**
   * Indicates if m_testSet is a training or a test set. 0 for test, >0 for training
   */
  protected int m_testOrTrain;

  /**
   * The last set number for this series
   */
  protected int m_maxSetNumber;
  
  public static int TEST = 0;
  public static int TRAINING = 1;

  /**
   * Creates a new <code>BatchClustererEvent</code> instance.
   *
   * @param source the source object
   * @param scheme a Clusterer
   * @param tstI the training/test instances
   * @param setNum the set number of the training/testinstances
   * @param maxSetNum the last set number in the series
   * @param testOrTrain 0 if the set is a test set, >0 if it is a training set
   */
  public BatchClustererEvent(Object source, Clusterer scheme, DataSetEvent tstI, int setNum, int maxSetNum, int testOrTrain) {
    super(source);
    //    m_trainingSet = trnI;
    m_clusterer = scheme;
    m_testSet = tstI;
    m_setNumber = setNum;
    m_maxSetNumber = maxSetNum;
    if(testOrTrain == 0)
        m_testOrTrain = TEST;
    else
        m_testOrTrain = TRAINING;
  }


  /**
   * Get the clusterer
   *
   * @return the clusterer
   */
  public Clusterer getClusterer() {
    return m_clusterer;
  }

  /**
   * Get the training/test set
   *
   * @return the training/testing instances
   */
  public DataSetEvent getTestSet() {
    return m_testSet;
  }

  /**
   * Get the set number (ie which fold this is)
   *
   * @return the set number for the training and testing data sets
   * encapsulated in this event
   */
  public int getSetNumber() {
    return m_setNumber;
  }

  /**
   * Get the maximum set number (ie the total number of training
   * and testing sets in the series).
   *
   * @return the maximum set number
   */
  public int getMaxSetNumber() {
    return m_maxSetNumber;
  }
  
  /**
   * Get whether the set of instances is a test or a training set
   *
   * @return 0 for test set, >0 fro training set
   */
  public int getTestOrTrain(){
    return m_testOrTrain;
  }
}


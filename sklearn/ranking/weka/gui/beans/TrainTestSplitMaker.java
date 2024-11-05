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
 *    TrainTestSplitMaker.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * Bean that accepts data sets, training sets, test sets and produces
 * both a training and test set by randomly spliting the data
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 4813 $
 */
public class TrainTestSplitMaker
  extends AbstractTrainAndTestSetProducer
  implements DataSourceListener, TrainingSetListener, TestSetListener,
	     UserRequestAcceptor, EventConstraints, Serializable {

  /** for serialization */
  private static final long serialVersionUID = 7390064039444605943L;

  private double m_trainPercentage = 66;
  private int m_randomSeed = 1;
  
  private Thread m_splitThread = null;

  public TrainTestSplitMaker() {
         m_visual.loadIcons(BeanVisual.ICON_PATH
		       +"TrainTestSplitMaker.gif",
		       BeanVisual.ICON_PATH
		       +"TrainTestSplittMaker_animated.gif");
    m_visual.setText("TrainTestSplitMaker");
  }

  /**
   * Set a custom (descriptive) name for this bean
   * 
   * @param name the name to use
   */
  public void setCustomName(String name) {
    m_visual.setText(name);
  }

  /**
   * Get the custom (descriptive) name for this bean (if one has been set)
   * 
   * @return the custom name (or the default name)
   */
  public String getCustomName() {
    return m_visual.getText();
  }

  /**
   * Global info for this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "Split an incoming data set into separate train and test sets." ;
  }

  /**
   * Tip text info for this property
   *
   * @return a <code>String</code> value
   */
  public String trainPercentTipText() {
    return "The percentage of data to go into the training set";
  }

  /**
   * Set the percentage of data to be in the training portion of the split
   *
   * @param newTrainPercent an <code>int</code> value
   */
  public void setTrainPercent(double newTrainPercent) {
    m_trainPercentage = newTrainPercent;
  }

  /**
   * Get the percentage of the data that will be in the training portion of
   * the split
   *
   * @return an <code>int</code> value
   */
  public double getTrainPercent() {
    return m_trainPercentage;
  }

  /**
   * Tip text for this property
   *
   * @return a <code>String</code> value
   */
  public String seedTipText() {
    return "The randomization seed";
  }

  /**
   * Set the random seed
   *
   * @param newSeed an <code>int</code> value
   */
  public void setSeed(int newSeed) {
    m_randomSeed = newSeed;
  }

  /**
   * Get the value of the random seed
   *
   * @return an <code>int</code> value
   */
  public int getSeed() {
    return m_randomSeed;
  }

  /**
   * Accept a training set
   *
   * @param e a <code>TrainingSetEvent</code> value
   */
  public void acceptTrainingSet(TrainingSetEvent e) {
    Instances trainingSet = e.getTrainingSet();
    DataSetEvent dse = new DataSetEvent(this, trainingSet);
    acceptDataSet(dse);
  }

  /**
   * Accept a test set
   *
   * @param e a <code>TestSetEvent</code> value
   */
  public void acceptTestSet(TestSetEvent e) {
    Instances testSet = e.getTestSet();
    DataSetEvent dse = new DataSetEvent(this, testSet);
    acceptDataSet(dse);
  }

  /**
   * Accept a data set
   *
   * @param e a <code>DataSetEvent</code> value
   */
  public void acceptDataSet(DataSetEvent e) {
    if (m_splitThread == null) {
      final Instances dataSet = new Instances(e.getDataSet());
      m_splitThread = new Thread() {
	  public void run() {
	    try {
	      dataSet.randomize(new Random(m_randomSeed));
	      int trainSize = 
                (int)Math.round(dataSet.numInstances() * m_trainPercentage / 100);
	      int testSize = dataSet.numInstances() - trainSize;
      
	      Instances train = new Instances(dataSet, 0, trainSize);
	      Instances test = new Instances(dataSet, trainSize, testSize);
      
	      TrainingSetEvent tse =
		new TrainingSetEvent(TrainTestSplitMaker.this, train);
	      tse.m_setNumber = 1; tse.m_maxSetNumber = 1;
	      if (m_splitThread != null) {
		notifyTrainingSetProduced(tse);
	      }
    
	      // inform all test set listeners
	      TestSetEvent teste = 
		new TestSetEvent(TrainTestSplitMaker.this, test);
	      teste.m_setNumber = 1; teste.m_maxSetNumber = 1;
	      if (m_splitThread != null) {
		notifyTestSetProduced(teste);
	      } else {
		if (m_logger != null) {
		  m_logger.logMessage("[TrainTestSplitMaker] "
		      + statusMessagePrefix() + " Split has been canceled!");
		  m_logger.statusMessage(statusMessagePrefix()
		      + "INTERRUPTED");
		}
	      }
	    } catch (Exception ex) {
	      stop(); // stop all processing
	      if (m_logger != null) {
	          m_logger.statusMessage(statusMessagePrefix()
	              + "ERROR (See log for details)");
	          m_logger.logMessage("[TrainTestSplitMaker] " 
	              + statusMessagePrefix()
	              + " problem during split creation. " 
	              + ex.getMessage());
	      }
	      ex.printStackTrace();
	    } finally {
	      if (isInterrupted()) {
	        if (m_logger != null) {
	          m_logger.logMessage("[TrainTestSplitMaker] "
	              + statusMessagePrefix() + " Split has been canceled!");
	          m_logger.statusMessage(statusMessagePrefix()
	              + "INTERRUPTED");
	        }
	      }
	      block(false);
	    }
	  }
	};
      m_splitThread.setPriority(Thread.MIN_PRIORITY);
      m_splitThread.start();

      //      if (m_splitThread.isAlive()) {
      block(true);
      //      }
      m_splitThread = null;
    }
  }

  /**
   * Notify test set listeners that a test set is available
   *
   * @param tse a <code>TestSetEvent</code> value
   */
  protected void notifyTestSetProduced(TestSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_testListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
        if (m_splitThread == null) {
          break;
        }
        //	System.err.println("Notifying test listeners "
        //			   +"(Train - test split maker)");
	((TestSetListener)l.elementAt(i)).acceptTestSet(tse);
      }
    }
  }

  /**
   * Notify training set listeners that a training set is available
   *
   * @param tse a <code>TrainingSetEvent</code> value
   */
  protected void notifyTrainingSetProduced(TrainingSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_trainingListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
        if (m_splitThread == null) {
          break;
        }
        //	System.err.println("Notifying training listeners "
        //			   +"(Train - test split fold maker)");
	((TrainingSetListener)l.elementAt(i)).acceptTrainingSet(tse);
      }
    }
  }

  /**
   * Function used to stop code that calls acceptDataSet. This is 
   * needed as split is performed inside a separate
   * thread of execution.
   *
   * @param tf a <code>boolean</code> value
   */
  private synchronized void block(boolean tf) {
    if (tf) {
      try {
	// make sure that the thread is still alive before blocking
	if (m_splitThread.isAlive()) {
	  wait();
	}
      } catch (InterruptedException ex) {
      }
    } else {
      notifyAll();
    }
  }

  /**
   * Stop processing
   */
  public void stop() {
    // tell the listenee (upstream bean) to stop
    if (m_listenee instanceof BeanCommon) {
      //      System.err.println("Listener is BeanCommon");
      ((BeanCommon)m_listenee).stop();
    }

    // stop the split thread
    if (m_splitThread != null) {
      Thread temp = m_splitThread;
      m_splitThread = null;
      temp.interrupt();
      temp.stop();
    }
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return (m_splitThread != null);
  }

  /**
   * Get list of user requests
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_splitThread != null) {
      newVector.addElement("Stop");
    }
    return newVector.elements();
  }

  /**
   * Perform the named request
   *
   * @param request a <code>String</code> value
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Stop") == 0) {
      stop();
    } else {
      throw new IllegalArgumentException(request
			 + " not supported (TrainTestSplitMaker)");
    }
  }

  /**
   * Returns true, if at the current time, the named event could
   * be generated. Assumes that the supplied event name is
   * an event that could be generated by this bean
   *
   * @param eventName the name of the event in question
   * @return true if the named event could be generated at this point in
   * time
   */
  public boolean eventGeneratable(String eventName) {
    if (m_listenee == null) {
      return false;
    }
    
    if (m_listenee instanceof EventConstraints) {
      if (((EventConstraints)m_listenee).eventGeneratable("dataSet") ||
	  ((EventConstraints)m_listenee).eventGeneratable("trainingSet") ||
	  ((EventConstraints)m_listenee).eventGeneratable("testSet")) {
	return true;
      } else {
	return false;
      }
    }
    return true;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|";
  }
}

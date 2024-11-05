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
 *    CrossValidationFoldMaker.java
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
 * Bean for splitting instances into training ant test sets according to
 * a cross validation
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 6003 $
 */
public class CrossValidationFoldMaker 
  extends AbstractTrainAndTestSetProducer
  implements DataSourceListener, TrainingSetListener, TestSetListener, 
	     UserRequestAcceptor, EventConstraints, Serializable {

  /** for serialization */
  private static final long serialVersionUID = -6350179298851891512L;

  private int m_numFolds = 10;
  private int m_randomSeed = 1;
  
  private boolean m_preserveOrder = false;

  private transient Thread m_foldThread = null;

  public CrossValidationFoldMaker() {
    m_visual.loadIcons(BeanVisual.ICON_PATH
		       +"CrossValidationFoldMaker.gif",
		       BeanVisual.ICON_PATH
		       +"CrossValidationFoldMaker_animated.gif");
    m_visual.setText("CrossValidationFoldMaker");
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
    return "Split an incoming data set into cross validation folds. "
      +"Separate train and test sets are produced for each of the k folds.";
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
    if (e.isStructureOnly()) {
      // Pass on structure to training and test set listeners
      TrainingSetEvent tse = new TrainingSetEvent(this, e.getDataSet());
      TestSetEvent tsee = new TestSetEvent(this, e.getDataSet());
      notifyTrainingSetProduced(tse);
      notifyTestSetProduced(tsee);
      return;
    }
    if (m_foldThread == null) {
      final Instances dataSet = new Instances(e.getDataSet());
      m_foldThread = new Thread() {
	  public void run() {
	    boolean errorOccurred = false;
	    try {
	      Random random = new Random(getSeed());
	      if (!m_preserveOrder) {
	        dataSet.randomize(random);
	      }
	      if (dataSet.classIndex() >= 0 && 
		  (dataSet.attribute(dataSet.classIndex()).isNominal()|| dataSet.attribute(dataSet.classIndex()).isNominal()) &&
		  !m_preserveOrder) {
		dataSet.stratify(getFolds());
		if (m_logger != null) {
		  m_logger.logMessage("[" + getCustomName() + "] "
				      +"stratifying data");
		}
	      }
	      
	      for (int i = 0; i < getFolds(); i++) {
		if (m_foldThread == null) {
		  if (m_logger != null) {
		    m_logger.logMessage("[" + getCustomName() + "] Cross validation has been canceled!");
		  }
		  // exit gracefully
		  break;
		}
		Instances train = (!m_preserveOrder) 
		  ? dataSet.trainCV(getFolds(), i, random)
		  : dataSet.trainCV(getFolds(), i); 
		Instances test  = dataSet.testCV(getFolds(), i);

		// inform all training set listeners
		TrainingSetEvent tse = new TrainingSetEvent(this, train);
		tse.m_setNumber = i+1; tse.m_maxSetNumber = getFolds();
		String msg = getCustomName() + "$" 
		  + CrossValidationFoldMaker.this.hashCode() + "|";
		if (m_logger != null) {
		  m_logger.statusMessage(msg + "seed: " + getSeed() + " folds: "
		      + getFolds() + "|Training fold " + (i+1));
		}
		if (m_foldThread != null) {
		  //		  System.err.println("--Just before notify training set");
		  notifyTrainingSetProduced(tse);
		  //		  System.err.println("---Just after notify");
		}
	      
		// inform all test set listeners
		TestSetEvent teste = new TestSetEvent(this, test);
		teste.m_setNumber = i+1; teste.m_maxSetNumber = getFolds();
		
		if (m_logger != null) {
		  m_logger.statusMessage(msg + "seed: " + getSeed() + " folds: "
		      + getFolds() + "|Test fold " + (i+1));
		}
		if (m_foldThread != null) {
		  notifyTestSetProduced(teste);
		}
	      }
	    } catch (Exception ex) {
	      // stop all processing
	      errorOccurred = true;
	      if (m_logger != null) {
	        m_logger.logMessage("[" + getCustomName() 
	            + "] problem during fold creation. "
	            + ex.getMessage());
	      }
	      ex.printStackTrace();
	      CrossValidationFoldMaker.this.stop();
	    } finally {
	      m_foldThread = null;
	      
	      if (errorOccurred) {
	        if (m_logger != null) {
	          m_logger.statusMessage(getCustomName() 
	              + "$" + CrossValidationFoldMaker.this.hashCode()
	              + "|"
	              + "ERROR (See log for details).");
	        }
	      } else if (isInterrupted()) {
	        String msg = "[" + getCustomName() + "] Cross validation interrupted";
	        if (m_logger != null) {
	          m_logger.logMessage("[" + getCustomName() + "] Cross validation interrupted");
	          m_logger.statusMessage(getCustomName() + "$"
	              + CrossValidationFoldMaker.this.hashCode() + "|"
	              + "INTERRUPTED");
	        } else {
	          System.err.println(msg);
	        }
	      } else {
	        String msg = getCustomName() + "$" 
	        + CrossValidationFoldMaker.this.hashCode() + "|";
	        if (m_logger != null) {
	          m_logger.statusMessage(msg + "Finished.");
	        }
	      }
	      block(false);
	    }
	  }
	};
      m_foldThread.setPriority(Thread.MIN_PRIORITY);
      m_foldThread.start();

      //      if (m_foldThread.isAlive()) {
      block(true);
	//      }
      m_foldThread = null;
    }
  }


  /**
   * Notify all test set listeners of a TestSet event
   *
   * @param tse a <code>TestSetEvent</code> value
   */
  private void notifyTestSetProduced(TestSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_testListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
        if (m_foldThread == null) {
          break;
        }
	//	System.err.println("Notifying test listeners "
	//			   +"(cross validation fold maker)");
	((TestSetListener)l.elementAt(i)).acceptTestSet(tse);
      }
    }
  }

  /**
   * Notify all listeners of a TrainingSet event
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
        if (m_foldThread == null) {
          break;
        }
	//	System.err.println("Notifying training listeners "
	//			   +"(cross validation fold maker)");
	((TrainingSetListener)l.elementAt(i)).acceptTrainingSet(tse);
      }
    }
  }

  /**
   * Set the number of folds for the cross validation
   *
   * @param numFolds an <code>int</code> value
   */
  public void setFolds(int numFolds) {
    m_numFolds = numFolds;
  }
  
  /**
   * Get the currently set number of folds
   *
   * @return an <code>int</code> value
   */
  public int getFolds() {
    return m_numFolds;
  }

  /**
   * Tip text for this property
   *
   * @return a <code>String</code> value
   */
  public String foldsTipText() {
    return "The number of train and test splits to produce";
  }
    
  /**
   * Set the seed
   *
   * @param randomSeed an <code>int</code> value
   */
  public void setSeed(int randomSeed) {
    m_randomSeed = randomSeed;
  }
  
  /**
   * Get the currently set seed
   *
   * @return an <code>int</code> value
   */
  public int getSeed() {
    return m_randomSeed;
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
   * Returns true if the order of the incoming instances is to
   * be preserved under cross-validation (no randomization or 
   * stratification is done in this case).
   * 
   * @return true if the order of the incoming instances is to
   * be preserved.
   */
  public boolean getPreserveOrder() {
    return m_preserveOrder;
  }
  
  /**
   * Sets whether the order of the incoming instances is to be
   * preserved under cross-validation (no randomization or 
   * stratification is done in this case).
   *  
   * @param p true if the order is to be preserved.
   */
  public void setPreserveOrder(boolean p) {
    m_preserveOrder = p;
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return (m_foldThread != null);
  }

  /**
   * Stop any action
   */
  public void stop() {
    // tell the listenee (upstream bean) to stop
    if (m_listenee instanceof BeanCommon) {
      //      System.err.println("Listener is BeanCommon");
      ((BeanCommon)m_listenee).stop();
    }

    // stop the fold thread
    if (m_foldThread != null) {
      Thread temp = m_foldThread;
      m_foldThread = null;
      temp.interrupt();
      temp.stop();
    }
  }

  /**
   * Function used to stop code that calls acceptDataSet. This is 
   * needed as cross validation is performed inside a separate
   * thread of execution.
   *
   * @param tf a <code>boolean</code> value
   */
  private synchronized void block(boolean tf) {
    if (tf) {
      try {
	// make sure the thread is still running before we block
	if (m_foldThread != null && m_foldThread.isAlive()) {
	  wait();
	}
      } catch (InterruptedException ex) {
      }
    } else {
      notifyAll();
    }
  }

  /**
   * Return an enumeration of user requests
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_foldThread != null) {
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
					 + " not supported (CrossValidation)");
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
}

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
 *    ClassifierPerformanceEvaluator.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.gui.explorer.ClassifierErrorsPlotInstances;
import weka.gui.explorer.ExplorerDefaults;
import weka.gui.visualize.PlotData2D;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

/**
 * A bean that evaluates the performance of batch trained classifiers
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 6804 $
 */
public class ClassifierPerformanceEvaluator 
  extends AbstractEvaluator
  implements BatchClassifierListener, 
	     Serializable, UserRequestAcceptor, EventConstraints {

  /** for serialization */
  private static final long serialVersionUID = -3511801418192148690L;

  /**
   * Evaluation object used for evaluating a classifier
   */
  private transient Evaluation m_eval;

  private transient Thread m_evaluateThread = null;
  
  private transient long m_currentBatchIdentifier;
  private transient int m_setsComplete;
  
  private Vector m_textListeners = new Vector();
  private Vector m_thresholdListeners = new Vector();
  private Vector m_visualizableErrorListeners = new Vector();

  public ClassifierPerformanceEvaluator() {
    m_visual.loadIcons(BeanVisual.ICON_PATH
		       +"ClassifierPerformanceEvaluator.gif",
		       BeanVisual.ICON_PATH
		       +"ClassifierPerformanceEvaluator_animated.gif");
    m_visual.setText("ClassifierPerformanceEvaluator");
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
    return "Evaluate the performance of batch trained classifiers.";
  }

  // ----- Stuff for ROC curves
  private boolean m_rocListenersConnected = false;
  
  /** for generating plottable instance with predictions appended. */
  private transient ClassifierErrorsPlotInstances m_PlotInstances = null;
  
  protected static Evaluation adjustForInputMappedClassifier(Evaluation eval, 
      weka.classifiers.Classifier classifier,
      Instances inst, ClassifierErrorsPlotInstances plotInstances) throws Exception {

    if (classifier instanceof weka.classifiers.misc.InputMappedClassifier) {
      Instances mappedClassifierHeader = 
        ((weka.classifiers.misc.InputMappedClassifier)classifier).
        getModelHeader(new Instances(inst, 0));

      eval = new Evaluation(new Instances(mappedClassifierHeader, 0));

      if (!eval.getHeader().equalHeaders(inst)) {
        // When the InputMappedClassifier is loading a model, 
        // we need to make a new dataset that maps the test instances to
        // the structure expected by the mapped classifier - this is only
        // to ensure that the ClassifierPlotInstances object is configured
        // in accordance with what the embeded classifier was trained with
        Instances mappedClassifierDataset = 
          ((weka.classifiers.misc.InputMappedClassifier)classifier).
          getModelHeader(new Instances(mappedClassifierHeader, 0));
        for (int zz = 0; zz < inst.numInstances(); zz++) {
          Instance mapped = ((weka.classifiers.misc.InputMappedClassifier)classifier).
          constructMappedInstance(inst.instance(zz));
          mappedClassifierDataset.add(mapped);
        }
        
        eval.setPriors(mappedClassifierDataset);
        plotInstances.setInstances(mappedClassifierDataset);
        plotInstances.setClassifier(classifier);
        plotInstances.setClassIndex(mappedClassifierDataset.classIndex());
        plotInstances.setEvaluation(eval);
      }
    }
    
    return eval;
  }

  /**
   * Accept a classifier to be evaluated.
   *
   * @param ce a <code>BatchClassifierEvent</code> value
   */
  public void acceptClassifier(final BatchClassifierEvent ce) {
    if (ce.getTestSet() == null || ce.getTestSet().isStructureOnly()) {
      return; // cant evaluate empty/non-existent test instances
    }
    try {
      if (m_evaluateThread == null) {
	m_evaluateThread = new Thread() {
	    public void run() {
	      boolean errorOccurred = false;
//	      final String oldText = m_visual.getText();
	      Classifier classifier = ce.getClassifier();
	      try {
		// if (ce.getSetNumber() == 1) {
	        if (ce.getGroupIdentifier() != m_currentBatchIdentifier) {
		  
		  if (ce.getTrainSet().getDataSet() == null ||
		      ce.getTrainSet().getDataSet().numInstances() == 0) {
		    // we have no training set to estimate majority class
		    // or mean of target from
		    m_eval = new Evaluation(ce.getTestSet().getDataSet());
		    m_PlotInstances = ExplorerDefaults.getClassifierErrorsPlotInstances();
		    m_PlotInstances.setInstances(ce.getTestSet().getDataSet());
		    m_PlotInstances.setClassifier(ce.getClassifier());
		    m_PlotInstances.setClassIndex(ce.getTestSet().getDataSet().classIndex());
		    m_PlotInstances.setEvaluation(m_eval);

		    m_eval = adjustForInputMappedClassifier(m_eval, ce.getClassifier(),
		        ce.getTestSet().getDataSet(), m_PlotInstances);
		    m_eval.useNoPriors();
		  } else {
		    // we can set up with the training set here
		    m_eval = new Evaluation(ce.getTrainSet().getDataSet());
		    m_PlotInstances = ExplorerDefaults.getClassifierErrorsPlotInstances();
		    m_PlotInstances.setInstances(ce.getTrainSet().getDataSet());
		    m_PlotInstances.setClassifier(ce.getClassifier());
		    m_PlotInstances.setClassIndex(ce.getTestSet().getDataSet().classIndex());
		    m_PlotInstances.setEvaluation(m_eval);
		    
		    m_eval = adjustForInputMappedClassifier(m_eval, ce.getClassifier(),
                        ce.getTrainSet().getDataSet(), m_PlotInstances);
		  }
//		  m_classifier = ce.getClassifier();

		  m_PlotInstances.setUp();
		  
		  m_currentBatchIdentifier = ce.getGroupIdentifier();
		  m_setsComplete = 0;
		}
//		if (ce.getSetNumber() <= ce.getMaxSetNumber()) {
	        if (m_setsComplete < ce.getMaxSetNumber()) {
		  
		  /*if (ce.getTrainSet().getDataSet() != null &&
		      ce.getTrainSet().getDataSet().numInstances() > 0) {
		    // set the priors
		    m_eval.setPriors(ce.getTrainSet().getDataSet());
		  } */
		  
//		  m_visual.setText("Evaluating ("+ce.getSetNumber()+")...");
		  if (m_logger != null) {
		    m_logger.statusMessage(statusMessagePrefix()
					   +"Evaluating ("+ce.getSetNumber()
					   +")...");
		  }
		  m_visual.setAnimated();
		  /*
		  m_eval.evaluateModel(ce.getClassifier(), 
		  ce.getTestSet().getDataSet()); */
		  for (int i = 0; i < ce.getTestSet().getDataSet().numInstances(); i++) {
		    Instance temp = ce.getTestSet().getDataSet().instance(i);
		    m_PlotInstances.process(temp, ce.getClassifier(), m_eval);
		  }
		  
		  m_setsComplete++;
		}
		
//		if (ce.getSetNumber() == ce.getMaxSetNumber()) {
	        if (m_setsComplete == ce.getMaxSetNumber()) {
                  //		  System.err.println(m_eval.toSummaryString());
		  // m_resultsString.append(m_eval.toSummaryString());
		  // m_outText.setText(m_resultsString.toString());
		  String textTitle = classifier.getClass().getName();
		  String textOptions = "";
		  if (classifier instanceof OptionHandler) {
	             textOptions = 
	               Utils.joinOptions(((OptionHandler)classifier).getOptions()); 
		  }
		  textTitle = 
		    textTitle.substring(textTitle.lastIndexOf('.')+1,
					textTitle.length());
		  String resultT = "=== Evaluation result ===\n\n"
		    + "Scheme: " + textTitle + "\n"
		    + ((textOptions.length() > 0) ? "Options: " + textOptions + "\n": "")
		    + "Relation: " + ce.getTestSet().getDataSet().relationName()
		    + "\n\n" + m_eval.toSummaryString();
                  
                  if (ce.getTestSet().getDataSet().
                      classAttribute().isNominal() || ce.getTestSet().getDataSet().classAttribute().isRanking()) {
                    resultT += "\n" + m_eval.toClassDetailsString()
                      + "\n" + m_eval.toMatrixString();
                  }
                  
		  TextEvent te = 
		    new TextEvent(ClassifierPerformanceEvaluator.this, 
				  resultT,
				  textTitle);
		  notifyTextListeners(te);

                  // set up visualizable errors
                  if (m_visualizableErrorListeners.size() > 0) {
                    PlotData2D errorD = m_PlotInstances.getPlotData(
                	textTitle + " " + textOptions);
                    VisualizableErrorEvent vel = 
                      new VisualizableErrorEvent(ClassifierPerformanceEvaluator.this, errorD);
                    notifyVisualizableErrorListeners(vel);
                    m_PlotInstances.cleanUp();
                  }
                  

		  if ((ce.getTestSet().getDataSet().classAttribute().isNominal() || ce.getTestSet().getDataSet().classAttribute().isRanking()) &&
		      m_thresholdListeners.size() > 0) {
		    ThresholdCurve tc = new ThresholdCurve();
		    Instances result = tc.getCurve(m_eval.predictions(), 0);
		    result.
		      setRelationName(ce.getTestSet().getDataSet().relationName());
		    PlotData2D pd = new PlotData2D(result);
		    String htmlTitle = "<html><font size=-2>"
		      + textTitle;
		    String newOptions = "";
		    if (classifier instanceof OptionHandler) {
		      String[] options = 
		        ((OptionHandler) classifier).getOptions();
		      if (options.length > 0) {
		        for (int ii = 0; ii < options.length; ii++) {
		          if (options[ii].length() == 0) {
		            continue;
		          }
		          if (options[ii].charAt(0) == '-' && 
		              !(options[ii].charAt(1) >= '0' &&
		                  options[ii].charAt(1)<= '9')) {
		            newOptions += "<br>";
		          }
		          newOptions += options[ii];
		        }
		      }
		    }
		    
		   htmlTitle += " " + newOptions + "<br>" 
		      + " (class: "
                      +ce.getTestSet().getDataSet().
                        classAttribute().value(0) + ")" 
                      + "</font></html>";
		    pd.setPlotName(textTitle + " (class: "
	                      +ce.getTestSet().getDataSet().
	                        classAttribute().value(0) + ")");
		    pd.setPlotNameHTML(htmlTitle);
		    boolean [] connectPoints = 
		      new boolean [result.numInstances()];
		    for (int jj = 1; jj < connectPoints.length; jj++) {
		      connectPoints[jj] = true;
		    }
		    pd.setConnectPoints(connectPoints);
		    ThresholdDataEvent rde = 
		      new ThresholdDataEvent(ClassifierPerformanceEvaluator.this,
				       pd, ce.getTestSet().getDataSet().classAttribute());
		    notifyThresholdListeners(rde);
		    /*te = new TextEvent(ClassifierPerformanceEvaluator.this,
				       result.toString(),
				       "ThresholdCurveInst");
				       notifyTextListeners(te); */
		  }
		  if (m_logger != null) {
		    m_logger.statusMessage(statusMessagePrefix() + "Finished.");
		  }

		  // save memory
		  m_PlotInstances = null;
		}
	      } catch (Exception ex) {
	        errorOccurred = true;
	        ClassifierPerformanceEvaluator.this.stop(); // stop all processing
	        if (m_logger != null) {
	          m_logger.logMessage("[ClassifierPerformanceEvaluator] "
	              + statusMessagePrefix() 
	              + " problem evaluating classifier. " 
	              + ex.getMessage());
	        }
		ex.printStackTrace();
	      } finally {
//		m_visual.setText(oldText);
		m_visual.setStatic();
		m_evaluateThread = null;
						
		if (m_logger != null) {
		  if (errorOccurred) {
		    m_logger.statusMessage(statusMessagePrefix() 
		        + "ERROR (See log for details)");
		  } else if (isInterrupted()) {
		    m_logger.logMessage("[" + getCustomName() +"] Evaluation interrupted!");
		    m_logger.statusMessage(statusMessagePrefix() 
		        + "INTERRUPTED");
		  }
		}
		block(false);
	      }
	    }
	  };
	m_evaluateThread.setPriority(Thread.MIN_PRIORITY);
	m_evaluateThread.start();

	// make sure the thread is still running before we block
	//	if (m_evaluateThread.isAlive()) {
	block(true);
	  //	}
	m_evaluateThread = null;
      }
    }  catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return (m_evaluateThread != null);
  }
    
  /**
   * Try and stop any action
   */
  public void stop() {
    // tell the listenee (upstream bean) to stop
    if (m_listenee instanceof BeanCommon) {
      //      System.err.println("Listener is BeanCommon");
      ((BeanCommon)m_listenee).stop();
    }

    // stop the evaluate thread
    if (m_evaluateThread != null) {
      m_evaluateThread.interrupt();
      m_evaluateThread.stop();
      m_evaluateThread = null;
      m_visual.setStatic();
    }
  }
  
  /**
   * Function used to stop code that calls acceptClassifier. This is 
   * needed as classifier evaluation is performed inside a separate
   * thread of execution.
   *
   * @param tf a <code>boolean</code> value
   */
  private synchronized void block(boolean tf) {
    if (tf) {
      try {
	// only block if thread is still doing something useful!
	if (m_evaluateThread != null && m_evaluateThread.isAlive()) {
	  wait();
	}
      } catch (InterruptedException ex) {
      }
    } else {
      notifyAll();
    }
  }

  /**
   * Return an enumeration of user activated requests for this bean
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_evaluateThread != null) {
      newVector.addElement("Stop");
    }
    return newVector.elements();
  }

  /**
   * Perform the named request
   *
   * @param request the request to perform
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Stop") == 0) {
      stop();
    } else {
      throw new 
	IllegalArgumentException(request

		    + " not supported (ClassifierPerformanceEvaluator)");
    }
  }

  /**
   * Add a text listener
   *
   * @param cl a <code>TextListener</code> value
   */
  public synchronized void addTextListener(TextListener cl) {
    m_textListeners.addElement(cl);
  }

  /**
   * Remove a text listener
   *
   * @param cl a <code>TextListener</code> value
   */
  public synchronized void removeTextListener(TextListener cl) {
    m_textListeners.remove(cl);
  }
  
  /**
   * Add a threshold data listener
   *
   * @param cl a <code>ThresholdDataListener</code> value
   */
  public synchronized void addThresholdDataListener(ThresholdDataListener cl) {
    m_thresholdListeners.addElement(cl);
  }

  /**
   * Remove a Threshold data listener
   *
   * @param cl a <code>ThresholdDataListener</code> value
   */
  public synchronized void removeThresholdDataListener(ThresholdDataListener cl) {
    m_thresholdListeners.remove(cl);
  }

  /**
   * Add a visualizable error listener
   *
   * @param vel a <code>VisualizableErrorListener</code> value
   */
  public synchronized void addVisualizableErrorListener(VisualizableErrorListener vel) {
    m_visualizableErrorListeners.add(vel);
  }

  /**
   * Remove a visualizable error listener
   *
   * @param vel a <code>VisualizableErrorListener</code> value
   */
  public synchronized void removeVisualizableErrorListener(VisualizableErrorListener vel) {
    m_visualizableErrorListeners.remove(vel);
  }

  /**
   * Notify all text listeners of a TextEvent
   *
   * @param te a <code>TextEvent</code> value
   */
  private void notifyTextListeners(TextEvent te) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_textListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	//	System.err.println("Notifying text listeners "
	//			   +"(ClassifierPerformanceEvaluator)");
	((TextListener)l.elementAt(i)).acceptText(te);
      }
    }
  }

  /**
   * Notify all ThresholdDataListeners of a ThresholdDataEvent
   *
   * @param te a <code>ThresholdDataEvent</code> value
   */
  private void notifyThresholdListeners(ThresholdDataEvent re) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_thresholdListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	//	System.err.println("Notifying text listeners "
	//			   +"(ClassifierPerformanceEvaluator)");
	((ThresholdDataListener)l.elementAt(i)).acceptDataSet(re);
      }
    }
  }

  /**
   * Notify all VisualizableErrorListeners of a VisualizableErrorEvent
   *
   * @param te a <code>VisualizableErrorEvent</code> value
   */
  private void notifyVisualizableErrorListeners(VisualizableErrorEvent re) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_visualizableErrorListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	//	System.err.println("Notifying text listeners "
	//			   +"(ClassifierPerformanceEvaluator)");
	((VisualizableErrorListener)l.elementAt(i)).acceptDataSet(re);
      }
    }
  }

  /**
   * Returns true, if at the current time, the named event could
   * be generated. Assumes that supplied event names are names of
   * events that could be generated by this bean.
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
      if (!((EventConstraints)m_listenee).
	  eventGeneratable("batchClassifier")) {
	return false;
      }
    }
    return true;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|";
  }
}


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
 *    Classifier.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileFilter;

import weka.classifiers.rules.ZeroR;
import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.xml.KOML;
import weka.core.xml.XStream;
import weka.experiment.Task;
import weka.experiment.TaskStatusInfo;
import weka.gui.ExtensionFileFilter;
import weka.gui.Logger;

/**
 * Bean that wraps around weka.classifiers
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 6804 $
 * @since 1.0
 * @see JPanel
 * @see BeanCommon
 * @see Visible
 * @see WekaWrapper
 * @see Serializable
 * @see UserRequestAcceptor
 * @see TrainingSetListener
 * @see TestSetListener
 * @see EnvironmentHandler
 */
public class Classifier
  extends JPanel
  implements BeanCommon, Visible, 
	     WekaWrapper, EventConstraints,
	     Serializable, UserRequestAcceptor,
	     TrainingSetListener, TestSetListener,
	     InstanceListener, ConfigurationProducer,
	     EnvironmentHandler {

  /** for serialization */
  private static final long serialVersionUID = 659603893917736008L;

  protected BeanVisual m_visual = 
    new BeanVisual("Classifier",
		   BeanVisual.ICON_PATH+"DefaultClassifier.gif",
		   BeanVisual.ICON_PATH+"DefaultClassifier_animated.gif");

  private static int IDLE = 0;
  private static int BUILDING_MODEL = 1;
  private static int CLASSIFYING = 2;

  private int m_state = IDLE;

  //private Thread m_buildThread = null;

  /**
   * Global info for the wrapped classifier (if it exists).
   */
  protected String m_globalInfo;

  /**
   * Objects talking to us
   */
  private Hashtable m_listenees = new Hashtable();

  /**
   * Objects listening for batch classifier events
   */
  private Vector m_batchClassifierListeners = new Vector();

  /**
   * Objects listening for incremental classifier events
   */
  private Vector m_incrementalClassifierListeners = new Vector();

  /**
   * Objects listening for graph events
   */
  private Vector m_graphListeners = new Vector();

  /**
   * Objects listening for text events
   */
  private Vector m_textListeners = new Vector();

  /**
   * Holds training instances for batch training. Not transient because
   * header is retained for validating any instance events that this
   * classifier might be asked to predict in the future.
   */
  private Instances m_trainingSet;
  private transient Instances m_testingSet;
  private weka.classifiers.Classifier m_Classifier = new ZeroR();
  /** Template used for creating copies when building in parallel */
  private weka.classifiers.Classifier m_ClassifierTemplate = m_Classifier;
  
  private IncrementalClassifierEvent m_ie = 
    new IncrementalClassifierEvent(this);

  /** the extension for serialized models (binary Java serialization) */
  public final static String FILE_EXTENSION = "model";

  private transient JFileChooser m_fileChooser = null; 

  protected FileFilter m_binaryFilter =
    new ExtensionFileFilter("."+FILE_EXTENSION, "Binary serialized model file (*"
                            + FILE_EXTENSION + ")");

  protected FileFilter m_KOMLFilter =
    new ExtensionFileFilter(KOML.FILE_EXTENSION + FILE_EXTENSION,
                            "XML serialized model file (*"
                            + KOML.FILE_EXTENSION + FILE_EXTENSION + ")");

  protected FileFilter m_XStreamFilter =
    new ExtensionFileFilter(XStream.FILE_EXTENSION + FILE_EXTENSION,
                            "XML serialized model file (*"
                            + XStream.FILE_EXTENSION + FILE_EXTENSION + ")");
  
  protected transient Environment m_env;

  /**
   * If the classifier is an incremental classifier, should we
   * update it (ie train it on incoming instances). This makes it
   * possible incrementally test on a separate stream of instances
   * without updating the classifier, or mix batch training/testing
   * with incremental training/testing
   */
  private boolean m_updateIncrementalClassifier = true;

  private transient Logger m_log = null;

  /**
   * Event to handle when processing incremental updates
   */
  private InstanceEvent m_incrementalEvent;
  
  /**
   * Number of threads to use to train models with
   */
  protected int m_executionSlots = 2;
  
//  protected int m_queueSize = 5;
  
  /**
   * Pool of threads to train models on incoming data 
   */
  protected transient ThreadPoolExecutor m_executorPool;  
  
  /**
   * Stores completed models and associated data sets.
   */
  protected transient BatchClassifierEvent[][] m_outputQueues;

  /**
   * Stores which sets from which runs have been completed.
   */
  protected transient boolean[][] m_completedSets;
  
  /**
   * Identifier for the current batch. A batch is a group
   * of related runs/sets.
   */
  protected transient Date m_currentBatchIdentifier;
  
  /**
   * Holds original icon label text
   */
  protected String m_oldText = "";
  
  /**
   * true if we should reject any further training 
   * data sets, until all processing has been finished,
   *  once we've received the last fold of
   * the last run.
   */
  protected boolean m_reject = false;
  
  /** 
   * True if we should block rather reject until
   * all processing has been completed.
   */
  protected boolean m_block = false;

  /**
   * Global info (if it exists) for the wrapped classifier
   *
   * @return the global info
   */
  public String globalInfo() {
    return m_globalInfo;
  }

  /**
   * Creates a new <code>Classifier</code> instance.
   */
  public Classifier() {
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
    setClassifierTemplate(m_ClassifierTemplate);
    
    //setupFileChooser();
  }
  
  private void startExecutorPool() {
    
    if (m_executorPool != null) {
      m_executorPool.shutdownNow();
    }
    
    m_executorPool = new ThreadPoolExecutor(m_executionSlots, m_executionSlots,
        120, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
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

  protected void setupFileChooser() {
    if (m_fileChooser == null) {
      m_fileChooser = 
        new JFileChooser(new File(System.getProperty("user.dir")));
    }

    m_fileChooser.addChoosableFileFilter(m_binaryFilter);
    if (KOML.isPresent()) {
      m_fileChooser.addChoosableFileFilter(m_KOMLFilter);
    }
    if (XStream.isPresent()) {
      m_fileChooser.addChoosableFileFilter(m_XStreamFilter);
    }
    m_fileChooser.setFileFilter(m_binaryFilter);
  }
  
  /**
   * Get the number of execution slots (threads) used
   * to train models.
   * 
   * @return the number of execution slots.
   */
  public int getExecutionSlots() {
    return m_executionSlots;
  }
  
  /**
   * Set the number of execution slots (threads) to use to
   * train models with.
   * 
   * @param slots the number of execution slots to use.
   */
  public void setExecutionSlots(int slots) {
    m_executionSlots = slots;
  }
  
  /**
   * Set whether to block on receiving the last fold
   * of the last run rather than rejecting any further
   * data until all processing is complete.
   * 
   * @param block true if we should block on the
   * last fold of the last run.
   */
  public void setBlockOnLastFold(boolean block) {
    m_block = block;
  }
  
  /**
   * Gets whether we are blocking on the last fold of the
   * last run rather than rejecting any further data until
   * all processing has been completed.
   * 
   * @return true if we are blocking on the last fold
   * of the last run
   */
  public boolean getBlockOnLastFold() {
    return m_block;
  }

  /**
   * Set the template classifier for this wrapper
   *
   * @param c a <code>weka.classifiers.Classifier</code> value
   */
  public void setClassifierTemplate(weka.classifiers.Classifier c) {
    boolean loadImages = true;
    if (c.getClass().getName().
	compareTo(m_ClassifierTemplate.getClass().getName()) == 0) {
      loadImages = false;
    } else {
      // classifier has changed so any batch training status is now
      // invalid
      m_trainingSet = null;
    }
    m_ClassifierTemplate = c;
    String classifierName = c.getClass().toString();
    classifierName = classifierName.substring(classifierName.
					      lastIndexOf('.')+1, 
					      classifierName.length());
    if (loadImages) {
      if (!m_visual.loadIcons(BeanVisual.ICON_PATH+classifierName+".gif",
		       BeanVisual.ICON_PATH+classifierName+"_animated.gif")) {
	useDefaultVisual();
      }
    }
    m_visual.setText(classifierName);

    if (!(m_ClassifierTemplate instanceof weka.classifiers.UpdateableClassifier) &&
	(m_listenees.containsKey("instance"))) {
      if (m_log != null) {
	m_log.logMessage("[Classifier] " + statusMessagePrefix() + " WARNING : "
	    + getCustomName() +" is not an incremental classifier");
      }
    }
    // get global info
    m_globalInfo = KnowledgeFlowApp.getGlobalInfo(m_ClassifierTemplate);
    
    try {
      m_Classifier = weka.classifiers.AbstractClassifier.makeCopy(m_ClassifierTemplate);
    } catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  /**
   * Return the classifier template currently in use.
   * 
   * @return the classifier template currently in use.
   */
  public weka.classifiers.Classifier getClassifierTemplate() {
    return m_ClassifierTemplate;
  }
  
  private void setTrainedClassifier(weka.classifiers.Classifier tc) {
    m_Classifier = tc;
    
    // set the template
    weka.classifiers.Classifier newTemplate = null;
    try {
      String[] options = ((OptionHandler)tc).getOptions();
      newTemplate = weka.classifiers.AbstractClassifier.forName(tc.getClass().getName(), options);
      setClassifierTemplate(newTemplate);
    } catch (Exception ex) {
      if (m_log != null) {
        m_log.logMessage("[Classifier] " + statusMessagePrefix() + ex.getMessage());
        String errorMessage = statusMessagePrefix()
        + "ERROR: see log for details.";
        m_log.statusMessage(errorMessage);        
      } else {
        ex.printStackTrace();
      }
    }    
  }

  /**
   * Returns true if this classifier has an incoming connection that is
   * an instance stream
   *
   * @return true if has an incoming connection that is an instance stream
   */
  public boolean hasIncomingStreamInstances() {
    if (m_listenees.size() == 0) {
      return false;
    }
    if (m_listenees.containsKey("instance")) {
      return true;
    }
    return false;
  }

  /**
   * Returns true if this classifier has an incoming connection that is
   * a batch set of instances
   *
   * @return a <code>boolean</code> value
   */
  public boolean hasIncomingBatchInstances() {
    if (m_listenees.size() == 0) {
      return false;
    }
    if (m_listenees.containsKey("trainingSet") ||
	m_listenees.containsKey("testSet")) {
      return true;
    }
    return false;
  }

  /**
   * Get the currently trained classifier.
   *
   * @return a <code>weka.classifiers.Classifier</code> value
   */
  public weka.classifiers.Classifier getClassifier() {
    return m_Classifier;
  }

  /**
   * Sets the algorithm (classifier) for this bean
   *
   * @param algorithm an <code>Object</code> value
   * @exception IllegalArgumentException if an error occurs
   */
  public void setWrappedAlgorithm(Object algorithm) 
    {

    if (!(algorithm instanceof weka.classifiers.Classifier)) { 
      throw new IllegalArgumentException(algorithm.getClass()+" : incorrect "
					 +"type of algorithm (Classifier)");
    }
    setClassifierTemplate((weka.classifiers.Classifier)algorithm);
  }

  /**
   * Returns the wrapped classifier
   *
   * @return an <code>Object</code> value
   */
  public Object getWrappedAlgorithm() {
    return getClassifierTemplate();
  }

  /**
   * Get whether an incremental classifier will be updated on the
   * incoming instance stream.
   * 
   * @return true if an incremental classifier is to be updated.
   */
  public boolean getUpdateIncrementalClassifier() {
    return m_updateIncrementalClassifier;
  }

  /**
   * Set whether an incremental classifier will be updated on the
   * incoming instance stream.
   * 
   * @param update true if an incremental classifier is to be updated.
   */
  public void setUpdateIncrementalClassifier(boolean update) {
    m_updateIncrementalClassifier = update;
  }

  /**
   * Accepts an instance for incremental processing.
   *
   * @param e an <code>InstanceEvent</code> value
   */
  public void acceptInstance(InstanceEvent e) {
    m_incrementalEvent = e;
    handleIncrementalEvent();
  }

  /**
   * Handles initializing and updating an incremental classifier
   */
  private void handleIncrementalEvent() {
    if (m_executorPool != null && 
        (m_executorPool.getQueue().size() > 0 || 
            m_executorPool.getActiveCount() > 0)) {
      
      String messg = "[Classifier] " + statusMessagePrefix() 
        + " is currently batch training!";
      if (m_log != null) {
	m_log.logMessage(messg);
	m_log.statusMessage(statusMessagePrefix() + "WARNING: "
	    + "Can't accept instance - batch training in progress.");
      } else {
	System.err.println(messg);
      }
      return;
    }

    if (m_incrementalEvent.getStatus() == InstanceEvent.FORMAT_AVAILABLE) {
      // clear any warnings/errors from the log
      if (m_log != null) {
        m_log.statusMessage(statusMessagePrefix() + "remove");
      }
      
      //      Instances dataset = m_incrementalEvent.getInstance().dataset();
      Instances dataset = m_incrementalEvent.getStructure();
      // default to the last column if no class is set
      if (dataset.classIndex() < 0) {
        stop();
        String errorMessage = statusMessagePrefix()
            + "ERROR: no class attribute set in incoming stream!";
        if (m_log != null) {
          m_log.statusMessage(errorMessage);
          m_log.logMessage("[" + getCustomName() + "] " + errorMessage);
        } else {
          System.err.println("[" + getCustomName() + "] " + errorMessage);
        }
        return;
        
	// System.err.println("Classifier : setting class index...");
	//dataset.setClassIndex(dataset.numAttributes()-1);
      }
      try {
	// initialize classifier if m_trainingSet is null
	// otherwise assume that classifier has been pre-trained in batch
	// mode, *if* headers match
	if (m_trainingSet == null || !m_trainingSet.equalHeaders(dataset)) {
	  if (!(m_ClassifierTemplate instanceof 
		weka.classifiers.UpdateableClassifier)) {
	    stop(); // stop all processing
	    if (m_log != null) {
	      String msg = (m_trainingSet == null)
		? statusMessagePrefix()
		+ "ERROR: classifier has not been batch "
		+"trained; can't process instance events."
		: statusMessagePrefix() 
		  + "ERROR: instance event's structure is different from "
		  +"the data that "
		  + "was used to batch train this classifier; can't continue.";
	      m_log.logMessage("[Classifier] " + msg);
	      m_log.statusMessage(msg);
	    }
	    return;
	  }
	  
	  if (m_trainingSet != null && 
	      (!dataset.equalHeaders(m_trainingSet))) {
	    if (m_log != null) {
	      String msg = statusMessagePrefix() 
              + " WARNING : structure of instance events differ "
              +"from data used in batch training this "
              +"classifier. Resetting classifier...";
	      m_log.logMessage("[Classifier] " + msg);
	      m_log.statusMessage(msg);
	    }
	    m_trainingSet = null;
	  }
	  if (m_trainingSet == null) {
	    // initialize the classifier if it hasn't been trained yet
	    m_trainingSet = new Instances(dataset, 0);
	    m_Classifier = weka.classifiers.AbstractClassifier.makeCopy(m_ClassifierTemplate);
	    if (m_Classifier instanceof EnvironmentHandler && m_env != null) {
	      ((EnvironmentHandler)m_Classifier).setEnvironment(m_env);
	    }
	    m_Classifier.buildClassifier(m_trainingSet);
	  }
	}
      } catch (Exception ex) {
        stop();
        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix()
              + "ERROR (See log for details)");
          m_log.logMessage("[Classifier] " + statusMessagePrefix()
              + " problem during incremental processing. " 
              + ex.getMessage());
        }
	ex.printStackTrace();
      }
      // Notify incremental classifier listeners of new batch
      System.err.println("NOTIFYING NEW BATCH");
      m_ie.setStructure(dataset); 
      m_ie.setClassifier(m_Classifier);

      notifyIncrementalClassifierListeners(m_ie);
      return;
    } else {
      if (m_trainingSet == null) {
	// simply return. If the training set is still null after
	// the first instance then the classifier must not be updateable
	// and hasn't been previously batch trained - therefore we can't
	// do anything meaningful
	return;
      }
    }

    try {
      // test on this instance
      if (m_incrementalEvent.getInstance().dataset().classIndex() < 0) {
        // System.err.println("Classifier : setting class index...");
        m_incrementalEvent.getInstance().dataset().setClassIndex(
            m_incrementalEvent.getInstance().dataset().numAttributes()-1);
      }
      
      int status = IncrementalClassifierEvent.WITHIN_BATCH;
      /*      if (m_incrementalEvent.getStatus() == InstanceEvent.FORMAT_AVAILABLE) {
	      status = IncrementalClassifierEvent.NEW_BATCH; */
      /* } else */ if (m_incrementalEvent.getStatus() ==
		       InstanceEvent.BATCH_FINISHED) {
	status = IncrementalClassifierEvent.BATCH_FINISHED;
      }

      m_ie.setStatus(status); m_ie.setClassifier(m_Classifier);
      m_ie.setCurrentInstance(m_incrementalEvent.getInstance());

      notifyIncrementalClassifierListeners(m_ie);

      // now update on this instance (if class is not missing and classifier
      // is updateable and user has specified that classifier is to be
      // updated)
      if (m_ClassifierTemplate instanceof weka.classifiers.UpdateableClassifier &&
	  m_updateIncrementalClassifier == true &&
	  !(m_incrementalEvent.getInstance().
	    isMissing(m_incrementalEvent.getInstance().
		      dataset().classIndex()))) {
	((weka.classifiers.UpdateableClassifier)m_Classifier).
	  updateClassifier(m_incrementalEvent.getInstance());
      }
      if (m_incrementalEvent.getStatus() == 
	  InstanceEvent.BATCH_FINISHED) {
	if (m_textListeners.size() > 0) {
	  String modelString = m_Classifier.toString();
	  String titleString = m_Classifier.getClass().getName();

	  titleString = titleString.
	    substring(titleString.lastIndexOf('.') + 1,
		      titleString.length());
	  modelString = "=== Classifier model ===\n\n" +
	    "Scheme:   " +titleString+"\n" +
	    "Relation: "  + m_trainingSet.relationName() + "\n\n"
	    + modelString;
	  titleString = "Model: " + titleString;
	  TextEvent nt = new TextEvent(this,
				       modelString,
				       titleString);
	  notifyTextListeners(nt);
	}
      }
    } catch (Exception ex) {
      stop();
      if (m_log != null) {
	m_log.logMessage("[Classifier] " + statusMessagePrefix()
	    + ex.getMessage());
	m_log.statusMessage(statusMessagePrefix()
	    + "ERROR (see log for details)");
	ex.printStackTrace();
      } else {
        ex.printStackTrace();
      }
    }
  }
  
  protected class TrainingTask implements Runnable, Task {
    private int m_runNum;
    private int m_maxRunNum;
    private int m_setNum;
    private int m_maxSetNum;
    private Instances m_train = null;
    private TaskStatusInfo m_taskInfo = new TaskStatusInfo();

    public TrainingTask(int runNum, int maxRunNum, 
        int setNum, int maxSetNum, Instances train) {
      m_runNum = runNum;
      m_maxRunNum = maxRunNum;
      m_setNum = setNum;
      m_maxSetNum = maxSetNum;
      m_train = train;
      m_taskInfo.setExecutionStatus(TaskStatusInfo.TO_BE_RUN);
    }

    public void run() {
      execute();
    }

    public void execute() {
      try {
        if (m_train != null) {
          if (m_train.classIndex() < 0) {
            // stop all processing
            stop();
            String errorMessage = statusMessagePrefix()
                + "ERROR: no class attribute set in test data!";
            if (m_log != null) {
              m_log.statusMessage(errorMessage);
              m_log.logMessage("[Classifier] " + errorMessage);
            } else {
              System.err.println("[Classifier] " + errorMessage);
            }
            return;
            
            // assume last column is the class
/*            m_train.setClassIndex(m_train.numAttributes()-1);
            if (m_log != null) {
              m_log.logMessage("[Classifier] " + statusMessagePrefix() 
                  + " : assuming last "
                  +"column is the class");
            } */
          }
          if (m_runNum == 1 && m_setNum == 1) {
            // set this back to idle once the last fold
            // of the last run has completed
            m_state = BUILDING_MODEL; // global state
            
            // local status of this runnable
            m_taskInfo.setExecutionStatus(TaskStatusInfo.PROCESSING);
          }
          
          //m_visual.setAnimated();
          //m_visual.setText("Building model...");
          String msg = statusMessagePrefix()
            + "Building model for run " + m_runNum + " fold " + m_setNum;
          if (m_log != null) {
            m_log.statusMessage(msg);
          } else {
            System.err.println(msg);
          }
          // buildClassifier();
          
          // copy the classifier configuration
          weka.classifiers.Classifier classifierCopy = 
            weka.classifiers.AbstractClassifier.makeCopy(m_ClassifierTemplate);
          if (classifierCopy instanceof EnvironmentHandler && m_env != null) {
            ((EnvironmentHandler)classifierCopy).setEnvironment(m_env);
          }
          
          // build this model
          classifierCopy.buildClassifier(m_train);
          if (m_runNum == m_maxRunNum && m_setNum == m_maxSetNum) {
            // Save the last classifier (might be used later on for
            // classifying further test sets.
            m_Classifier = classifierCopy;
            m_trainingSet = m_train;
          }
                              
          //if (m_batchClassifierListeners.size() > 0) {
            // notify anyone who might be interested in just the model
            // and training set.
            BatchClassifierEvent ce = 
              new BatchClassifierEvent(Classifier.this, classifierCopy, 
                  new DataSetEvent(this, m_train),
                  null, // no test set (yet)
                  m_setNum, m_maxSetNum);
            ce.setGroupIdentifier(m_currentBatchIdentifier.getTime());
            notifyBatchClassifierListeners(ce);
                        
            // store in the output queue (if we have incoming test set events)
            ce = 
              new BatchClassifierEvent(Classifier.this, classifierCopy, 
                  new DataSetEvent(this, m_train),
                  null, // no test set (yet)
                  m_setNum, m_maxSetNum);
            ce.setGroupIdentifier(m_currentBatchIdentifier.getTime());
            classifierTrainingComplete(ce);
          //}

          if (classifierCopy instanceof weka.core.Drawable && 
              m_graphListeners.size() > 0) {
            String grphString = 
              ((weka.core.Drawable)classifierCopy).graph();
            int grphType = ((weka.core.Drawable)classifierCopy).graphType();
            String grphTitle = classifierCopy.getClass().getName();
            grphTitle = grphTitle.substring(grphTitle.
                lastIndexOf('.')+1, 
                grphTitle.length());
            grphTitle = "Set " + m_setNum + " ("
            + m_train.relationName() + ") "
            + grphTitle;

            GraphEvent ge = new GraphEvent(Classifier.this, 
                grphString, 
                grphTitle,
                grphType);
            notifyGraphListeners(ge);
          }

          if (m_textListeners.size() > 0) {
            String modelString = classifierCopy.toString();
            String titleString = classifierCopy.getClass().getName();

            titleString = titleString.
            substring(titleString.lastIndexOf('.') + 1,
                titleString.length());
            modelString = "=== Classifier model ===\n\n" +
            "Scheme:   " +titleString+"\n" +
            "Relation: "  + m_train.relationName() + 
            ((m_maxSetNum > 1) 
                ? "\nTraining Fold: " + m_setNum
                    :"")
                    + "\n\n"
                    + modelString;
            titleString = "Model: " + titleString;

            TextEvent nt = new TextEvent(Classifier.this,
                modelString,
                titleString);
            notifyTextListeners(nt);
          }
        }
      } catch (Exception ex) {
        ex.printStackTrace();
        if (m_log != null) {
          String titleString = "[Classifier] " + statusMessagePrefix();

          titleString += " run " + m_runNum + " fold " + m_setNum
          + " failed to complete.";
          m_log.logMessage(titleString + " (build classifier). " 
              + ex.getMessage());
          m_log.statusMessage(statusMessagePrefix() 
              + "ERROR (see log for details)");
          ex.printStackTrace();
        }
        m_taskInfo.setExecutionStatus(TaskStatusInfo.FAILED);
        // Stop all processing
        stop();
      } finally {
        m_visual.setStatic();
        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix() + "Finished.");
        }
        m_state = IDLE;
        if (Thread.currentThread().isInterrupted()) {
          // prevent any classifier events from being fired
          m_trainingSet = null;
          if (m_log != null) {
            String titleString = "[Classifier] " + statusMessagePrefix();                 
         
            m_log.logMessage(titleString + " ("
               + " run " + m_runNum + " fold " + m_setNum + ") interrupted!");
            m_log.statusMessage(statusMessagePrefix() + "INTERRUPTED");
            
            /* // are we the last active thread?
            if (m_executorPool.getActiveCount() == 1) {
              String msg = "[Classifier] " + statusMessagePrefix() 
              + " last classifier unblocking...";
              System.err.println(msg + " (interrupted)");
              m_log.logMessage(msg + " (interrupted)");
//              m_log.statusMessage(statusMessagePrefix() + "finished.");
              m_block = false;
              m_state = IDLE;
              block(false);
            } */
          }
          /*System.err.println("Queue size: " + m_executorPool.getQueue().size() +
              " Active count: " + m_executorPool.getActiveCount()); */
        } /* else {
          // check to see if we are the last active thread
          if (m_executorPool == null || 
              (m_executorPool.getQueue().size() == 0 && 
                  m_executorPool.getActiveCount() == 1)) {

            String msg = "[Classifier] " + statusMessagePrefix() 
            + " last classifier unblocking...";
            System.err.println(msg);
            if (m_log != null) {
              m_log.logMessage(msg);
            } else {
              System.err.println(msg);
            }
            //m_visual.setText(m_oldText);

            if (m_log != null) {
              m_log.statusMessage(statusMessagePrefix() + "Finished.");
            }
            // m_outputQueues = null; // free memory
            m_block = false;
            block(false);
          }
        } */
      }
    }
  
    public TaskStatusInfo getTaskStatus() {
      // TODO
      return null;
    }     
  }     
  
  /**
   * Accepts a training set and builds batch classifier
   *
   * @param e a <code>TrainingSetEvent</code> value
   */
  public void acceptTrainingSet(final TrainingSetEvent e) {
    
    if (e.isStructureOnly()) {
      // no need to build a classifier, instead just generate a dummy
      // BatchClassifierEvent in order to pass on instance structure to
      // any listeners - eg. PredictionAppender can use it to determine
      // the final structure of instances with predictions appended
      BatchClassifierEvent ce = 
        new BatchClassifierEvent(this, m_Classifier, 
                                 new DataSetEvent(this, e.getTrainingSet()),
                                 new DataSetEvent(this, e.getTrainingSet()),
                                 e.getSetNumber(), e.getMaxSetNumber());

      notifyBatchClassifierListeners(ce);
      return;
    }
    
    if (m_reject) {
      //block(true);
      if (m_log != null) {
        m_log.statusMessage(statusMessagePrefix() + "BUSY. Can't accept data "
            + "at this time.");
        m_log.logMessage("[Classifier] " + statusMessagePrefix()
            + " BUSY. Can't accept data at this time.");
      }
      return;
    }
    
    // Do some initialization if this is the first set of the first run
    if (e.getRunNumber() == 1 && e.getSetNumber() == 1) {
//      m_oldText = m_visual.getText();
      // store the training header
      m_trainingSet = new Instances(e.getTrainingSet(), 0);
      m_state = BUILDING_MODEL;
      
      String msg = "[Classifier] " + statusMessagePrefix() 
        + " starting executor pool ("
        + getExecutionSlots() + " slots)...";
      if (m_log != null) {
        m_log.logMessage(msg);
      } else {
        System.err.println(msg);
      }
      // start the execution pool
      if (m_executorPool == null) {
        startExecutorPool();
      }
            
      // setup output queues
      msg = "[Classifier] " + statusMessagePrefix() + " setup output queues.";
      if (m_log != null) {
        m_log.logMessage(msg);
      } else {
        System.err.println(msg);
      }
      
      m_outputQueues = 
        new BatchClassifierEvent[e.getMaxRunNumber()][e.getMaxSetNumber()];
      m_completedSets = new boolean[e.getMaxRunNumber()][e.getMaxSetNumber()];
      m_currentBatchIdentifier = new Date();
    }
    
    // create a new task and schedule for execution
    TrainingTask newTask = new TrainingTask(e.getRunNumber(), e.getMaxRunNumber(),
        e.getSetNumber(), e.getMaxSetNumber(), e.getTrainingSet());
    String msg = "[Classifier] " + statusMessagePrefix() + " scheduling run " 
    + e.getRunNumber() +" fold " + e.getSetNumber() + " for execution...";
    if (m_log != null) {
      m_log.logMessage(msg);
    } else {
      System.err.println(msg);
    }
    
    // delay just a little bit
    /*try {
      Thread.sleep(10);
    } catch (Exception ex){} */
    m_executorPool.execute(newTask);
  }

  /**
   * Accepts a test set for a batch trained classifier
   *
   * @param e a <code>TestSetEvent</code> value
   */    
  public synchronized void acceptTestSet(TestSetEvent e) {
    if (m_reject) {
      if (m_log != null) {
        m_log.statusMessage(statusMessagePrefix() + "BUSY. Can't accept data "
            + "at this time.");
        m_log.logMessage("[Classifier] " + statusMessagePrefix()
            + " BUSY. Can't accept data at this time.");
      }
      return;
    }
    
    
    weka.classifiers.Classifier classifierToUse = m_Classifier;
    
    Instances testSet = e.getTestSet();
    if (testSet != null) {
      if (testSet.classIndex() < 0) {
  //        testSet.setClassIndex(testSet.numAttributes() - 1);
        // stop all processing
        stop();
        String errorMessage = statusMessagePrefix()
            + "ERROR: no class attribute set in test data!";
        if (m_log != null) {
          m_log.statusMessage(errorMessage);
          m_log.logMessage("[Classifier] " + errorMessage);
        } else {
          System.err.println("[Classifier] " + errorMessage);
        }
        return;
      }
    }

    // If we just have a test set connection or
    // there is just one run involving one set (and we are not
    // currently building a model), then use the
    // last saved model
    if (classifierToUse != null && m_state == IDLE && 
        (!m_listenees.containsKey("trainingSet") || 
        (e.getMaxRunNumber() == 1 && e.getMaxSetNumber() == 1))) {
      // if this is structure only then just return at this point
      if (e.getTestSet() != null && e.isStructureOnly()) {
        return;
      }
      
      if (classifierToUse instanceof EnvironmentHandler && m_env != null) {
        ((EnvironmentHandler)classifierToUse).setEnvironment(m_env);
      }
      
      if (classifierToUse instanceof weka.classifiers.misc.InputMappedClassifier) {
        // make sure that we have the correct training header (if InputMappedClassifier
        // is loading a model from a file).
        try {
          m_trainingSet = 
            ((weka.classifiers.misc.InputMappedClassifier)classifierToUse).
              getModelHeader(m_trainingSet); // this returns the argument if a model is not being loaded
        } catch (Exception e1) {
          // TODO Auto-generated catch block
          e1.printStackTrace();
        }
      }
      
      // check that we have a training set/header (if we don't,
      // then it means that no model has been loaded
      if (m_trainingSet == null) {
        stop();
        String errorMessage = statusMessagePrefix()
            + "ERROR: no trained/loaded classifier to use for prediction!";
        if (m_log != null) {
          m_log.statusMessage(errorMessage);
          m_log.logMessage("[Classifier] " + errorMessage);
        } else {
          System.err.println("[Classifier] " + errorMessage);
        }
        return;
      }
      
      testSet = e.getTestSet();
      if (e.getRunNumber() == 1 && e.getSetNumber() == 1) {
        m_currentBatchIdentifier = new Date();
      }
      
      if (testSet != null) {        
        if (!m_trainingSet.equalHeaders(testSet) && 
            !(classifierToUse instanceof weka.classifiers.misc.InputMappedClassifier)) {
          boolean wrapClassifier = false;
          if (!Utils.
              getDontShowDialog("weka.gui.beans.Classifier.AutoWrapInInputMappedClassifier")) {
            
            java.awt.GraphicsEnvironment ge = 
              java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
            if (!ge.isHeadless()) {
              JCheckBox dontShow = new JCheckBox("Do not show this message again");
              Object[] stuff = new Object[2];
              stuff[0] = "Data used to train model and test set are not compatible.\n" +
              "Would you like to automatically wrap the classifier in\n" + 
              "an \"InputMappedClassifier\" before proceeding?.\n";
              stuff[1] = dontShow;

              int result = JOptionPane.showConfirmDialog(this, stuff, 
                  "KnowledgeFlow:Classifier", JOptionPane.YES_OPTION);

              if (result == JOptionPane.YES_OPTION) {
                wrapClassifier = true;
              }

              if (dontShow.isSelected()) {
                String response = (wrapClassifier) ? "yes" : "no";
                try {
                  Utils.
                  setDontShowDialogResponse("weka.gui.explorer.ClassifierPanel.AutoWrapInInputMappedClassifier", 
                      response);
                } catch (Exception e1) {
                  // TODO Auto-generated catch block
                  e1.printStackTrace();
                }
              }
            } else {
              // running headless, so just go ahead and wrap anyway
              wrapClassifier = true;
            }
          } else {
            // What did the user say - do they want to autowrap or not?
            String response;
            try {
              response = Utils.getDontShowDialogResponse("weka.gui.explorer.ClassifierPanel.AutoWrapInInputMappedClassifier");
              if (response != null && response.equalsIgnoreCase("yes")) {
                wrapClassifier = true;
              }
            } catch (Exception e1) {
              // TODO Auto-generated catch block
              e1.printStackTrace();
            }
          }
          
          if (wrapClassifier) {
            weka.classifiers.misc.InputMappedClassifier temp =
              new weka.classifiers.misc.InputMappedClassifier();

            temp.setClassifier(classifierToUse);
            temp.setModelHeader(new Instances(m_trainingSet, 0));
            classifierToUse = temp;
          }          
        }         
        
        if (m_trainingSet.equalHeaders(testSet) || 
            (classifierToUse instanceof weka.classifiers.misc.InputMappedClassifier)) {
          BatchClassifierEvent ce =
            new BatchClassifierEvent(this, classifierToUse,                                       
                new DataSetEvent(this, m_trainingSet),
                new DataSetEvent(this, e.getTestSet()),
           e.getRunNumber(), e.getMaxRunNumber(), 
           e.getSetNumber(), e.getMaxSetNumber());
          ce.setGroupIdentifier(m_currentBatchIdentifier.getTime());
          
          if (m_log != null && !e.isStructureOnly()) {
            m_log.statusMessage(statusMessagePrefix() + "Finished.");
          }
          notifyBatchClassifierListeners(ce);          
        } else {
          // if headers do not match check to see if it's
          // just the class that is different and that
          // all class values are missing
          if (testSet.numInstances() > 0) {
            if (testSet.classIndex() == m_trainingSet.classIndex() && 
                testSet.attributeStats(testSet.classIndex()).missingCount ==
                testSet.numInstances()) {
              // now check the other attributes against the training
              // structure
              boolean ok = true;
              for (int i = 0; i < testSet.numAttributes(); i++) {
                if (i != testSet.classIndex()) {
                  ok = testSet.attribute(i).equals(m_trainingSet.attribute(i));
                  if (!ok) {
                    break;
                  }
                }
              }
              
              if (ok) {
                BatchClassifierEvent ce =
                  new BatchClassifierEvent(this, classifierToUse,                                       
                      new DataSetEvent(this, m_trainingSet),
                      new DataSetEvent(this, e.getTestSet()),
                 e.getRunNumber(), e.getMaxRunNumber(), 
                 e.getSetNumber(), e.getMaxSetNumber());
                ce.setGroupIdentifier(m_currentBatchIdentifier.getTime());
                
                if (m_log != null && !e.isStructureOnly()) {
                  m_log.statusMessage(statusMessagePrefix() + "Finished.");
                }
                notifyBatchClassifierListeners(ce);
              } else {
                stop();
                String errorMessage = statusMessagePrefix()
                + "ERROR: structure of training and test sets is not compatible!";
                if (m_log != null) {
                  m_log.statusMessage(errorMessage);
                  m_log.logMessage("[Classifier] " + errorMessage);
                } else {
                  System.err.println("[Classifier] " + errorMessage);
                }
              }
            }
          }
        }
      }
    } else {
/*      System.err.println("[Classifier] accepting test set: run " 
          + e.getRunNumber() + " fold " + e.getSetNumber()); */
      
      if (m_outputQueues[e.getRunNumber() - 1][e.getSetNumber() - 1] == null) {
        // store an event with a null model and training set (to be filled in later)
        m_outputQueues[e.getRunNumber() - 1][e.getSetNumber() - 1] =
          new BatchClassifierEvent(this, null, null, 
              new DataSetEvent(this, e.getTestSet()),
              e.getRunNumber(), e.getMaxRunNumber(),
              e.getSetNumber(), e.getMaxSetNumber());
        if (e.getRunNumber() == e.getMaxRunNumber() && 
            e.getSetNumber() == e.getMaxSetNumber()) {
          
          // block on the last fold of the last run
          /* System.err.println("[Classifier] blocking on last fold of last run...");
          block(true); */
          m_reject = true;
          if (m_block) {
            block(true);
          }
        }
      } else {
        // Otherwise, there is a model here waiting for a test set...
        m_outputQueues[e.getRunNumber() - 1][e.getSetNumber() - 1].
          setTestSet(new DataSetEvent(this, e.getTestSet()));
        checkCompletedRun(e.getRunNumber(), e.getMaxRunNumber(), e.getMaxSetNumber());
      }
    }
  }
  
  private synchronized void classifierTrainingComplete(BatchClassifierEvent ce) {
    // check the output queues if we have an incoming test set connection
    if (m_listenees.containsKey("testSet")) {
      String msg = "[Classifier] " + statusMessagePrefix() 
      + " storing model for run " + ce.getRunNumber() 
      + " fold " + ce.getSetNumber();
      if (m_log != null) {
        m_log.logMessage(msg);
      } else {
        System.err.println(msg);
      }
      
      if (m_outputQueues[ce.getRunNumber() - 1][ce.getSetNumber() - 1] == null) {
        // store the event - test data filled in later
        m_outputQueues[ce.getRunNumber() - 1][ce.getSetNumber() - 1] = ce;
      } else {
        // there is a test set here waiting for a model and training set
        m_outputQueues[ce.getRunNumber() - 1][ce.getSetNumber() - 1].
          setClassifier(ce.getClassifier());
        m_outputQueues[ce.getRunNumber() - 1][ce.getSetNumber() - 1].
        setTrainSet(ce.getTrainSet());
        
      }
      checkCompletedRun(ce.getRunNumber(), ce.getMaxRunNumber(), ce.getMaxSetNumber());
    }
  }
  
  private synchronized void checkCompletedRun(int runNum, int maxRunNum, int maxSets) {
    // look to see if there are any completed classifiers that we can pass
    // on for evaluation
    for (int i = 0; i < maxSets; i++) {
      if (m_outputQueues[runNum - 1][i] != null) {
        if (m_outputQueues[runNum - 1][i].getClassifier() != null &&
            m_outputQueues[runNum - 1][i].getTestSet() != null) {
          String msg = "[Classifier] " + statusMessagePrefix() 
          + " dispatching run/set " + runNum + "/" + (i+1) + " to listeners.";
          if (m_log != null) {
            m_log.logMessage(msg);
          } else {
            System.err.println(msg);
          }
          
          // dispatch this one
          m_outputQueues[runNum - 1][i].setGroupIdentifier(m_currentBatchIdentifier.getTime());
          notifyBatchClassifierListeners(m_outputQueues[runNum - 1][i]);
          // save memory
          m_outputQueues[runNum - 1][i] = null;
          // mark as done
          m_completedSets[runNum - 1][i] = true;
        }
      }
    }
    
    // scan for completion
    boolean done = true;
    for (int i = 0; i < maxRunNum; i++) {
      for (int j = 0; j < maxSets; j++) {
        if (!m_completedSets[i][j]) {
          done = false;
          break;
        }
      }
      if (!done) {
        break;
      }
    }
    
    if (done) {
      String msg = "[Classifier] " + statusMessagePrefix() 
      + " last classifier unblocking...";

      if (m_log != null) {
        m_log.logMessage(msg);
      } else {
        System.err.println(msg);
      }
      //m_visual.setText(m_oldText);

      if (m_log != null) {
        m_log.statusMessage(statusMessagePrefix() + "Finished.");
      }
      // m_outputQueues = null; // free memory
      m_reject = false;
      block(false);
      m_state = IDLE;
    }
  }

  /*private synchronized void checkCompletedRun(int runNum, int maxRunNum, int  maxSets) {
    boolean runOK = true;
    for (int i = 0; i < maxSets; i++) {
      if (m_outputQueues[runNum - 1][i] == null) {
        runOK = false;
        break;
      } else if (m_outputQueues[runNum - 1][i].getClassifier() == null ||
          m_outputQueues[runNum - 1][i].getTestSet() == null) {
        runOK = false;
        break;       
      }
    }
    
    if (runOK) {
      String msg = "[Classifier] " + statusMessagePrefix() 
        + " dispatching run " + runNum + " to listeners.";
      if (m_log != null) {
        m_log.logMessage(msg);
      } else {
        System.err.println(msg);
      }
      // dispatch this run to listeners
      for (int i = 0; i < maxSets; i++) {
        notifyBatchClassifierListeners(m_outputQueues[runNum - 1][i]);
        // save memory
        m_outputQueues[runNum - 1][i] = null;
      }
      
      if (runNum == maxRunNum) {
        // unblock
        msg = "[Classifier] " + statusMessagePrefix() 
        + " last classifier unblocking...";

        if (m_log != null) {
          m_log.logMessage(msg);
        } else {
          System.err.println(msg);
        }
        //m_visual.setText(m_oldText);

        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix() + "Finished.");
        }
        // m_outputQueues = null; // free memory
        m_reject = false;
        block(false);
        m_state = IDLE;
      }
    }
  } */

  /**
   * Sets the visual appearance of this wrapper bean
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Gets the visual appearance of this wrapper bean
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Use the default visual appearance for this bean
   */
  public void useDefaultVisual() {
    // try to get a default for this package of classifiers
    String name = m_ClassifierTemplate.getClass().toString();
    String packageName = name.substring(0, name.lastIndexOf('.'));
    packageName = 
      packageName.substring(packageName.lastIndexOf('.')+1,
                            packageName.length());
    if (!m_visual.loadIcons(BeanVisual.ICON_PATH+"Default_"+packageName
                            +"Classifier.gif",
                            BeanVisual.ICON_PATH+"Default_"+packageName
                            +"Classifier_animated.gif")) {
      m_visual.loadIcons(BeanVisual.
                         ICON_PATH+"DefaultClassifier.gif",
                         BeanVisual.
                         ICON_PATH+"DefaultClassifier_animated.gif");
    }
  }

  /**
   * Add a batch classifier listener
   *
   * @param cl a <code>BatchClassifierListener</code> value
   */
  public synchronized void 
    addBatchClassifierListener(BatchClassifierListener cl) {
    m_batchClassifierListeners.addElement(cl);
  }

  /**
   * Remove a batch classifier listener
   *
   * @param cl a <code>BatchClassifierListener</code> value
   */
  public synchronized void 
    removeBatchClassifierListener(BatchClassifierListener cl) {
    m_batchClassifierListeners.remove(cl);
  }

  /**
   * Notify all batch classifier listeners of a batch classifier event
   *
   * @param ce a <code>BatchClassifierEvent</code> value
   */
  private synchronized void notifyBatchClassifierListeners(BatchClassifierEvent ce) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_batchClassifierListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((BatchClassifierListener)l.elementAt(i)).acceptClassifier(ce);
      }
    }
  }

  /**
   * Add a graph listener
   *
   * @param cl a <code>GraphListener</code> value
   */
  public synchronized void addGraphListener(GraphListener cl) {
    m_graphListeners.addElement(cl);
  }

  /**
   * Remove a graph listener
   *
   * @param cl a <code>GraphListener</code> value
   */
  public synchronized void removeGraphListener(GraphListener cl) {
    m_graphListeners.remove(cl);
  }

  /**
   * Notify all graph listeners of a graph event
   *
   * @param ge a <code>GraphEvent</code> value
   */
  private void notifyGraphListeners(GraphEvent ge) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_graphListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((GraphListener)l.elementAt(i)).acceptGraph(ge);
      }
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
   * We don't have to keep track of configuration listeners (see the
   * documentation for ConfigurationListener/ConfigurationEvent).
   * 
   * @param cl a ConfigurationListener.
   */
  public synchronized void addConfigurationListener(ConfigurationListener cl) {
    
  }
  
  /**
   * We don't have to keep track of configuration listeners (see the
   * documentation for ConfigurationListener/ConfigurationEvent).
   * 
   * @param cl a ConfigurationListener.
   */
  public synchronized void removeConfigurationListener(ConfigurationListener cl) {
    
  }

  /**
   * Notify all text listeners of a text event
   *
   * @param ge a <code>TextEvent</code> value
   */
  private void notifyTextListeners(TextEvent ge) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_textListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((TextListener)l.elementAt(i)).acceptText(ge);
      }
    }
  }

  /**
   * Add an incremental classifier listener
   *
   * @param cl an <code>IncrementalClassifierListener</code> value
   */
  public synchronized void 
    addIncrementalClassifierListener(IncrementalClassifierListener cl) {
    m_incrementalClassifierListeners.add(cl);
  }

  /**
   * Remove an incremental classifier listener
   *
   * @param cl an <code>IncrementalClassifierListener</code> value
   */
  public synchronized void 
    removeIncrementalClassifierListener(IncrementalClassifierListener cl) {
    m_incrementalClassifierListeners.remove(cl);
  }

  /**
   * Notify all incremental classifier listeners of an incremental classifier
   * event
   *
   * @param ce an <code>IncrementalClassifierEvent</code> value
   */
  private void 
    notifyIncrementalClassifierListeners(IncrementalClassifierEvent ce) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_incrementalClassifierListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	((IncrementalClassifierListener)l.elementAt(i)).acceptClassifier(ce);
      }
    }
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection with respect to the named event
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    /*    if (eventName.compareTo("instance") == 0) {
      if (!(m_Classifier instanceof weka.classifiers.UpdateableClassifier)) {
	return false;
      }
      } */
    if (m_listenees.containsKey(eventName)) {
      return false;
    }
    return true;
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * EventSetDescriptor
   *
   * @param esd the EventSetDescriptor
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(EventSetDescriptor esd) {
    return connectionAllowed(esd.getName());
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source with respect to the named event
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void connectionNotification(String eventName,
						  Object source) {
    if (eventName.compareTo("instance") == 0) {
      if (!(m_ClassifierTemplate instanceof weka.classifiers.UpdateableClassifier)) {
	if (m_log != null) {
	  String msg = statusMessagePrefix() + "WARNING: "
          + m_ClassifierTemplate.getClass().getName() 
          + " Is not an updateable classifier. This "
          +"classifier will only be evaluated on incoming "
          +"instance events and not trained on them.";
	  m_log.logMessage("[Classifier] " + msg);
	  m_log.statusMessage(msg);
	}
      }
    }

    if (connectionAllowed(eventName)) {
      m_listenees.put(eventName, source);
      /*      if (eventName.compareTo("instance") == 0) {
	startIncrementalHandler();
	} */
    }
  }

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source with respect to the supplied event name
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void disconnectionNotification(String eventName,
						     Object source) {
    m_listenees.remove(eventName);
    if (eventName.compareTo("instance") == 0) {
      stop(); // kill the incremental handler thread if it is running
    }
  }

  /**
   * Function used to stop code that calls acceptTrainingSet. This is 
   * needed as classifier construction is performed inside a separate
   * thread of execution.
   *
   * @param tf a <code>boolean</code> value
   */
  private synchronized void block(boolean tf) {

    if (tf) {
      try {
	  // only block if thread is still doing something useful!
//	if (m_state != IDLE) {
	  wait();
	  //}
      } catch (InterruptedException ex) {
      }
    } else {
      notifyAll();
    }
  }


  /**
   * Stop any classifier action
   */
  public void stop() {
    // tell all listenees (upstream beans) to stop
    Enumeration en = m_listenees.keys();
    while (en.hasMoreElements()) {
      Object tempO = m_listenees.get(en.nextElement());
      if (tempO instanceof BeanCommon) {
	((BeanCommon)tempO).stop();
      }
    }
    
    // shutdown the executor pool and reclaim storage
    if (m_executorPool != null) {
      m_executorPool.shutdownNow();
      m_executorPool.purge();
      m_executorPool = null;
    }
    m_reject = false;
    block(false);
    m_visual.setStatic();
    if (m_oldText.length() > 0) {
      //m_visual.setText(m_oldText);
    }

    // stop the build thread
    /*if (m_buildThread != null) {
      m_buildThread.interrupt();
      m_buildThread.stop();
      m_buildThread = null;
      m_visual.setStatic();
    } */
  }

  public void loadModel() {
    try {
      if (m_fileChooser == null) {
        // i.e. after de-serialization
        setupFileChooser();
      }
      int returnVal = m_fileChooser.showOpenDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File loadFrom = m_fileChooser.getSelectedFile();

        // add extension if necessary
        if (m_fileChooser.getFileFilter() == m_binaryFilter) {
          if (!loadFrom.getName().toLowerCase().endsWith("." + FILE_EXTENSION)) {
            loadFrom = new File(loadFrom.getParent(),
                                loadFrom.getName() + "." + FILE_EXTENSION);
          }
        } else if (m_fileChooser.getFileFilter() == m_KOMLFilter) {
          if (!loadFrom.getName().toLowerCase().endsWith(KOML.FILE_EXTENSION 
                                                         + FILE_EXTENSION)) {
            loadFrom = new File(loadFrom.getParent(),
                                loadFrom.getName() + KOML.FILE_EXTENSION 
                                + FILE_EXTENSION);
          }
        } else if (m_fileChooser.getFileFilter() == m_XStreamFilter) {
          if (!loadFrom.getName().toLowerCase().endsWith(XStream.FILE_EXTENSION 
                                                        + FILE_EXTENSION)) {
            loadFrom = new File(loadFrom.getParent(),
                                loadFrom.getName() + XStream.FILE_EXTENSION 
                                + FILE_EXTENSION);
          }
        }

        weka.classifiers.Classifier temp = null;
        Instances tempHeader = null;
        // KOML ?
        if ((KOML.isPresent()) &&
            (loadFrom.getAbsolutePath().toLowerCase().
             endsWith(KOML.FILE_EXTENSION + FILE_EXTENSION))) {
          Vector v = (Vector) KOML.read(loadFrom.getAbsolutePath());
          temp = (weka.classifiers.Classifier) v.elementAt(0);
          if (v.size() == 2) {
            // try and grab the header
            tempHeader = (Instances) v.elementAt(1);
          }
        } /* XStream */ else if ((XStream.isPresent()) &&
                                 (loadFrom.getAbsolutePath().toLowerCase().
                                  endsWith(XStream.FILE_EXTENSION + FILE_EXTENSION))) {
          Vector v = (Vector) XStream.read(loadFrom.getAbsolutePath());
          temp = (weka.classifiers.Classifier) v.elementAt(0);
          if (v.size() == 2) {
            // try and grab the header
            tempHeader = (Instances) v.elementAt(1);
          } 
        } /* binary */ else {

          ObjectInputStream is = 
            new ObjectInputStream(new BufferedInputStream(
                                                          new FileInputStream(loadFrom)));
          // try and read the model
          temp = (weka.classifiers.Classifier)is.readObject();
          // try and read the header (if present)
          try {
            tempHeader = (Instances)is.readObject();
          } catch (Exception ex) {
            //            System.err.println("No header...");
            // quietly ignore
          }
          is.close();
        }

        // Update name and icon
        setTrainedClassifier(temp);
        // restore header
        m_trainingSet = tempHeader;

        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix() + "Loaded model.");
          m_log.logMessage("[Classifier] " + statusMessagePrefix() 
              + "Loaded classifier: "
              + m_Classifier.getClass().toString());
        }
      }
    } catch (Exception ex) {
      JOptionPane.showMessageDialog(Classifier.this,
                                    "Problem loading classifier.\n",
                                    "Load Model",
                                    JOptionPane.ERROR_MESSAGE);
      if (m_log != null) {
        m_log.statusMessage(statusMessagePrefix() + "ERROR: unable to load " +
        		"model (see log).");
        m_log.logMessage("[Classifier] " + statusMessagePrefix() 
            + "Problem loading classifier. " 
            + ex.getMessage());
      }
    }
  }

  public void saveModel() {
    try {
      if (m_fileChooser == null) {
        // i.e. after de-serialization
        setupFileChooser();
      }
      int returnVal = m_fileChooser.showSaveDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File saveTo = m_fileChooser.getSelectedFile();
        String fn = saveTo.getAbsolutePath();
        if (m_fileChooser.getFileFilter() == m_binaryFilter) {
          if (!fn.toLowerCase().endsWith("." + FILE_EXTENSION)) {
            fn += "." + FILE_EXTENSION;
          }
        } else if (m_fileChooser.getFileFilter() == m_KOMLFilter) {
          if (!fn.toLowerCase().endsWith(KOML.FILE_EXTENSION + FILE_EXTENSION)) {
            fn += KOML.FILE_EXTENSION + FILE_EXTENSION;
          }
        } else if (m_fileChooser.getFileFilter() == m_XStreamFilter) {
          if (!fn.toLowerCase().endsWith(XStream.FILE_EXTENSION + FILE_EXTENSION)) {
            fn += XStream.FILE_EXTENSION + FILE_EXTENSION;
          }
        }
        saveTo = new File(fn);

        // now serialize model
        // KOML?
        if ((KOML.isPresent()) &&
            saveTo.getAbsolutePath().toLowerCase().
            endsWith(KOML.FILE_EXTENSION + FILE_EXTENSION)) {
          SerializedModelSaver.saveKOML(saveTo,
                                        m_Classifier,
                                        (m_trainingSet != null)
                                        ? new Instances(m_trainingSet, 0)
                                        : null);
          /*          Vector v = new Vector();
          v.add(m_Classifier);
          if (m_trainingSet != null) {
            v.add(new Instances(m_trainingSet, 0));
          }
          v.trimToSize();
          KOML.write(saveTo.getAbsolutePath(), v); */
        } /* XStream */ else if ((XStream.isPresent()) &&
                                 saveTo.getAbsolutePath().toLowerCase().
            endsWith(XStream.FILE_EXTENSION + FILE_EXTENSION)) {

          SerializedModelSaver.saveXStream(saveTo,
                                           m_Classifier,
                                           (m_trainingSet != null)
                                           ? new Instances(m_trainingSet, 0)
                                           : null);
          /*          Vector v = new Vector();
          v.add(m_Classifier);
          if (m_trainingSet != null) {
            v.add(new Instances(m_trainingSet, 0));
          }
          v.trimToSize();
          XStream.write(saveTo.getAbsolutePath(), v); */
        } else /* binary */ {
          ObjectOutputStream os = 
            new ObjectOutputStream(new BufferedOutputStream(
                                   new FileOutputStream(saveTo)));
          os.writeObject(m_Classifier);
          if (m_trainingSet != null) {
            Instances header = new Instances(m_trainingSet, 0);
            os.writeObject(header);
          }
          os.close();
        }
        if (m_log != null) {
          m_log.statusMessage(statusMessagePrefix() + "Model saved.");
          m_log.logMessage("[Classifier] " + statusMessagePrefix() 
              + " Saved classifier " + getCustomName());
        }
      }
    } catch (Exception ex) {
      JOptionPane.showMessageDialog(Classifier.this,
                                    "Problem saving classifier.\n",
                                    "Save Model",
                                    JOptionPane.ERROR_MESSAGE);
      if (m_log != null) {
        m_log.statusMessage(statusMessagePrefix() + "ERROR: unable to" +
        		" save model (see log).");
        m_log.logMessage("[Classifier] " + statusMessagePrefix() 
            + " Problem saving classifier " + getCustomName() 
            + ex.getMessage());
      }
    }
  }

  /**
   * Set a logger
   *
   * @param logger a <code>Logger</code> value
   */
  public void setLog(Logger logger) {
    m_log = logger;
  }

  /**
   * Return an enumeration of requests that can be made by the user
   *
   * @return an <code>Enumeration</code> value
   */
  public Enumeration enumerateRequests() {
    Vector newVector = new Vector(0);
    if (m_executorPool != null && 
        (m_executorPool.getQueue().size() > 0 || 
            m_executorPool.getActiveCount() > 0)) {
      newVector.addElement("Stop");
    }

    if ((m_executorPool == null || 
        (m_executorPool.getQueue().size() == 0 && 
            m_executorPool.getActiveCount() == 0)) && 
        m_Classifier != null) {
      newVector.addElement("Save model");
    }

    if (m_executorPool == null || 
        (m_executorPool.getQueue().size() == 0 && 
            m_executorPool.getActiveCount() == 0)) {
      newVector.addElement("Load model");
    }
    return newVector.elements();
  }

  /**
   * Perform a particular request
   *
   * @param request the request to perform
   * @exception IllegalArgumentException if an error occurs
   */
  public void performRequest(String request) {
    if (request.compareTo("Stop") == 0) {
      stop();
    } else if (request.compareTo("Save model") == 0) {
      saveModel();
    } else if (request.compareTo("Load model") == 0) {
      loadModel();
    } else {
      throw new IllegalArgumentException(request
					 + " not supported (Classifier)");
    }
  }

  /**
   * Returns true, if at the current time, the event described by the
   * supplied event descriptor could be generated.
   *
   * @param esd an <code>EventSetDescriptor</code> value
   * @return a <code>boolean</code> value
   */
  public boolean eventGeneratable(EventSetDescriptor esd) {
    String eventName = esd.getName();
    return eventGeneratable(eventName);
  }
  
  /**
   * @param name of the event to check
   * @return true if eventName is one of the possible events
   * that this component can generate
   */
  private boolean generatableEvent(String eventName) {
    if (eventName.compareTo("graph") == 0
	|| eventName.compareTo("text") == 0
	|| eventName.compareTo("batchClassifier") == 0
	|| eventName.compareTo("incrementalClassifier") == 0
	|| eventName.compareTo("configuration") == 0) {
      return true;
    }
    return false;
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
    if (!generatableEvent(eventName)) {
      return false;
    }
    if (eventName.compareTo("graph") == 0) {
      // can't generate a GraphEvent if classifier is not drawable
      if (!(m_ClassifierTemplate instanceof weka.core.Drawable)) {
	return false;
      }
      // need to have a training set before the classifier
      // can generate a graph!
      if (!m_listenees.containsKey("trainingSet")) {
	return false;
      }
      // Source needs to be able to generate a trainingSet
      // before we can generate a graph
      Object source = m_listenees.get("trainingSet");
       if (source instanceof EventConstraints) {
	if (!((EventConstraints)source).eventGeneratable("trainingSet")) {
	  return false;
	}
      }
    }

    if (eventName.compareTo("batchClassifier") == 0) {
      /*      if (!m_listenees.containsKey("testSet")) {
        return false;
      }
      if (!m_listenees.containsKey("trainingSet") && 
          m_trainingSet == null) {
	return false;
        } */
      if (!m_listenees.containsKey("testSet") && 
          !m_listenees.containsKey("trainingSet")) {
        return false;
      }
      Object source = m_listenees.get("testSet");
      if (source instanceof EventConstraints) {
	if (!((EventConstraints)source).eventGeneratable("testSet")) {
	  return false;
	}
      }
      /*      source = m_listenees.get("trainingSet");
      if (source instanceof EventConstraints) {
	if (!((EventConstraints)source).eventGeneratable("trainingSet")) {
	  return false;
	}
        } */
    }

    if (eventName.compareTo("text") == 0) {
      if (!m_listenees.containsKey("trainingSet") &&
	  !m_listenees.containsKey("instance")) {
	return false;
      }
      Object source = m_listenees.get("trainingSet");
      if (source != null && source instanceof EventConstraints) {
	if (!((EventConstraints)source).eventGeneratable("trainingSet")) {
	  return false;
	}
      }
      source = m_listenees.get("instance");
      if (source != null && source instanceof EventConstraints) {
	if (!((EventConstraints)source).eventGeneratable("instance")) {
	  return false;
	}
      }
    }

    if (eventName.compareTo("incrementalClassifier") == 0) {
      /*      if (!(m_Classifier instanceof weka.classifiers.UpdateableClassifier)) {
	return false;
	} */
      if (!m_listenees.containsKey("instance")) {
	return false;
      }
      Object source = m_listenees.get("instance");
      if (source instanceof EventConstraints) {
	if (!((EventConstraints)source).eventGeneratable("instance")) {
	  return false;
	}
      }
    }
    
    if (eventName.equals("configuration") && m_Classifier == null) {
      return false;
    }
    
    return true;
  }
    
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    if (m_executorPool == null || 
        (m_executorPool.getQueue().size() == 0 && 
            m_executorPool.getActiveCount() == 0) && m_state == IDLE) {
      return false;
    }
    /* System.err.println("isBusy() Q:" + m_executorPool.getQueue().size()
        +" A:" + m_executorPool.getActiveCount()); */
    return true;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|"
    + ((m_Classifier instanceof OptionHandler &&
        Utils.joinOptions(((OptionHandler)m_Classifier).getOptions()).length() > 0) 
        ? Utils.joinOptions(((OptionHandler)m_Classifier).getOptions()) + "|"
            : "");
  }

  /**
   * Set environment variables to pass on to the classifier (if
   * if is an EnvironmentHandler)
   */
  public void setEnvironment(Environment env) {
    m_env = env;
  }
}

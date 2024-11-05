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
 *    SerializedModelSaver.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.io.ObjectInputStream;
import java.io.Serializable;
import java.io.File;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.util.ArrayList;
import java.util.Vector;
import javax.swing.JPanel;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Instances;
import weka.core.Environment;
import weka.core.EnvironmentHandler;
import weka.core.xml.KOML;
import weka.core.xml.XStream;
import weka.core.Tag;
import weka.core.Utils;

/**
 * A bean that saves serialized models
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org
 * @version $Revision: 6804 $
 */
public class SerializedModelSaver
  extends JPanel
  implements BeanCommon, Visible, BatchClassifierListener, 
             IncrementalClassifierListener, BatchClustererListener,
	     EnvironmentHandler, Serializable {

  /** for serialization */
  private static final long serialVersionUID = 3956528599473814287L;

  /**
   * Default visual for data sources
   */
  protected BeanVisual m_visual = 
    new BeanVisual("AbstractDataSink", 
		   BeanVisual.ICON_PATH+"SerializedModelSaver.gif",
		   BeanVisual.ICON_PATH+"SerializedModelSaver_animated.gif");

  /**
   * Non null if this object is a target for any events.
   * Provides for the simplest case when only one incomming connection
   * is allowed.
   */
  protected Object m_listenee = null;

  /**
   * The log for this bean
   */
  protected transient weka.gui.Logger m_logger = null;

  /**
   * The prefix for the file name (model + training set info will be appended)
   */
  private String m_filenamePrefix = "";

  /**
   * The directory to hold the saved model(s)
   */
  private File m_directory = new File(System.getProperty("user.dir"));

  /**
   * File format stuff
   */
  private Tag m_fileFormat;

  public final static int BINARY = 0;
  public final static int KOMLV = 1;
  public final static int XSTREAM = 2;

  /** the extension for serialized models (binary Java serialization) */
  public final static String FILE_EXTENSION = "model";

  /** relative path for the directory (relative to the user.dir (startup directory))? */
  private boolean m_useRelativePath = false;
  
  /** include relation name in filename */
  private boolean m_includeRelationName = false;

  /**
   * Available file formats. Reflection is used to check if classes
   * are available for deep object serialization to XML
   */
  public static ArrayList<Tag> s_fileFormatsAvailable;
  static {
    s_fileFormatsAvailable = new ArrayList<Tag>();
    s_fileFormatsAvailable.add(new Tag(BINARY, "Binary serialized model file (*"
                                       + FILE_EXTENSION + ")", "", false));
    if (KOML.isPresent()) {
      s_fileFormatsAvailable.add(new Tag(KOMLV,
                                         "XML serialized model file (*"
                                         + KOML.FILE_EXTENSION + FILE_EXTENSION + ")", "", false));
    }

    if (XStream.isPresent()) {
      s_fileFormatsAvailable.add(new Tag(XSTREAM,
                                         "XML serialized model file (*"
                                         + XStream.FILE_EXTENSION + FILE_EXTENSION + ")", "", false));
    }
  }
  
  /**
   * The environment variables.
   */
  protected transient Environment m_env;

  /**
   * Constructor.
   */
  public SerializedModelSaver() {
    useDefaultVisual();
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
    m_fileFormat = s_fileFormatsAvailable.get(0);
    
    m_env = Environment.getSystemWide();
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
   * Use the default images for this bean.
   *
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"SerializedModelSaver.gif",
		       BeanVisual.ICON_PATH+"SerializedModelSaver_animated.gif");
    m_visual.setText("SerializedModelSaver");
  }

  /**
   * Set the visual for this data source.
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Get the visual being used by this data source.
   *
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * EventSetDescriptor.
   *
   * @param esd the EventSetDescriptor
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(EventSetDescriptor esd) {
    return connectionAllowed(esd.getName());
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * event name.
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    return (m_listenee == null);
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source with respect to the supplied event name.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void connectionNotification(String eventName,
						  Object source) {
    if (connectionAllowed(eventName)) {
      m_listenee = source;
    }
  }

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source with respect to the supplied event name.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void disconnectionNotification(String eventName,
						     Object source) {
    if (m_listenee == source) {
      m_listenee = null;
    }
  }
  
  /**
   * Set a log for this bean.
   *
   * @param logger a <code>weka.gui.Logger</code> value
   */
  public void setLog(weka.gui.Logger logger) {
    m_logger = logger;
  }

  /**
   * Stop any processing that the bean might be doing.
   */
  public void stop() {
    // tell the listenee (upstream bean) to stop
    if (m_listenee instanceof BeanCommon) {
      ((BeanCommon)m_listenee).stop();
    }
  }
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return false;
  }

  /**
   * makes sure that the filename is valid, i.e., replaces slashes,
   * backslashes and colons with underscores ("_").
   * 
   * @param filename	the filename to cleanse
   * @return		the cleansed filename
   */
  protected String sanitizeFilename(String filename) {
    return filename.replaceAll("\\\\", "_").replaceAll(":", "_").replaceAll("/", "_");
  }

  /**
   * Accept and save a batch trained clusterer.
   *
   * @param ce a <code>BatchClassifierEvent</code> value
   */
  public void acceptClusterer(BatchClustererEvent ce) {
    if (ce.getTestSet() == null || 
        ce.getTestOrTrain() == BatchClustererEvent.TEST ||
        ce.getTestSet().isStructureOnly()) {
      return;
    }

    Instances trainHeader = new Instances(ce.getTestSet().getDataSet(), 0);
    String titleString = ce.getClusterer().getClass().getName();		      
    titleString = titleString.
      substring(titleString.lastIndexOf('.') + 1,
                titleString.length());

    String prefix = "";
    String relationName = (m_includeRelationName)
    ? trainHeader.relationName()
    : "";
    try {
      prefix = m_env.substitute(m_filenamePrefix);
    } catch (Exception ex) {
      stop(); // stop all processing
      String message = "[SerializedModelSaver] " 
        + statusMessagePrefix() 
        + " Can't save model. Reason: " 
        + ex.getMessage();
      if (m_logger != null) {
        m_logger.logMessage(message);
        m_logger.statusMessage(statusMessagePrefix()
            + "ERROR (See log for details)");
      } else {
        System.err.println(message);
      }
      return;
    }
    String fileName = "" 
      + prefix
      + relationName
      + titleString
      + "_"
      + ce.getSetNumber() 
      + "_" + ce.getMaxSetNumber();
    fileName = sanitizeFilename(fileName);
    
    String dirName = m_directory.getPath();
    try {
      dirName = m_env.substitute(dirName);
    } catch (Exception ex) {
      stop(); // stop all processing
      String message = "[SerializedModelSaver] "
        + statusMessagePrefix() + " Can't save model. Reason: " 
                           + ex.getMessage();
      if (m_logger != null) {
        m_logger.logMessage(message);
        m_logger.statusMessage(statusMessagePrefix()
            + "ERROR (See log for details)");
      } else {
        System.err.println(message);
      }
      return;
    }
    File tempFile = new File(dirName);
    fileName = tempFile.getAbsolutePath() 
      + File.separator
      + fileName;

    saveModel(fileName, trainHeader, ce.getClusterer());
  }

  /**
   * Accept and save an incrementally trained classifier.
   *
   * @param ce the BatchClassifierEvent containing the classifier
   */
  public void acceptClassifier(final IncrementalClassifierEvent ce) {
    if (ce.getStatus() == IncrementalClassifierEvent.BATCH_FINISHED) {
      // Only save model when the end of the stream is reached
      Instances header = ce.getStructure();
      String titleString = ce.getClassifier().getClass().getName();		      
      titleString = titleString.
        substring(titleString.lastIndexOf('.') + 1,
                  titleString.length());

      String prefix = "";
      String relationName = (m_includeRelationName)
        ? header.relationName()
        : "";
        
      try {
        prefix = m_env.substitute(m_filenamePrefix);
      } catch (Exception ex) {
        stop(); // stop processing
        String message = "[SerializedModelSaver] "
          + statusMessagePrefix() + " Can't save model. Reason: " 
          + ex.getMessage();
        if (m_logger != null) {
          m_logger.logMessage(message);
          m_logger.statusMessage(statusMessagePrefix()
              + "ERROR (See log for details)");
        } else {
          System.err.println(message);
        }
        return;
      }
      
      String fileName = "" + prefix + relationName + titleString;
      fileName = sanitizeFilename(fileName);

      String dirName = m_directory.getPath();
      try {
        dirName = m_env.substitute(dirName);
      } catch (Exception ex) {
        stop(); // stop processing
        String message = "[SerializedModelSaver] "
          + statusMessagePrefix() + " Can't save model. Reason: " 
          + ex.getMessage();
        if (m_logger != null) {
          m_logger.logMessage(message);
          m_logger.statusMessage(statusMessagePrefix()
              + "ERROR (See log for details)");
        } else {
          System.err.println(message);
        }
        return;
      }
      File tempFile = new File(dirName);

      fileName = tempFile.getAbsolutePath() 
        + File.separator
        + fileName;
      
      saveModel(fileName, header, ce.getClassifier());
    }
  }
  
  /**
   * Accept and save a batch trained classifier.
   *
   * @param ce the BatchClassifierEvent containing the classifier
   */
  public void acceptClassifier(final BatchClassifierEvent ce) {
    if (ce.getTrainSet() == null || 
        ce.getTrainSet().isStructureOnly()) {
      return;
    }
    Instances trainHeader = new Instances(ce.getTrainSet().getDataSet(), 0);
    
    // adjust for InputMappedClassifier (if necessary)
    if (ce.getClassifier() instanceof weka.classifiers.misc.InputMappedClassifier) {
      try {
        trainHeader = 
          ((weka.classifiers.misc.InputMappedClassifier)ce.getClassifier()).
            getModelHeader(trainHeader);
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
    String titleString = ce.getClassifier().getClass().getName();		      
    titleString = titleString.
      substring(titleString.lastIndexOf('.') + 1,
                titleString.length());

    String prefix = "";
    String relationName = (m_includeRelationName)
    ? trainHeader.relationName()
    : "";
    try {
      prefix = m_env.substitute(m_filenamePrefix);
    } catch (Exception ex) {
      stop(); // stop processing
      String message = "[SerializedModelSaver] "
        + statusMessagePrefix() + " Can't save model. Reason: " 
        + ex.getMessage();
      if (m_logger != null) {
        m_logger.logMessage(message);
        m_logger.statusMessage(statusMessagePrefix()
            + "ERROR (See log for details)");
      } else {
        System.err.println(message);
      }
      return;
    }

    String fileName = "" 
      + prefix
      + relationName
      + titleString
      + "_"
      + ce.getSetNumber() 
      + "_" + ce.getMaxSetNumber();
    fileName = sanitizeFilename(fileName);
    
    String dirName = m_directory.getPath();
    try {
      dirName = m_env.substitute(dirName);
    } catch (Exception ex) {
      stop(); // stop processing
      String message = "[SerializedModelSaver] "
        + statusMessagePrefix() + " Can't save model. Reason: " 
                           + ex.getMessage();
      if (m_logger != null) {
        m_logger.logMessage(message);
        m_logger.statusMessage(statusMessagePrefix()
            + "ERROR (See log for details)");
      } else {
        System.err.println(message);
      }
      return;
    }
    File tempFile = new File(dirName);

    fileName = tempFile.getAbsolutePath() 
      + File.separator
      + fileName;

    saveModel(fileName, trainHeader, ce.getClassifier());
  }

  /**
   * Helper routine to actually save the models.
   */
  private void saveModel(String fileName, Instances trainHeader, Object model) {
    m_fileFormat = validateFileFormat(m_fileFormat);
    if (m_fileFormat == null) {
      // default to binary if validation fails
      m_fileFormat = s_fileFormatsAvailable.get(0);
    }
    try {
      switch (m_fileFormat.getID()) {
      case KOMLV:
        fileName = fileName + KOML.FILE_EXTENSION + FILE_EXTENSION;
        saveKOML(new File(fileName), model, trainHeader);
        break;
      case XSTREAM:
        fileName = fileName + XStream.FILE_EXTENSION + FILE_EXTENSION;
        saveXStream(new File(fileName), model, trainHeader);
        break;
      default:
        fileName = fileName + "." + FILE_EXTENSION;
        saveBinary(new File(fileName), model, trainHeader);
        break;
      }        
    } catch (Exception ex) {
      stop(); // stop all processing
      System.err.println("[SerializedModelSaver] Problem saving model");
      if (m_logger != null) {
        m_logger.logMessage("[SerializedModelSaver] "
            + statusMessagePrefix() + " Problem saving model");
        m_logger.statusMessage(statusMessagePrefix()
            + "ERROR (See log for details)");
      }
    }
  }

  /**
   * Save a model in binary form.
   *
   * @param saveTo the file name to save to
   * @param model the model to save
   * @param header the header of the data that was used to train the model (optional)
   */
  public static void saveBinary(File saveTo, Object model, Instances header) throws IOException {
    ObjectOutputStream os =
      new ObjectOutputStream(new BufferedOutputStream(
                             new FileOutputStream(saveTo)));
    os.writeObject(model);
    // now the header
    if (header != null) {
      os.writeObject(header);
    }
    os.close();
  }

  /**
   * Save a model in KOML deep object serialized XML form.
   *
   * @param saveTo the file name to save to
   * @param model the model to save
   * @param header the header of the data that was used to train the model (optional)
   */
  public static void saveKOML(File saveTo, Object model, Instances header) throws Exception {
    Vector v = new Vector();
    v.add(model);
    if (header != null) {
      v.add(header);
    }
    v.trimToSize();
    KOML.write(saveTo.getAbsolutePath(), v);
  }

  /**
   * Save a model in XStream deep object serialized XML form.
   *
   * @param saveTo the file name to save to
   * @param model the model to save
   * @param header the header of the data that was used to train the model (optional)
   */
  public static void saveXStream(File saveTo, Object model, Instances header) throws Exception {
    Vector v = new Vector();
    v.add(model);
    if (header != null) {
      v.add(header);
    }
    v.trimToSize();
    XStream.write(saveTo.getAbsolutePath(), v);
  }

  /**
   * Get the directory that the model(s) will be saved into
   *
   * @return the directory to save to
   */
  public File getDirectory() {
    return m_directory;
  }
  
  /**
   * Set the directory that the model(s) will be saved into.
   *
   * @param d the directory to save to
   */
  public void setDirectory(File d) {
    m_directory = d;
    if (m_useRelativePath) {
      try {
        m_directory = Utils.convertToRelativePath(m_directory);
      } catch (Exception ex) {
      }
    }
  }

  /**
   * Set whether to use relative paths for the directory.
   * I.e. relative to the startup (user.dir) directory
   *
   * @param rp true if relative paths are to be used
   */
  public void setUseRelativePath(boolean rp) {
    m_useRelativePath = rp;
  }
  
  /**
   * Get whether to use relative paths for the directory.
   * I.e. relative to the startup (user.dir) directory
   *
   * @return true if relative paths are to be used
   */
  public boolean getUseRelativePath() {
    return m_useRelativePath;
  }
  
  /**
   * Set whether the relation name of the training data
   * used to create the model should be included as part
   * of the filename for the serialized model.
   * 
   * @param rn true if the relation name should be included
   * in the file name
   */
  public void setIncludeRelationName(boolean rn) {
    m_includeRelationName = rn;
  }
  
  /**
   * Get whether the relation name of the training
   * data used to create the model is to be included
   * in the filename of the serialized model.
   * 
   * @return true if the relation name is to be included
   * in the file name
   */
  public boolean getIncludeRelationName() {
    return m_includeRelationName;
  }

  /**
   * Get the prefix to prepend to the model file names.
   *
   * @return the prefix to prepend
   */
  public String getPrefix() {
    return m_filenamePrefix;
  }

  /**
   * Set the prefix to prepend to the model file names.
   *
   * @param p the prefix to prepend
   */
  public void setPrefix(String p) {
    m_filenamePrefix = p;
  }

  /**
   * Global info for this bean. Gets displayed in the GUI.
   *
   * @return information about this bean.
   */
  public String globalInfo() {
    return "Save trained models to serialized object files.";
  }

  /**
   * Set the file format to use for saving.
   *
   * @param ff the file format to use
   */
  public void setFileFormat(Tag ff) {
    m_fileFormat = ff;
  }

  /**
   * Get the file format to use for saving.
   *
   * @return the file format to use
   */
  public Tag getFileFormat() {
    return m_fileFormat;
  }

  /**
   * Validate the file format. After this bean is deserialized, classes for
   * XML serialization may not be in the classpath any more.
   *
   * @param ff the current file format to validate
   */
  public Tag validateFileFormat(Tag ff) {
    Tag r = ff;
    if (ff.getID() == BINARY) {
      return ff;
    }

    if (ff.getID() == KOMLV && !KOML.isPresent()) {
      r = null;
    }

    if (ff.getID() == XSTREAM && !XStream.isPresent()) {
      r = null;
    }

    return r;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|";
  }
  
  /**
   * Set environment variables to use.
   * 
   * @param env the environment variables to
   * use
   */
  public void setEnvironment(Environment env) {
    m_env = env;
  }
  
  // Custom de-serialization in order to set default
  // environment variables on de-serialization
  private void readObject(ObjectInputStream aStream) 
    throws IOException, ClassNotFoundException {
    aStream.defaultReadObject();
    
    // set a default environment to use
    m_env = Environment.getSystemWide();
  }
}

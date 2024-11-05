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
 *    ClassAssigner.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.Serializable;
import java.util.Vector;

import javax.swing.JPanel;

/**
 * Bean that assigns a class attribute to a data set.
 *
 * @author Mark Hall
 * @version $Revision: 5667 $
 */
public class ClassAssigner
  extends JPanel
  implements Visible, DataSourceListener, TrainingSetListener, TestSetListener,
	     DataSource, TrainingSetProducer, TestSetProducer,
	     BeanCommon, EventConstraints, Serializable,
	     InstanceListener {

  /** for serialization */
  private static final long serialVersionUID = 4011131665025817924L;
  
  private String m_classColumn = "last";

  /** format of instances for current incoming connection (if any) */
  private Instances m_connectedFormat;

  private Object m_trainingProvider;
  private Object m_testProvider;
  private Object m_dataProvider;
  private Object m_instanceProvider;

  private Vector m_trainingListeners = new Vector();
  private Vector m_testListeners = new Vector();
  private Vector m_dataListeners = new Vector();
  private Vector m_instanceListeners = new Vector();

  private Vector m_dataFormatListeners = new Vector();

  protected transient weka.gui.Logger m_logger = null;

  protected BeanVisual m_visual = 
    new BeanVisual("ClassAssigner", 
		   BeanVisual.ICON_PATH+"ClassAssigner.gif",
		   BeanVisual.ICON_PATH+"ClassAssigner_animated.gif");

  /**
   * Global info for this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "Designate which column is to be considered the class column "
      +"in incoming data.";
  }

  public ClassAssigner() {
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);    
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
   * Tool tip text for this property
   *
   * @return a <code>String</code> value
   */
  public String classColumnTipText() {
    return "Specify the number of the column that contains the class attribute";
  }

  /**
   * Returns the structure of the incoming instances (if any)
   *
   * @return an <code>Instances</code> value
   */
  public Instances getConnectedFormat() {
    return m_connectedFormat;
  }

  public void setClassColumn(String col) {
    m_classColumn = col;
    if (m_connectedFormat != null) {
      assignClass(m_connectedFormat);
    }
  }

  public String getClassColumn() {
    return m_classColumn;
  }

  public void acceptDataSet(DataSetEvent e) {
    Instances dataSet = e.getDataSet();
    assignClass(dataSet);
    notifyDataListeners(e);
    if (e.isStructureOnly()) {
      m_connectedFormat = e.getDataSet();
      // tell any listening customizers (or other
      notifyDataFormatListeners();
    }
  }

  public void acceptTrainingSet(TrainingSetEvent e) {
    Instances trainingSet = e.getTrainingSet();
    assignClass(trainingSet);
    notifyTrainingListeners(e);
    
    if (e.isStructureOnly()) {
      m_connectedFormat = e.getTrainingSet();
      // tell any listening customizers (or other
      notifyDataFormatListeners();
    }
  }

  public void acceptTestSet(TestSetEvent e) {
    Instances testSet = e.getTestSet();
    assignClass(testSet);
    notifyTestListeners(e);
    if (e.isStructureOnly()) {
      m_connectedFormat = e.getTestSet();
      // tell any listening customizers (or other
      notifyDataFormatListeners();
    }
  }

  public void acceptInstance(InstanceEvent e) {
    if (e.getStatus() == InstanceEvent.FORMAT_AVAILABLE) {
      //      Instances dataSet = e.getInstance().dataset();
      m_connectedFormat = e.getStructure();
      
      //      System.err.println("Assigning class column...");
      assignClass(m_connectedFormat);
      notifyInstanceListeners(e);

      // tell any listening customizers (or other interested parties)
      System.err.println("Notifying customizer...");
      notifyDataFormatListeners();
    } else {
      //      Instances dataSet = e.getInstance().dataset();
      //      assignClass(dataSet);
      notifyInstanceListeners(e);
    }
  }

  private void assignClass(Instances dataSet) {
    int classCol = -1;
    if (m_classColumn.toLowerCase().compareTo("last") == 0) {
      dataSet.setClassIndex(dataSet.numAttributes()-1);
    } else if (m_classColumn.toLowerCase().compareTo("first") == 0) {
      dataSet.setClassIndex(0);
    } else {
      classCol = Integer.parseInt(m_classColumn) - 1;
      if (/*classCol < 0 ||*/ classCol > dataSet.numAttributes()-1) {
	if (m_logger != null) {
	  m_logger.logMessage("Class column outside range of data "
			      +"(ClassAssigner)");
	}
      } else {
	dataSet.setClassIndex(classCol);
      }
    }
  }

  protected void notifyTestListeners(TestSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_testListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	System.err.println("Notifying test listeners "
			   +"(ClassAssigner)");
	((TestSetListener)l.elementAt(i)).acceptTestSet(tse);
      }
    }
  }

  protected void notifyTrainingListeners(TrainingSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_trainingListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	System.err.println("Notifying training listeners "
			   +"(ClassAssigner)");
	((TrainingSetListener)l.elementAt(i)).acceptTrainingSet(tse);
      }
    }
  }

  protected void notifyDataListeners(DataSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_dataListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	System.err.println("Notifying data listeners "
			   +"(ClassAssigner)");
	((DataSourceListener)l.elementAt(i)).acceptDataSet(tse);
      }
    }
  }

  protected void notifyInstanceListeners(InstanceEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_instanceListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	//	System.err.println("Notifying instance listeners "
	//			   +"(ClassAssigner)");
	
	((InstanceListener)l.elementAt(i)).acceptInstance(tse);
      }
    }
  }

  protected void notifyDataFormatListeners() {
    Vector l;
    synchronized (this) {
      l = (Vector)m_dataFormatListeners.clone();
    }
    if (l.size() > 0) {
      DataSetEvent dse = new DataSetEvent(this, m_connectedFormat);
      for(int i = 0; i < l.size(); i++) {
	//	System.err.println("Notifying instance listeners "
	//			   +"(ClassAssigner)");
	((DataFormatListener)l.elementAt(i)).newDataFormat(dse);
      }
    }
  }

  public synchronized void addInstanceListener(InstanceListener tsl) {
    m_instanceListeners.addElement(tsl);
    if (m_connectedFormat != null) {
      InstanceEvent e = new InstanceEvent(this, m_connectedFormat);
      tsl.acceptInstance(e);
    }
  }

  public synchronized void removeInstanceListener(InstanceListener tsl) {
    m_instanceListeners.removeElement(tsl);
  }

  public synchronized void addDataSourceListener(DataSourceListener tsl) {
    m_dataListeners.addElement(tsl);
    // pass on any format that we might know about
    if (m_connectedFormat != null) {
      DataSetEvent e = new DataSetEvent(this, m_connectedFormat);
      tsl.acceptDataSet(e);
    }
  }

  public synchronized void removeDataSourceListener(DataSourceListener tsl) {
    m_dataListeners.removeElement(tsl);
  }

  public synchronized void addTrainingSetListener(TrainingSetListener tsl) {
    m_trainingListeners.addElement(tsl);
    // pass on any format that we might know about
    if (m_connectedFormat != null) {
      TrainingSetEvent e = new TrainingSetEvent(this, m_connectedFormat);
      tsl.acceptTrainingSet(e);
    }
  }

  public synchronized void removeTrainingSetListener(TrainingSetListener tsl) {
    m_trainingListeners.removeElement(tsl);
  }

  public synchronized void addTestSetListener(TestSetListener tsl) {
    m_testListeners.addElement(tsl);
    // pass on any format that we might know about
    if (m_connectedFormat != null) {
      TestSetEvent e = new TestSetEvent(this, m_connectedFormat);
      tsl.acceptTestSet(e);
    }
  }

  public synchronized void removeTestSetListener(TestSetListener tsl) {
    m_testListeners.removeElement(tsl);
  }

  public synchronized void addDataFormatListener(DataFormatListener dfl) {
    m_dataFormatListeners.addElement(dfl);
  }

  public synchronized void removeDataFormatListener(DataFormatListener dfl) {
    m_dataFormatListeners.removeElement(dfl);
  }

  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  public BeanVisual getVisual() {
    return m_visual;
  }
  
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"ClassAssigner.gif",
		       BeanVisual.ICON_PATH+"ClassAssigner_animated.gif");
  }

  /**
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * event name
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    if (eventName.compareTo("trainingSet") == 0 && 
	(m_trainingProvider != null || m_dataProvider != null ||
	 m_instanceProvider != null)) { 
      return false;
    }
    
    if (eventName.compareTo("testSet") == 0 && 
	m_testProvider != null) { 
      return false;
    }

     if (eventName.compareTo("instance") == 0 &&
	m_instanceProvider != null || m_trainingProvider != null ||
	 m_dataProvider != null) {
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
   * a source with respect to the supplied event name
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public synchronized void connectionNotification(String eventName,
						  Object source) {
    if (connectionAllowed(eventName)) {
      if (eventName.compareTo("trainingSet") == 0) {
	m_trainingProvider = source;
      } else if (eventName.compareTo("testSet") == 0) {
	m_testProvider = source;
      } else if (eventName.compareTo("dataSet") == 0) {
	m_dataProvider = source;
      } else if (eventName.compareTo("instance") == 0) {
	m_instanceProvider = source;
      }
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

    if (eventName.compareTo("trainingSet") == 0) {
      if (m_trainingProvider == source) {
	m_trainingProvider = null;
      }
    }
    if (eventName.compareTo("testSet") == 0) {
      if (m_testProvider == source) {
	m_testProvider = null;
      }
    }
    if (eventName.compareTo("dataSet") == 0) {
      if (m_dataProvider == source) {
	m_dataProvider = null;
      }
    }

    if (eventName.compareTo("instance") == 0) {
      if (m_instanceProvider == source) {
	m_instanceProvider = null;
      }
    }
  }
  
  public void setLog(weka.gui.Logger logger) {
    m_logger = logger;
  }

  public void stop() {
    // Pass on to upstream beans
    if (m_trainingProvider != null && m_trainingProvider instanceof BeanCommon) {
      ((BeanCommon)m_trainingProvider).stop();
    }
    
    if (m_testProvider != null && m_testProvider instanceof BeanCommon) {
      ((BeanCommon)m_testProvider).stop();
    }
    
    if (m_dataProvider != null && m_dataProvider instanceof BeanCommon) {
      ((BeanCommon)m_dataProvider).stop();
    }
    
    if (m_instanceProvider != null && m_instanceProvider instanceof BeanCommon) {
      ((BeanCommon)m_instanceProvider).stop();
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
   * Returns true, if at the current time, the named event could
   * be generated. Assumes that the supplied event name is
   * an event that could be generated by this bean
   *
   * @param eventName the name of the event in question
   * @return true if the named event could be generated at this point in
   * time
   */
  public boolean eventGeneratable(String eventName) {
    if (eventName.compareTo("trainingSet") == 0) { 
      if (m_trainingProvider == null) {
	return false;
      } else {
	if (m_trainingProvider instanceof EventConstraints) {
	  if (!((EventConstraints)m_trainingProvider).
	      eventGeneratable("trainingSet")) {
	    return false;
	  }
	}
      }
    }

    if (eventName.compareTo("dataSet") == 0) { 
      if (m_dataProvider == null) {
	if (m_instanceProvider == null) {
	  m_connectedFormat = null;
	  notifyDataFormatListeners();
	}
	return false;
      } else {
	if (m_dataProvider instanceof EventConstraints) {
	  if (!((EventConstraints)m_dataProvider).
	      eventGeneratable("dataSet")) {
	    m_connectedFormat = null;
	    notifyDataFormatListeners();
	    return false;
	  }
	}
      }
    }

    if (eventName.compareTo("instance") == 0) { 
      if (m_instanceProvider == null) {
	if (m_dataProvider == null) {
	  m_connectedFormat = null;
	  notifyDataFormatListeners();
	}
	return false;
      } else {
	if (m_instanceProvider instanceof EventConstraints) {
	  if (!((EventConstraints)m_instanceProvider).
	      eventGeneratable("instance")) {
	    m_connectedFormat = null;
	    notifyDataFormatListeners();
	    return false;
	  }
	}
      }
    }

    if (eventName.compareTo("testSet") == 0) {
      if (m_testProvider == null) {
	return false;
      } else {
	if (m_testProvider instanceof EventConstraints) {
	  if (!((EventConstraints)m_testProvider).
	      eventGeneratable("testSet")) {
	    return false;
	  }
	}
      }
    }
    return true;
  }
}

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
 *    ClassValuePicker.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instances;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.SwapValues;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.Serializable;
import java.util.Vector;

import javax.swing.JPanel;

/**
 * @author Mark Hall
 * @version $Revision: 5731 $
 */
public class ClassValuePicker
  extends JPanel
  implements Visible, DataSourceListener, BeanCommon,
	     EventConstraints, Serializable {

  /** for serialization */
  private static final long serialVersionUID = -1196143276710882989L;

  /** the class value index considered to be the positive class */
  private int m_classValueIndex = 0;

  /** format of instances for the current incoming connection (if any) */
  private Instances m_connectedFormat;

  private Object m_dataProvider;

  private Vector m_dataListeners = new Vector();
  private Vector m_dataFormatListeners = new Vector();

  protected transient weka.gui.Logger m_logger = null;
  
  protected BeanVisual m_visual = 
    new BeanVisual("ClassValuePicker", 
		   BeanVisual.ICON_PATH+"ClassValuePicker.gif",
		   BeanVisual.ICON_PATH+"ClassValuePicker_animated.gif");

  /**
   * Global info for this bean
   *
   * @return a <code>String</code> value
   */
  public String globalInfo() {
    return "Designate which class value is to be considered the \"positive\" "
      +"class value (useful for ROC style curves).";
  }

  public ClassValuePicker() {
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
   * Returns the structure of the incoming instances (if any)
   *
   * @return an <code>Instances</code> value
   */
  public Instances getConnectedFormat() {
    if (m_connectedFormat ==null) {
      System.err.println("Is null!!!!!!");
    }
    return m_connectedFormat;
  }

  /**
   * Set the class value index considered to be the "positive"
   * class value.
   *
   * @param index the class value index to use
   */
  public void setClassValueIndex(int index) {
    m_classValueIndex = index;
    if (m_connectedFormat != null) {
      notifyDataFormatListeners();
    }
  }

  /**
   * Gets the class value index considered to be the "positive"
   * class value.
   *
   * @return the class value index
   */
  public int getClassValueIndex() {
    return m_classValueIndex;
  }

  public void acceptDataSet(DataSetEvent e) {
    if (e.isStructureOnly()) {
      if (m_connectedFormat == null ||
	  !m_connectedFormat.equalHeaders(e.getDataSet())) { 
	m_connectedFormat = new Instances(e.getDataSet(), 0);
	// tell any listening customizers (or other
	notifyDataFormatListeners();
      }
    }
    Instances dataSet = e.getDataSet();
    Instances newDataSet = assignClassValue(dataSet);
    e = new DataSetEvent(this, newDataSet);
    notifyDataListeners(e);
  }

  private Instances assignClassValue(Instances dataSet) {
    if (dataSet.classIndex() < 0) {
      if (m_logger != null) {
	m_logger.
	  logMessage("[ClassValuePicker] " 
	      + statusMessagePrefix() 
	      + " No class attribute defined in data set.");
	m_logger.statusMessage(statusMessagePrefix()
	    + "WARNING: No class attribute defined in data set.");
      }
      return dataSet;
    }
    
    if (dataSet.classAttribute().isNumeric()) {
      if (m_logger != null) {
	m_logger.
	  logMessage("[ClassValuePicker] "
	      + statusMessagePrefix()
	      + " Class attribute must be nominal (ClassValuePicker)");
	m_logger.statusMessage(statusMessagePrefix()
	    + "WARNING: Class attribute must be nominal.");
      }
      return dataSet;
    } else {
      if (m_logger != null) {
        m_logger.statusMessage(statusMessagePrefix() + "remove");
      }
    }

    if (m_classValueIndex != 0) { // nothing to do if == 0
      // swap selected index with index 0
      try {
	SwapValues sv = new SwapValues();
	sv.setAttributeIndex(""+(dataSet.classIndex()+1));
	sv.setFirstValueIndex("first");
	sv.setSecondValueIndex(""+(m_classValueIndex+1));
	sv.setInputFormat(dataSet);
	Instances newDataSet = Filter.useFilter(dataSet, sv);
	newDataSet.setRelationName(dataSet.relationName());
	return newDataSet;
      } catch (Exception ex) {
	if (m_logger != null) {
	  m_logger.
	    logMessage("[ClassValuePicker] "
	        +statusMessagePrefix()
	        + " Unable to swap class attibute values.");
	  m_logger.statusMessage(statusMessagePrefix()
	      + "ERROR (See log for details)");
	}
      }
    }
    return dataSet;
  }

  protected void notifyDataListeners(DataSetEvent tse) {
    Vector l;
    synchronized (this) {
      l = (Vector)m_dataListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
	System.err.println("Notifying data listeners "
			   +"(ClassValuePicker)");
	((DataSourceListener)l.elementAt(i)).acceptDataSet(tse);
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
	((DataFormatListener)l.elementAt(i)).newDataFormat(dse);
      }
    }
  }

  public synchronized void addDataSourceListener(DataSourceListener tsl) {
    m_dataListeners.addElement(tsl);
  }

  public synchronized void removeDataSourceListener(DataSourceListener tsl) {
    m_dataListeners.removeElement(tsl);
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
    m_visual.loadIcons(BeanVisual.ICON_PATH+"ClassValuePicker.gif",
		       BeanVisual.ICON_PATH+"ClassValuePicker_animated.gif");
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
    if (eventName.compareTo("dataSet") == 0 && 
	(m_dataProvider != null)) { 
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
      if (eventName.compareTo("dataSet") == 0) {
	m_dataProvider = source;
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

    if (eventName.compareTo("dataSet") == 0) {
      if (m_dataProvider == source) {
	m_dataProvider = null;
      }
    }
  }

  public void setLog(weka.gui.Logger logger) {
    m_logger = logger;
  }

  public void stop() {
    // nothing to do
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
    if (eventName.compareTo("dataSet") != 0) {
      return false;
    }

    if (eventName.compareTo("dataSet") == 0) { 
      if (m_dataProvider == null) {
	m_connectedFormat = null;
	notifyDataFormatListeners();
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
    return true;
  }
  
  private String statusMessagePrefix() {
    return getCustomName() + "$" + hashCode() + "|";
  }
}

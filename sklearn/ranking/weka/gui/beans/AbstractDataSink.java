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
 *    AbstractDataSink.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.Serializable;

import javax.swing.JPanel;

/**
 * Abstract class for objects that store instances to some destination.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5753 $
 * @since 1.0
 * @see JPanel
 * @see Serializable
 */

public abstract class AbstractDataSink
  extends JPanel
  implements DataSink, BeanCommon, Visible, 
	     DataSourceListener, TrainingSetListener,
	     TestSetListener, InstanceListener, ThresholdDataListener,
	     Serializable {

  /** for serialization */
  private static final long serialVersionUID = 3956528599473814287L;

  /**
   * Default visual for data sources
   */
  protected BeanVisual m_visual = 
    new BeanVisual("AbstractDataSink", 
		   BeanVisual.ICON_PATH+"DefaultDataSink.gif",
		   BeanVisual.ICON_PATH+"DefaultDataSink_animated.gif");

  /**
   * Non null if this object is a target for any events.
   * Provides for the simplest case when only one incomming connection
   * is allowed. Subclasses can overide the appropriate BeanCommon methods
   * to change this behaviour and allow multiple connections if desired
   */
  protected Object m_listenee = null;

  protected transient weka.gui.Logger m_logger = null;

  public AbstractDataSink() {
    useDefaultVisual();
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }

  /**
   * Accept a training set
   *
   * @param e a <code>TrainingSetEvent</code> value
   */
  public abstract void acceptTrainingSet(TrainingSetEvent e);

  /**
   * Accept a test set
   *
   * @param e a <code>TestSetEvent</code> value
   */
  public abstract void acceptTestSet(TestSetEvent e);

  /**
   * Accept a data set
   *
   * @param e a <code>DataSetEvent</code> value
   */
  public abstract void acceptDataSet(DataSetEvent e);
  
  /**
   * Accept a threshold data set
   *
   * @param e a <code>ThresholdDataEvent</code> value
   */
  public abstract void acceptDataSet(ThresholdDataEvent e);
  
  /**
   * Accept an instance
   *
   * @param e an <code>InstanceEvent</code> value
   */
  public abstract void acceptInstance(InstanceEvent e);

  /**
   * Set the visual for this data source
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
   * Use the default images for a data source
   *
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultDataSink.gif",
		       BeanVisual.ICON_PATH+"DefaultDataSink_animated.gif");
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
   * Returns true if, at this time, 
   * the object will accept a connection according to the supplied
   * event name
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    return (m_listenee == null);
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
      m_listenee = source;
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
    if (m_listenee == source) {
      m_listenee = null;
    }
  }
  
  /**
   * Set a log for this bean
   *
   * @param logger a <code>weka.gui.Logger</code> value
   */
  public void setLog(weka.gui.Logger logger) {
    m_logger = logger;
  }

  /**
   * Stop any processing that the bean might be doing.
   * Subclass must implement
   */
  public abstract void stop();
}

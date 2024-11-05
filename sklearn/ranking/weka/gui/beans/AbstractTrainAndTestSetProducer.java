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
 *    AbstractTrainAndTestSetProducer.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.io.Serializable;
import java.util.Vector;

import javax.swing.JPanel;

/**
 * Abstract base class for TrainAndTestSetProducers that contains default
 * implementations of add/remove listener methods and defualt
 * visual representation.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.5 $
 */
public abstract class AbstractTrainAndTestSetProducer 
  extends JPanel
  implements Visible, TrainingSetProducer, TestSetProducer, 
	     BeanCommon, Serializable, DataSourceListener {

  /** for serialization */
  private static final long serialVersionUID = -1809339823613492037L;

  /**
   * Objects listening for trainin set events
   */
  protected Vector m_trainingListeners = new Vector();

  /**
   * Objects listening for test set events
   */
  protected Vector m_testListeners = new Vector();

  protected BeanVisual m_visual = 
    new BeanVisual("AbstractTrainingSetProducer", 
		   BeanVisual.ICON_PATH+"DefaultTrainTest.gif",
		   BeanVisual.ICON_PATH+"DefaultTrainTest_animated.gif");

  /**
   * non null if this object is a target for any events.
   */
  protected Object m_listenee = null;

  protected transient weka.gui.Logger m_logger = null;

  /**
   * Creates a new <code>AbstractTrainAndTestSetProducer</code> instance.
   */
  public AbstractTrainAndTestSetProducer() {
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }

  /**
   * Subclass must implement
   *
   * @param e a <code>DataSetEvent</code> value
   */
  public abstract void acceptDataSet(DataSetEvent e);

  /**
   * Add a training set listener
   *
   * @param tsl a <code>TrainingSetListener</code> value
   */
  public synchronized void addTrainingSetListener(TrainingSetListener tsl) {
    m_trainingListeners.addElement(tsl);
  }

  /**
   * Remove a training set listener
   *
   * @param tsl a <code>TrainingSetListener</code> value
   */
  public synchronized void removeTrainingSetListener(TrainingSetListener tsl) {
    m_trainingListeners.removeElement(tsl);
  }

  /**
   * Add a test set listener
   *
   * @param tsl a <code>TestSetListener</code> value
   */
  public synchronized void addTestSetListener(TestSetListener tsl) {
    m_testListeners.addElement(tsl);
  }

  /**
   * Remove a test set listener
   *
   * @param tsl a <code>TestSetListener</code> value
   */
  public synchronized void removeTestSetListener(TestSetListener tsl) {
    m_testListeners.removeElement(tsl);
  }

  /**
   * Set the visual for this bean
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Get the visual for this bean
   *
   * @return a <code>BeanVisual</code> value
   */
  public BeanVisual getVisual() {
    return m_visual;
  }
  
  /**
   * Use the default visual for this bean
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultTrainTest.gif",
		       BeanVisual.ICON_PATH+"DefaultTrainTest_animated.gif");
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


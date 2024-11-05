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
 *    InstanceStreamToBatchMaker.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.beans.EventSetDescriptor;
import java.util.ArrayList;

import javax.swing.JPanel;

import weka.core.Instance;
import weka.core.Instances;
import weka.gui.Logger;

/**
 * Bean that converts an instance stream into a (batch) data set.
 * Useful in conjunction with the Reservoir sampling filter.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5563 $
 */
public class InstanceStreamToBatchMaker extends JPanel 
  implements BeanCommon, Visible, InstanceListener, 
  EventConstraints, DataSource {
  
  /**
   * For serialization
   */
  private static final long serialVersionUID = -7037141087208627799L;

  protected BeanVisual m_visual = 
    new BeanVisual("InstanceStreamToBatchMaker",
                   BeanVisual.ICON_PATH+"InstanceStreamToBatchMaker.gif",
                   BeanVisual.ICON_PATH+"InstanceStreamToBatchMaker_animated.gif");
  
  /**
   * The log.
   */
  private transient Logger m_log;
  
  /**
   * Component connected to us.
   */
  private Object m_listenee;
  
  private ArrayList<DataSourceListener> m_dataListeners = 
    new ArrayList<DataSourceListener>();
  
  /**
   * Collects up the instances. 
   */
  private ArrayList<Instance> m_batch;
  
  private Instances m_structure;
  
  public InstanceStreamToBatchMaker() {
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
  }
  
  /**
   * Accept an instance to add to the batch.
   *
   * @param e an <code>InstanceEvent</code> value
   */
  public void acceptInstance(InstanceEvent e) {
    if (e.getStatus() == InstanceEvent.FORMAT_AVAILABLE) {
      m_batch = new ArrayList<Instance>();
      m_structure = e.getStructure();
      
      // notify dataset listeners of structure available
      if (m_log != null) {
        m_log.logMessage("[InstanceStreamToBatch] passing on structure.");
      }
      DataSetEvent dse = new DataSetEvent(this, m_structure);
      notifyDataListeners(dse);
    } else if (e.getStatus() == InstanceEvent.INSTANCE_AVAILABLE) {
      m_batch.add(e.getInstance());
    } else {
      // batch finished
      
      // add the last instance
      m_batch.add(e.getInstance());
      
      // create the new Instances
      Instances dataSet = new Instances(m_structure, m_batch.size());
      for (Instance i : m_batch) {
        dataSet.add(i);
      }
      dataSet.compactify();
      
      // save memory
      m_batch = null;
      
      if (m_log != null) {
        m_log.logMessage("[InstanceStreamToBatch] sending batch to listeners.");
      }
      
      // notify dataset listeners
      DataSetEvent dse = new DataSetEvent(this, dataSet);
      notifyDataListeners(dse);
    }
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
   * the object will accept a connection with respect to the named event
   *
   * @param eventName the event
   * @return true if the object will accept a connection
   */
  public boolean connectionAllowed(String eventName) {
    if (m_listenee != null || !eventName.equals("instance")) {
      return false;
    }
    return true;
  }

  /**
   * Notify this object that it has been registered as a listener with
   * a source with respect to the named event
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  public void connectionNotification(String eventName, Object source) {
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
  public void disconnectionNotification(String eventName, Object source) {
    m_listenee = null;
  }
  
  /**
   * Returns true if, at the current time, the named event could be
   * generated.
   *
   * @param eventName the name of the event in question
   * @return true if the named event could be generated
   */
  public boolean eventGeneratable(String eventName) {
    if (!eventName.equals("dataSet")) {
      return false;
    }
    
    if (m_listenee == null) {
      return false;
    }
    
    if (m_listenee instanceof EventConstraints) {
      if (!((EventConstraints)m_listenee).eventGeneratable("instance")) {
        return false;
      }
    }
    
    return true;
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
   * Set a custom (descriptive) name for this bean
   * 
   * @param name the name to use
   */
  public void setCustomName(String name) {
    m_visual.setText(name);
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
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  public boolean isBusy() {
    return false;
  }

  /**
   * Stop any action (if possible).
   */
  public void stop() {
    // not much we can do. Stopping depends on upstream components.
  }

  /**
   * Gets the visual appearance of this wrapper bean
   */
  public BeanVisual getVisual() {
    return m_visual;
  }

  /**
   * Sets the visual appearance of this wrapper bean
   *
   * @param newVisual a <code>BeanVisual</code> value
   */
  public void setVisual(BeanVisual newVisual) {
    m_visual = newVisual;
  }

  /**
   * Use the default visual appearance for this bean
   */
  public void useDefaultVisual() {
    m_visual.loadIcons(BeanVisual.ICON_PATH+"InstanceStreamToBatchMaker.gif",
        BeanVisual.ICON_PATH+"InstanceStreamToBatchMaker_animated.gif");
  }
  
  /**
   * Notify all data source listeners.
   * 
   * @param tse the DataSetEvent to pass on.
   */
  protected void notifyDataListeners(DataSetEvent tse) {
    ArrayList<DataSourceListener> l;
    synchronized (this) {
      l = (ArrayList<DataSourceListener>)m_dataListeners.clone();
    }
    if (l.size() > 0) {
      for(int i = 0; i < l.size(); i++) {
        l.get(i).acceptDataSet(tse);
      }
    }
  }

  public synchronized void addDataSourceListener(DataSourceListener tsl) {
    m_dataListeners.add(tsl);
    // pass on any format that we might know about
    if (m_structure != null) {
      DataSetEvent e = new DataSetEvent(this, m_structure);
      tsl.acceptDataSet(e);
    }
  }

  public synchronized void removeDataSourceListener(DataSourceListener tsl) {
    m_dataListeners.remove(tsl);
  }
  
  public synchronized void addInstanceListener(InstanceListener il) {
    // we don't produce instance events
  }
  
  public synchronized void removeInstanceListener(InstanceListener il) {
    // we don't produce instance events
  }
}

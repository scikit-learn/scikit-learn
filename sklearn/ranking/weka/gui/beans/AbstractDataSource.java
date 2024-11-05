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
 *    AbstractDataSource.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.awt.BorderLayout;
import java.beans.PropertyChangeListener;
import java.beans.VetoableChangeListener;
import java.beans.beancontext.BeanContext;
import java.beans.beancontext.BeanContextChild;
import java.beans.beancontext.BeanContextChildSupport;
import java.io.Serializable;
import java.util.Vector;

import javax.swing.JPanel;

/**
 * Abstract class for objects that can provide instances from some source
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.4 $
 * @since 1.0
 * @see JPanel
 * @see DataSource
 * @see Serializable
 */
public abstract class AbstractDataSource
  extends JPanel
  implements DataSource, Visible, Serializable, BeanContextChild {

  /** for serialization */
  private static final long serialVersionUID = -4127257701890044793L;

  /**
   * True if this bean's appearance is the design mode appearance
   */
  protected boolean m_design;
  
  /**
   * BeanContex that this bean might be contained within
   */
  protected transient BeanContext m_beanContext = null;

  /**
   * BeanContextChild support
   */
  protected BeanContextChildSupport m_bcSupport = 
    new BeanContextChildSupport(this);

  /**
   * Default visual for data sources
   */
  protected BeanVisual m_visual = 
    new BeanVisual("AbstractDataSource", 
		   BeanVisual.ICON_PATH+"DefaultDataSource.gif",
		   BeanVisual.ICON_PATH+"DefaultDataSource_animated.gif");
  
  /**
   * Objects listening for events from data sources
   */
  protected Vector m_listeners;

  /**
   * Creates a new <code>AbstractDataSource</code> instance.
   *
   */
  public AbstractDataSource() {    
    useDefaultVisual();
    setLayout(new BorderLayout());
    add(m_visual, BorderLayout.CENTER);
    m_listeners = new Vector();
  }

  /**
   * Add a listener
   *
   * @param dsl a <code>DataSourceListener</code> value
   */
  public synchronized void addDataSourceListener(DataSourceListener dsl) {
    m_listeners.addElement(dsl);
  }
  
  /**
   * Remove a listener
   *
   * @param dsl a <code>DataSourceListener</code> value
   */
  public synchronized void removeDataSourceListener(DataSourceListener dsl) {
    m_listeners.remove(dsl);
  }

  /**
   * Add an instance listener
   *
   * @param dsl a <code>InstanceListener</code> value
   */
  public synchronized void addInstanceListener(InstanceListener dsl) {
    m_listeners.addElement(dsl);
  }
  
  /**
   * Remove an instance listener
   *
   * @param dsl a <code>InstanceListener</code> value
   */
  public synchronized void removeInstanceListener(InstanceListener dsl) {
    m_listeners.remove(dsl);
  }

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
    m_visual.loadIcons(BeanVisual.ICON_PATH+"DefaultDataSource.gif",
		       BeanVisual.ICON_PATH+"DefaultDataSource_animated.gif");
  }

  /**
   * Set a bean context for this bean
   *
   * @param bc a <code>BeanContext</code> value
   */
  public void setBeanContext(BeanContext bc) {
    m_beanContext = bc;
    m_design = m_beanContext.isDesignTime();
  }

  /**
   * Return the bean context (if any) that this bean is embedded in
   *
   * @return a <code>BeanContext</code> value
   */
  public BeanContext getBeanContext() {
    return m_beanContext;
  }

  /**
   * Add a property change listener to this bean
   *
   * @param name the name of the property of interest
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void addPropertyChangeListener(String name,
					PropertyChangeListener pcl) {
    m_bcSupport.addPropertyChangeListener(name, pcl);
  }

  /**
   * Remove a property change listener from this bean
   *
   * @param name the name of the property of interest
   * @param pcl a <code>PropertyChangeListener</code> value
   */
  public void removePropertyChangeListener(String name,
					   PropertyChangeListener pcl) {
    m_bcSupport.removePropertyChangeListener(name, pcl);
  }

  /**
   * Add a vetoable change listener to this bean
   *
   * @param name the name of the property of interest
   * @param vcl a <code>VetoableChangeListener</code> value
   */
  public void addVetoableChangeListener(String name,
				       VetoableChangeListener vcl) {
    m_bcSupport.addVetoableChangeListener(name, vcl);
  }
  
  /**
   * Remove a vetoable change listener from this bean
   *
   * @param name the name of the property of interest
   * @param vcl a <code>VetoableChangeListener</code> value
   */
  public void removeVetoableChangeListener(String name,
					   VetoableChangeListener vcl) {
    m_bcSupport.removeVetoableChangeListener(name, vcl);
  }
}








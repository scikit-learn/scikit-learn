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
 *    InstanceEvent.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import weka.core.Instance;
import weka.core.Instances;

import java.util.EventObject;

/**
 * Event that encapsulates a single instance or header information only
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.5 $
 * @see EventObject
 */
public class InstanceEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 6104920894559423946L;
  
  public static final int FORMAT_AVAILABLE = 0;
  public static final int INSTANCE_AVAILABLE = 1;
  public static final int BATCH_FINISHED = 2;
  
  private Instances m_structure;
  private Instance m_instance;
  private int m_status;

  /**
   * Creates a new <code>InstanceEvent</code> instance that encapsulates
   * a single instance only.
   *
   * @param source the source of the event
   * @param instance the instance
   * @param status status code (either INSTANCE_AVAILABLE or BATCH_FINISHED)
   */
  public InstanceEvent(Object source, Instance instance, int status) {
    super(source);
    m_instance = instance;
    m_status = status;
  }

  /**
   * Creates a new <code>InstanceEvent</code> instance which encapsulates
   * header information only.
   *
   * @param source an <code>Object</code> value
   * @param structure an <code>Instances</code> value
   */
  public InstanceEvent(Object source, Instances structure) {
    super(source);
    m_structure = structure;
    m_status = FORMAT_AVAILABLE;
  }

  public InstanceEvent(Object source) {
    super(source);
  }
  
  /**
   * Get the instance
   *
   * @return an <code>Instance</code> value
   */
  public Instance getInstance() {
    return m_instance;
  }
  
  /**
   * Set the instance
   *
   * @param i an <code>Instance</code> value
   */
  public void setInstance(Instance i) {
    m_instance = i;
  }

  /**
   * Get the status
   *
   * @return an <code>int</code> value
   */
  public int getStatus() {
    return m_status;
  }

  /**
   * Set the status
   *
   * @param s an <code>int</code> value
   */
  public void setStatus(int s) {
    m_status = s;
  }

  /**
   * Set the instances structure
   *
   * @param structure an <code>Instances</code> value
   */
  public void setStructure(Instances structure) {
    m_structure = structure;
    m_instance = null;
    m_status = FORMAT_AVAILABLE;
  }

  /**
   * Get the instances structure (may be null if this is not
   * a FORMAT_AVAILABLE event)
   *
   * @return an <code>Instances</code> value
   */
  public Instances getStructure() {
    return m_structure;
  }
}

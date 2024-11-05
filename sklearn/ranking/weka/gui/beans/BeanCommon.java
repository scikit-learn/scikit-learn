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
 *    BeanCommon.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.beans.EventSetDescriptor;

/**
 * Interface specifying routines that all weka beans should implement
 * in order to allow the bean environment to exercise some control over the
 * bean and also to allow the bean to exercise some control over connections.
 *
 * Beans may want to  impose constraints in terms of
 * the number of connections they will allow via a particular 
 * listener interface. Some beans may only want to be registered
 * as a listener for a particular event type with only one source, or
 * perhaps a limited number of sources.
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 5820 $
 * @since 1.0
 */
public interface BeanCommon {

  /**
   * Set a custom (descriptive) name for this bean
   * 
   * @param name the name to use
   */
  void setCustomName(String name);
  
  /**
   * Get the custom (descriptive) name for this bean (if one has been set)
   * 
   * @return the custom name (or the default name)
   */
  String getCustomName();

  /**
   * Stop any processing that the bean might be doing.
   */
  void stop();
  
  /**
   * Returns true if. at this time, the bean is busy with some
   * (i.e. perhaps a worker thread is performing some calculation).
   * 
   * @return true if the bean is busy.
   */
  boolean isBusy();

  /**
   * Set a logger
   *
   * @param logger a <code>weka.gui.Logger</code> value
   */
  void setLog(weka.gui.Logger logger);


  /**
   * Returns true if, at this time, 
   * the object will accept a connection via the named event
   *
   * @param esd the EventSetDescriptor for the event in question
   * @return true if the object will accept a connection
   */
  boolean connectionAllowed(EventSetDescriptor esd);

  /**
   * Returns true if, at this time, 
   * the object will accept a connection via the named event
   *
   * @param eventName the name of the event
   * @return true if the object will accept a connection
   */
  boolean connectionAllowed(String eventName);

  /**
   * Notify this object that it has been registered as a listener with
   * a source for recieving events described by the named event
   * This object is responsible for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  void connectionNotification(String eventName, Object source);

  /**
   * Notify this object that it has been deregistered as a listener with
   * a source for named event. This object is responsible
   * for recording this fact.
   *
   * @param eventName the event
   * @param source the source with which this object has been registered as
   * a listener
   */
  void disconnectionNotification(String eventName, Object source);
}

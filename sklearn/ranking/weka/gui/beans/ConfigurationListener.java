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
 *    ConfigurationListener.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

/**
 * Matching listener for ConfigurationEvent. Implementers of
 * ConfigurationProducer do not actually have to generate the
 * event (nor will listeners ever need to process ConfigurationEvent).
 * Configurations will be pulled (rather than pushed) by
 * ConfigurationListeners. It is a listener's responsibility (if
 * they are interested in utilizing configurations) to implement
 * BeanCommon and store/delete reference(s) to ConfigurationProducers 
 * when connectionNotification() and disconnectionNotification() are
 * called on them.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org)
 * @version $Revision $
 */
public interface ConfigurationListener {
  
  /**
   * Implementers do not have to do anything in this
   * method (see the above documentation).
   * 
   * @param e a ConfigurationEvent
   */
  void acceptConfiguration(ConfigurationEvent e);
}

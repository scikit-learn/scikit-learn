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
 *    ConfigurationProducer.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

/**
 * Marker interface for components that can share their configuration.
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org)
 * @version $Revision $
 */
public interface ConfigurationProducer {
  
  /**
   * We don't have to keep track of configuration listeners (see the
   * documentation for ConfigurationListener/ConfigurationEvent).
   * 
   * @param cl a ConfigurationListener.
   */
  void addConfigurationListener(ConfigurationListener cl);
  
  /**
   * We don't have to keep track of configuration listeners (see the
   * documentation for ConfigurationListener/ConfigurationEvent).
   * 
   * @param cl a ConfigurationListener.
   */
  void removeConfigurationListener(ConfigurationListener cl);
}

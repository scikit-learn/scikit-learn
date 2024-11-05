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
 * EnvironmentHandler.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import weka.core.Environment;

/**
 * Interface for something that can utilize environment
 * variables. NOTE: since environment variables should
 * be transient, the implementer needs to be careful
 * of state after de-serialization. Default system-wide
 * environment variables can be got via a call to
 * <code>weka.core.Environment.getSystemWide()</code>
 * 
 * @author mhall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 5563 $
 */
public interface EnvironmentHandler {
  
  /**
   * Set environment variables to use.
   * 
   * @param env the environment variables to
   * use
   */
  void setEnvironment(Environment env);
}

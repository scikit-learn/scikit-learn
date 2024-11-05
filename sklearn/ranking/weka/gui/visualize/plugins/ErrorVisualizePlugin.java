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
 * ErrorVisualizePlugin.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize.plugins;

import weka.core.Instances;

import javax.swing.JMenuItem;

/**
 * Interface implemented by classes loaded dynamically to
 * visualize classifier errors in the explorer.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5012 $
 */
public interface ErrorVisualizePlugin {

  /**
   * Get a JMenu or JMenuItem which contain action listeners
   * that perform the visualization of the classifier errors.  
   * <p/>
   * The actual class is the attribute declared as class attribute, the
   * predicted class values is found in the attribute prior to the class
   * attribute's position. In other words, if the <code>classIndex()</code>
   * method returns 10, then the attribute position for the predicted class 
   * values is 9.
   * <p/>
   * Exceptions thrown because of changes in Weka since compilation need to 
   * be caught by the implementer.
   *
   * @see NoClassDefFoundError
   * @see IncompatibleClassChangeError
   *
   * @param predInst 	the instances with the actual and predicted class values
   * @return menuitem 	for opening visualization(s), or null
   *         		to indicate no visualization is applicable for the input
   */
  public JMenuItem getVisualizeMenuItem(Instances predInst);

  /**
   * Get the minimum version of Weka, inclusive, the class
   * is designed to work with.  eg: <code>3.5.0</code>
   * 
   * @return		the minimum version
   */
  public String getMinVersion();

  /**
   * Get the maximum version of Weka, exclusive, the class
   * is designed to work with.  eg: <code>3.6.0</code>
   * 
   * @return		the maximum version
   */
  public String getMaxVersion();

  /**
   * Get the specific version of Weka the class is designed for.
   * eg: <code>3.5.1</code>
   * 
   * @return		the version the plugin was designed for
   */
  public String getDesignVersion();
}




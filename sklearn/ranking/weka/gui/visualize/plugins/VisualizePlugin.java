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
 * VisualizePlugin.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 * Written by Jeffery Grajkowski of the AICML
 *
 */

package weka.gui.visualize.plugins;

import weka.core.FastVector;
import weka.core.Attribute;
import javax.swing.JMenuItem;

/**
 * Interface implemented by classes loaded dynamically to
 * visualize classifier results in the explorer.
 *
 * @author Jeffery Grajkowski (grajkows@cs.ualberta.ca)
 * @version $Revision: 1.1 $
 */
public interface VisualizePlugin {

  /**
   * Get a JMenu or JMenuItem which contain action listeners
   * that perform the visualization, using some but not
   * necessarily all of the data.  Exceptions thrown because of
   * changes in Weka since compilation need to be caught by
   * the implementer.
   *
   * @see NoClassDefFoundError
   * @see IncompatibleClassChangeError
   *
   * @param  preds predictions
   * @param  classAtt class attribute
   * @return menuitem for opening visualization(s), or null
   *         to indicate no visualization is applicable for the input
   */
  public JMenuItem getVisualizeMenuItem(FastVector preds, Attribute classAtt);

  /**
   * Get the minimum version of Weka, inclusive, the class
   * is designed to work with.  eg: <code>3.5.0</code>
   */
  public String getMinVersion();

  /**
   * Get the maximum version of Weka, exclusive, the class
   * is designed to work with.  eg: <code>3.6.0</code>
   */
  public String getMaxVersion();

  /**
   * Get the specific version of Weka the class is designed for.
   * eg: <code>3.5.1</code>
   */
  public String getDesignVersion();
}




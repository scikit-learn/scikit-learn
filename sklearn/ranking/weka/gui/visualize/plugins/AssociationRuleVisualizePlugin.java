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
 *    AssociationRulesVisualizePlugin.java
 *    Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.visualize.plugins;

import java.util.List;

import javax.swing.JMenuItem;

import weka.associations.AssociationRules;

/**
 * Interface implemented by classes loaded dynamically to
 * visualize association results in the explorer.
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6499 $
 */
public interface AssociationRuleVisualizePlugin {

  /**
   * Get a JMenu or JMenuItem which contain action listeners
   * that perform the visualization of the association rules.
   *
   * @see NoClassDefFoundError
   * @see IncompatibleClassChangeError
   *
   * @param rules       the association rules
   * @param name        the name of the item (in the Explorer's history list)
   * @return menuitem   for opening visualization(s), or null
   *                    to indicate no visualization is applicable for the input
   */
  public JMenuItem getVisualizeMenuItem(AssociationRules rules, String name);
  
}

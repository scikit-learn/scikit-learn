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
 *    AttributePanelListener.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.visualize;

/**
 * Interface for classes that want to listen for Attribute selection
 * changes in the attribute panel
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public interface AttributePanelListener {

  /**
   * Called when the user clicks on an attribute bar
   * @param e the event encapsulating what happened
   */
  void attributeSelectionChange(AttributePanelEvent e);

}

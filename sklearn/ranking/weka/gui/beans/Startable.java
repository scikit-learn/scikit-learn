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
 *    Startable.java
 *    Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

/**
 * Interface to something that is a start point for a flow and
 * can be launched programatically.
 *
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}org)
 * @version $Revision: 5928 $
 * @since 1.0
 */
public interface Startable {
  
  /**
   * Start the flow running
   *
   * @exception Exception if something goes wrong
   */
  void start() throws Exception;
  
  /**
   * Gets a string that describes the start action. The
   * KnowledgeFlow uses this in the popup contextual menu
   * for the component. The string can be proceeded by
   * a '$' character to indicate that the component can't
   * be started at present.
   * 
   * @return a string describing the start action.
   */
  String getStartMessage();
}

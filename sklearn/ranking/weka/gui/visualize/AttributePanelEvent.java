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
 *    AttributePanelEvent.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.visualize;

/**
 * Class encapsulating a change in the AttributePanel's selected x and y
 * attributes.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class AttributePanelEvent {

  /** True if the x selection changed */
  public boolean m_xChange;

  /** True if the y selection changed */
  public boolean m_yChange;

  /** The index for the new attribute */
  public int m_indexVal;

  /**
   * Constructor
   * @param xChange true if a change occured to the x selection
   * @param yChange true if a change occured to the y selection
   * @param indexVal the index of the new attribute
   */
  public AttributePanelEvent(boolean xChange, boolean yChange, int indexVal) {
    m_xChange = xChange;
    m_yChange = yChange;
    m_indexVal = indexVal;
  }
}

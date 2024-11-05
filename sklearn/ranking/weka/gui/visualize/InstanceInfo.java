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
 * InstanceInfo.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.visualize;

import weka.core.Instances;

import java.util.Vector;

/**
 * Interface for JFrames that display instance info.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5014 $
 */
public interface InstanceInfo {
  
  /**
   * Sets the text to display.
   * 
   * @param text	the text to display
   */
  public void setInfoText(String text);
  
  /**
   * Returns the currently displayed info text.
   * 
   * @return		the info text
   */
  public String getInfoText();
  
  /**
   * Sets the underlying data.
   * 
   * @param data	the data of the info text
   */
  public void setInfoData(Vector<Instances> data);
  
  /**
   * Returns the underlying data.
   * 
   * @return		the data of the info text, can be null
   */
  public Vector<Instances> getInfoData();
}

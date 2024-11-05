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
 *    NamedColor.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.treevisualizer;

import java.awt.*;

/**
 * This class contains a color name and the rgb values of that color
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class NamedColor {

  /** The name of the color */
  public String m_name;

  /** The actual color object */
  public Color m_col;
  
  /**
   * @param n The name of the color.
   * @param r The red component of the color.
   * @param g The green component of the color.
   * @param b The blue component of the color.
   */   
  public NamedColor(String n,int r,int g,int b) {
    m_name = n;
    m_col = new Color(r,g,b);
  }
}








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
 *    VisualizePanelEvent.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.visualize;

import weka.core.*;


/** 
 * This event Is fired to a listeners 'userDataEvent' function when
 * The user on the VisualizePanel clicks submit. It contains the attributes
 * selected at the time and a FastVector containing the various shapes
 * that had been drawn into the panel.
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class VisualizePanelEvent {
  
  /** No longer used */
  public static int NONE = 0;
  public static int RECTANGLE = 1;
  public static int OVAL = 2;
  public static int POLYGON = 3;
  public static int LINE = 4;
  public static int VLINE = 5;
  public static int HLINE = 6;
 
  
  /** Contains FastVectors, each one containing the points for an object. */
  private FastVector m_values;
  /** The instances that fall inside the shapes described in m_values. */
  private Instances m_inst;
  /** The instances that fall outside the shapes described in m_values. */
  private Instances m_inst2;
  /** The attribute along the x axis. */
  private int m_attrib1;
  /** The attribute along the y axis. */
  private int m_attrib2;
  

  /**
   * This constructor creates the event with all the parameters set.
   * @param ar The list of shapes.
   * @param i The instances that lie in these shapes.
   * @param i2 The instances that lie outside these shapes.
   * @param at1 The attribute that was along the x axis.
   * @param at2 The attribute that was along the y axis.
   */
  public VisualizePanelEvent(FastVector ar, Instances i, Instances i2, int at1,
			     int at2) {
    m_values = ar;
    m_inst = i;
    m_inst2 = i2;
    m_attrib1 = at1;
    m_attrib2 = at2;
    
    
  }

  /**
   * @return The list of shapes.
   */
  public FastVector getValues() {
    return m_values;
  }
  
  /**
   * @return The instances that lie in the shapes.
   */
  public Instances getInstances1() {
    return m_inst;
  }
  
  /**
   * @return The instances that lie outside the shapes.
   */
  public Instances getInstances2() {
    return m_inst2;
  }

  /**
   * @return The x axis attribute.
   */
  public int getAttribute1() {
    return m_attrib1;
  }
  
  /**
   * @return The y axis attribute.
   */
  public int getAttribute2() {
    return m_attrib2;
  }





}

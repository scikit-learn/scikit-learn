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
 *    Edge.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;

import java.util.*;
import java.awt.*;


/**
 * This class is used in conjunction with the Node class to form a tree 
 * structure.
 * This in particular contains information about an edges in the tree.
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $ 
 */
public class Edge {
  /** The text caption for the edge. */
  private String m_label;

  /** The ID string of the parent Node of this edge (used for consrtuction 
   * purposes). */
  private String m_rsource;

  /** The ID string of the child Node of this edge (used for construction 
   * purposes). */
  private String m_rtarget;

  /** The parent Node of this edge. */
  private Node m_source;

  /** The child Node of this edge. */
  private Node m_target;

  /** The label broken up into lines. */
  private Vector m_lines;

  /** 
   * This constructs an Edge with the specified label 
   * and parent , child serial tags.
   *
   * @param label The text caption for the edge.
   * @param source The ID string for this edges parent.
   * @param target The ID string for this edges child.
   */
  public Edge(String label,String source,String target) {
    m_label = label;
    m_rsource = source;
    m_rtarget = target;
    m_lines = new Vector(3,2);
    breakupLabel();
  }

  
  /**
   * Get the value of label.
   *
   * @return Value of label.
   */
  public String getLabel() {
    
    return m_label;
  }
  
  /**
   * This function is called to break the label of the edge up in to 
   * seperate lines
   */
  private void breakupLabel() {
    int prev = 0,noa;
    for (noa = 0;noa < m_label.length();noa++) {
      if (m_label.charAt(noa) == '\n') {
	m_lines.addElement(m_label.substring(prev,noa));
	prev = noa+1;
      }
    }
    m_lines.addElement(m_label.substring(prev,noa));
  }
  
  /**
   * This will calculate how large a rectangle using the <i>FontMetrics</i>
   * passed that the lines of the label will take up
   *
   * @param f The size information for a particular Font
   * @return A Dimension containing the size and width of the text
   */
  public Dimension stringSize(FontMetrics f) {
    Dimension d = new Dimension();
    int old = 0;
    String s;
    int noa = 0;
    while ((s = getLine(noa)) != null) {
      noa++;
      old = f.stringWidth(s);
      
      if (old > d.width) {
	d.width = old;
      }
    }
    d.height = noa * f.getHeight();
    return d;
  }
 
  /**
   * Returns line number <i>n</i>
   *
   * @param n The number of the line requested
   * @return The string for the line number or NULL if it didn't exist
   */ 
  public String getLine(int n) {
    if (n < m_lines.size()) {
      return (String)m_lines.elementAt(n);
    }
    else {
      return null;
    }
  }
  
  
  /**
   * Get the value of rsource.
   *
   * @return Value of rsource.
   */
  public String getRsource() {
    
    return m_rsource;
  }
  
  /**
   * Set the value of rsource.
   *
   * @param v  Value to assign to rsource.
   */
  public void setRsource(String v) {
    
    m_rsource = v;
  }
  
  
  
  /**
   * Get the value of rtarget.
   *
   * @return Value of rtarget.
   */
  public String getRtarget() {
    
    return m_rtarget;
  }
  
  /**
   * Set the value of rtarget.
   *
   * @param v Value to assign to rtarget.
   */
  public void setRtarget(String v) {
    
    m_rtarget = v;
  }
  
  /**
   * Get the value of source.
   *
   * @return Value of source.
   */
  public Node getSource() {
    
    return m_source;
  }
  
  /**
   * Set the value of source. And then call v.addChild to add the edge to 
   * the Node.
   *
   * @param v  Value to assign to source.
   */
  public void setSource(Node v) {
    
    m_source = v;
    v.addChild(this);
  }
  
  /**
   * Get the value of target.
   *
   * @return Value of target.
   */
  public Node getTarget() {
    
    return m_target;
  }
  
  /**
   * Set the value of target. And then call v.addParent to add the edge to 
   * the Node.
   *
   * @param v Value to assign to target.
   */
  public void setTarget(Node v) {
    
    m_target = v;
    v.setParent(this);
  }
}






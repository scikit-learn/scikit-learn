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
 *    PlaceNode1.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.treevisualizer;

/**
 * This class will place the Nodes of a tree. <p>
 * 
 * It will place these nodes so that they symetrically fill each row. 
 * This is simple to calculate but is not visually nice for most trees.<p>
 *
 * @author Malcolm F Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class PlaceNode1 implements NodePlace {
  /** An array containing the spacing value for each level */
  private double[] m_levels; //contains num of nodes one each level

  /** The number of levels in the tree */ 
  private int m_noLevels;//contains num of levels

  /** An array containing the current node place for each level to place 
   * each node accordingly. */
  private int[] m_levelNode; //contains num of node upto on particular level

  /** The distance between each level. */
  private double m_yRatio; //for quicker running y_ratio is a constant after 

                         //being calculated
  /**
   * Call this function to have each node in the tree starting at 'r' placed 
   * in a visual
   * (not logical, they already are) tree position.
   *
   * @param r The top of the tree.
   */
  public void place(Node r) {
    /* this is the first and most basic algorithm to write
       I will use this as a reference to test the classes 

       this works by counting up the nodes on each level and spacing the
       level evenly so that it is all used
    */

    /* this loop will work by starting at the first node
       and systematically going through all their children from left
       to right.but first it will do a quick pass to find out the number
       of levels there are*/

    //+ 1 so that no nodes are on edge of screen
    m_noLevels = r.getHeight(r,0)+1;
    
    m_yRatio = 1 / (double) m_noLevels;
    
    m_levels = new double[m_noLevels];
    m_levelNode = new int[m_noLevels];
    for (int noa = 0;noa < m_noLevels;noa++) {
      m_levels[noa] = 1;
      m_levelNode[noa] = 0;
    }
    
    setNumOfNodes(r,0);
    
    for (int noa = 0;noa < m_noLevels;noa++) {
      m_levels[noa] = 1 / m_levels[noa];
    }
    
    placer(r,0);
  }

  /**
   * This function finds the number of nodes on each level recursively.
   *
   * @param r The current Node upto.
   * @param l The current level upto.
   */
  private void setNumOfNodes(Node r,int l) {
    Edge e;
    l++;
    
    m_levels[l]++;
    for (int noa = 0;(e = r.getChild(noa)) != null && r.getCVisible();noa++) {
      setNumOfNodes(e.getTarget(),l);
    }
  }
  
  /**
   * This function goes through and sets the position of each node
   *
   * @param r The current node upto.
   * @param l the current level upto.
   */
  private void placer(Node r,int l) {
    Edge e;
    l++;
    m_levelNode[l]++;
    r.setCenter(m_levelNode[l] * m_levels[l]);
    r.setTop(l * m_yRatio);
    for (int noa = 0;(e = r.getChild(noa)) != null && r.getCVisible();noa++) {
      placer(e.getTarget(),l);
    }
  }
}

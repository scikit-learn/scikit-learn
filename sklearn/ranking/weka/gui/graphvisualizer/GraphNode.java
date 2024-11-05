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
 *    GraphNode.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.graphvisualizer;

/**
 * This class represents a node in the Graph.
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 4995 $ - 23 Apr 2003 - Initial version (Ashraf M. Kibriya)
 */

public class GraphNode extends Object implements GraphConstants {
  /** ID and label for the node */
  public String ID, lbl;
  /** The outcomes for the given node */
  public String [] outcomes;
  /** probability table for each outcome given outcomes of parents, if any */
  public double [][] probs;   //probabilities
  /** The x and y position of the node */
  public int x=0, y=0;
  /** The indices of parent nodes */
  public int [] prnts;       //parent nodes
  /** The indices of nodes to which there are edges from this
   * node, plus the type of edge */
  public int [][] edges;
  /**  Type of node. Default is Normal node type */
  public int nodeType=NORMAL;
  
  /**
   *  Constructor
   *
   */
  public GraphNode(String id, String label) {
    ID = id; lbl = label; nodeType=NORMAL;
  }
  
  /**
   *  Constructor
   *
   */
  public GraphNode(String id, String label, int type ) {
    ID = id; lbl = label; nodeType = type;
  }
  
  /**
   *  Returns true if passed in argument is an instance
   *  of GraphNode and is equal to this node.
   *  Implemented to enable the use of contains method
   *  in Vector/FastVector class.
   */
  public boolean equals(Object n) {
    if(n instanceof GraphNode && ((GraphNode) n).ID.equalsIgnoreCase(this.ID)) {
      //System.out.println("returning true, n.ID >"+((GraphNode)n).ID+
      //                   "< this.ID >"+this.ID+"<");
      return true;
    }
    else
      return false;
  }
} // GraphNode

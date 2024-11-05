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
 * KDTreeNode.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.kdtrees;

import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/**
 * A class representing a KDTree node. A node does not explicitly
 * store the instances that it contains. Instead, it only stores 
 * the start and end index of a portion in a master index array. Each
 * node is assigned a portion in the master index array that stores 
 * the indices of the instances that the node contains. Every time a 
 * node is split by the KDTree's contruction method, the instances of 
 * its left child are moved to the left and the instances of its 
 * right child are moved to the right, in the portion of the master 
 * index array belonging to the node. The start and end index in each
 * of its children are then set accordingly within that portion so 
 * that each have their own portion which contains their instances.   
 * P.S.: The master index array is only stored in KDTree class.
 * 
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public class KDTreeNode
  implements Serializable, RevisionHandler {
   
  /** for serialization. */
  private static final long serialVersionUID = -3660396067582792648L;

  /** node number (only for debug). */
  public int m_NodeNumber;

  /** left subtree; contains instances with smaller or equal to split value. */
  public KDTreeNode m_Left = null;

  /** right subtree; contains instances with larger than split value. */
  public KDTreeNode m_Right = null;

  /** value to split on. */
  public double m_SplitValue;

  /** attribute to split on. */
  public int m_SplitDim;

  /**
   * lowest and highest value and width (= high - low) for each
   * dimension.
   */
  public double[][] m_NodeRanges;

  /** 
   * The lo and high bounds of the hyper rectangle described by the
   * node.
   */
  public double[][] m_NodesRectBounds;

  /**
   * The start index of the portion of the master index array, 
   * which stores the indices of the instances/points the node 
   * contains.
   */
  public int m_Start = 0;
  
  /**
   * The end index of the portion of the master index array, 
   * which stores indices of the instances/points the node 
   * contains.
   */
  public int m_End = 0;

  /**
   * Constructor.
   */
  public KDTreeNode() {}

  /**
   * Constructor.
   * 
   * @param nodeNum The node number/id.
   * @param startidx The start index of node's portion 
   * in master index array.
   * @param endidx The start index of node's portion 
   * in master index array.
   * @param nodeRanges The attribute ranges of the 
   * Instances/points contained in this node.
   */
  public KDTreeNode(int nodeNum, int startidx, int endidx, double[][] nodeRanges) {
    m_NodeNumber = nodeNum;
    m_Start = startidx; m_End = endidx;
    m_NodeRanges = nodeRanges;
  }

  /**
   * 
   * @param nodeNum The node number/id.
   * @param startidx The start index of node's portion 
   * in master index array.
   * @param endidx The start index of node's portion 
   * in master index array.
   * @param nodeRanges The attribute ranges of the 
   * Instances/points contained in this node.
   * @param rectBounds The range of the rectangular 
   * region in the point space that this node 
   * represents (points inside this rectangular
   * region can have different range).
   */
  public KDTreeNode(int nodeNum, int startidx, int endidx, double[][] nodeRanges, double[][] rectBounds) {
    m_NodeNumber = nodeNum;
    m_Start = startidx; m_End = endidx;
    m_NodeRanges = nodeRanges;
    m_NodesRectBounds = rectBounds;
  }

  /**
   * Gets the splitting dimension.
   * 
   * @return 		splitting dimension
   */
  public int getSplitDim() {
    return m_SplitDim;
  }

  /**
   * Gets the splitting value.
   * 
   * @return 		splitting value
   */
  public double getSplitValue() {
    return m_SplitValue;
  }

  /**
   * Checks if node is a leaf.
   * 
   * @return 		true if it is a leaf
   */
  public boolean isALeaf() {
    return (m_Left == null);
  }         

  /**
   * Returns the number of Instances 
   * in the rectangular region defined 
   * by this node.
   * @return The number of instances in
   * this KDTreeNode.
   */
  public int numInstances() {
    return (m_End-m_Start+1);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
}

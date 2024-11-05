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
 * BallNode.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.balltrees;

import weka.core.DistanceFunction;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;

/**
 * Class representing a node of a BallTree.
 * 
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5987 $
 */
public class BallNode
  implements Serializable, RevisionHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = -8289151861759883510L;
  
  /**
   * The start index of the portion of the master index array, 
   * which stores the indices of the instances/points the node 
   * contains.
   */
  public int m_Start;
  
  /**
   * The end index of the portion of the master index array, 
   * which stores indices of the instances/points the node 
   * contains.
   */
  public int m_End;
  
  /** The number of instances/points in the node. */
  public int m_NumInstances;
  
  /** The node number/id. */
  public int m_NodeNumber;
  
  /** The attribute that splits this node (not 
   * always used). */
  public int m_SplitAttrib = -1;
  
  /** The value of m_SpiltAttrib that splits this
   * node (not always used).
   */
  public double m_SplitVal = -1;
  
  /** The left child of the node. */
  public BallNode m_Left = null;
  
  /** The right child of the node. */
  public BallNode m_Right = null;
  
  /** 
   * The pivot/centre of the ball. 
   */
  protected Instance m_Pivot;
  
  /** The radius of this ball (hyper sphere). */
  protected double m_Radius;
  
  /**
   * Constructor.
   * @param nodeNumber The node's number/id.
   */
  public BallNode(int nodeNumber) {
    m_NodeNumber = nodeNumber;
  }
  
  /**
   * Creates a new instance of BallNode.
   * @param start The begining index of the portion of
   * the master index array belonging to this node.
   * @param end The end index of the portion of the 
   * master index array belonging to this node. 
   * @param nodeNumber The node's number/id.
   */
  public BallNode(int start, int end, int nodeNumber) {
    m_Start = start;
    m_End = end;
    m_NodeNumber = nodeNumber;
    m_NumInstances = end - start + 1;
  }
  
  /**
   * Creates a new instance of BallNode.
   * @param start The begining index of the portion of
   * the master index array belonging to this node.
   * @param end The end index of the portion of the 
   * master index array belonging to this node. 
   * @param nodeNumber The node's number/id.
   * @param pivot The pivot/centre of the node's ball.
   * @param radius The radius of the node's ball.
   */
  public BallNode(int start, int end, int nodeNumber, Instance pivot, double radius) {
    m_Start = start;
    m_End = end;
    m_NodeNumber = nodeNumber; 
    m_Pivot = pivot;
    m_Radius = radius;
    m_NumInstances = end - start + 1;
  }
  
  /** 
   * Returns true if the node is a leaf node (if
   * both its left and right child are null).
   * @return true if the node is a leaf node.
   */
  public boolean isALeaf() {
    return (m_Left==null && m_Right==null);
  }
  
  /** 
   * Sets the the start and end index of the
   * portion of the master index array that is
   * assigned to this node.  
   * @param start The start index of the 
   * master index array. 
   * @param end The end index of the master
   * indext array. 
   */
  public void setStartEndIndices(int start, int end) {
    m_Start = start;
    m_End = end;
    m_NumInstances = end - start + 1;    
  }

  /**
   * Sets the pivot/centre of this nodes
   * ball.
   * @param pivot The centre/pivot.
   */
  public void setPivot(Instance pivot) {
    m_Pivot = pivot;
  }
  
  /**
   * Returns the pivot/centre of the
   * node's ball.
   * @return The ball pivot/centre.
   */
  public Instance getPivot() {
    return m_Pivot;
  }
  
  /** 
   * Sets the radius of the node's 
   * ball.
   * @param radius The radius of the nodes ball.
   */
  public void setRadius(double radius) {
    m_Radius = radius;
  }
  
  /**
   * Returns the radius of the node's ball.
   * @return Radius of node's ball.
   */
  public double getRadius() {
    return m_Radius;
  }
  
  /** 
   * Returns the number of instances in the
   * hyper-spherical region of this node. 
   * @return The number of instances in the
   * node. 
   */
  public int numInstances() {
    return (m_End-m_Start+1);
  }
  
  /**
   * Calculates the centroid pivot of a node. The node is given
   * in the form of an indices array that contains the 
   * indices of the points inside the node.   
   * @param instList The indices array pointing to the 
   * instances in the node.
   * @param insts The actual instances. The instList
   * points to instances in this object.  
   * @return The calculated centre/pivot of the node.  
   */
  public static Instance calcCentroidPivot(int[] instList, Instances insts) {
    double[] attrVals = new double[insts.numAttributes()];
    
    Instance temp;
    for(int i=0; i<instList.length; i++) {
      temp = insts.instance(instList[i]);
      for(int j=0; j<temp.numValues(); j++) {
        attrVals[j] += temp.valueSparse(j);
      }
    }
    for(int j=0, numInsts=instList.length; j<attrVals.length; j++) {
      attrVals[j] /= numInsts;
    }
    temp = new DenseInstance(1.0, attrVals);
    return temp;
  }
  
  /**
   * Calculates the centroid pivot of a node. The node is given
   * in the form of the portion of an indices array that 
   * contains the indices of the points inside the node.
   * @param start The start index marking the start of 
   * the portion belonging to the node.
   * @param end The end index marking the end of the
   * portion in the indices array that belongs to the node.    
   * @param instList The indices array pointing to the 
   * instances in the node.
   * @param insts The actual instances. The instList
   * points to instances in this object.  
   * @return The calculated centre/pivot of the node.  
   */
  public static Instance calcCentroidPivot(int start, int end, int[] instList, 
                                          Instances insts) {
    double[] attrVals = new double[insts.numAttributes()];
    Instance temp;
    for(int i=start; i<=end; i++) {
      temp = insts.instance(instList[i]);
      for(int j=0; j<temp.numValues(); j++) {
        attrVals[j] += temp.valueSparse(j);
      }
    }
    for(int j=0, numInsts=end-start+1; j<attrVals.length; j++) {
      attrVals[j] /= numInsts;
    }
    
    temp = new DenseInstance(1.0, attrVals);    
    return temp;
  }
  
  /**
   * Calculates the radius of node.
   *  
   * @param instList The indices array containing the indices of the 
   * instances inside the node. 
   * @param insts The actual instances object. instList points to 
   * instances in this object.
   * @param pivot The centre/pivot of the node.
   * @param distanceFunction The distance fuction to use to calculate 
   * the radius. 
   * @return The radius of the node. 
   * @throws Exception If there is some problem in calculating the 
   * radius. 
   */
  public static double calcRadius(int[] instList, Instances insts,Instance pivot, 
                                 DistanceFunction distanceFunction) 
                                                  throws Exception {
    return calcRadius(0, instList.length-1, instList, insts, 
                      pivot, distanceFunction);
  }
  
  /**
   * Calculates the radius of a node.
   * 
   * @param start The start index of the portion in indices array 
   * that belongs to the node.
   * @param end The end index of the portion in indices array 
   * that belongs to the node. 
   * @param instList The indices array holding indices of 
   * instances. 
   * @param insts The actual instances. instList points to 
   * instances in this object. 
   * @param pivot The centre/pivot of the node. 
   * @param distanceFunction The distance function to use to 
   * calculate the radius. 
   * @return The radius of the node. 
   * @throws Exception If there is some problem calculating the 
   * radius. 
   */
  public static double calcRadius(int start, int end, int[] instList, 
                                 Instances insts, Instance pivot, 
                                 DistanceFunction distanceFunction) 
                                                             throws Exception {
    double radius = Double.NEGATIVE_INFINITY;
    
    for(int i=start; i<=end; i++) {
      double dist = distanceFunction.distance(pivot, 
                                              insts.instance(instList[i]), Double.POSITIVE_INFINITY);
      
      if(dist>radius)
        radius = dist;
    }
    return Math.sqrt(radius);
  }
 
  /**
   * Calculates the centroid pivot of a node based on its
   * two child nodes (if merging two nodes).
   * @param child1 The first child of the node.
   * @param child2 The second child of the node.
   * @param insts The set of instances on which 
   * the tree is (or is to be) built.
   * @return The centre/pivot of the node.
   * @throws Exception If there is some problem calculating
   * the pivot.
   */
  public static Instance calcPivot(BallNode child1, BallNode child2, 
                                         Instances insts)  throws Exception {
    Instance p1 = child1.getPivot(), p2 = child2.getPivot();
    double[] attrVals = new double[p1.numAttributes()];
    
    for(int j=0; j<attrVals.length; j++) {
      attrVals[j] += p1.value(j);
      attrVals[j] += p2.value(j);
      attrVals[j] /= 2D;
    }
    
    p1 = new DenseInstance(1.0, attrVals);
    return p1;
  }

  /**
   * Calculates the radius of a node based on its two 
   * child nodes (if merging two nodes).
   * @param child1 The first child of the node.
   * @param child2 The second child of the node.
   * @param pivot The centre/pivot of the node. 
   * @param distanceFunction The distance function to 
   * use to calculate the radius
   * @return The radius of the node. 
   * @throws Exception If there is some problem 
   * in calculating the radius.
   */
  public static double calcRadius(BallNode child1, BallNode child2, 
                                  Instance pivot, 
                                  DistanceFunction distanceFunction) 
                                                             throws Exception {
    Instance p1 = child1.getPivot(), p2 = child2.getPivot();                                                               
    
    double radius = child1.getRadius() + distanceFunction.distance(p1, p2) + 
                    child2.getRadius();
    
    return radius/2;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
}

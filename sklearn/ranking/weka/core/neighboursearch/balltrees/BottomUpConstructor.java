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
 * BottomUpConstructor.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.balltrees;

import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.ArrayList;

/**
 <!-- globalinfo-start -->
 * The class that constructs a ball tree bottom up.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;techreport{Omohundro1989,
 *    author = {Stephen M. Omohundro},
 *    institution = {International Computer Science Institute},
 *    month = {December},
 *    number = {TR-89-063},
 *    title = {Five Balltree Construction Algorithms},
 *    year = {1989}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N &lt;value&gt;
 *  Set maximum number of instances in a leaf node
 *  (default: 40)</pre>
 * 
 * <pre> -R
 *  Set internal nodes' radius to the sum 
 *  of the child balls radii. So that it 
 * contains the child balls.</pre>
 * 
 <!-- options-end --> 
 *
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5987 $
 */
public class BottomUpConstructor
  extends BallTreeConstructor 
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = 5864250777657707687L;

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return "The class that constructs a ball tree bottom up.";
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing detailed
   * information about the technical background of this class, e.g., paper
   * reference or book this class is based on.
   * 
   * @return 		the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation result;

    result = new TechnicalInformation(Type.TECHREPORT);
    result.setValue(Field.AUTHOR, "Stephen M. Omohundro");
    result.setValue(Field.YEAR, "1989");
    result.setValue(Field.TITLE, "Five Balltree Construction Algorithms");
    result.setValue(Field.MONTH, "December");
    result.setValue(Field.NUMBER, "TR-89-063");
    result.setValue(Field.INSTITUTION, "International Computer Science Institute");

    return result;
  }

  /**
   * Creates a new instance of BottomUpConstructor.
   */
  public BottomUpConstructor() {
  }

  /**
   * Builds the ball tree bottom up. 
   * @return The root node of the tree. 
   * @throws Exception If there is problem building
   * the tree.
   */
  public BallNode buildTree() throws Exception {
    ArrayList<TempNode> list = new ArrayList<TempNode>();
    
    for(int i=0; i<m_InstList.length; i++) {
      TempNode n = new TempNode();
      n.points = new int[1]; n.points[0] = m_InstList[i];
      n.anchor = m_Instances.instance(m_InstList[i]);
      n.radius = 0.0;
      list.add(n);
    }
    
    return mergeNodes(list, 0, m_InstList.length-1, m_InstList);
  }

  /**
   * Merges nodes into one top node.
   *  
   * @param list List of bottom most nodes (the actual
   * instances).
   * @param startIdx The index marking the start of 
   * the portion of master index array containing 
   * instances that need to be merged. 
   * @param endIdx The index marking the end of 
   * the portion of master index array containing 
   * instances that need to be merged.
   * @param instList The master index array.
   * @return The root node of the tree resulting
   * from merging of bottom most nodes.
   * @throws Exception If there is some problem
   * merging the nodes. 
   */
  protected BallNode mergeNodes(ArrayList<TempNode> list, int startIdx, int endIdx, 
                                int[] instList) throws Exception {
    double minRadius=Double.POSITIVE_INFINITY, tmpRadius;
    Instance pivot, minPivot=null; int min1=-1, min2=-1;
    int [] minInstList=null; int merge=1;
    TempNode parent;
    
    while(list.size() > 1) { //main merging loop
      System.err.print("merge step: "+merge+++"               \r");
      minRadius = Double.POSITIVE_INFINITY;
      min1 = -1; min2 = -1; 
   
      for(int i=0; i<list.size(); i++) {
        TempNode first = (TempNode) list.get(i);
        for(int j=i+1; j<list.size(); j++) {
          TempNode second = (TempNode) list.get(j);
          pivot = calcPivot(first, second, m_Instances);
          tmpRadius = calcRadius(first, second); 
          if(tmpRadius < minRadius) {
            minRadius = tmpRadius; 
            min1=i; min2=j;
            minPivot = pivot;
          }
        }//end for(j)
      }//end for(i)
      parent = new TempNode();
      parent.left  = (TempNode) list.get(min1);
      parent.right = (TempNode) list.get(min2);
      minInstList = new int[parent.left.points.length+parent.right.points.length]; 
      System.arraycopy(parent.left.points, 0, minInstList, 0, parent.left.points.length);
      System.arraycopy(parent.right.points, 0, minInstList, parent.left.points.length, 
          	       parent.right.points.length);
      parent.points = minInstList;
      parent.anchor = minPivot;
      parent.radius = BallNode.calcRadius(parent.points, m_Instances, minPivot, m_DistanceFunction);
      list.remove(min1); list.remove(min2-1);
      list.add(parent);
    }//end while
    System.err.println("");
    TempNode tmpRoot = (TempNode)list.get(0);
    
    if(m_InstList.length != tmpRoot.points.length)
      throw new Exception("Root nodes instance list is of irregular length. " +
                          "Please check code.");
    System.arraycopy(tmpRoot.points, 0, m_InstList, 0, tmpRoot.points.length);

    m_NumNodes = m_MaxDepth = m_NumLeaves = 0;
    tmpRadius = BallNode.calcRadius(instList, m_Instances, tmpRoot.anchor, m_DistanceFunction);    
    BallNode node = makeBallTree(tmpRoot, startIdx, endIdx, instList, 0, tmpRadius); 
    
    return node;    
  }
  
  /**
   * Makes ball tree nodes of temp nodes that were used
   * in the merging process. 
   * @param node The temp root node.
   * @param startidx The index marking the start of the 
   * portion of master index array containing instances 
   * to be merged. 
   * @param endidx The index marking the end of the 
   * portion of master index array containing instances 
   * to be merged. 
   * @param instList The master index array.
   * @param depth The depth of the provided temp node.
   * @param rootRadius The smallest ball enclosing all
   * data points.
   * @return The proper top BallTreeNode. 
   * @throws Exception If there is some problem.
   */
  protected BallNode makeBallTree(TempNode node, int startidx, int endidx, 
                                int[] instList, int depth, final double rootRadius) throws Exception {
    BallNode ball=null;
    Instance pivot;
    
    if(m_MaxDepth < depth)
      m_MaxDepth = depth;
    
    if(node.points.length > m_MaxInstancesInLeaf && 
       (rootRadius==0 ? false : node.radius/rootRadius >= m_MaxRelLeafRadius) && 
       node.left!=null && node.right!=null) { //make an internal node
      ball = new BallNode(
      startidx, endidx, m_NumNodes, 
      (pivot=BallNode.calcCentroidPivot(startidx, endidx, instList, m_Instances)),
      BallNode.calcRadius(startidx, endidx, instList, m_Instances, pivot, 
                          m_DistanceFunction)
      );
      m_NumNodes += 1;
      ball.m_Left = makeBallTree(node.left, startidx, startidx+node.left.points.length-1, instList, depth+1, rootRadius);
      ball.m_Right= makeBallTree(node.right, startidx+node.left.points.length, endidx, instList, depth+1, rootRadius);
    }
    else { //make a leaf node
      ball = new BallNode(startidx, endidx, m_NumNodes,       
      (pivot=BallNode.calcCentroidPivot(startidx, endidx, instList, m_Instances)),
      BallNode.calcRadius(startidx, endidx, instList, m_Instances, pivot, 
                          m_DistanceFunction)
                         );
      m_NumNodes += 1;
      m_NumLeaves++;
    }
    return ball;
  }
  
  /**
   * Adds an instance to the ball tree. 
   * @param node The root node of the tree.
   * @param inst The instance to add to the tree.
   * @return The new master index array after adding the 
   * instance. 
   * @throws Exception Always as BottomUpConstructor
   * does not allow addition of instances after batch 
   * construction. 
   */
  
  public int[] addInstance(BallNode node, Instance inst) throws Exception {
    throw new Exception("BottomUpConstruction method does not allow addition " +
                        "of new Instances.");
  }

  /**
   * Calculates the centroid pivot of a node based on its
   * two child nodes. 
   * @param node1 The first child node.
   * @param node2 The second child node.
   * @param insts The instance on which the tree is to be
   * built.
   * @return The centre/pivot of the node. 
   * @throws Exception If there is some problem calculating 
   * the centre/pivot of the node.
   */
  public Instance calcPivot(TempNode node1, TempNode node2, Instances insts) 
  throws Exception {
    int classIdx = m_Instances.classIndex();
    double[] attrVals = new double[insts.numAttributes()];
    Instance temp;
    double anchr1Ratio = (double)node1.points.length / 
    (node1.points.length+node2.points.length),
    anchr2Ratio = (double)node2.points.length / 
    (node1.points.length+node2.points.length);                         
    for(int k=0; k<node1.anchor.numValues(); k++) {
      if(node1.anchor.index(k)==classIdx)
	continue;
      attrVals[k] += node1.anchor.valueSparse(k)*anchr1Ratio;
    }
    for(int k=0; k<node2.anchor.numValues(); k++) {
      if(node2.anchor.index(k)==classIdx)
	continue;
      attrVals[k] += node2.anchor.valueSparse(k)*anchr2Ratio;
    }
    temp = new DenseInstance(1.0, attrVals);
    return temp;
  }
  
  /**
   * Calculates the radius of a node based on its two
   * child nodes. 
   * @param n1 The first child node. 
   * @param n2 The second child node.
   * @return The calculated radius of the the node. 
   * @throws Exception If there is some problem 
   * in calculating the radius. 
   */
  public double calcRadius(TempNode n1, TempNode n2) throws Exception {
    Instance a1 = n1.anchor, a2 = n2.anchor;
    double radius = n1.radius + m_DistanceFunction.distance(a1, a2) + n2.radius;
    return radius/2;
  }

  /** 
   * Temp class to represent either a leaf node or an internal node. Should only 
   * have two children (could be the case one child is an instance and the 
   * other another node).
   *
   * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
   * @version $Revision: 5987 $
   */
  protected class TempNode
    implements RevisionHandler {
    
    /** The centre/pivot of the node. */
    Instance anchor;
    /** The radius of the node. */
    double radius;
    /** Indices of the points in the node. */
    int [] points;
    /** The node's left child. */
    TempNode left = null;
    /** The node's right child. */
    TempNode right = null;
    
    /** 
     * Prints the node.
     * @return The node as a string.
     */
    public String toString() {
      StringBuffer bf = new StringBuffer();
      bf.append("p: ");
      for(int i=0; i<points.length; i++) 
        if(i!=0)
          bf.append(", "+points[i]);
        else
          bf.append(""+points[i]);
      return bf.toString();
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
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
}

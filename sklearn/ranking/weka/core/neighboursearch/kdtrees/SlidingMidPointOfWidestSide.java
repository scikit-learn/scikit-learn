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
 * SlidingMidPointOfWidestSide.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.kdtrees;

import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

/**
 <!-- globalinfo-start -->
 * The class that splits a node into two based on the midpoint value of the dimension in which the node's rectangle is widest. If after splitting one side is empty then it is slided towards the non-empty side until there is at least one point on the empty side.<br/>
 * <br/>
 * For more information see also:<br/>
 * <br/>
 * David M. Mount (2006). ANN Programming Manual. College Park, MD, USA.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;manual{Mount2006,
 *    address = {College Park, MD, USA},
 *    author = {David M. Mount},
 *    organization = {Department of Computer Science, University of Maryland},
 *    title = {ANN Programming Manual},
 *    year = {2006},
 *    HTTP = {Available from http://www.cs.umd.edu/\~mount/ANN/}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 <!-- options-end -->
 *
 * @author  Ashraf M. Kibriya (amk14@waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class SlidingMidPointOfWidestSide
  extends KDTreeNodeSplitter 
  implements TechnicalInformationHandler {

  /** for serialization. */
  private static final long serialVersionUID = 852857628205680562L;

  /** The floating point error to tolerate in finding the widest 
   * rectangular side. */
  protected static double ERR = 0.001;

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "The class that splits a node into two based on the midpoint value of "
      + "the dimension in which the node's rectangle is widest. If after "
      + "splitting one side is empty then it is slided towards the non-empty "
      + "side until there is at least one point on the empty side.\n\n"
      + "For more information see also:\n\n"
      + getTechnicalInformation().toString();
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

    result = new TechnicalInformation(Type.MANUAL);
    result.setValue(Field.AUTHOR, "David M. Mount");
    result.setValue(Field.YEAR, "2006");
    result.setValue(Field.TITLE, "ANN Programming Manual");
    result.setValue(Field.ORGANIZATION, "Department of Computer Science, University of Maryland");
    result.setValue(Field.ADDRESS,
        "College Park, MD, USA");
    result.setValue(Field.HTTP,
        "Available from http://www.cs.umd.edu/~mount/ANN/");

    return result;
  }

  /** 
   * Splits a node into two based on the midpoint value of the dimension 
   * in which the node's rectangle is widest. If after splitting one side
   * is empty then it is slided towards the non-empty side until there is 
   * at least one point on the empty side. The two nodes created after the 
   * whole splitting are correctly initialised. And, node.left and 
   * node.right are set appropriately.  
   * @param node The node to split.
   * @param numNodesCreated The number of nodes that so far have been
   * created for the tree, so that the newly created nodes are 
   * assigned correct/meaningful node numbers/ids.
   * @param nodeRanges The attributes' range for the points inside
   * the node that is to be split.
   * @param universe The attributes' range for the whole 
   * point-space.
   * @throws Exception If there is some problem in splitting the
   * given node.
   */
  public void splitNode(KDTreeNode node, int numNodesCreated,
      double[][] nodeRanges, double[][] universe) throws Exception {

    correctlyInitialized();

    if (node.m_NodesRectBounds == null) {
      node.m_NodesRectBounds = new double[2][node.m_NodeRanges.length];
      for (int i = 0; i < node.m_NodeRanges.length; i++) {
        node.m_NodesRectBounds[MIN][i] = node.m_NodeRanges[i][MIN];
        node.m_NodesRectBounds[MAX][i] = node.m_NodeRanges[i][MAX];
      }
    }

    // finding widest side of the hyper rectangle
    double maxRectWidth = Double.NEGATIVE_INFINITY, maxPtWidth = Double.NEGATIVE_INFINITY, tempval;
    int splitDim = -1, classIdx = m_Instances.classIndex();

    for (int i = 0; i < node.m_NodesRectBounds[0].length; i++) {
      if (i == classIdx)
        continue;
      tempval = node.m_NodesRectBounds[MAX][i] - node.m_NodesRectBounds[MIN][i];
      if (m_NormalizeNodeWidth) {
        tempval = tempval / universe[i][WIDTH];
      }
      if (tempval > maxRectWidth && node.m_NodeRanges[i][WIDTH] > 0.0)
        maxRectWidth = tempval;
    }

    for (int i = 0; i < node.m_NodesRectBounds[0].length; i++) {
      if (i == classIdx)
        continue;
      tempval = node.m_NodesRectBounds[MAX][i] - node.m_NodesRectBounds[MIN][i];
      if (m_NormalizeNodeWidth) {
        tempval = tempval / universe[i][WIDTH];
      }
      if (tempval >= maxRectWidth * (1 - ERR)
          && node.m_NodeRanges[i][WIDTH] > 0.0) {
        if (node.m_NodeRanges[i][WIDTH] > maxPtWidth) {
          maxPtWidth = node.m_NodeRanges[i][WIDTH];
          if (m_NormalizeNodeWidth)
            maxPtWidth = maxPtWidth / universe[i][WIDTH];
          splitDim = i;
        }
      }
    }

    double splitVal = node.m_NodesRectBounds[MIN][splitDim]
        + (node.m_NodesRectBounds[MAX][splitDim] - node.m_NodesRectBounds[MIN][splitDim])
        * 0.5;
    // might want to try to slide it further to contain more than one point on
    // the
    // side that is resulting empty
    if (splitVal < node.m_NodeRanges[splitDim][MIN])
      splitVal = node.m_NodeRanges[splitDim][MIN];
    else if (splitVal >= node.m_NodeRanges[splitDim][MAX])
      splitVal = node.m_NodeRanges[splitDim][MAX]
          - node.m_NodeRanges[splitDim][WIDTH] * 0.001;

    int rightStart = rearrangePoints(m_InstList, node.m_Start, node.m_End,
        splitDim, splitVal);

    if (rightStart == node.m_Start || rightStart > node.m_End) {
      if (rightStart == node.m_Start)
        throw new Exception("Left child is empty in node " + node.m_NodeNumber
            + ". Not possible with "
            + "SlidingMidPointofWidestSide splitting method. Please "
            + "check code.");
      else
        throw new Exception("Right child is empty in node " + node.m_NodeNumber
            + ". Not possible with "
            + "SlidingMidPointofWidestSide splitting method. Please "
            + "check code.");
    }

    node.m_SplitDim = splitDim;
    node.m_SplitValue = splitVal;

    double[][] widths = new double[2][node.m_NodesRectBounds[0].length];

    System.arraycopy(node.m_NodesRectBounds[MIN], 0, widths[MIN], 0,
        node.m_NodesRectBounds[MIN].length);
    System.arraycopy(node.m_NodesRectBounds[MAX], 0, widths[MAX], 0,
        node.m_NodesRectBounds[MAX].length);
    widths[MAX][splitDim] = splitVal;

    node.m_Left = new KDTreeNode(numNodesCreated + 1, node.m_Start,
        rightStart - 1, m_EuclideanDistance.initializeRanges(m_InstList,
            node.m_Start, rightStart - 1), widths);

    widths = new double[2][node.m_NodesRectBounds[0].length];
    System.arraycopy(node.m_NodesRectBounds[MIN], 0, widths[MIN], 0,
        node.m_NodesRectBounds[MIN].length);
    System.arraycopy(node.m_NodesRectBounds[MAX], 0, widths[MAX], 0,
        node.m_NodesRectBounds[MAX].length);
    widths[MIN][splitDim] = splitVal;

    node.m_Right = new KDTreeNode(numNodesCreated + 2, rightStart, node.m_End,
        m_EuclideanDistance.initializeRanges(m_InstList, rightStart, node.m_End), widths);
  }
  
  /** 
   * Re-arranges the indices array such that the points <= to the splitVal 
   * are on the left of the array and those > the splitVal are on the right.
   * 
   * @param indices The master index array.
   * @param startidx The begining index of portion of indices that needs 
   * re-arranging. 
   * @param endidx The end index of portion of indices that needs 
   * re-arranging. 
   * @param splitDim The split dimension/attribute.
   * @param splitVal The split value.
   * @return The startIdx of the points > the splitVal (the points 
   * belonging to the right child of the node).
   */
  protected int rearrangePoints(int[] indices, final int startidx,
      final int endidx, final int splitDim, final double splitVal) {

    int tmp, left = startidx - 1;
    for (int i = startidx; i <= endidx; i++) {
      if (m_EuclideanDistance.valueIsSmallerEqual(m_Instances
          .instance(indices[i]), splitDim, splitVal)) {
        left++;
        tmp = indices[left];
        indices[left] = indices[i];
        indices[i] = tmp;
      }// end valueIsSmallerEqual
    }// endfor
    return left + 1;
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

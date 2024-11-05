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
 * MidPointOfWidestDimension.java
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
 * The class that splits a KDTree node based on the midpoint value of a dimension in which the node's points have the widest spread.<br/>
 * <br/>
 * For more information see also:<br/>
 * <br/>
 * Andrew Moore (1991). A tutorial on kd-trees.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;techreport{Moore1991,
 *    author = {Andrew Moore},
 *    booktitle = {University of Cambridge Computer Laboratory Technical Report No. 209},
 *    howpublished = {Extract from PhD Thesis},
 *    title = {A tutorial on kd-trees},
 *    year = {1991},
 *    HTTP = {http://www.autonlab.org/autonweb/14665.html}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 <!-- options-end -->
 *
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public class MidPointOfWidestDimension
  extends KDTreeNodeSplitter 
  implements TechnicalInformationHandler {

  /** for serialization. */
  private static final long serialVersionUID = -7617277960046591906L;

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "The class that splits a KDTree node based on the midpoint value of "
      + "a dimension in which the node's points have the widest spread.\n\n"
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

    result = new TechnicalInformation(Type.TECHREPORT);
    result.setValue(Field.AUTHOR, "Andrew Moore");
    result.setValue(Field.YEAR, "1991");
    result.setValue(Field.TITLE, "A tutorial on kd-trees");
    result.setValue(Field.HOWPUBLISHED, "Extract from PhD Thesis");
    result.setValue(Field.BOOKTITLE, "University of Cambridge Computer Laboratory Technical Report No. 209");
    result.setValue(Field.HTTP, "http://www.autonlab.org/autonweb/14665.html");

    return result;
  }
  
  /** 
   * Splits a node into two based on the midpoint value of the dimension 
   * in which the points have the widest spread. After splitting two 
   * new nodes are created and correctly initialised. And, node.left 
   * and node.right are set appropriately.  
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

    int splitDim = widestDim(nodeRanges, universe);

    double splitVal = m_EuclideanDistance.getMiddle(nodeRanges[splitDim]);

    int rightStart = rearrangePoints(m_InstList, node.m_Start, node.m_End,
        splitDim, splitVal);

    if (rightStart == node.m_Start || rightStart > node.m_End) {
      if (rightStart == node.m_Start)
        throw new Exception("Left child is empty in node " 
                            + node.m_NodeNumber + 
                            ". Not possible with " + 
                            "MidPointofWidestDim splitting method. Please " + 
                            "check code.");
      else
        throw new Exception("Right child is empty in node " + node.m_NodeNumber + 
                            ". Not possible with " + 
                            "MidPointofWidestDim splitting method. Please " + 
                            "check code.");
    }
    
    node.m_SplitDim = splitDim;
    node.m_SplitValue = splitVal;
    node.m_Left = new KDTreeNode(numNodesCreated + 1, node.m_Start,
        rightStart - 1, m_EuclideanDistance.initializeRanges(m_InstList,
            node.m_Start, rightStart - 1));
    node.m_Right = new KDTreeNode(numNodesCreated + 2, rightStart, node.m_End,
        m_EuclideanDistance
            .initializeRanges(m_InstList, rightStart, node.m_End));	
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
  protected int rearrangePoints(int[] indices, final int startidx, final int endidx,
      			      final int splitDim, final double splitVal) {
    
    int tmp, left = startidx - 1;
    for (int i = startidx; i <= endidx; i++) {
      if (m_EuclideanDistance.valueIsSmallerEqual(m_Instances
          .instance(indices[i]), splitDim, splitVal)) {
        left++;
        tmp = indices[left];
        indices[left] = indices[i];
        indices[i] = tmp;
      }//end if valueIsSmallerEqual
    }//end for
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

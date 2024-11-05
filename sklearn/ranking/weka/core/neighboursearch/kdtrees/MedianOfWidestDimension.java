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
 * MedianOfWidestDimension.java
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
 * The class that splits a KDTree node based on the median value of a dimension in which the node's points have the widest spread.<br/>
 * <br/>
 * For more information see also:<br/>
 * <br/>
 * Jerome H. Friedman, Jon Luis Bentley, Raphael Ari Finkel (1977). An Algorithm for Finding Best Matches in Logarithmic Expected Time. ACM Transactions on Mathematics Software. 3(3):209-226.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Friedman1977,
 *    author = {Jerome H. Friedman and Jon Luis Bentley and Raphael Ari Finkel},
 *    journal = {ACM Transactions on Mathematics Software},
 *    month = {September},
 *    number = {3},
 *    pages = {209-226},
 *    title = {An Algorithm for Finding Best Matches in Logarithmic Expected Time},
 *    volume = {3},
 *    year = {1977}
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
public class MedianOfWidestDimension
  extends KDTreeNodeSplitter 
  implements TechnicalInformationHandler {

  /** for serialization. */
  private static final long serialVersionUID = 1383443320160540663L;

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "The class that splits a KDTree node based on the median value of "
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

    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Jerome H. Friedman and Jon Luis Bentley and Raphael Ari Finkel");
    result.setValue(Field.YEAR, "1977");
    result.setValue(Field.TITLE, "An Algorithm for Finding Best Matches in Logarithmic Expected Time");
    result.setValue(Field.JOURNAL, "ACM Transactions on Mathematics Software");
    result.setValue(Field.PAGES, "209-226");
    result.setValue(Field.MONTH, "September");
    result.setValue(Field.VOLUME, "3");
    result.setValue(Field.NUMBER, "3");

    return result;
  }
  
  /** 
   * Splits a node into two based on the median value of the dimension 
   * in which the points have the widest spread. After splitting two 
   * new nodes are created and correctly initialised. And, node.left 
   * and node.right are set appropriately.
   * 
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
    
    //In this case median is defined to be either the middle value (in case of
    //odd number of values) or the left of the two middle values (in case of 
    //even number of values).
    int medianIdxIdx = node.m_Start + (node.m_End-node.m_Start)/2;
    //the following finds the median and also re-arranges the array so all 
    //elements to the left are < median and those to the right are > median.
    int medianIdx = select(splitDim, m_InstList, node.m_Start, node.m_End, (node.m_End-node.m_Start)/2+1);
    
    
    node.m_SplitDim = splitDim;
    node.m_SplitValue = m_Instances.instance(m_InstList[medianIdx]).value(splitDim);
    
    node.m_Left  = new KDTreeNode(numNodesCreated+1, node.m_Start, medianIdxIdx,
	m_EuclideanDistance.initializeRanges(m_InstList, node.m_Start, medianIdxIdx));
    node.m_Right = new KDTreeNode(numNodesCreated+2, medianIdxIdx+1, node.m_End,
	m_EuclideanDistance.initializeRanges(m_InstList, medianIdxIdx+1, node.m_End));	
  }
  
  /**
   * Partitions the instances around a pivot. Used by quicksort and
   * kthSmallestValue.
   *
   * @param attIdx The attribution/dimension based on which the 
   * instances should be partitioned.
   * @param index The master index array containing indices of the 
   * instances.
   * @param l The begining index of the portion of master index 
   * array that should be partitioned. 
   * @param r The end index of the portion of master index array 
   * that should be partitioned.
   * @return the index of the middle element
   */
  protected int partition(int attIdx, int[] index, int l, int r) {
    
    double pivot = m_Instances.instance(index[(l + r) / 2]).value(attIdx);
    int help;

    while (l < r) {
      while ((m_Instances.instance(index[l]).value(attIdx) < pivot) && (l < r)) {
        l++;
      }
      while ((m_Instances.instance(index[r]).value(attIdx) > pivot) && (l < r)) {
        r--;
      }
      if (l < r) {
        help = index[l];
        index[l] = index[r];
        index[r] = help;
        l++;
        r--;
      }
    }
    if ((l == r) && (m_Instances.instance(index[r]).value(attIdx) > pivot)) {
      r--;
    } 

    return r;
  }

  /**
   * Implements computation of the kth-smallest element according
   * to Manber's "Introduction to Algorithms".
   *
   * @param attIdx The dimension/attribute of the instances in 
   * which to find the kth-smallest element.
   * @param indices The master index array containing indices of 
   * the instances.
   * @param left The begining index of the portion of the master 
   * index array in which to find the kth-smallest element.
   * @param right The end index of the portion of the master index 
   * array in which to find the kth-smallest element.
   * @param k The value of k
   * @return The index of the kth-smallest element
   */
  public int select(int attIdx, int[] indices, int left, int right, int k) {
    
    if (left == right) {
      return left;
    } else {
      int middle = partition(attIdx, indices, left, right);
      if ((middle - left + 1) >= k) {
        return select(attIdx, indices, left, middle, k);
      } else {
        return select(attIdx, indices, middle + 1, right, k - (middle - left + 1));
      }
    }
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

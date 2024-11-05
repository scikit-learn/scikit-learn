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
 * KMeansInpiredMethod.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.kdtrees;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

/**
 <!-- globalinfo-start -->
 * The class that splits a node into two such that the overall sum of squared distances of points to their centres on both sides of the (axis-parallel) splitting plane is minimum.<br/>
 * <br/>
 * For more information see also:<br/>
 * <br/>
 * Ashraf Masood Kibriya (2007). Fast Algorithms for Nearest Neighbour Search. Hamilton, New Zealand.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;mastersthesis{Kibriya2007,
 *    address = {Hamilton, New Zealand},
 *    author = {Ashraf Masood Kibriya},
 *    school = {Department of Computer Science, School of Computing and Mathematical Sciences, University of Waikato},
 *    title = {Fast Algorithms for Nearest Neighbour Search},
 *    year = {2007}
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
public class KMeansInpiredMethod
  extends KDTreeNodeSplitter
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = -866783749124714304L;

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "The class that splits a node into two such that the overall sum "
      + "of squared distances of points to their centres on both sides " 
      + "of the (axis-parallel) splitting plane is minimum.\n\n"
      + "For more information see also:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing detailed
   * information about the technical background of this class, e.g., paper
   * reference or book this class is based on.
   * 
   * @return		the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation result;

    result = new TechnicalInformation(Type.MASTERSTHESIS);
    result.setValue(Field.AUTHOR, "Ashraf Masood Kibriya");
    result.setValue(Field.TITLE, "Fast Algorithms for Nearest Neighbour Search");
    result.setValue(Field.YEAR, "2007");
    result.setValue(Field.SCHOOL, "Department of Computer Science, School of Computing and Mathematical Sciences, University of Waikato");
    result.setValue(Field.ADDRESS, "Hamilton, New Zealand");

    return result;
  }

  /** 
   * Splits a node into two such that the overall sum of squared distances 
   * of points to their centres on both sides of the (axis-parallel) 
   * splitting plane is minimum. The two nodes created after the whole 
   * splitting are correctly initialised. And, node.left and node.right 
   * are set appropriately.
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

    int splitDim = -1;
    double splitVal = Double.NEGATIVE_INFINITY;

    double leftAttSum[] = new double[m_Instances.numAttributes()], 
           rightAttSum[] = new double[m_Instances.numAttributes()], 
           leftAttSqSum[] = new double[m_Instances.numAttributes()], 
           rightAttSqSum[] = new double[m_Instances.numAttributes()], 
           rightSqMean, leftSqMean, leftSqSum, rightSqSum, 
           minSum = Double.POSITIVE_INFINITY, val;

    for (int dim = 0; dim < m_Instances.numAttributes(); dim++) {
      // m_MaxRelativeWidth in KDTree ensure there'll be atleast one dim with
      // width > 0.0
      if (node.m_NodeRanges[dim][WIDTH] == 0.0
          || dim == m_Instances.classIndex())
        continue;

      quickSort(m_Instances, m_InstList, dim, node.m_Start, node.m_End);

      for (int i = node.m_Start; i <= node.m_End; i++) {
        for (int j = 0; j < m_Instances.numAttributes(); j++) {
          if (j == m_Instances.classIndex())
            continue;
          val = m_Instances.instance(m_InstList[i]).value(j);
          if (m_NormalizeNodeWidth) {
            if (Double.isNaN(universe[j][MIN])
                || universe[j][MIN] == universe[j][MAX])
              val = 0.0;
            else
              val = ((val - universe[j][MIN]) / universe[j][WIDTH]); // normalizing
                                                                      // value
          }
          if (i == node.m_Start) {
            leftAttSum[j] = rightAttSum[j] = leftAttSqSum[j] = rightAttSqSum[j] = 0.0;
          }
          rightAttSum[j] += val;
          rightAttSqSum[j] += val * val;
        }
      }

      for (int i = node.m_Start; i <= node.m_End - 1; i++) {
        Instance inst = m_Instances.instance(m_InstList[i]);
        leftSqSum = rightSqSum = 0.0;
        for (int j = 0; j < m_Instances.numAttributes(); j++) {
          if (j == m_Instances.classIndex())
            continue;
          val = inst.value(j);

          if (m_NormalizeNodeWidth) {
            if (Double.isNaN(universe[j][MIN])
                || universe[j][MIN] == universe[j][MAX])
              val = 0.0;
            else
              val = ((val - universe[j][MIN]) / universe[j][WIDTH]); // normalizing
                                                                      // value
          }

          leftAttSum[j] += val;
          rightAttSum[j] -= val;
          leftAttSqSum[j] += val * val;
          rightAttSqSum[j] -= val * val;
          leftSqMean = leftAttSum[j] / (i - node.m_Start + 1);
          leftSqMean *= leftSqMean;
          rightSqMean = rightAttSum[j] / (node.m_End - i);
          rightSqMean *= rightSqMean;

          leftSqSum += leftAttSqSum[j] - (i - node.m_Start + 1) * leftSqMean;
          rightSqSum += rightAttSqSum[j] - (node.m_End - i) * rightSqMean;
        }

        if (minSum > (leftSqSum + rightSqSum)) {
          minSum = leftSqSum + rightSqSum;

          if (i < node.m_End)
            splitVal = (m_Instances.instance(m_InstList[i]).value(dim) + m_Instances
                .instance(m_InstList[i + 1]).value(dim)) / 2;
          else
            splitVal = m_Instances.instance(m_InstList[i]).value(dim);

          splitDim = dim;
        }
      }// end for instance i
    }// end for attribute dim

    int rightStart = rearrangePoints(m_InstList, node.m_Start, node.m_End,
        splitDim, splitVal);

    if (rightStart == node.m_Start || rightStart > node.m_End) {
      System.out.println("node.m_Start: " + node.m_Start + " node.m_End: "
          + node.m_End + " splitDim: " + splitDim + " splitVal: " + splitVal
          + " node.min: " + node.m_NodeRanges[splitDim][MIN] + " node.max: "
          + node.m_NodeRanges[splitDim][MAX] + " node.numInstances: "
          + node.numInstances());

      if (rightStart == node.m_Start)
        throw new Exception("Left child is empty in node " + node.m_NodeNumber
            + ". Not possible with "
            + "KMeanInspiredMethod splitting method. Please " + "check code.");
      else
        throw new Exception("Right child is empty in node " + node.m_NodeNumber
            + ". Not possible with "
            + "KMeansInspiredMethod splitting method. Please " + "check code.");
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
   * Partitions the instances around a pivot. Used by quicksort and
   * kthSmallestValue.
   *
   * @param insts	The instances on which the tree is (or is 
   * to be) built.
   * @param index The master index array containing indices 
   * of the instances.
   * @param attidx The attribution/dimension based on which
   * the instances should be partitioned.
   * @param l	The begining index of the portion of master index 
   * array that should be partitioned. 
   * @param r	The end index of the portion of master index array 
   * that should be partitioned.
   * @return the index of the middle element
   */
  protected static int partition(Instances insts, int[] index, int attidx, int l, int r) {
    
    double pivot = insts.instance(index[(l + r) / 2]).value(attidx);
    int help;

    while (l < r) {
      while ((insts.instance(index[l]).value(attidx) < pivot) && (l < r)) {
        l++;
      }
      while ((insts.instance(index[r]).value(attidx) > pivot) && (l < r)) {
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
    if ((l == r) && (insts.instance(index[r]).value(attidx) > pivot)) {
      r--;
    } 

    return r;
  }
  
  /**
   * Sorts the instances according to the given attribute/dimension.
   * The sorting is done on the master index array and not on the
   * actual instances object.
   * 
   * @param insts The instances on which the tree is (or is 
   * to be) built.
   * @param indices The master index array containing indices 
   * of the instances.
   * @param attidx The dimension/attribute based on which 
   * the instances should be sorted.
   * @param left The begining index of the portion of the master 
   * index array that needs to be sorted.
   * @param right The end index of the portion of the master index 
   * array that needs to be sorted.
   */
  protected static void quickSort(Instances insts, int[] indices, int attidx, int left, int right) {

    if (left < right) {
      int middle = partition(insts, indices, attidx, left, right);
      quickSort(insts, indices, attidx, left, middle);
      quickSort(insts, indices, attidx, middle + 1, right);
    }
  }  

  /**
   * Method to validate the sorting done by quickSort().
   * 
   * @param insts The instances on which the tree is (or is 
   * to be) built.
   * @param indices The master index array containing indices 
   * of the instances.
   * @param attidx The dimension/attribute based on which 
   * the instances should be sorted.
   * @param start The start of the portion in master index
   * array that needs to be sorted.
   * @param end The end of the portion in master index 
   * array that needs to be sorted.
   * @throws Exception If the indices of the instances 
   * are not in sorted order.
   */
  private static void checkSort(Instances insts, int[] indices, int attidx, 
                               int start, int end) throws Exception {
    for(int i=start+1; i<=end; i++) {
      if( insts.instance(indices[i-1]).value(attidx) > 
          insts.instance(indices[i]).value(attidx) ) {
        System.out.println("value[i-1]: "+insts.instance(indices[i-1]).value(attidx));
        System.out.println("value[i]: "+insts.instance(indices[i]).value(attidx));
        System.out.println("indices[i-1]: "+indices[i-1]);
        System.out.println("indices[i]: "+indices[i]);
        System.out.println("i: "+i);
        if(insts.instance(indices[i-1]).value(attidx) > insts.instance(indices[i]).value(attidx))
          System.out.println("value[i-1] > value[i]");
        
        throw new Exception("Indices not sorted correctly.");
      }//end if
    }
  }
  
  /** 
   * Re-arranges the indices array so that in the portion of the array
   * belonging to the node to be split, the points <= to the splitVal 
   * are on the left of the portion and those > the splitVal are on the right.
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

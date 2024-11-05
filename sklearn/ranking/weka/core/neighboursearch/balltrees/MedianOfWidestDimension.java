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

package weka.core.neighboursearch.balltrees;

import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class that splits a BallNode of a ball tree based on the median value of the widest dimension of the points in the ball. It essentially implements Omohundro's  KD construction algorithm.
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
 * <pre> -N
 *  Normalize dimensions' widths.</pre>
 * 
 <!-- options-end --> 
 * 
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public class MedianOfWidestDimension
  extends BallSplitter 
  implements OptionHandler, TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = 3054842574468790421L;
  
  /** 
   * Should we normalize the widths(ranges) of the dimensions (attributes) 
   * before selecting the widest one. 
   */
  protected boolean m_NormalizeDimWidths = true;
  
  /**
   * Constructor.
   */
  public MedianOfWidestDimension() {
  }
  
  /**
   * Constructor. 
   * @param instList The master index array.
   * @param insts The instances on which the tree
   * is (or is to be) built.
   * @param e The Euclidean distance function to 
   * use for splitting.
   */
  public MedianOfWidestDimension(int[] instList, Instances insts, 
                                 EuclideanDistance e) {
    super(instList, insts, e);
  }
  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Class that splits a BallNode of a ball tree based on the "
      + "median value of the widest dimension of the points in the ball. "
      + "It essentially implements Omohundro's  KD construction algorithm.";
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
   * Splits a ball into two. 
   * @param node The node to split.
   * @param numNodesCreated The number of nodes that so far have been
   * created for the tree, so that the newly created nodes are 
   * assigned correct/meaningful node numbers/ids.
   * @throws Exception If there is some problem in splitting the
   * given node.
   */
  public void splitNode(BallNode node, int numNodesCreated) throws Exception {
    correctlyInitialized();
    //int[] instList = getNodesInstsList(node); 
    double[][] ranges = m_DistanceFunction.initializeRanges(m_Instlist, 
                                                            node.m_Start, 
                                                            node.m_End);
    
    int splitAttrib = widestDim(ranges, m_DistanceFunction.getRanges());
    
    //In this case median is defined to be either the middle value (in case of
    //odd number of values) or the left of the two middle values (in case of 
    //even number of values).
    int medianIdxIdx = node.m_Start + (node.m_End-node.m_Start)/2;
    //the following finds the median and also re-arranges the array so all 
    //elements to the left are < median and those to the right are > median.
    int medianIdx = select(splitAttrib, m_Instlist, node.m_Start, node.m_End, (node.m_End-node.m_Start)/2+1); //Utils.select(array, indices, node.m_Start, node.m_End, (node.m_End-node.m_Start)/2+1); //(int) (node.m_NumInstances/2D+0.5D);

    Instance pivot;
    
    node.m_SplitAttrib = splitAttrib;
    node.m_SplitVal = m_Instances.instance(m_Instlist[medianIdx])
                                            .value(splitAttrib);

    node.m_Left = new BallNode(node.m_Start, medianIdxIdx, numNodesCreated+1,
                              (pivot=BallNode.calcCentroidPivot(node.m_Start,
                                                       medianIdxIdx, m_Instlist, 
                                                       m_Instances)), 
                              BallNode.calcRadius(node.m_Start, medianIdxIdx, 
                                                  m_Instlist, m_Instances, 
                                                  pivot, m_DistanceFunction)
                              );
    node.m_Right = new BallNode(medianIdxIdx+1, node.m_End, numNodesCreated+2,
                              (pivot=BallNode.calcCentroidPivot(medianIdxIdx+1,
                                                       node.m_End, m_Instlist, 
                                                       m_Instances)), 
                              BallNode.calcRadius(medianIdxIdx+1, node.m_End, 
                                                  m_Instlist, m_Instances, 
                                                  pivot, m_DistanceFunction)
                              );
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
   * @return the index of the middle element (in the master 
   * index array, i.e. index of the index of middle element).
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
  public int select(int attIdx, int[] indices, 
                            int left, int right, int k) {
    
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
   * Returns the widest dimension. The width of each 
   * dimension (for the points inside the node) is 
   * normalized, if m_NormalizeNodeWidth is set to 
   * true.
   * @param nodeRanges The attributes' range of the 
   * points inside the node that is to be split.
   * @param universe The attributes' range for the
   * whole point-space.
   * @return The index of the attribute/dimension
   * in which the points of the node have widest
   * spread.
   */
  protected int widestDim(double[][] nodeRanges, double[][] universe) {
    final int classIdx = m_Instances.classIndex();
    double widest = 0.0;
    int w = -1;
    if (m_NormalizeDimWidths) {
	for (int i = 0; i < nodeRanges.length; i++) {
	  double newWidest = nodeRanges[i][m_DistanceFunction.R_WIDTH] / 
	  		     universe[i][m_DistanceFunction.R_WIDTH];
	  if (newWidest > widest) {
	    if(i == classIdx) continue;
	    widest = newWidest;
	    w = i;
	  }
	}
    }
    else {
	for (int i = 0; i < nodeRanges.length; i++) {
	  if (nodeRanges[i][m_DistanceFunction.R_WIDTH] > widest) {
	    if(i == classIdx) continue;
	    widest = nodeRanges[i][m_DistanceFunction.R_WIDTH];
	    w = i;
	  }
	}
    }
    return w;
  }  
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String normalizeDimWidthsTipText() {
    return 
        "Whether to normalize the widths(ranges) of the dimensions "
      + "(attributes) before selecting the widest one.";
  }
  
  /** 
   * Should we normalize the widths(ranges) of the dimensions (attributes) 
   * before selecting the widest one. 
   * @param normalize Should be true if the widths are to be
   * normalized. 
   */
  public void setNormalizeDimWidths(boolean normalize) {
    m_NormalizeDimWidths = normalize;
  }
  
  /** 
   * Whether we are normalizing the widths(ranges) of the dimensions (attributes) 
   * or not.
   * @return true if widths are being normalized.
   */
  public boolean getNormalizeDimWidths() {
    return m_NormalizeDimWidths;
  }
  
  /**
   * Returns an enumeration describing the available options.
   * 
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> newVector = new Vector<Option>();
    
    newVector.addElement(new Option(
	"\tNormalize dimensions' widths.",
	"N", 0, "-N"));
    
    return newVector.elements();
  }

  /** 
   * Parses a given list of options.
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -N
   *  Normalize dimensions' widths.</pre>
   * 
   <!-- options-end --> 
   * 
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options)
    throws Exception {
    
    setNormalizeDimWidths(Utils.getFlag('N', options));
  }

  /** 
   * Gets the current settings.
   * @return An array of strings suitable for passing to 
   * setOptions or to be displayed by a 
   * GenericObjectEditor.
   */
  public String[] getOptions() {
    Vector<String>	result;
    
    result = new Vector<String>();

    if (getNormalizeDimWidths())
      result.add("-N");
    
    return result.toArray(new String[result.size()]);
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

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
 * MedianDistanceFromArbitraryPoint.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.balltrees;

import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class that splits a BallNode of a ball tree using Uhlmann's described method.<br/>
 * <br/>
 * For information see:<br/>
 * <br/>
 * Jeffrey K. Uhlmann (1991). Satisfying general proximity/similarity queries with metric trees. Information Processing Letters. 40(4):175-179.<br/>
 * <br/>
 * Ashraf Masood Kibriya (2007). Fast Algorithms for Nearest Neighbour Search. Hamilton, New Zealand.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Uhlmann1991,
 *    author = {Jeffrey K. Uhlmann},
 *    journal = {Information Processing Letters},
 *    month = {November},
 *    number = {4},
 *    pages = {175-179},
 *    title = {Satisfying general proximity/similarity queries with metric trees},
 *    volume = {40},
 *    year = {1991}
 * }
 * 
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
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;num&gt;
 *  The seed value for the random number generator.
 *  (default: 17)</pre>
 * 
 <!-- options-end -->
 *
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public class MedianDistanceFromArbitraryPoint
  extends BallSplitter 
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = 5617378551363700558L;
  
  /** Seed for random number generator. */
  protected int m_RandSeed = 17;

  /** 
   * Random number generator for selecting
   * an abitrary (random) point. 
   */
  protected Random m_Rand;
  
  /** Constructor. */
  public MedianDistanceFromArbitraryPoint() {
  }
  
  /**
   * Constructor. 
   * @param instList The master index array.
   * @param insts The instances on which the tree
   * is (or is to be) built.
   * @param e The Euclidean distance function to 
   * use for splitting.
   */
  public MedianDistanceFromArbitraryPoint(int[] instList, Instances insts, EuclideanDistance e) {
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
        "Class that splits a BallNode of a ball tree using Uhlmann's "
      + "described method.\n\n"
      + "For information see:\n\n"
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
    TechnicalInformation additional;

    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Jeffrey K. Uhlmann");
    result.setValue(Field.TITLE, "Satisfying general proximity/similarity queries with metric trees");
    result.setValue(Field.JOURNAL, "Information Processing Letters");
    result.setValue(Field.MONTH, "November");
    result.setValue(Field.YEAR, "1991");
    result.setValue(Field.NUMBER, "4");
    result.setValue(Field.VOLUME, "40");
    result.setValue(Field.PAGES, "175-179");

    additional = result.add(Type.MASTERSTHESIS);
    additional.setValue(Field.AUTHOR, "Ashraf Masood Kibriya");
    additional.setValue(Field.TITLE, "Fast Algorithms for Nearest Neighbour Search");
    additional.setValue(Field.YEAR, "2007");
    additional.setValue(Field.SCHOOL, "Department of Computer Science, School of Computing and Mathematical Sciences, University of Waikato");
    additional.setValue(Field.ADDRESS, "Hamilton, New Zealand");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();
    
    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement((Option)enm.nextElement());
      
    result.addElement(new Option(
        "\tThe seed value for the random number generator.\n"
        + "\t(default: 17)",
        "S", 1, "-S <num>"));
    
    return result.elements();
  }

  /**
   * Parses a given list of options.
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;num&gt;
   *  The seed value for the random number generator.
   *  (default: 17)</pre>
   * 
   <!-- options-end -->
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    super.setOptions(options);

    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() > 0)
      setRandomSeed(Integer.parseInt(tmpStr));
    else
      setRandomSeed(17);
  }

  /**
   * Gets the current settings of the object.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    String[]		options;
    int			i;
    
    result  = new Vector<String>();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-S");
    result.add("" + getRandomSeed());
    
    return result.toArray(new String[result.size()]);
  }
  
  /**
   * Sets the seed for random number generator.
   * @param seed The seed value to set.
   */
  public void setRandomSeed(int seed) {
    m_RandSeed = seed;
  }
  
  /**
   * Returns the seed value of random 
   * number generator.
   * @return The random seed currently in use.
   */
  public int getRandomSeed() {
    return m_RandSeed;
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui.
   */
  public String randomSeedTipText() {
    return "The seed value for the random number generator.";
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

    m_Rand = new Random(m_RandSeed);
    
    int ridx = node.m_Start+m_Rand.nextInt(node.m_NumInstances);
    Instance randomInst = (Instance)
                            m_Instances.instance( m_Instlist[ridx] ).copy();
    double [] distList = new double[node.m_NumInstances-1];
    Instance temp;
    for(int i=node.m_Start, j=0; i<node.m_End; i++, j++) {
      temp = m_Instances.instance( m_Instlist[i] );
      distList[j] = m_DistanceFunction.distance(randomInst, temp, 
	  Double.POSITIVE_INFINITY);
    }
    
    int medianIdx = select(distList, m_Instlist, 0, distList.length-1, 
                           node.m_Start, (node.m_End-node.m_Start)/2+1) + 
                    node.m_Start;
    
    Instance pivot;
    node.m_Left = new BallNode(node.m_Start, medianIdx, numNodesCreated+1,
                              (pivot=BallNode.calcCentroidPivot(node.m_Start,
                                                       medianIdx, m_Instlist, 
                                                       m_Instances)), 
                              BallNode.calcRadius(node.m_Start, medianIdx, 
                                                  m_Instlist, m_Instances, 
                                                  pivot, m_DistanceFunction)
                              );
    
    node.m_Right = new BallNode(medianIdx+1, node.m_End, numNodesCreated+2,
                              (pivot=BallNode.calcCentroidPivot(medianIdx+1,
                                                       node.m_End, m_Instlist, 
                                                       m_Instances)), 
                              BallNode.calcRadius(medianIdx+1, node.m_End, 
                                                  m_Instlist, m_Instances, 
                                                  pivot, m_DistanceFunction)
                              );
  }
  
  /**
   * Partitions the instances around a pivot. Used by quicksort and
   * kthSmallestValue.
   *
   * @param array The array of distances of the points to the
   * arbitrarily selected point.
   * @param index The master index array containing indices of the 
   * instances.
   * @param l The relative begining index of the portion of master 
   * index array that should be partitioned. 
   * @param r The relative end index of the portion of master index 
   * array that should be partitioned.
   * @param indexStart The absolute begining index of the portion 
   * of master index array that should be partitioned. 
   * @return the index of the middle element (in the master 
   * index array, i.e. index of the index of middle element).
   */
  protected int partition(double[] array, int[] index, int l, int r, 
                        final int indexStart) {
    
    double pivot = array[(l + r) / 2];
    int help;

    while (l < r) {
      while ((array[l] < pivot) && (l < r)) {
        l++;
      }
      while ((array[r] > pivot) && (l < r)) {
        r--;
      }
      if (l < r) {
        help = index[indexStart+l];
        index[indexStart+l] = index[indexStart+r];
        index[indexStart+r] = help;
        l++;
        r--;
      }
    }
    if ((l == r) && (array[r] > pivot)) {
      r--;
    } 

    return r;
  }
  
  /**
   * Implements computation of the kth-smallest element according
   * to Manber's "Introduction to Algorithms".
   *
   * @param array Array containing the distances of points from
   * the arbitrarily selected.
   * @param indices The master index array containing indices of 
   * the instances.
   * @param left The relative begining index of the portion of the 
   * master index array in which to find the kth-smallest element.
   * @param right The relative end index of the portion of the 
   * master index array in which to find the kth-smallest element.
   * @param indexStart The absolute begining index of the portion 
   * of the master index array in which to find the kth-smallest 
   * element. 
   * @param k The value of k
   * @return The index of the kth-smallest element
   */
  protected int select(double[] array, int[] indices, 
                            int left, int right, final int indexStart, int k) {
    
    if (left == right) {
      return left;
    } else {
      int middle = partition(array, indices, left, right, indexStart);
      if ((middle - left + 1) >= k) {
        return select(array, indices, left, middle, indexStart, k);
      } else {
        return select(array, indices, middle + 1, right, 
                      indexStart, k - (middle - left + 1));
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

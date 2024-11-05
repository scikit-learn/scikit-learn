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
 *    KDTree.java
 *    Copyright (C) 2000-2007 University of Waikato
 *    
 */

package weka.core.neighboursearch;

import weka.core.DistanceFunction;
import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.neighboursearch.kdtrees.KDTreeNode;
import weka.core.neighboursearch.kdtrees.KDTreeNodeSplitter;
import weka.core.neighboursearch.kdtrees.SlidingMidPointOfWidestSide;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class implementing the KDTree search algorithm for nearest neighbour search.<br/>
 * The connection to dataset is only a reference. For the tree structure the indexes are stored in an array. <br/>
 * Building the tree:<br/>
 * If a node has &lt;maximal-inst-number&gt; (option -L) instances no further splitting is done. Also if the split would leave one side empty, the branch is not split any further even if the instances in the resulting node are more than &lt;maximal-inst-number&gt; instances.<br/>
 * **PLEASE NOTE:** The algorithm can not handle missing values, so it is advisable to run ReplaceMissingValues filter if there are any missing values in the dataset.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * Jerome H. Friedman, Jon Luis Bentley, Raphael Ari Finkel (1977). An Algorithm for Finding Best Matches in Logarithmic Expected Time. ACM Transactions on Mathematics Software. 3(3):209-226.<br/>
 * <br/>
 * Andrew Moore (1991). A tutorial on kd-trees.
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
 * 
 * &#64;techreport{Moore1991,
 *    author = {Andrew Moore},
 *    booktitle = {University of Cambridge Computer Laboratory Technical Report No. 209},
 *    howpublished = {Extract from PhD Thesis},
 *    title = {A tutorial on kd-trees},
 *    year = {1991},
 *    HTTP = {Available from http://www.autonlab.org/autonweb/14665.html}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;classname and options&gt;
 *  Node splitting method to use.
 *  (default: weka.core.neighboursearch.kdtrees.SlidingMidPointOfWidestSide)</pre>
 * 
 * <pre> -W &lt;value&gt;
 *  Set minimal width of a box
 *  (default: 1.0E-2).</pre>
 * 
 * <pre> -L
 *  Maximal number of instances in a leaf
 *  (default: 40).</pre>
 * 
 * <pre> -N
 *  Normalizing will be done
 *  (Select dimension for split, with normalising to universe).</pre>
 * 
 <!-- options-end --> 
 * 
 * @author Gabi Schmidberger (gabi[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @author Malcolm Ware (mfw4[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5987 $
 */
public class KDTree
  extends NearestNeighbourSearch
  implements TechnicalInformationHandler {

  /** For serialization. */
  private static final long serialVersionUID = 1505717283763272533L;

  /**
   * Array holding the distances of the nearest neighbours. It is filled up both
   * by nearestNeighbour() and kNearestNeighbours().
   */
  protected double[] m_DistanceList;

  /**
   * Indexlist of the instances of this kdtree. Instances get sorted according
   * to the splits. the nodes of the KDTree just hold their start and end
   * indices
   */
  protected int[] m_InstList;

  /** The root node of the tree. */
  protected KDTreeNode m_Root;

  /** The node splitter. */
  protected KDTreeNodeSplitter m_Splitter = new SlidingMidPointOfWidestSide();

  /** Tree stats. */
  protected int m_NumNodes, m_NumLeaves, m_MaxDepth;

  /** Tree Stats variables. */
  protected TreePerformanceStats m_TreeStats = null;

  // Constants
  /** The index of MIN value in attributes' range array. */
  public static final int MIN = EuclideanDistance.R_MIN;
  
  /** The index of MAX value in attributes' range array. */
  public static final int MAX = EuclideanDistance.R_MAX; 
  
  /** The index of WIDTH (MAX-MIN) value in attributes' range array. */
  public static final int WIDTH = EuclideanDistance.R_WIDTH;

  /**
   * Returns an instance of a TechnicalInformation object, containing detailed
   * information about the technical background of this class, e.g., paper
   * reference or book this class is based on.
   * 
   * @return		the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation result;
    TechnicalInformation additional;

    result = new TechnicalInformation(Type.ARTICLE);
    result.setValue(Field.AUTHOR, "Jerome H. Friedman and Jon Luis Bentley and Raphael Ari Finkel");
    result.setValue(Field.YEAR, "1977");
    result.setValue(Field.TITLE, "An Algorithm for Finding Best Matches in Logarithmic Expected Time");
    result.setValue(Field.JOURNAL, "ACM Transactions on Mathematics Software");
    result.setValue(Field.PAGES, "209-226");
    result.setValue(Field.MONTH, "September");
    result.setValue(Field.VOLUME, "3");
    result.setValue(Field.NUMBER, "3");

    additional = result.add(Type.TECHREPORT);
    additional.setValue(Field.AUTHOR, "Andrew Moore");
    additional.setValue(Field.YEAR, "1991");
    additional.setValue(Field.TITLE, "A tutorial on kd-trees");
    additional.setValue(Field.HOWPUBLISHED, "Extract from PhD Thesis");
    additional.setValue(Field.BOOKTITLE, "University of Cambridge Computer Laboratory Technical Report No. 209");
    additional.setValue(Field.HTTP, "Available from http://www.autonlab.org/autonweb/14665.html");

    return result;
  }

  /**
   * Creates a new instance of KDTree.
   */
  public KDTree() {
    super();
    if (getMeasurePerformance())
      m_Stats = m_TreeStats = new TreePerformanceStats();
  }

  /**
   * Creates a new instance of KDTree. 
   * It also builds the tree on supplied set of Instances.
   * @param insts The instances/points on which the BallTree 
   * should be built on.
   */
  public KDTree(Instances insts) {
    super(insts);
    if (getMeasurePerformance())
      m_Stats = m_TreeStats = new TreePerformanceStats();
  }

  /**
   * Builds the KDTree on the supplied set of instances/points. It 
   * is adviseable to run the replace missing attributes filter 
   * on the passed instances first.
   * NOTE: This method should not be called from outside this 
   * class. Outside classes should call setInstances(Instances)
   * instead.
   * 
   * @param instances	The instances to build the tree on
   * @throws Exception	if something goes wrong
   */
  protected void buildKDTree(Instances instances) throws Exception {

    checkMissing(instances);
    if (m_EuclideanDistance == null)
      m_DistanceFunction = m_EuclideanDistance = new EuclideanDistance(
          instances);
    else
      m_EuclideanDistance.setInstances(instances);

    m_Instances = instances;
    int numInst = m_Instances.numInstances();

    // Make the global index list
    m_InstList = new int[numInst];

    for (int i = 0; i < numInst; i++) {
      m_InstList[i] = i;
    }

    double[][] universe = m_EuclideanDistance.getRanges();

    // initializing internal fields of KDTreeSplitter
    m_Splitter.setInstances(m_Instances);
    m_Splitter.setInstanceList(m_InstList);
    m_Splitter.setEuclideanDistanceFunction(m_EuclideanDistance);
    m_Splitter.setNodeWidthNormalization(m_NormalizeNodeWidth);

    // building tree
    m_NumNodes = m_NumLeaves = 1;
    m_MaxDepth = 0;
    m_Root = new KDTreeNode(m_NumNodes, 0, m_Instances.numInstances() - 1,
        universe);

    splitNodes(m_Root, universe, m_MaxDepth + 1);
  }

  /** 
   * Recursively splits nodes of a tree starting from the supplied node.
   * The splitting stops for any node for which the number of instances/points
   * falls below a given threshold (given by m_MaxInstInLeaf), or if the 
   * maximum relative width/range of the instances/points 
   * (i.e. max_i(max(att_i) - min(att_i)) ) falls below a given threshold 
   * (given by m_MinBoxRelWidth). 
   * 
   * @param node The node to start splitting from.
   * @param universe The attribute ranges of the whole dataset.
   * @param depth The depth of the supplied node.  
   * @throws Exception If there is some problem 
   * splitting.
   */
  protected void splitNodes(KDTreeNode node, double[][] universe,
      int depth) throws Exception {
    double[][] nodeRanges = m_EuclideanDistance.initializeRanges(m_InstList,
                                                 node.m_Start, node.m_End);
    if (node.numInstances() <= m_MaxInstInLeaf
        || getMaxRelativeNodeWidth(nodeRanges, universe) <= m_MinBoxRelWidth)
      return;

    // splitting a node so it is no longer a leaf
    m_NumLeaves--;

    if (depth > m_MaxDepth)
      m_MaxDepth = depth;

    m_Splitter.splitNode(node, m_NumNodes, nodeRanges, universe);
    m_NumNodes += 2;
    m_NumLeaves += 2;

    splitNodes(node.m_Left, universe, depth + 1);
    splitNodes(node.m_Right, universe, depth + 1);
  }

  /**
   * Returns (in the supplied heap object) the k nearest 
   * neighbours of the given instance starting from the give 
   * tree node. &gt;k neighbours are returned if there are more than 
   * one neighbours at the kth boundary. NOTE: This method should 
   * not be used from outside this class. Outside classes should 
   * call kNearestNeighbours(Instance, int).
   * 
   * @param target  The instance to find the nearest neighbours for.
   * @param node The KDTreeNode to start the search from.
   * @param k    The number of neighbours to find.
   * @param heap The MyHeap object to store/update the kNNs found
   * during the search.
   * @param distanceToParents The distance of the supplied target 
   * to the parents of the supplied tree node. 
   * @throws Exception  if the nearest neighbour could not be found.
   */
  protected void findNearestNeighbours(Instance target, KDTreeNode node, int k,
      MyHeap heap, double distanceToParents) throws Exception {
    if (node.isALeaf()) {
      if (m_TreeStats != null) {
        m_TreeStats.updatePointCount(node.numInstances());
        m_TreeStats.incrLeafCount();
      }
      double distance;
      // look at all the instances in this leaf
      for (int idx = node.m_Start; idx <= node.m_End; idx++) {
        if (target == m_Instances.instance(m_InstList[idx])) // for
                                                              // hold-one-out
                                                              // cross-validation
          continue;
        if (heap.size() < k) {
          distance = m_EuclideanDistance.distance(target, m_Instances
              .instance(m_InstList[idx]), Double.POSITIVE_INFINITY, m_Stats);
          heap.put(m_InstList[idx], distance);
        } else {
          MyHeapElement temp = heap.peek();
          distance = m_EuclideanDistance.distance(target, m_Instances
              .instance(m_InstList[idx]), temp.distance, m_Stats);
          if (distance < temp.distance) {
            heap.putBySubstitute(m_InstList[idx], distance);
          } else if (distance == temp.distance) {
            heap.putKthNearest(m_InstList[idx], distance);
          }
        }// end else heap.size==k
      }// end for

    } else {
      if (m_TreeStats != null) {
        m_TreeStats.incrIntNodeCount();
      }
      KDTreeNode nearer, further;
      boolean targetInLeft = m_EuclideanDistance.valueIsSmallerEqual(target,
          node.m_SplitDim, node.m_SplitValue);

      if (targetInLeft) {
        nearer = node.m_Left;
        further = node.m_Right;
      } else {
        nearer = node.m_Right;
        further = node.m_Left;
      }
      findNearestNeighbours(target, nearer, k, heap, distanceToParents);

      // ... now look in further half if maxDist reaches into it
      if (heap.size() < k) { // if haven't found the first k
        double distanceToSplitPlane = distanceToParents
            + m_EuclideanDistance.sqDifference(node.m_SplitDim, target
                .value(node.m_SplitDim), node.m_SplitValue);
        findNearestNeighbours(target, further, k, heap, distanceToSplitPlane);
        return;
      } else { // else see if ball centered at query intersects with the other
                // side.
        double distanceToSplitPlane = distanceToParents
            + m_EuclideanDistance.sqDifference(node.m_SplitDim, target
                .value(node.m_SplitDim), node.m_SplitValue);
        if (heap.peek().distance >= distanceToSplitPlane) {
          findNearestNeighbours(target, further, k, heap, distanceToSplitPlane);
        }
      }// end else
    }// end else_if an internal node
  }

  /**
   * Returns the k nearest neighbours of the supplied instance.
   * &gt;k neighbours are returned if there are more than one 
   * neighbours at the kth boundary. 
   * 
   * @param target	The instance to find the nearest neighbours for.
   * @param k 		The number of neighbours to find.
   * @return The k nearest neighbours (or &gt;k if more there are than
   * one neighbours at the kth boundary). 
   * @throws Exception 	if the nearest neighbour could not be found.
   */
  public Instances kNearestNeighbours(Instance target, int k) throws Exception {
    checkMissing(target);

    if (m_Stats != null)
      m_Stats.searchStart();

    MyHeap heap = new MyHeap(k);
    findNearestNeighbours(target, m_Root, k, heap, 0.0);

    if (m_Stats != null)
      m_Stats.searchFinish();

    Instances neighbours = new Instances(m_Instances, (heap.size() + heap
        .noOfKthNearest()));
    m_DistanceList = new double[heap.size() + heap.noOfKthNearest()];
    int[] indices = new int[heap.size() + heap.noOfKthNearest()];
    int i = indices.length - 1;
    MyHeapElement h;
    while (heap.noOfKthNearest() > 0) {
      h = heap.getKthNearest();
      indices[i] = h.index;
      m_DistanceList[i] = h.distance;
      i--;
    }
    while (heap.size() > 0) {
      h = heap.get();
      indices[i] = h.index;
      m_DistanceList[i] = h.distance;
      i--;
    }
    m_DistanceFunction.postProcessDistances(m_DistanceList);

    for (int idx = 0; idx < indices.length; idx++) {
      neighbours.add(m_Instances.instance(indices[idx]));
    }

    return neighbours;
  }
  

  /**
   * Returns the nearest neighbour of the supplied target 
   * instance. 
   *  
   * @param target	The instance to find the nearest neighbour for.
   * @return The nearest neighbour from among the previously 
   * supplied training instances.
   * @throws Exception 	if the neighbours could not be found.
   */
  public Instance nearestNeighbour(Instance target) throws Exception {
    return (kNearestNeighbours(target, 1)).instance(0);
  }
  
  /**
   * Returns the distances to the kNearest or 1 nearest neighbour currently
   * found with either the kNearestNeighbours or the nearestNeighbour method.
   * 
   * @return array containing the distances of the
   *         nearestNeighbours. The length and ordering of the array 
   *         is the same as that of the instances returned by 
   *         nearestNeighbour functions.
   * @throws Exception 	if called before calling kNearestNeighbours or
   * 			nearestNeighbours.
   */
  public double[] getDistances() throws Exception {
    if (m_Instances == null || m_DistanceList == null)
      throw new Exception("The tree has not been supplied with a set of "
          + "instances or getDistances() has been called "
          + "before calling kNearestNeighbours().");
    return m_DistanceList;
  }
  

  /**
   * Builds the KDTree on the given set of instances.
   * @param instances The insts on which the KDTree is to be 
   * built. 
   * @throws Exception If some error occurs while 
   * building the KDTree
   */
  public void setInstances(Instances instances) throws Exception {
    super.setInstances(instances);
    buildKDTree(instances);
  }
  

  /**
   * Adds one instance to the KDTree. This updates the KDTree structure to take
   * into account the newly added training instance.
   * 
   * @param instance 	the instance to be added. Usually the newly added instance in the
   *          		training set.
   * @throws Exception If the instance cannot be added.
   */
  public void update(Instance instance) throws Exception { // better to change
                                                            // to addInstance
    if (m_Instances == null)
      throw new Exception("No instances supplied yet. Have to call "
          + "setInstances(instances) with a set of Instances " + "first.");

    addInstanceInfo(instance);
    addInstanceToTree(instance, m_Root);
  }

  /**
   * Recursively adds an instance to the tree starting from
   * the supplied KDTreeNode.
   * NOTE: This should not be called by outside classes,
   * outside classes should instead call update(Instance)
   * method. 
   *  
   * @param inst The instance to add to the tree
   * @param node The node to start the recursive search 
   * from, for the leaf node where the supplied instance 
   * would go.
   * @throws Exception If some error occurs while adding
   * the instance.
   */
  protected void addInstanceToTree(Instance inst, KDTreeNode node)
      throws Exception {
    if (node.isALeaf()) {
      int instList[] = new int[m_Instances.numInstances()];
      try {
        System.arraycopy(m_InstList, 0, instList, 0, node.m_End + 1); // m_InstList.squeezeIn(m_End,
                                                                      // index);
        if (node.m_End < m_InstList.length - 1)
          System.arraycopy(m_InstList, node.m_End + 1, instList,
              node.m_End + 2, m_InstList.length - node.m_End - 1);
        instList[node.m_End + 1] = m_Instances.numInstances() - 1;
      } catch (ArrayIndexOutOfBoundsException ex) {
        System.err.println("m_InstList.length: " + m_InstList.length
            + " instList.length: " + instList.length + "node.m_End+1: "
            + (node.m_End + 1) + "m_InstList.length-node.m_End+1: "
            + (m_InstList.length - node.m_End - 1));
        throw ex;
      }
      m_InstList = instList;

      node.m_End++;
      node.m_NodeRanges = m_EuclideanDistance.updateRanges(inst,
          node.m_NodeRanges);

      m_Splitter.setInstanceList(m_InstList);

      // split this leaf node if necessary
      double[][] universe = m_EuclideanDistance.getRanges();
      if (node.numInstances() > m_MaxInstInLeaf
          && getMaxRelativeNodeWidth(node.m_NodeRanges, universe) > m_MinBoxRelWidth) {
        m_Splitter.splitNode(node, m_NumNodes, node.m_NodeRanges, universe);
        m_NumNodes += 2;
      }
    }// end if node is a leaf
    else {
      if (m_EuclideanDistance.valueIsSmallerEqual(inst, node.m_SplitDim,
          node.m_SplitValue)) {
        addInstanceToTree(inst, node.m_Left);
        afterAddInstance(node.m_Right);
      } else
        addInstanceToTree(inst, node.m_Right);

      node.m_End++;
      node.m_NodeRanges = m_EuclideanDistance.updateRanges(inst,
          node.m_NodeRanges);
    }
  }

  /**
   * Corrects the start and end indices of a 
   * KDTreeNode after an instance is added to
   * the tree. The start and end indices for
   * the master index array (m_InstList) 
   * stored in the nodes need to be updated
   * for all nodes in the subtree on the 
   * right of a node where the instance 
   * was added. 
   * NOTE: No outside class should call this
   * method.
   * 
   * @param node KDTreeNode whose start and end indices 
   * need to be updated.
   */
  protected void afterAddInstance(KDTreeNode node) {
    node.m_Start++;
    node.m_End++;
    if (!node.isALeaf()) {
      afterAddInstance(node.m_Left);
      afterAddInstance(node.m_Right);
    }
  }

  /**
   * Adds one instance to KDTree loosly. It only changes the ranges in
   * EuclideanDistance, and does not affect the structure of the KDTree.
   * 
   * @param instance	the new instance. Usually this is the test instance 
   * 			supplied to update the range of attributes in the distance function.
   */
  public void addInstanceInfo(Instance instance) {
    m_EuclideanDistance.updateRanges(instance);
  }

  /**
   * Checks if there is any instance with missing values. Throws an exception if
   * there is, as KDTree does not handle missing values.
   * 
   * @param instances	the instances to check
   * @throws Exception	if missing values are encountered
   */
  protected void checkMissing(Instances instances) throws Exception {
    for (int i = 0; i < instances.numInstances(); i++) {
      Instance ins = instances.instance(i);
      for (int j = 0; j < ins.numValues(); j++) {
        if (ins.index(j) != ins.classIndex())
          if (ins.isMissingSparse(j)) {
            throw new Exception("ERROR: KDTree can not deal with missing "
                + "values. Please run ReplaceMissingValues filter "
                + "on the dataset before passing it on to the KDTree.");
          }
      }
    }
  }

  /**
   * Checks if there is any missing value in the given 
   * instance.
   * @param ins The instance to check missing values in.
   * @throws Exception If there is a missing value in the 
   * instance.
   */
  protected void checkMissing(Instance ins) throws Exception {
    for (int j = 0; j < ins.numValues(); j++) {
      if (ins.index(j) != ins.classIndex())
        if (ins.isMissingSparse(j)) {
          throw new Exception("ERROR: KDTree can not deal with missing "
              + "values. Please run ReplaceMissingValues filter "
              + "on the dataset before passing it on to the KDTree.");
        }
    }
  }
  
  /** 
   * Returns the maximum attribute width of instances/points 
   * in a KDTreeNode relative to the whole dataset. 
   * 
   * @param nodeRanges The attribute ranges of the 
   * KDTreeNode whose maximum relative width is to be 
   * determined.
   * @param universe The attribute ranges of the whole
   * dataset (training instances + test instances so 
   * far encountered).
   * @return The maximum relative width
   */
  protected double getMaxRelativeNodeWidth(double[][] nodeRanges,
      double[][] universe) {
    int widest = widestDim(nodeRanges, universe);
    if(widest < 0)
    	return 0.0;
    else
    	return nodeRanges[widest][WIDTH] / universe[widest][WIDTH];
  }

  /**
   * Returns the widest dimension/attribute in a 
   * KDTreeNode (widest after normalizing).
   * @param nodeRanges The attribute ranges of 
   * the KDTreeNode.
   * @param universe The attribute ranges of the 
   * whole dataset (training instances + test 
   * instances so far encountered).
   * @return The index of the widest 
   * dimension/attribute.
   */
  protected int widestDim(double[][] nodeRanges, double[][] universe) {
    final int classIdx = m_Instances.classIndex();
    double widest = 0.0;
    int w = -1;
    if (m_NormalizeNodeWidth) {
      for (int i = 0; i < nodeRanges.length; i++) {
        double newWidest = nodeRanges[i][WIDTH] / universe[i][WIDTH];
        if (newWidest > widest) {
          if (i == classIdx)
            continue;
          widest = newWidest;
          w = i;
        }
      }
    } else {
      for (int i = 0; i < nodeRanges.length; i++) {
        if (nodeRanges[i][WIDTH] > widest) {
          if (i == classIdx)
            continue;
          widest = nodeRanges[i][WIDTH];
          w = i;
        }
      }
    }
    return w;
  }

  /**
   * Returns the size of the tree.
   * 
   * @return 		the size of the tree
   */
  public double measureTreeSize() {
    return m_NumNodes;
  }

  /**
   * Returns the number of leaves.
   * 
   * @return 		the number of leaves
   */
  public double measureNumLeaves() {
    return m_NumLeaves;
  }

  /**
   * Returns the depth of the tree.
   * 
   * @return The depth of the tree
   */
  public double measureMaxDepth() {
    return m_MaxDepth;
  }

  /**
   * Returns an enumeration of the additional measure names.
   * 
   * @return 		an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector<String> newVector = new Vector<String>();
    newVector.addElement("measureTreeSize");
    newVector.addElement("measureNumLeaves");
    newVector.addElement("measureMaxDepth");
    if (m_Stats != null) {
      for (Enumeration e = m_Stats.enumerateMeasures(); e.hasMoreElements();) {
        newVector.addElement((String)e.nextElement());
      }
    }
    return newVector.elements();
  }

  /**
   * Returns the value of the named measure.
   * 
   * @param additionalMeasureName	the name of 
   * the measure to query for its value.
   * @return The value of the named measure
   * @throws IllegalArgumentException	If the named measure 
   * is not supported.
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.compareToIgnoreCase("measureMaxDepth") == 0) {
      return measureMaxDepth();
    } else if (additionalMeasureName.compareToIgnoreCase("measureTreeSize") == 0) {
      return measureTreeSize();
    } else if (additionalMeasureName.compareToIgnoreCase("measureNumLeaves") == 0) {
      return measureNumLeaves();
    } else if (m_Stats != null) {
      return m_Stats.getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException(additionalMeasureName
          + " not supported (KDTree)");
    }
  }

  /**
   * Sets whether to calculate the performance statistics or not.
   * @param measurePerformance Should be true if performance 
   * statistics are to be measured.
   */
  public void setMeasurePerformance(boolean measurePerformance) {
    m_MeasurePerformance = measurePerformance;
    if (m_MeasurePerformance) {
      if (m_Stats == null)
        m_Stats = m_TreeStats = new TreePerformanceStats();
    } else
      m_Stats = m_TreeStats = null;
  }

  /**
   * Assigns instances to centers using KDTree.
   * 
   * @param centers	the current centers
   * @param assignments	the centerindex for each instance
   * @param pc		the threshold value for pruning.
   * @throws Exception If there is some problem 
   * assigning instances to centers.
   */
  public void centerInstances(Instances centers, int[] assignments, double pc)
      throws Exception {

    int[] centList = new int[centers.numInstances()];
    for (int i = 0; i < centers.numInstances(); i++)
      centList[i] = i;

    determineAssignments(m_Root, centers, centList, assignments, pc);
  }

  /**
   * Assigns instances to the current centers called candidates.
   * 
   * @param node The node to start assigning the instances from.
   * @param centers	all the current centers.
   * @param candidates	the current centers the method works on.
   * @param assignments	the center index for each instance.
   * @param pc the threshold value for pruning.
   * @throws Exception If there is some problem assigning 
   * instances to centers.
   */
  protected void determineAssignments(KDTreeNode node, Instances centers,
      int[] candidates, int[] assignments, double pc) throws Exception {

    // reduce number of owners for current hyper rectangle
    int[] owners = refineOwners(node, centers, candidates);

    // only one owner
    if (owners.length == 1) {
      // all instances of this node are owned by one center
      for (int i = node.m_Start; i <= node.m_End; i++) {
        assignments[m_InstList[i]] // the assignment of this instance
        = owners[0]; // is the current owner
      }
    } else if (!node.isALeaf()) {
      // more than one owner and it is not a leaf
      determineAssignments(node.m_Left, centers, owners, assignments, pc);
      determineAssignments(node.m_Right, centers, owners, assignments, pc);
    } else {
      // this is a leaf and there are more than 1 owner
      // XMeans.
      assignSubToCenters(node, centers, owners, assignments);
    }
  }

  /**
   * Refines the ownerlist.
   * 
   * @param node The current tree node.
   * @param centers	all centers
   * @param candidates	the indexes of those centers that are candidates.
   * @return list of owners
   * @throws Exception If some problem occurs in refining.
   */
  protected int[] refineOwners(KDTreeNode node, Instances centers,
      int[] candidates) throws Exception {

    int[] owners = new int[candidates.length];
    double minDistance = Double.POSITIVE_INFINITY;
    int ownerIndex = -1;
    Instance owner;
    int numCand = candidates.length;
    double[] distance = new double[numCand];
    boolean[] inside = new boolean[numCand];
    for (int i = 0; i < numCand; i++) {
      distance[i] = distanceToHrect(node, centers.instance(candidates[i]));
      inside[i] = (distance[i] == 0.0);
      if (distance[i] < minDistance) {
        minDistance = distance[i];
        ownerIndex = i;
      }
    }
    owner = (Instance)centers.instance(candidates[ownerIndex]).copy();

    // are there other owners
    // loop also goes over already found owner, keeps order
    // in owner list
    int index = 0;
    for (int i = 0; i < numCand; i++) {
      // 1. all centers that are points within rectangle are owners
      if ((inside[i])

      // 2. take all points with same distance to the rect. as the owner
          || (distance[i] == distance[ownerIndex])) {

        // add competitor to owners list
        owners[index++] = candidates[i];
      } else {

        Instance competitor = (Instance)centers.instance(candidates[i]).copy();
        if

        // 3. point has larger distance to rectangle but still can compete
        // with owner for some points in the rectangle
        (!candidateIsFullOwner(node, owner, competitor))

        {
          // also add competitor to owners list
          owners[index++] = candidates[i];
        }
      }
    }
    int[] result = new int[index];
    for (int i = 0; i < index; i++)
      result[i] = owners[i];
    return result;
  }

  /**
   * Returns the distance between a point and an hyperrectangle.
   * 
   * @param node The current node from whose hyperrectangle 
   * the distance is to be measured.
   * @param x		the point
   * @return 		the distance
   * @throws Exception If some problem occurs in determining 
   * the distance to the hyperrectangle.
   */
  protected double distanceToHrect(KDTreeNode node, Instance x) throws Exception {
    double distance = 0.0;

    Instance closestPoint = (Instance)x.copy();
    boolean inside;
    inside = clipToInsideHrect(node, closestPoint);
    if (!inside)
      distance = m_EuclideanDistance.distance(closestPoint, x);
    return distance;
  }

  /**
   * Finds the closest point in the hyper rectangle to a given point. Change the
   * given point to this closest point by clipping of at all the dimensions to
   * be clipped of. If the point is inside the rectangle it stays unchanged. The
   * return value is true if the point was not changed, so the the return value
   * is true if the point was inside the rectangle.
   * 
   * @param node The current KDTreeNode in whose hyperrectangle the closest 
   * point is to be found.
   * @param x		a point
   * @return 		true if the input point stayed unchanged.
   */
  protected boolean clipToInsideHrect(KDTreeNode node, Instance x) {
    boolean inside = true;
    for (int i = 0; i < m_Instances.numAttributes(); i++) {
      // TODO treat nominals differently!??

      if (x.value(i) < node.m_NodeRanges[i][MIN]) {
        x.setValue(i, node.m_NodeRanges[i][MIN]);
        inside = false;
      } else if (x.value(i) > node.m_NodeRanges[i][MAX]) {
        x.setValue(i, node.m_NodeRanges[i][MAX]);
        inside = false;
      }
    }
    return inside;
  }

  /**
   * Returns true if candidate is a full owner in respect to a competitor.
   * <p>
   * 
   * The candidate has been the closer point to the current rectangle or even
   * has been a point within the rectangle. The competitor is competing with the
   * candidate for a few points out of the rectangle although it is a point
   * further away from the rectangle then the candidate. The extrem point is the
   * corner of the rectangle that is furthest away from the candidate towards
   * the direction of the competitor.
   * 
   * If the distance candidate to this extreme point is smaller then the
   * distance competitor to this extreme point, then it is proven that none of
   * the points in the rectangle can be owned be the competitor and the
   * candidate is full owner of the rectangle in respect to this competitor. See
   * also D. Pelleg and A. Moore's paper 'Accelerating exact k-means Algorithms
   * with Geometric Reasoning'.
   * <p>
   * 
   * @param node The current KDTreeNode / hyperrectangle.
   * @param candidate	instance that is candidate to be owner
   * @param competitor	instance that competes against the candidate
   * @return 		true if candidate is full owner
   * @throws Exception If some problem occurs.
   */
  protected boolean candidateIsFullOwner(KDTreeNode node, Instance candidate,
      Instance competitor) throws Exception {
    // get extreme point
    Instance extreme = (Instance)candidate.copy();
    for (int i = 0; i < m_Instances.numAttributes(); i++) {
      if ((competitor.value(i) - candidate.value(i)) > 0) {
        extreme.setValue(i, node.m_NodeRanges[i][MAX]);
      } else {
        extreme.setValue(i, node.m_NodeRanges[i][MIN]);
      }
    }
    boolean isFullOwner = m_EuclideanDistance.distance(extreme, candidate) < m_EuclideanDistance
        .distance(extreme, competitor);

    return isFullOwner;
  }

  /**
   * Assigns instances of this node to center. Center to be assign to is decided
   * by the distance function.
   * 
   * @param node The KDTreeNode whose instances are to be assigned.
   * @param centers	all the input centers
   * @param centList	the list of centers to work with
   * @param assignments	index list of last assignments
   * @throws Exception If there is error assigning the instances.
   */
  public void assignSubToCenters(KDTreeNode node, Instances centers,
      int[] centList, int[] assignments) throws Exception {
    // todo: undecided situations
    int numCent = centList.length;

    // WARNING: assignments is "input/output-parameter"
    // should not be null and the following should not happen
    if (assignments == null) {
      assignments = new int[m_Instances.numInstances()];
      for (int i = 0; i < assignments.length; i++) {
        assignments[i] = -1;
      }
    }

    // set assignments for all instances of this node
    for (int i = node.m_Start; i <= node.m_End; i++) {
      int instIndex = m_InstList[i];
      Instance inst = m_Instances.instance(instIndex);
      // if (instList[i] == 664) System.out.println("664***");
      int newC = m_EuclideanDistance.closestPoint(inst, centers, centList);
      // int newC = clusterProcessedInstance(inst, centers);
      assignments[instIndex] = newC;
    }
  }

  /**
   * Properties' variables =====================================================
   */

  /** flag for normalizing. */
  boolean m_NormalizeNodeWidth = true;

  /** The euclidean distance function to use. */
  protected EuclideanDistance m_EuclideanDistance;
  { // to make sure we have only one object of EuclideanDistance
    if (m_DistanceFunction instanceof EuclideanDistance)
      m_EuclideanDistance = (EuclideanDistance) m_DistanceFunction;
    else
      m_DistanceFunction = m_EuclideanDistance = new EuclideanDistance();
  }

  /** minimal relative width of a KDTree rectangle. */
  protected double m_MinBoxRelWidth = 1.0E-2;

  /** maximal number of instances in a leaf. */
  protected int m_MaxInstInLeaf = 40;

  /**
   * the GET and SET - functions ===============================================
   */

  /**
   * Tip text for this property.
   * 
   * @return 		the tip text for this property
   */
  public String minBoxRelWidthTipText() {
    return "The minimum relative width of the box. A node is only made a leaf "
        + "if the width of the split dimension of the instances in a node "
        + "normalized over the width of the split dimension of all the "
        + "instances is less than or equal to this minimum relative width.";
  }

  /**
   * Sets the minimum relative box width.
   * 
   * @param i		the minimum relative box width
   */
  public void setMinBoxRelWidth(double i) {
    m_MinBoxRelWidth = i;
  }

  /**
   * Gets the minimum relative box width.
   * 
   * @return 		the minimum relative box width
   */
  public double getMinBoxRelWidth() {
    return m_MinBoxRelWidth;
  }

  /**
   * Tip text for this property.
   * 
   * @return 		the tip text for this property
   */
  public String maxInstInLeafTipText() {
    return "The max number of instances in a leaf.";
  }

  /**
   * Sets the maximum number of instances in a leaf.
   * 
   * @param i		the maximum number of instances in a leaf
   */
  public void setMaxInstInLeaf(int i) {
    m_MaxInstInLeaf = i;
  }

  /**
   * Get the maximum number of instances in a leaf.
   * 
   * @return 		the maximum number of instances in a leaf
   */
  public int getMaxInstInLeaf() {
    return m_MaxInstInLeaf;
  }

  /**
   * Tip text for this property.
   * 
   * @return 		the tip text for this property
   */
  public String normalizeNodeWidthTipText() {
    return "Whether if the widths of the KDTree node should be normalized "
        + "by the width of the universe or not. "
        + "Where, width of the node is the range of the split attribute "
        + "based on the instances in that node, and width of the "
        + "universe is the range of the split attribute based on all the "
        + "instances (default: false).";
  }

  /**
   * Sets the flag for normalizing the widths of a KDTree Node by the width of
   * the dimension in the universe.
   * 
   * @param n		true to use normalizing.
   */
  public void setNormalizeNodeWidth(boolean n) {
    m_NormalizeNodeWidth = n;
  }

  /**
   * Gets the normalize flag.
   * 
   * @return 		True if normalizing
   */
  public boolean getNormalizeNodeWidth() {
    return m_NormalizeNodeWidth;
  }

  /**
   * returns the distance function currently in use.
   * 
   * @return 		the distance function
   */
  public DistanceFunction getDistanceFunction() {
    return (DistanceFunction) m_EuclideanDistance;
  }

  /**
   * sets the distance function to use for nearest neighbour search.
   * 
   * @param df		the distance function to use
   * @throws Exception	if not EuclideanDistance
   */
  public void setDistanceFunction(DistanceFunction df) throws Exception {
    if (!(df instanceof EuclideanDistance))
      throw new Exception("KDTree currently only works with "
          + "EuclideanDistanceFunction.");
    m_DistanceFunction = m_EuclideanDistance = (EuclideanDistance) df;
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String nodeSplitterTipText() {
    return "The the splitting method to split the nodes of the KDTree.";
  }

  /**
   * Returns the splitting method currently in use to split the nodes of the
   * KDTree.
   * 
   * @return The KDTreeNodeSplitter currently in use.
   */
  public KDTreeNodeSplitter getNodeSplitter() {
    return m_Splitter;
  }

  /**
   * Sets the splitting method to use to split the nodes of the KDTree.
   * 
   * @param splitter The KDTreeNodeSplitter to use.
   */
  public void setNodeSplitter(KDTreeNodeSplitter splitter) {
    m_Splitter = splitter;
  }

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Class implementing the KDTree search algorithm for nearest "
      + "neighbour search.\n"
      + "The connection to dataset is only a reference. For the tree "
      + "structure the indexes are stored in an array. \n"
      + "Building the tree:\n"
      + "If a node has <maximal-inst-number> (option -L) instances no "
      + "further splitting is done. Also if the split would leave one "
      + "side empty, the branch is not split any further even if the "
      + "instances in the resulting node are more than "
      + "<maximal-inst-number> instances.\n"
      + "**PLEASE NOTE:** The algorithm can not handle missing values, so it "
      + "is advisable to run ReplaceMissingValues filter if there are any "
      + "missing values in the dataset.\n\n"
      + "For more information see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an enumeration describing the available options.
   * 
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> newVector = new Vector<Option>();
    
    newVector.add(new Option(
	"\tNode splitting method to use.\n"
	+ "\t(default: weka.core.neighboursearch.kdtrees.SlidingMidPointOfWidestSide)",
	"S", 1, "-S <classname and options>"));
    
    newVector.addElement(new Option(
	"\tSet minimal width of a box\n"
        + "\t(default: 1.0E-2).", 
        "W", 0, "-W <value>"));
    
    newVector.addElement(new Option(
	"\tMaximal number of instances in a leaf\n"
        + "\t(default: 40).",
        "L", 0, "-L"));
    
    newVector.addElement(new Option(
	"\tNormalizing will be done\n"
        + "\t(Select dimension for split, with normalising to universe).",
        "N", 0, "-N"));
    
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;classname and options&gt;
   *  Node splitting method to use.
   *  (default: weka.core.neighboursearch.kdtrees.SlidingMidPointOfWidestSide)</pre>
   * 
   * <pre> -W &lt;value&gt;
   *  Set minimal width of a box
   *  (default: 1.0E-2).</pre>
   * 
   * <pre> -L
   *  Maximal number of instances in a leaf
   *  (default: 40).</pre>
   * 
   * <pre> -N
   *  Normalizing will be done
   *  (Select dimension for split, with normalising to universe).</pre>
   * 
   <!-- options-end -->
   * 
   * @param options	the list of options as an array of strings
   * @throws Exception	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    super.setOptions(options);

    String optionString = Utils.getOption('S', options);
    if (optionString.length() != 0) {
      String splitMethodSpec[] = Utils.splitOptions(optionString);
      if (splitMethodSpec.length == 0) {
        throw new Exception("Invalid DistanceFunction specification string.");
      }
      String className = splitMethodSpec[0];
      splitMethodSpec[0] = "";

      setNodeSplitter((KDTreeNodeSplitter) Utils.forName(
          KDTreeNodeSplitter.class, className, splitMethodSpec));
    }
    else {
      setNodeSplitter(new SlidingMidPointOfWidestSide());
    }

    optionString = Utils.getOption('W', options);
    if (optionString.length() != 0)
      setMinBoxRelWidth(Double.parseDouble(optionString));
    else
      setMinBoxRelWidth(1.0E-2);

    optionString = Utils.getOption('L', options);
    if (optionString.length() != 0)
      setMaxInstInLeaf(Integer.parseInt(optionString));
    else
      setMaxInstInLeaf(40);

    setNormalizeNodeWidth(Utils.getFlag('N', options));
  }

  /**
   * Gets the current settings of KDtree.
   * 
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    String[]		options;
    int			i;
    
    result = new Vector<String>();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-S");
    result.add(
	(m_Splitter.getClass().getName() + " " +
	 Utils.joinOptions(m_Splitter.getOptions())).trim());

    result.add("-W");
    result.add("" + getMinBoxRelWidth());

    result.add("-L");
    result.add("" + getMaxInstInLeaf());

    if (getNormalizeNodeWidth())
      result.add("-N");

    return result.toArray(new String[result.size()]);
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

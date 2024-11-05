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
 * BallTree.java
 * Copyright (C) 2007 University of Waikato
 */

package weka.core.neighboursearch;

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
import weka.core.neighboursearch.balltrees.BallNode;
import weka.core.neighboursearch.balltrees.BallTreeConstructor;
import weka.core.neighboursearch.balltrees.TopDownConstructor;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class implementing the BallTree/Metric Tree algorithm for nearest neighbour search.<br/>
 * The connection to dataset is only a reference. For the tree structure the indexes are stored in an array.<br/>
 * See the implementing classes of different construction methods of the trees for details on its construction.<br/>
 * <br/>
 * For more information see also:<br/>
 * <br/>
 * Stephen M. Omohundro (1989). Five Balltree Construction Algorithms.<br/>
 * <br/>
 * Jeffrey K. Uhlmann (1991). Satisfying general proximity/similarity queries with metric trees. Information Processing Letters. 40(4):175-179.
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
 * 
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
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;classname and options&gt;
 *  The construction method to employ. Either TopDown or BottomUp
 *  (default: weka.core.TopDownConstructor)</pre>
 * 
 <!-- options-end --> 
 *
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public class BallTree
  extends NearestNeighbourSearch 
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = 728763855952698328L;

  /** 
   * The instances indices sorted inorder of appearence in the tree from left
   * most leaf node to the right most leaf node.
   */
  protected int[] m_InstList;
  
  /** 
   * The maximum number of instances in a leaf. A node is made into a leaf if 
   * the number of instances in it become less than or equal to this value.
   */
  protected int m_MaxInstancesInLeaf = 40;
  
  /** Tree Stats variables. */
  protected TreePerformanceStats m_TreeStats = null;

  /** The root node of the BallTree. */
  protected BallNode m_Root;
  
  /** The constructor method to use to build the tree. */
  protected BallTreeConstructor m_TreeConstructor = new TopDownConstructor();
  
  /** Array holding the distances of the nearest neighbours. It is filled up
   *  both by nearestNeighbour() and kNearestNeighbours(). 
   */
  protected double[] m_Distances;

  /**
   * Creates a new instance of BallTree.
   */
  public BallTree() {
    super();
    if(getMeasurePerformance())
      m_Stats = m_TreeStats = new TreePerformanceStats();
  }
  
  /**
   * Creates a new instance of BallTree. 
   * It also builds the tree on supplied set of Instances.
   * @param insts The instances/points on which the BallTree 
   * should be built on.
   */
  public BallTree(Instances insts) {
    super(insts);
    if(getMeasurePerformance())
      m_Stats = m_TreeStats = new TreePerformanceStats();
  }

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Class implementing the BallTree/Metric Tree algorithm for "
      + "nearest neighbour search.\n"
      + "The connection to dataset is only a reference. For the tree "
      + "structure the indexes are stored in an array.\n"
      + "See the implementing classes of different construction methods of "
      + "the trees for details on its construction.\n\n"
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
    TechnicalInformation additional;

    result = new TechnicalInformation(Type.TECHREPORT);
    result.setValue(Field.AUTHOR, "Stephen M. Omohundro");
    result.setValue(Field.YEAR, "1989");
    result.setValue(Field.TITLE, "Five Balltree Construction Algorithms");
    result.setValue(Field.MONTH, "December");
    result.setValue(Field.NUMBER, "TR-89-063");
    result.setValue(Field.INSTITUTION, "International Computer Science Institute");

    additional = result.add(Type.ARTICLE);
    additional.setValue(Field.AUTHOR, "Jeffrey K. Uhlmann");
    additional.setValue(Field.TITLE, "Satisfying general proximity/similarity queries with metric trees");
    additional.setValue(Field.JOURNAL, "Information Processing Letters");
    additional.setValue(Field.MONTH, "November");
    additional.setValue(Field.YEAR, "1991");
    additional.setValue(Field.NUMBER, "4");
    additional.setValue(Field.VOLUME, "40");
    additional.setValue(Field.PAGES, "175-179");
    
    return result;
  }
  
  /**
   * Builds the BallTree on the supplied set of 
   * instances/points (supplied with setInstances(Instances) 
   * method and referenced by the m_Instances field). This
   * method should not be called by outside classes. They
   * should only use setInstances(Instances) method.
   * 
   * @throws Exception If no instances are supplied 
   * (m_Instances is null), or if some other error in the 
   * supplied BallTreeConstructor occurs while building 
   * the tree.  
   */
  protected void buildTree() throws Exception {
    if(m_Instances==null)
      throw new Exception("No instances supplied yet. Have to call " +
                          "setInstances(instances) with a set of Instances " +
                          "first.");
    m_InstList = new int[m_Instances.numInstances()];
    
    for(int i=0; i<m_InstList.length; i++) {
      m_InstList[i] = i;
    } //end for
    
    m_DistanceFunction.setInstances(m_Instances);
    m_TreeConstructor.setInstances(m_Instances);
    m_TreeConstructor.setInstanceList(m_InstList);
    m_TreeConstructor.setEuclideanDistanceFunction(
                      (EuclideanDistance)m_DistanceFunction);
    
    m_Root = m_TreeConstructor.buildTree();
  }
   
  /**
   * Returns k nearest instances in the current neighbourhood to the supplied
   * instance. &gt;k instances can be returned if there is more than one instance 
   * at the kth boundary (i.e. if there are more than 1 instance with the 
   * same distance as the kth nearest neighbour).
   * 
   * @param target 	The instance to find the k nearest neighbours for.
   * @param k		The number of nearest neighbours to find.
   * @throws Exception 	If the neighbours could not be found.
   * @return The k nearest neighbours of the given target instance 
   * (&gt;k nearest neighbours, if there are more instances that have same 
   *  distance as the kth nearest neighbour).
   */
  public Instances kNearestNeighbours(Instance target, int k) throws Exception {
    MyHeap heap = new MyHeap(k);

    if(m_Stats!=null)
      m_Stats.searchStart();
    
    nearestNeighbours(heap, m_Root, target, k);
    
    if(m_Stats!=null)
      m_Stats.searchFinish();

    Instances neighbours = new Instances(m_Instances, heap.totalSize());
    m_Distances = new double[heap.totalSize()];
    int [] indices = new int[heap.totalSize()];
    int i=1; MyHeapElement h;
    while(heap.noOfKthNearest()>0) {
      h = heap.getKthNearest();
      indices[indices.length-i] = h.index;
      m_Distances[indices.length-i] = h.distance;
      i++;
    }
    while(heap.size()>0) {
      h = heap.get();
      indices[indices.length-i] = h.index;
      m_Distances[indices.length-i] = h.distance;
      i++;
    }
    
    m_DistanceFunction.postProcessDistances(m_Distances);
    
    for(i=0; i<indices.length; i++)
      neighbours.add(m_Instances.instance(indices[i]));
    
    return neighbours;  // <---Check this statement
  }

  /** 
   * Does NN search according to Moore's method. 
   * Should not be used by outside classes. They should instead
   * use kNearestNeighbours(Instance, int).
   * P.S.: The distance returned are squared. Need to post process the 
   * distances. 
   * @param heap MyHeap object to store/update NNs found during the search.
   * @param node The BallNode to do the NN search on.
   * @param target The target instance for which the NNs are required.
   * @param k The number of NNs to find.
   * @throws Exception If the structure of the BallTree is not correct, 
   * or if there is some problem putting NNs in the heap.
   */
  protected void nearestNeighbours(MyHeap heap, BallNode node, Instance target,
                                   int k) throws Exception{
    double distance = Double.NEGATIVE_INFINITY;

    if (heap.totalSize() >= k)
      distance = m_DistanceFunction.distance(target, node.getPivot());

    // The radius is not squared so need to take sqrt before comparison
    if (distance > -0.000001
        && Math.sqrt(heap.peek().distance) < distance - node.getRadius()) {
      return;
    } else if (node.m_Left != null && node.m_Right != null) { // if node is not
                                                              // a leaf
      if (m_TreeStats != null) {
        m_TreeStats.incrIntNodeCount();
      }
      double leftPivotDist = Math.sqrt(m_DistanceFunction.distance(target,
          node.m_Left.getPivot(), Double.POSITIVE_INFINITY));
      double rightPivotDist = Math.sqrt(m_DistanceFunction.distance(target,
          node.m_Right.getPivot(), Double.POSITIVE_INFINITY));
      double leftBallDist = leftPivotDist - node.m_Left.getRadius();
      double rightBallDist = rightPivotDist - node.m_Right.getRadius();
      // if target is inside both balls then see which center is closer
      if (leftBallDist < 0 && rightBallDist < 0) {
        if (leftPivotDist < rightPivotDist) {
          nearestNeighbours(heap, node.m_Left, target, k);
          nearestNeighbours(heap, node.m_Right, target, k);
        } else {
          nearestNeighbours(heap, node.m_Right, target, k);
          nearestNeighbours(heap, node.m_Left, target, k);
        }
      }
      // else see which ball is closer (if dist < 0 target is inside a ball, and
      // hence the ball is closer).
      else {
        if (leftBallDist < rightBallDist) {
          nearestNeighbours(heap, node.m_Left, target, k);
          nearestNeighbours(heap, node.m_Right, target, k);
        } else {
          nearestNeighbours(heap, node.m_Right, target, k);
          nearestNeighbours(heap, node.m_Left, target, k);
        }
      }
    } else if (node.m_Left != null || node.m_Right != null) { // invalid leaves
                                                              // assignment
      throw new Exception("Error: Only one leaf of the built ball tree is " + 
                          "assigned. Please check code.");
    } else if (node.m_Left == null && node.m_Right == null) { // if node is a
                                                              // leaf
      if (m_TreeStats != null) {
        m_TreeStats.updatePointCount(node.numInstances());
        m_TreeStats.incrLeafCount();
      }
      for (int i = node.m_Start; i <= node.m_End; i++) {
        if (target == m_Instances.instance(m_InstList[i])) //for hold-one-out cross-validation
          continue;
        if (heap.totalSize() < k) {
          distance = m_DistanceFunction.distance(target, m_Instances
              .instance(m_InstList[i]), Double.POSITIVE_INFINITY, m_Stats);
          heap.put(m_InstList[i], distance);
        } else {
          MyHeapElement head = heap.peek();
          distance = m_DistanceFunction.distance(target, 
              m_Instances.instance(m_InstList[i]), head.distance, m_Stats);
          if (distance < head.distance) {
            heap.putBySubstitute(m_InstList[i], distance);
          } else if (distance == head.distance) {
            heap.putKthNearest(m_InstList[i], distance);
          }
        }//end else(heap.totalSize())
      }
    }//end else if node is a leaf
  }
  
  /**
   * Returns the nearest instance in the current neighbourhood to the supplied
   * instance.
   * 
   * @param target 	The instance to find the nearest neighbour for.
   * @throws Exception 	if the nearest neighbour could not be found.
   * @return The nearest neighbour of the given target instance.
   */
  public Instance nearestNeighbour(Instance target) throws Exception {
    return kNearestNeighbours(target, 1).instance(0);
  }

 /** 
   * Returns the distances of the k nearest neighbours. The kNearestNeighbours
   * or nearestNeighbour must always be called before calling this function. If
   * this function is called before calling either the kNearestNeighbours or 
   * the nearestNeighbour, then it throws an exception. If, however, any 
   * one of the two functions is called at any point in the past, then no   
   * exception is thrown and the distances of NN(s) from the training set for 
   * the last supplied target instance (to either one of the nearestNeighbour 
   * functions) is/are returned.
   *
   * @return 		array containing the distances of the 
   *            	nearestNeighbours. The length and ordering of the 
   *            	array is the same as that of the instances returned 
   *            	by nearestNeighbour functions.
   * @throws Exception 	if called before calling kNearestNeighbours
   *            	or nearestNeighbours.
   */
  public double[] getDistances() throws Exception {
    if(m_Distances==null)
      throw new Exception("No distances available. Please call either "+
                          "kNearestNeighbours or nearestNeighbours first.");
    return m_Distances;    
  }

  /**
   * Adds one instance to the BallTree. This involves creating/updating the 
   * structure to reflect the newly added training instance 
   * 
   * @param ins The instance to be added. Usually the newly added instance in the
   *            training set.
   * @throws Exception If the instance cannot be added to the tree.
   */
  public void update(Instance ins) throws Exception {
    addInstanceInfo(ins);
    m_InstList = m_TreeConstructor.addInstance(m_Root, ins);    
  }
  
  /** 
   * Adds the given instance's info. This implementation updates the attributes' 
   * range datastructures of EuclideanDistance class.
   * 
   * @param ins		The instance to add the information of. Usually this is
   * 			the test instance supplied to update the range of 
   * 			attributes in the distance function.
   */
  public void addInstanceInfo(Instance ins) {
    if(m_Instances!=null)
      m_DistanceFunction.update(ins);
  }  

  /**
   * Builds the BallTree based on the given set of instances.
   * @param insts The insts for which the BallTree is to be 
   * built. 
   * @throws Exception If some error occurs while 
   * building the BallTree
   */
  public void setInstances(Instances insts) throws Exception {
    super.setInstances(insts);
    buildTree();
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String ballTreeConstructorTipText() {
    return "The tree constructor being used.";
  }
      
  /**
   * Returns the BallTreeConstructor currently in use.
   * @return The BallTreeConstructor currently in use.
   */
  public BallTreeConstructor getBallTreeConstructor() {
    return m_TreeConstructor;
  }
  
  /**
   * Sets the BallTreeConstructor for building the BallTree 
   * (default TopDownConstructor).
   * @param constructor The new BallTreeConstructor. 
   */
  public void setBallTreeConstructor(BallTreeConstructor constructor) {
    m_TreeConstructor = constructor;
  }
  
  /**
   * Returns the size of the tree.
   * 
   * @return 		the size of the tree
   */
  public double measureTreeSize() {
    return m_TreeConstructor.getNumNodes();
  }
  
  /**
   * Returns the number of leaves.
   * 
   * @return 		the number of leaves
   */
  public double measureNumLeaves() {
    return m_TreeConstructor.getNumLeaves();
  }
  
  /**
   * Returns the depth of the tree. 
   * 
   * @return 		the number of rules
   */
  public double measureMaxDepth() {
    return m_TreeConstructor.getMaxDepth();
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
   * @param additionalMeasureName 	the name of the measure to query for 
   * 					its value.
   * @return 				the value of the named measure.
   * @throws IllegalArgumentException 	if the named measure is not supported.
   */
  public double getMeasure(String additionalMeasureName) {
    if (additionalMeasureName.compareToIgnoreCase("measureMaxDepth") == 0) {
      return measureMaxDepth();
    } else if (additionalMeasureName.compareToIgnoreCase("measureTreeSize") == 0) {
      return measureTreeSize();
    } else if (additionalMeasureName.compareToIgnoreCase("measureNumLeaves") == 0) {
      return measureNumLeaves();
    } else if(m_Stats!=null) {
      return m_Stats.getMeasure(additionalMeasureName);
    } else {
      throw new IllegalArgumentException(additionalMeasureName 
			  + " not supported (BallTree)");
    }
  }
  	
  /**
   * Sets whether to calculate the performance statistics or not.
   * @param measurePerformance This should be true if performance
   * statistics are to be calculated.
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
   * Returns an enumeration describing the available options.
   * 
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> newVector = new Vector<Option>();
    
    newVector.addElement(new Option(
	"\tThe construction method to employ. Either TopDown or BottomUp\n"
	+ "\t(default: weka.core.TopDownConstructor)",
	"C", 1, "-C <classname and options>"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options.
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;classname and options&gt;
   *  The construction method to employ. Either TopDown or BottomUp
   *  (default: weka.core.TopDownConstructor)</pre>
   * 
   <!-- options-end --> 
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options)
    throws Exception {

    super.setOptions(options);

    String optionString = Utils.getOption('C', options);
    if(optionString.length() != 0) {
      String constructorSpec[] = Utils.splitOptions(optionString);
      if(constructorSpec.length == 0) { 
        throw new Exception("Invalid BallTreeConstructor specification string."); 
      }
      String className = constructorSpec[0];
      constructorSpec[0] = "";

      setBallTreeConstructor( (BallTreeConstructor)
                            Utils.forName( BallTreeConstructor.class, 
                                           className, constructorSpec) );
    }
    else {
      setBallTreeConstructor(new TopDownConstructor());  
    }
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
    
    result.add("-C");
    result.add(
	(m_TreeConstructor.getClass().getName() + " " +
	 Utils.joinOptions(m_TreeConstructor.getOptions())).trim());

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

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
 * TopDownConstructor.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.balltrees;

import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.neighboursearch.balltrees.BallNode;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * The class implementing the TopDown construction method of ball trees. It further uses one of a number of different splitting methods to split a ball while constructing the tree top down.<br/>
 * <br/>
 * For more information see also:<br/>
 * <br/>
 * Stephen M. Omohundro (1989). Five Balltree Construction Algorithms.
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
 * <pre> -S &lt;classname and options&gt;
 *  Ball splitting algorithm to use.</pre>
 * 
 <!-- options-end --> 
 * 
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
public class TopDownConstructor
  extends BallTreeConstructor 
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = -5150140645091889979L;
  
  /** 
   * The BallSplitter algorithm used by the TopDown BallTree constructor, if it 
   * is selected. 
   */
  protected BallSplitter m_Splitter = new PointsClosestToFurthestChildren();

  /** 
   * Creates a new instance of TopDownConstructor.
   */
  public TopDownConstructor() {
  }

  /**
   * Returns a string describing this nearest neighbour search algorithm.
   * 
   * @return 		a description of the algorithm for displaying in the
   *         		explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "The class implementing the TopDown construction method of "
      + "ball trees. It further uses one of a number of different splitting "
      + "methods to split a ball while constructing the tree top down.\n\n"
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
    result.setValue(Field.AUTHOR, "Stephen M. Omohundro");
    result.setValue(Field.YEAR, "1989");
    result.setValue(Field.TITLE, "Five Balltree Construction Algorithms");
    result.setValue(Field.MONTH, "December");
    result.setValue(Field.NUMBER, "TR-89-063");
    result.setValue(Field.INSTITUTION, "International Computer Science Institute");

    return result;
  }

  /**
   * Builds the ball tree top down. 
   * @return The root node of the tree. 
   * @throws Exception If there is problem building
   * the tree.
   */
  public BallNode buildTree() throws Exception {
    BallNode root;
    
    m_NumNodes = m_MaxDepth = 0;
    m_NumLeaves = 1;
    
    m_Splitter.setInstances(m_Instances);
    m_Splitter.setInstanceList(m_InstList);
    m_Splitter.
    setEuclideanDistanceFunction((EuclideanDistance)m_DistanceFunction);
    
    root = new BallNode(0, m_InstList.length-1, 0);
    root.setPivot(BallNode.calcCentroidPivot(m_InstList, m_Instances));
    root.setRadius(BallNode.calcRadius(m_InstList, m_Instances, root.getPivot(), m_DistanceFunction));
    
    splitNodes(root, m_MaxDepth+1, root.m_Radius);
    
    return root; 
  }
    
  /**
   * Recursively splits nodes of a ball tree until 
   * <=m_MaxInstancesInLeaf instances remain in a node.
   * @param node The node to split.
   * @param depth The depth of this node in the tree, 
   * so that m_MaxDepth is correctly updated.
   * @param rootRadius The smallest ball enclosing all
   * the data points.
   * @throws Exception If there is some problem in 
   * splitting.
   */
  protected void splitNodes(BallNode node, int depth, final double rootRadius) throws Exception {
    
    if(node.m_NumInstances <= m_MaxInstancesInLeaf || 
       (rootRadius==0 ? true : node.m_Radius/rootRadius < m_MaxRelLeafRadius))
      return;
    
    m_NumLeaves--;
    m_Splitter.splitNode(node, m_NumNodes);
    m_NumNodes += 2;
    m_NumLeaves += 2;
    
    if(m_MaxDepth < depth)
      m_MaxDepth = depth;
  
    splitNodes(node.m_Left, depth+1, rootRadius);
    splitNodes(node.m_Right, depth+1, rootRadius);
    
    if(m_FullyContainChildBalls) {
      double radius = BallNode.calcRadius(node.m_Left, node.m_Right, 
                                         node.getPivot(), m_DistanceFunction);
      Instance pivot = BallNode.calcPivot(node.m_Left, node.m_Right, m_Instances);
//      System.err.println("Left Radius: "+node.m_Left.getRadius()+
//                         " Right Radius: "+node.m_Right.getRadius()+
//                         " d(p1,p2): "+
//                         m_DistanceFunction.distance(node.m_Left.getPivot(), node.m_Right.getPivot())+
//                         " node's old radius: "+node.getRadius()+
//                         " node's new Radius: "+radius+
//                         " node;s old pivot: "+node.getPivot()+
//                         " node's new pivot: "+pivot);
      node.setRadius(radius);
    }    
  }
    
  /**
   * Adds an instance to the ball tree. 
   * @param node The root node of the tree.
   * @param inst The instance to add to the tree.
   * @return The new master index array after adding the 
   * instance. 
   * @throws Exception If there is some problem adding the 
   * given instance to the tree. 
   */
  public int[] addInstance(BallNode node, Instance inst) throws Exception {
    
    double leftDist, rightDist;
    
    if (node.m_Left!=null && node.m_Right!=null) { //if node is not a leaf      
      // go further down the tree to look for the leaf the instance should be in
      
      leftDist = m_DistanceFunction.distance(inst, node.m_Left.getPivot(), 
                                    Double.POSITIVE_INFINITY); //instance.value(m_SplitDim);
      rightDist = m_DistanceFunction.distance(inst, node.m_Right.getPivot(), 
                                    Double.POSITIVE_INFINITY); 
      if (leftDist < rightDist) {
        addInstance(node.m_Left, inst);
        // go into right branch to correct instance list boundaries
        processNodesAfterAddInstance(node.m_Right);
      }
      else {
        addInstance(node.m_Right, inst);
      }
      // correct end index of instance list of this node
      node.m_End++;
    }
    else if(node.m_Left!=null || node.m_Right!=null) {
      throw new Exception("Error: Only one leaf of the built ball tree is " +
                          "assigned. Please check code.");
    }
    else { // found the leaf to insert instance
           
      int index = m_Instances.numInstances() - 1;
      
      int instList[] = new int[m_Instances.numInstances()];
      System.arraycopy(m_InstList, 0, instList, 0, node.m_End+1);
      if(node.m_End < m_InstList.length-1)
        System.arraycopy(m_InstList, node.m_End+2, instList, node.m_End+2, m_InstList.length-node.m_End-1);      
      instList[node.m_End+1] = index;
      node.m_End++;
      node.m_NumInstances++;
      m_InstList = instList;
      
      m_Splitter.setInstanceList(m_InstList);
      
      if(node.m_NumInstances > m_MaxInstancesInLeaf) {
        m_Splitter.splitNode(node, m_NumNodes);
        m_NumNodes += 2;
      }
    }
    return m_InstList;
  }
  
  /**
   * Post process method to correct the start and end 
   * indices of BallNodes on the right of the 
   * node where the instance was added.  
   * @param node The node whose m_Start and m_End 
   * need to be updated.
   */
  protected void processNodesAfterAddInstance(BallNode node) {
    //updating start and end indices for the node
    node.m_Start++;
    node.m_End++;    
    //processing child nodes
    if(node.m_Left!=null && node.m_Right!=null) {
      processNodesAfterAddInstance(node.m_Left);
      processNodesAfterAddInstance(node.m_Right);
    }
  }
  
  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String ballSplitterTipText() {
    return 
        "The BallSplitter algorithm set that would be used by the TopDown "
      + "BallTree constructor.";
  }
  
  /** 
   * Returns the BallSplitter algorithm set that would be 
   * used by the TopDown BallTree constructor.
   * @return The BallSplitter currently in use.
   */
  public BallSplitter getBallSplitter() {
    return m_Splitter;
  }
  
  /**
   * Sets the ball splitting algorithm to be used by the 
   * TopDown constructor.
   * @param splitter The BallSplitter to use.
   */
  public void setBallSplitter(BallSplitter splitter) {
    m_Splitter = splitter;
  }
  
  /**
   * Returns an enumeration describing the available options.
   * 
   * @return 		an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> newVector = new Vector<Option>();

    newVector.addElement(new Option(
	"\tBall splitting algorithm to use.",
	"S", 1, "-S <classname and options>"));
    
    return newVector.elements();
  }

  /**
   * Parses a given list of options.
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;classname and options&gt;
   *  Ball splitting algorithm to use.</pre>
   * 
   <!-- options-end --> 
   * 
   * @param options 	the list of options as an array of strings
   * @throws Exception	if an option is not supported
   */
  public void setOptions(String[] options)
    throws Exception {

    super.setOptions(options);
    
    String optionString = Utils.getOption('S', options);
    if(optionString.length() != 0) {
      String nnSearchClassSpec[] = Utils.splitOptions(optionString);
      if(nnSearchClassSpec.length == 0) { 
        throw new Exception("Invalid BallSplitter specification string."); 
      }
      String className = nnSearchClassSpec[0];
      nnSearchClassSpec[0] = "";

      setBallSplitter( (BallSplitter)
                            Utils.forName( BallSplitter.class, 
                                           className, nnSearchClassSpec) );
    }
    else {
      setBallSplitter(new PointsClosestToFurthestChildren());  
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
    
    result.add("-S");
    result.add(m_Splitter.getClass().getName());

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

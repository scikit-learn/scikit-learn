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
 * PointsClosestToFurthestChildren.java
 * Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 */

package weka.core.neighboursearch.balltrees;

import weka.core.EuclideanDistance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

/**
 <!-- globalinfo-start -->
 * Implements the Moore's method to split a node of a ball tree.<br/>
 * <br/>
 * For more information please see section 2 of the 1st and 3.2.3 of the 2nd:<br/>
 * <br/>
 * Andrew W. Moore: The Anchors Hierarchy: Using the Triangle Inequality to Survive High Dimensional Data. In: UAI '00: Proceedings of the 16th Conference on Uncertainty in Artificial Intelligence, San Francisco, CA, USA, 397-405, 2000.<br/>
 * <br/>
 * Ashraf Masood Kibriya (2007). Fast Algorithms for Nearest Neighbour Search. Hamilton, New Zealand.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Moore2000,
 *    address = {San Francisco, CA, USA},
 *    author = {Andrew W. Moore},
 *    booktitle = {UAI '00: Proceedings of the 16th Conference on Uncertainty in Artificial Intelligence},
 *    pages = {397-405},
 *    publisher = {Morgan Kaufmann Publishers Inc.},
 *    title = {The Anchors Hierarchy: Using the Triangle Inequality to Survive High Dimensional Data},
 *    year = {2000}
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
 <!-- options-end -->
 *
 * @author Ashraf M. Kibriya (amk14[at-the-rate]cs[dot]waikato[dot]ac[dot]nz)
 * @version $Revision: 5953 $
 */
//better rename to MidPoint of Furthest Pair/Children
public class PointsClosestToFurthestChildren
  extends BallSplitter
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = -2947177543565818260L;

  /**
   * Returns a string describing this object.
   * 
   * @return A description of the algorithm for displaying in the
   * explorer/experimenter gui.
   */
  public String globalInfo() {
    return 
        "Implements the Moore's method to split a node of a ball tree.\n\n"
      + "For more information please see section 2 of the 1st and 3.2.3 of "
      + "the 2nd:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing detailed
   * information about the technical background of this class, e.g., paper
   * reference or book this class is based on.
   * 
   * @return The technical information about this class.
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation result;
    TechnicalInformation additional;

    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Andrew W. Moore");
    result.setValue(Field.TITLE, "The Anchors Hierarchy: Using the Triangle Inequality to Survive High Dimensional Data");
    result.setValue(Field.YEAR, "2000");
    result.setValue(Field.BOOKTITLE, "UAI '00: Proceedings of the 16th Conference on Uncertainty in Artificial Intelligence");
    result.setValue(Field.PAGES, "397-405");
    result.setValue(Field.PUBLISHER, "Morgan Kaufmann Publishers Inc.");
    result.setValue(Field.ADDRESS, "San Francisco, CA, USA");

    additional = result.add(Type.MASTERSTHESIS);
    additional.setValue(Field.AUTHOR, "Ashraf Masood Kibriya");
    additional.setValue(Field.TITLE, "Fast Algorithms for Nearest Neighbour Search");
    additional.setValue(Field.YEAR, "2007");
    additional.setValue(Field.SCHOOL, "Department of Computer Science, School of Computing and Mathematical Sciences, University of Waikato");
    additional.setValue(Field.ADDRESS, "Hamilton, New Zealand");
    
    return result;
  }

  /**  Constructor. */
  public PointsClosestToFurthestChildren() {
  }
  
  /**
   * Constructor. 
   * @param instList The master index array.
   * @param insts The instances on which the tree
   * is (or is to be) built.
   * @param e The Euclidean distance function to 
   * use for splitting.
   */
  public PointsClosestToFurthestChildren(int[] instList, Instances insts, 
                                         EuclideanDistance e) {
    super(instList, insts, e);
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
    
    double maxDist = Double.NEGATIVE_INFINITY, dist = 0.0;
    Instance furthest1=null, furthest2=null, pivot=node.getPivot(), temp;
    double distList[] = new double[node.m_NumInstances];
    for(int i=node.m_Start; i<=node.m_End; i++) {
      temp = m_Instances.instance(m_Instlist[i]);
      dist = m_DistanceFunction.distance(pivot, temp, Double.POSITIVE_INFINITY);
      if(dist > maxDist) {
        maxDist = dist; furthest1 = temp;
      }
    }
    maxDist = Double.NEGATIVE_INFINITY;
    furthest1 = (Instance)furthest1.copy();
    for(int i=0; i < node.m_NumInstances; i++) {
      temp = m_Instances.instance(m_Instlist[i+node.m_Start]);
      distList[i] = m_DistanceFunction.distance(furthest1, temp, 
                                                Double.POSITIVE_INFINITY);
      if(distList[i] > maxDist) {
        maxDist = distList[i]; furthest2 = temp; //tempidx = i+node.m_Start;
      }
    }
    furthest2 = (Instance) furthest2.copy();
    dist = 0.0; int numRight=0;
    //moving indices in the right branch to the right end of the array
    for(int i=0, j=0; i < node.m_NumInstances-numRight; i++, j++) {
      temp = m_Instances.instance(m_Instlist[i+node.m_Start]);
      dist = m_DistanceFunction.distance(furthest2, temp, Double.POSITIVE_INFINITY);
      if(dist < distList[i]) {
        int t = m_Instlist[node.m_End-numRight];
        m_Instlist[node.m_End-numRight] = m_Instlist[i+node.m_Start];
        m_Instlist[i+node.m_Start] = t;
        double d = distList[distList.length-1-numRight];
        distList[distList.length-1-numRight] = distList[i];
        distList[i] = d;
        numRight++;
        i--;
      }
    }
    
    if(!(numRight > 0 && numRight < node.m_NumInstances)) 
      throw new Exception("Illegal value for numRight: "+numRight);
    
    node.m_Left = new BallNode(node.m_Start, node.m_End-numRight, numNodesCreated+1,
                              (pivot=BallNode.calcCentroidPivot(node.m_Start,
                                                node.m_End-numRight, m_Instlist, 
                                                m_Instances)), 
                              BallNode.calcRadius(node.m_Start, 
                                                node.m_End-numRight, m_Instlist, 
                                                m_Instances, pivot, 
                                                m_DistanceFunction)
                              );
    
    node.m_Right = new BallNode(node.m_End-numRight+1, node.m_End, numNodesCreated+2,
                       (pivot=BallNode.calcCentroidPivot(node.m_End-numRight+1,
                                                         node.m_End, m_Instlist, 
                                                         m_Instances)), 
                          BallNode.calcRadius(node.m_End-numRight+1, node.m_End, 
                                              m_Instlist, m_Instances, pivot, 
                                              m_DistanceFunction)
                              );
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

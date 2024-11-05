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
 *    EuclideanDistance.java
 *    Copyright (C) 1999-2007 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.neighboursearch.PerformanceStats;

/**
 <!-- globalinfo-start -->
 * Implementing Euclidean distance (or similarity) function.<br/>
 * <br/>
 * One object defines not one distance but the data model in which the distances between objects of that data model can be computed.<br/>
 * <br/>
 * Attention: For efficiency reasons the use of consistency checks (like are the data models of the two instances exactly the same), is low.<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * Wikipedia. Euclidean distance. URL http://en.wikipedia.org/wiki/Euclidean_distance.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;misc{missing_id,
 *    author = {Wikipedia},
 *    title = {Euclidean distance},
 *    URL = {http://en.wikipedia.org/wiki/Euclidean_distance}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Turns off the normalization of attribute 
 *  values in distance calculation.</pre>
 * 
 * <pre> -R &lt;col1,col2-col4,...&gt;
 *  Specifies list of columns to used in the calculation of the 
 *  distance. 'first' and 'last' are valid indices.
 *  (default: first-last)</pre>
 * 
 * <pre> -V
 *  Invert matching sense of column indices.</pre>
 * 
 <!-- options-end --> 
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class EuclideanDistance
  extends NormalizableDistance
  implements Cloneable, TechnicalInformationHandler {

  /** for serialization. */
  private static final long serialVersionUID = 1068606253458807903L;

  /**
   * Constructs an Euclidean Distance object, Instances must be still set.
   */
  public EuclideanDistance() {
    super();
  }

  /**
   * Constructs an Euclidean Distance object and automatically initializes the
   * ranges.
   * 
   * @param data 	the instances the distance function should work on
   */
  public EuclideanDistance(Instances data) {
    super(data);
  }

  /**
   * Returns a string describing this object.
   * 
   * @return 		a description of the evaluator suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Implementing Euclidean distance (or similarity) function.\n\n"
      + "One object defines not one distance but the data model in which "
      + "the distances between objects of that data model can be computed.\n\n"
      + "Attention: For efficiency reasons the use of consistency checks "
      + "(like are the data models of the two instances exactly the same), "
      + "is low.\n\n"
      + "For more information, see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return 		the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.MISC);
    result.setValue(Field.AUTHOR, "Wikipedia");
    result.setValue(Field.TITLE, "Euclidean distance");
    result.setValue(Field.URL, "http://en.wikipedia.org/wiki/Euclidean_distance");

    return result;
  }
  
  /**
   * Calculates the distance between two instances.
   * 
   * @param first 	the first instance
   * @param second 	the second instance
   * @return 		the distance between the two given instances
   */
  public double distance(Instance first, Instance second) {
    return Math.sqrt(distance(first, second, Double.POSITIVE_INFINITY));
  }
  
  /**
   * Calculates the distance (or similarity) between two instances. Need to
   * pass this returned distance later on to postprocess method to set it on
   * correct scale. <br/>
   * P.S.: Please don't mix the use of this function with
   * distance(Instance first, Instance second), as that already does post
   * processing. Please consider passing Double.POSITIVE_INFINITY as the cutOffValue to
   * this function and then later on do the post processing on all the
   * distances.
   *
   * @param first 	the first instance
   * @param second 	the second instance
   * @param stats 	the structure for storing performance statistics.
   * @return 		the distance between the two given instances or 
   * 			Double.POSITIVE_INFINITY.
   */
  public double distance(Instance first, Instance second, PerformanceStats stats) { //debug method pls remove after use
    return Math.sqrt(distance(first, second, Double.POSITIVE_INFINITY, stats));
  }
  
  /**
   * Updates the current distance calculated so far with the new difference
   * between two attributes. The difference between the attributes was 
   * calculated with the difference(int,double,double) method.
   * 
   * @param currDist	the current distance calculated so far
   * @param diff	the difference between two new attributes
   * @return		the update distance
   * @see		#difference(int, double, double)
   */
  protected double updateDistance(double currDist, double diff) {
    double	result;
    
    result  = currDist;
    result += diff * diff;
    
    return result;
  }
  
  /**
   * Does post processing of the distances (if necessary) returned by
   * distance(distance(Instance first, Instance second, double cutOffValue). It
   * is necessary to do so to get the correct distances if
   * distance(distance(Instance first, Instance second, double cutOffValue) is
   * used. This is because that function actually returns the squared distance
   * to avoid inaccuracies arising from floating point comparison.
   * 
   * @param distances	the distances to post-process
   */
  public void postProcessDistances(double distances[]) {
    for(int i = 0; i < distances.length; i++) {
      distances[i] = Math.sqrt(distances[i]);
    }
  }
  
  /**
   * Returns the squared difference of two values of an attribute.
   * 
   * @param index	the attribute index
   * @param val1	the first value
   * @param val2	the second value
   * @return		the squared difference
   */
  public double sqDifference(int index, double val1, double val2) {
    double val = difference(index, val1, val2);
    return val*val;
  }
  
  /**
   * Returns value in the middle of the two parameter values.
   * 
   * @param ranges 	the ranges to this dimension
   * @return 		the middle value
   */
  public double getMiddle(double[] ranges) {

    double middle = ranges[R_MIN] + ranges[R_WIDTH] * 0.5;
    return middle;
  }
  
  /**
   * Returns the index of the closest point to the current instance.
   * Index is index in Instances object that is the second parameter.
   *
   * @param instance 	the instance to assign a cluster to
   * @param allPoints 	all points
   * @param pointList 	the list of points
   * @return 		the index of the closest point
   * @throws Exception	if something goes wrong
   */
  public int closestPoint(Instance instance, Instances allPoints,
      			  int[] pointList) throws Exception {
    double minDist = Integer.MAX_VALUE;
    int bestPoint = 0;
    for (int i = 0; i < pointList.length; i++) {
      double dist = distance(instance, allPoints.instance(pointList[i]), Double.POSITIVE_INFINITY);
      if (dist < minDist) {
        minDist = dist;
        bestPoint = i;
      }
    }
    return pointList[bestPoint];
  }
  
  /**
   * Returns true if the value of the given dimension is smaller or equal the
   * value to be compared with.
   * 
   * @param instance 	the instance where the value should be taken of
   * @param dim 	the dimension of the value
   * @param value 	the value to compare with
   * @return 		true if value of instance is smaller or equal value
   */
  public boolean valueIsSmallerEqual(Instance instance, int dim,
      				     double value) {  //This stays
    return instance.value(dim) <= value;
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

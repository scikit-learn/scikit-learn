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
 *    ChebyshevDistance.java
 *    Copyright (C) 2007 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

/**
 <!-- globalinfo-start -->
 * Implements the Chebyshev distance. The distance between two vectors is the greatest of their differences along any coordinate dimension.<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * Wikipedia. Chebyshev distance. URL http://en.wikipedia.org/wiki/Chebyshev_distance.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;misc{missing_id,
 *    author = {Wikipedia},
 *    title = {Chebyshev distance},
 *    URL = {http://en.wikipedia.org/wiki/Chebyshev_distance}
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
 * @author Fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public class ChebyshevDistance
  extends NormalizableDistance
  implements TechnicalInformationHandler {
  
  /** for serialization. */
  private static final long serialVersionUID = -7739904999895461429L;

  /**
   * Constructs an Chebyshev Distance object, Instances must be still set.
   */
  public ChebyshevDistance() {
    super();
  }

  /**
   * Constructs an Chebyshev Distance object and automatically initializes the
   * ranges.
   * 
   * @param data 	the instances the distance function should work on
   */
  public ChebyshevDistance(Instances data) {
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
        "Implements the Chebyshev distance. The distance between two vectors "
      + "is the greatest of their differences along any coordinate dimension.\n\n"
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
    result.setValue(Field.TITLE, "Chebyshev distance");
    result.setValue(Field.URL, "http://en.wikipedia.org/wiki/Chebyshev_distance");

    return result;
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
    
    result = currDist;
    
    diff = Math.abs(diff);
    if (diff > result)
      result = diff;
    
    return result;
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

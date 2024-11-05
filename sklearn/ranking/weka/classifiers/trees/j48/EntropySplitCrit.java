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
 *    EntropySplitCrit.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 * Class for computing the entropy for a given distribution.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public final class EntropySplitCrit
  extends EntropyBasedSplitCrit {

  /** for serialization */
  private static final long serialVersionUID = 5986252682266803935L;

  /**
   * Computes entropy for given distribution.
   */
  public final double splitCritValue(Distribution bags) {
    
    return newEnt(bags);
  }

  /**
   * Computes entropy of test distribution with respect to training distribution.
   */
  public final double splitCritValue(Distribution train, Distribution test) {

    double result = 0;
    int numClasses = 0;
    int i, j;
    
    // Find out relevant number of classes
    for (j = 0; j < test.numClasses(); j++)
      if (Utils.gr(train.perClass(j), 0) || Utils.gr(test.perClass(j), 0))
	numClasses++;

    // Compute entropy of test data with respect to training data
    for (i = 0; i < test.numBags(); i++)
      if (Utils.gr(test.perBag(i),0)) {
	for (j = 0; j < test.numClasses(); j++)
	  if (Utils.gr(test.perClassPerBag(i, j), 0))
	    result -= test.perClassPerBag(i, j)*
	      Math.log(train.perClassPerBag(i, j) + 1);
	result += test.perBag(i) * Math.log(train.perBag(i) + numClasses);
      }
  
    return result / log2;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.8 $");
  }
}

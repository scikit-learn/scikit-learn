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
 *    InfoGainSplitCrit.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.RevisionUtils;
import weka.core.Utils;

/**
 * Class for computing the information gain for a given distribution.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.10 $
 */
public final class InfoGainSplitCrit
  extends EntropyBasedSplitCrit{

  /** for serialization */
  private static final long serialVersionUID = 4892105020180728499L;

  /**
   * This method is a straightforward implementation of the information
   * gain criterion for the given distribution.
   */
  public final double splitCritValue(Distribution bags) {

    double numerator;
        
    numerator = oldEnt(bags)-newEnt(bags);

    // Splits with no gain are useless.
    if (Utils.eq(numerator,0))
      return Double.MAX_VALUE;
        
    // We take the reciprocal value because we want to minimize the
    // splitting criterion's value.
    return bags.total()/numerator;
  }

  /**
   * This method computes the information gain in the same way 
   * C4.5 does.
   *
   * @param bags the distribution
   * @param totalNoInst weight of ALL instances (including the
   * ones with missing values).
   */
  public final double splitCritValue(Distribution bags, double totalNoInst) {
    
    double numerator;
    double noUnknown;
    double unknownRate;
    int i;
    
    noUnknown = totalNoInst-bags.total();
    unknownRate = noUnknown/totalNoInst;
    numerator = (oldEnt(bags)-newEnt(bags));
    numerator = (1-unknownRate)*numerator;
    
    // Splits with no gain are useless.
    if (Utils.eq(numerator,0))
      return 0;
    
    return numerator/bags.total();
  }

  /**
   * This method computes the information gain in the same way 
   * C4.5 does.
   *
   * @param bags the distribution
   * @param totalNoInst weight of ALL instances 
   * @param oldEnt entropy with respect to "no-split"-model.
   */
  public final double splitCritValue(Distribution bags,double totalNoInst,
                                     double oldEnt) {
    
    double numerator;
    double noUnknown;
    double unknownRate;
    int i;
    
    noUnknown = totalNoInst-bags.total();
    unknownRate = noUnknown/totalNoInst;
    numerator = (oldEnt-newEnt(bags));
    numerator = (1-unknownRate)*numerator;
    
    // Splits with no gain are useless.
    if (Utils.eq(numerator,0))
      return 0;
    
    return numerator/bags.total();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.10 $");
  }
}

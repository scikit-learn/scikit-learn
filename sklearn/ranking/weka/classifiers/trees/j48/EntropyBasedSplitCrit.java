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
 *    EntropyBasedSplitCrit.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

/**
 * "Abstract" class for computing splitting criteria
 * based on the entropy of a class distribution.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public abstract class EntropyBasedSplitCrit
  extends SplitCriterion {

  /** for serialization */
  private static final long serialVersionUID = -2618691439791653056L;

  /** The log of 2. */
  protected static double log2 = Math.log(2);

  /**
   * Help method for computing entropy.
   */
  public final double logFunc(double num) {

    // Constant hard coded for efficiency reasons
    if (num < 1e-6)
      return 0;
    else
      return num*Math.log(num)/log2;
  }

  /**
   * Computes entropy of distribution before splitting.
   */
  public final double oldEnt(Distribution bags) {

    double returnValue = 0;
    int j;

    for (j=0;j<bags.numClasses();j++)
      returnValue = returnValue+logFunc(bags.perClass(j));
    return logFunc(bags.total())-returnValue; 
  }

  /**
   * Computes entropy of distribution after splitting.
   */
  public final double newEnt(Distribution bags) {
    
    double returnValue = 0;
    int i,j;

    for (i=0;i<bags.numBags();i++){
      for (j=0;j<bags.numClasses();j++)
	returnValue = returnValue+logFunc(bags.perClassPerBag(i,j));
      returnValue = returnValue-logFunc(bags.perBag(i));
    }
    return -returnValue;
  }

  /**
   * Computes entropy after splitting without considering the
   * class values.
   */
  public final double splitEnt(Distribution bags) {

    double returnValue = 0;
    int i;

    for (i=0;i<bags.numBags();i++)
      returnValue = returnValue+logFunc(bags.perBag(i));
    return logFunc(bags.total())-returnValue;
  }
}


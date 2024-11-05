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
 *    SplitCriterion.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.RevisionHandler;

import java.io.Serializable;

/**
 * Abstract class for computing splitting criteria
 * with respect to distributions of class values.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public abstract class SplitCriterion
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 5490996638027101259L;

  /**
   * Computes result of splitting criterion for given distribution.
   *
   * @return value of splitting criterion. 0 by default
   */
  public double splitCritValue(Distribution bags){

    return 0;
  }

  /**
   * Computes result of splitting criterion for given training and
   * test distributions.
   *
   * @return value of splitting criterion. 0 by default
   */
  public double splitCritValue(Distribution train, Distribution test){

    return 0;
  }

  /**
   * Computes result of splitting criterion for given training and
   * test distributions and given number of classes.
   *
   * @return value of splitting criterion. 0 by default
   */
  public double splitCritValue(Distribution train, Distribution test,
			       int noClassesDefault){

    return 0;
  }

  /**
   * Computes result of splitting criterion for given training and
   * test distributions and given default distribution.
   *
   * @return value of splitting criterion. 0 by default
   */
  public double splitCritValue(Distribution train, Distribution test,
			       Distribution defC){

    return 0;
  }
}



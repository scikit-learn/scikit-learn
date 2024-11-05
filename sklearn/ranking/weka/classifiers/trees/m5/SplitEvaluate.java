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
 *    SplitEvaluate.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *    
 */

package weka.classifiers.trees.m5;

import weka.core.Instances;

/**
 * Interface for objects that determine a split point on an attribute
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.3 $
 */
public interface SplitEvaluate {
  
  /**
   * makes a copy of the SplitEvaluate object
   * @return a copy of the object
   */
  SplitEvaluate copy () throws Exception;

  /** 
   * Finds the best splitting point for an attribute in the instances
   * @param attr the splitting attribute
   * @param inst the instances
   * @exception Exception if something goes wrong
   */
   void attrSplit (int attr, Instances inst) throws Exception;

  /**
   * Returns the impurity of this split
   *
   * @return the impurity of this split
   */
   double maxImpurity();

  /**
   * Returns the position of the split in the sorted values. -1 indicates that
   * a split could not be found.
   *
   * @return an <code>int</code> value
   */
   int position();
  
  /**
   * Returns the attribute used in this split
   *
   * @return the attribute used in this split
   */
   int splitAttr();

  /**
   * Returns the split value
   *
   * @return the split value
   */
   double splitValue();

}

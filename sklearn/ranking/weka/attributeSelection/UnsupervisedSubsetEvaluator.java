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
 *    UnsupervisedSubsetEvaluator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.clusterers.Clusterer;

/** 
 * Abstract unsupervised attribute subset evaluator.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.10 $
 */
public abstract class UnsupervisedSubsetEvaluator 
  extends ASEvaluation
  implements SubsetEvaluator {

  /** for serialization */
  static final long serialVersionUID = 627934376267488763L;
  
  /**
   * Return the number of clusters used by the subset evaluator
   *
   * @return the number of clusters used
   * @exception Exception if an error occurs
   */
  public abstract int getNumClusters() throws Exception;

  /**
   * Get the clusterer
   *
   * @return the clusterer
   */
  public abstract Clusterer getClusterer();

  /**
   * Set the clusterer to use
   *
   * @param d the clusterer to use
   */
  public abstract void setClusterer(Clusterer d);
}

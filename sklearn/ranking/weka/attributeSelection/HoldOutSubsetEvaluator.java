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
 *    HoldOutSubsetEvaluator.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.Instance;
import weka.core.Instances;

import java.util.BitSet;

/** 
 * Abstract attribute subset evaluator capable of evaluating subsets with
 * respect to a data set that is distinct from that used to initialize/
 * train the subset evaluator.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $
 */
public abstract class HoldOutSubsetEvaluator 
  extends ASEvaluation
  implements SubsetEvaluator {

  /** for serialization */
  private static final long serialVersionUID = 8280529785412054174L;
  
  /**
   * Evaluates a subset of attributes with respect to a set of instances.
   * @param subset a bitset representing the attribute subset to be
   * evaluated
   * @param holdOut a set of instances (possibly seperate and distinct
   * from those use to build/train the evaluator) with which to
   * evaluate the merit of the subset
   * @return the "merit" of the subset on the holdOut data
   * @exception Exception if the subset cannot be evaluated
   */
  public abstract double evaluateSubset(BitSet subset, Instances holdOut)
    throws Exception;

  /**
   * Evaluates a subset of attributes with respect to a single instance.
   * @param subset a bitset representing the attribute subset to be
   * evaluated
   * @param holdOut a single instance (possibly not one of those used to
   * build/train the evaluator) with which to evaluate the merit of the subset
   * @param retrain true if the classifier should be retrained with respect
   * to the new subset before testing on the holdOut instance.
   * @return the "merit" of the subset on the holdOut instance
   * @exception Exception if the subset cannot be evaluated
   */
  public abstract double evaluateSubset(BitSet subset, 
					Instance holdOut,
					boolean retrain)
    throws Exception;
}

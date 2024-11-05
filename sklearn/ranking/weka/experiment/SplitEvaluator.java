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
 *    SplitEvaluator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.Instances;
import java.io.Serializable;

/**
 * Interface to objects able to generate a fixed set of results for
 * a particular split of a dataset. The set of results should contain
 * fields related to any settings of the SplitEvaluator (not including
 * the dataset name. For example, one field for the classifier used to
 * get the results, another for the classifier options, etc). <p>
 *
 * Possible implementations of SplitEvaluator: <br>
 * <ul>
 *   <li>StdClassification results
 *   <li>StdRegression results
 * </ul>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $
 */
public interface SplitEvaluator extends Serializable {
  
  /**
   * Sets a list of method names for additional measures to look for
   * in SplitEvaluators.
   * @param additionalMeasures a list of method names
   */
  void setAdditionalMeasures(String [] additionalMeasures);

  /**
   * Gets the names of each of the key columns produced for a single run.
   * The names should not contain spaces (use '_' instead for easy 
   * translation.) The number of key fields must be constant for a given 
   * SplitEvaluator.
   *
   * @return an array containing the name of each key column
   */
  String [] getKeyNames();

  /**
   * Gets the data types of each of the key columns produced for a single run.
   * The number of key fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing objects of the type of each key column. The 
   * objects should be Strings, or Doubles.
   */
  Object [] getKeyTypes();

  /**
   * Gets the names of each of the result columns produced for a single run.
   * The names should not contain spaces (use '_' instead for easy 
   * translation.) The number of result fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing the name of each result column
   */
  String [] getResultNames();

  /**
   * Gets the data types of each of the result columns produced for a 
   * single run. The number of result fields must be constant
   * for a given SplitEvaluator.
   *
   * @return an array containing objects of the type of each result column. 
   * The objects should be Strings, or Doubles.
   */
  Object [] getResultTypes();

  /**
   * Gets the key describing the current SplitEvaluator. For example
   * This may contain the name of the classifier used for classifier
   * predictive evaluation. The number of key fields must be constant
   * for a given SplitEvaluator.
   *
   * @return a value of type 'Object'
   */
  Object [] getKey();

  /**
   * Gets the results for the supplied train and test datasets.
   *
   * @param train the training Instances.
   * @param test the testing Instances.
   * @return the results stored in an array. The objects stored in
   * the array may be Strings, Doubles, or null (for the missing value).
   * @exception Exception if a problem occurs while getting the results
   */
  Object [] getResult(Instances train, Instances test) throws Exception;

  /**
   * Returns the raw output for the most recent call to getResult. Useful
   * for debugging splitEvaluators.
   * 
   * @return the raw output corresponding to the most recent call
   * to getResut
   */
  String getRawResultOutput();

} // SplitEvaluator






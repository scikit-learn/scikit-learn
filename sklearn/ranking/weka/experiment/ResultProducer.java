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
 *    ResultProducer.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.experiment;

import weka.core.Instances;
import java.io.Serializable;

/**
 * This interface defines the methods required for an object 
 * that produces results for different randomizations of a dataset. <p>
 *
 * Possible implementations of ResultProducer: <br>
 * <ul>
 *   <li>Random test/train splits
 *   <li>CrossValidation splits
 *   <li>LearningCurve splits (multiple results per run?)
 *   <li>Averaging results of other result producers
* </ul>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $
 */

public interface ResultProducer extends Serializable {
  
  /**
   * Sets the dataset that results will be obtained for.
   *
   * @param instances a value of type 'Instances'.
   */
  void setInstances(Instances instances);

  /**
   * Sets the object to send results of each run to.
   *
   * @param listener a value of type 'ResultListener'
   */
  void setResultListener(ResultListener listener);

  /**
   * Sets a list of method names for additional measures to look for
   * in SplitEvaluators.
   * @param additionalMeasures a list of method names
   */
  void setAdditionalMeasures(String [] additionalMeasures);

  /**
   * Prepare to generate results. The ResultProducer should call
   * preProcess(this) on the ResultListener it is to send results to.
   *
   * @exception Exception if an error occurs during preprocessing.
   */
  void preProcess() throws Exception;
  
  /**
   * Perform any postprocessing. When this method is called, it indicates
   * that no more requests to generate results for the current experiment
   * will be sent. The ResultProducer should call
   * preProcess(this) on the ResultListener it is to send results to.
   *
   * @exception Exception if an error occurs
   */
  void postProcess() throws Exception;
  
  /**
   * Gets the results for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Results
   * produced should be sent to the current ResultListener, but only
   * if the ResultListener says the result is required (it may already
   * have that result). A single run may produce multiple results.
   *
   * @param run the run number to generate results for.
   * @exception Exception if a problem occurs while getting the results
   */
  void doRun(int run) throws Exception;
  
  /**
   * Gets the keys for a specified run number. Different run
   * numbers correspond to different randomizations of the data. Keys
   * produced should be sent to the current ResultListener
   *
   * @param run the run number to get keys for.
   * @exception Exception if a problem occurs while getting the keys
   */
  void doRunKeys(int run) throws Exception;

  /**
   * Gets the names of each of the key columns produced for a single run.
   * The names should not contain spaces (use '_' instead for easy 
   * translation.)
   *
   * @return an array containing the name of each key column
   * @exception Exception if the key names could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  String [] getKeyNames() throws Exception;

  /**
   * Gets the data types of each of the key columns produced for a single run.
   *
   * @return an array containing objects of the type of each key column. The 
   * objects should be Strings, or Doubles.
   * @exception Exception if the key types could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  Object [] getKeyTypes() throws Exception;

  /**
   * Gets the names of each of the result columns produced for a single run.
   * The names should not contain spaces (use '_' instead for easy 
   * translation.)
   *
   * @return an array containing the name of each result column
   * @exception Exception if the result names could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  String [] getResultNames() throws Exception;

  /**
   * Gets the data types of each of the result columns produced for a 
   * single run.
   *
   * @return an array containing objects of the type of each result column. 
   * The objects should be Strings, or Doubles.
   * @exception Exception if the result types could not be determined (perhaps
   * because of a problem from a nested sub-resultproducer)
   */
  Object [] getResultTypes() throws Exception;

  /**
   * Gets a description of the internal settings of the result
   * producer, sufficient for distinguishing a ResultProducer
   * instance from another with different settings (ignoring
   * those settings set through this interface). For example,
   * a cross-validation ResultProducer may have a setting for the
   * number of folds. For a given state, the results produced should
   * be compatible. Typically if a ResultProducer is an OptionHandler,
   * this string will represent those command line arguments required
   * to set the ResultProducer to that state.
   *
   * @return the description of the ResultProducer state, or null
   * if no state is defined
   */
  String getCompatibilityState();

} // ResultProducer

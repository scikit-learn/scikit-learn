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
 *    ASEvaluation.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.attributeSelection;

import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SerializedObject;
import weka.core.Utils;

import java.io.Serializable;

/** 
 * Abstract attribute selection evaluation class
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 6624 $
 */
public abstract class ASEvaluation
  implements Serializable, CapabilitiesHandler, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 2091705669885950849L;
  
  // ===============
  // Public methods.
  // ===============

  /**
   * Generates a attribute evaluator. Has to initialize all fields of the 
   * evaluator that are not being set via options.
   *
   * @param data set of instances serving as training data 
   * @exception Exception if the evaluator has not been 
   * generated successfully
   */
  public abstract void buildEvaluator(Instances data) throws Exception;

  /**
   * Provides a chance for a attribute evaluator to do any special
   * post processing of the selected attribute set.
   *
   * @param attributeSet the set of attributes found by the search
   * @return a possibly ranked list of postprocessed attributes
   * @exception Exception if postprocessing fails for some reason
   */
  public int [] postProcess(int [] attributeSet) 
    throws Exception {
    return attributeSet;
  }

  /**
   * Creates a new instance of an attribute/subset evaluator 
   * given it's class name and
   * (optional) arguments to pass to it's setOptions method. If the
   * evaluator implements OptionHandler and the options parameter is
   * non-null, the evaluator will have it's options set.
   *
   * @param evaluatorName the fully qualified class name of the evaluator
   * @param options an array of options suitable for passing to setOptions. May
   * be null.
   * @return the newly created evaluator, ready for use.
   * @exception Exception if the evaluator name is invalid, or the options
   * supplied are not acceptable to the evaluator
   */
  public static ASEvaluation forName(String evaluatorName,
				     String [] options) throws Exception {
    return (ASEvaluation)Utils.forName(ASEvaluation.class,
				       evaluatorName,
				       options);
  }

  /**
   * Creates copies of the current evaluator. Note that this method
   * now uses Serialization to perform a deep copy, so the evaluator
   * object must be fully Serializable. Any currently built model will
   * now be copied as well.
   *
   * @param model an example evaluator to copy
   * @param num the number of evaluator copies to create.
   * @return an array of evaluators.
   * @exception Exception if an error occurs 
   */
  public static ASEvaluation [] makeCopies(ASEvaluation model,
					 int num) throws Exception {

    if (model == null) {
      throw new Exception("No model evaluator set");
    }
    ASEvaluation [] evaluators = new ASEvaluation [num];
    SerializedObject so = new SerializedObject(model);
    for(int i = 0; i < evaluators.length; i++) {
      evaluators[i] = (ASEvaluation) so.getObject();
    }
    return evaluators;
  }

  /**
   * Returns the capabilities of this evaluator.
   *
   * @return            the capabilities of this evaluator
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = new Capabilities(this);
    result.enableAll();
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6624 $");
  }
  
  /**
   * runs the evaluator with the given commandline options
   * 
   * @param evaluator	the evaluator to run
   * @param options	the commandline options
   */
  public static void runEvaluator(ASEvaluation evaluator, String[] options) {
    try {
      System.out.println(
	  AttributeSelection.SelectAttributes(evaluator, options));
    }
    catch (Exception e) {
      String msg = e.toString().toLowerCase();
      if (    (msg.indexOf("help requested") == -1)
           && (msg.indexOf("no training file given") == -1) )
        e.printStackTrace();
      System.err.println(e.getMessage());
    }
  }
}

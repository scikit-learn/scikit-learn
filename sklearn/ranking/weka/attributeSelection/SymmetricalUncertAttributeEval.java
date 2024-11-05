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
 *    SymmetricalUncertAttributeEval.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.Capabilities;
import weka.core.ContingencyTables;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.supervised.attribute.Discretize;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * SymmetricalUncertAttributeEval :<br/>
 * <br/>
 * Evaluates the worth of an attribute by measuring the symmetrical uncertainty with respect to the class. <br/>
 * <br/>
 *  SymmU(Class, Attribute) = 2 * (H(Class) - H(Class | Attribute)) / H(Class) + H(Attribute).<br/>
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -M
 *  treat missing values as a seperate value.</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5447 $
 * @see Discretize
 */
public class SymmetricalUncertAttributeEval
  extends ASEvaluation
  implements AttributeEvaluator, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -8096505776132296416L;

  /** The training instances */
  private Instances m_trainInstances;

  /** The class index */
  private int m_classIndex;

  /** The number of attributes */
  private int m_numAttribs;

  /** The number of instances */
  private int m_numInstances;

  /** The number of classes */
  private int m_numClasses;

  /** Treat missing values as a seperate value */
  private boolean m_missing_merge;

  /**
   * Returns a string describing this attribute evaluator
   * @return a description of the evaluator suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "SymmetricalUncertAttributeEval :\n\nEvaluates the worth of an attribute "
      +"by measuring the symmetrical uncertainty with respect to the class. "
      +"\n\n SymmU(Class, Attribute) = 2 * (H(Class) - H(Class | Attribute)) "
      +"/ H(Class) + H(Attribute).\n";
  }
  
  /**
   * Constructor
   */
  public SymmetricalUncertAttributeEval () {
    resetOptions();
  }


  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   **/
  public Enumeration listOptions () {
    Vector newVector = new Vector(1);
    newVector.addElement(new Option("\ttreat missing values as a seperate " 
				    + "value.", "M", 0, "-M"));
    return  newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -M
   *  treat missing values as a seperate value.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   **/
  public void setOptions (String[] options)
    throws Exception {
    resetOptions();
    setMissingMerge(!(Utils.getFlag('M', options)));
  }

  /**
   * Returns the tip text for this property
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String missingMergeTipText() {
    return "Distribute counts for missing values. Counts are distributed "
      +"across other values in proportion to their frequency. Otherwise, "
      +"missing is treated as a separate value.";
  }

  /**
   * distribute the counts for missing values across observed values
   *
   * @param b true=distribute missing values.
   */
  public void setMissingMerge (boolean b) {
    m_missing_merge = b;
  }


  /**
   * get whether missing values are being distributed or not
   *
   * @return true if missing values are being distributed.
   */
  public boolean getMissingMerge () {
    return  m_missing_merge;
  }


  /**
   * Gets the current settings of WrapperSubsetEval.
   * @return an array of strings suitable for passing to setOptions()
   */
  public String[] getOptions () {
    String[] options = new String[1];
    int current = 0;

    if (!getMissingMerge()) {
      options[current++] = "-M";
    }

    while (current < options.length) {
      options[current++] = "";
    }

    return  options;
  }

  /**
   * Returns the capabilities of this evaluator.
   *
   * @return            the capabilities of this evaluator
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();
    
    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Initializes a symmetrical uncertainty attribute evaluator. 
   * Discretizes all attributes that are numeric.
   *
   * @param data set of instances serving as training data 
   * @throws Exception if the evaluator has not been 
   * generated successfully
   */
  public void buildEvaluator (Instances data)
    throws Exception {

    // can evaluator handle data?
    getCapabilities().testWithFail(data);

    m_trainInstances = data;
    m_classIndex = m_trainInstances.classIndex();
    m_numAttribs = m_trainInstances.numAttributes();
    m_numInstances = m_trainInstances.numInstances();
    Discretize disTransform = new Discretize();
    disTransform.setUseBetterEncoding(true);
    disTransform.setInputFormat(m_trainInstances);
    m_trainInstances = Filter.useFilter(m_trainInstances, disTransform);
    m_numClasses = m_trainInstances.attribute(m_classIndex).numValues();
  }


  /**
   * set options to default values
   */
  protected void resetOptions () {
    m_trainInstances = null;
    m_missing_merge = true;
  }


  /**
   * evaluates an individual attribute by measuring the symmetrical
   * uncertainty between it and the class.
   *
   * @param attribute the index of the attribute to be evaluated
   * @return the uncertainty
   * @throws Exception if the attribute could not be evaluated
   */
  public double evaluateAttribute (int attribute)
    throws Exception {
    int i, j, ii, jj;
    int nnj, nni, ni, nj;
    double sum = 0.0;
    ni = m_trainInstances.attribute(attribute).numValues() + 1;
    nj = m_numClasses + 1;
    double[] sumi, sumj;
    Instance inst;
    double temp = 0.0;
    sumi = new double[ni];
    sumj = new double[nj];
    double[][] counts = new double[ni][nj];
    sumi = new double[ni];
    sumj = new double[nj];

    for (i = 0; i < ni; i++) {
      sumi[i] = 0.0;

      for (j = 0; j < nj; j++) {
	sumj[j] = 0.0;
	counts[i][j] = 0.0;
      }
    }

    // Fill the contingency table
    for (i = 0; i < m_numInstances; i++) {
      inst = m_trainInstances.instance(i);

      if (inst.isMissing(attribute)) {
	ii = ni - 1;
      }
      else {
	ii = (int)inst.value(attribute);
      }

      if (inst.isMissing(m_classIndex)) {
	jj = nj - 1;
      }
      else {
	jj = (int)inst.value(m_classIndex);
      }

      counts[ii][jj]++;
    }

    // get the row totals
    for (i = 0; i < ni; i++) {
      sumi[i] = 0.0;

      for (j = 0; j < nj; j++) {
	sumi[i] += counts[i][j];
	sum += counts[i][j];
      }
    }

    // get the column totals
    for (j = 0; j < nj; j++) {
      sumj[j] = 0.0;

      for (i = 0; i < ni; i++) {
	sumj[j] += counts[i][j];
      }
    }

    // distribute missing counts
    if (m_missing_merge && 
	(sumi[ni-1] < m_numInstances) && 
	(sumj[nj-1] < m_numInstances)) {
      double[] i_copy = new double[sumi.length];
      double[] j_copy = new double[sumj.length];
      double[][] counts_copy = new double[sumi.length][sumj.length];

      for (i = 0; i < ni; i++) {
	System.arraycopy(counts[i], 0, counts_copy[i], 0, sumj.length);
      }

      System.arraycopy(sumi, 0, i_copy, 0, sumi.length);
      System.arraycopy(sumj, 0, j_copy, 0, sumj.length);
      double total_missing = (sumi[ni - 1] + sumj[nj - 1] 
			      - counts[ni - 1][nj - 1]);

      // do the missing i's
      if (sumi[ni - 1] > 0.0) {
	for (j = 0; j < nj - 1; j++) {
	  if (counts[ni - 1][j] > 0.0) {
	    for (i = 0; i < ni - 1; i++) {
	      temp = ((i_copy[i]/(sum - i_copy[ni - 1])) * 
		      counts[ni - 1][j]);
	      counts[i][j] += temp;
	      sumi[i] += temp;
	    }

	    counts[ni - 1][j] = 0.0;
	  }
	}
      }

      sumi[ni - 1] = 0.0;

      // do the missing j's
      if (sumj[nj - 1] > 0.0) {
	for (i = 0; i < ni - 1; i++) {
	  if (counts[i][nj - 1] > 0.0) {
	    for (j = 0; j < nj - 1; j++) {
	      temp = ((j_copy[j]/(sum - j_copy[nj - 1]))*counts[i][nj - 1]);
	      counts[i][j] += temp;
	      sumj[j] += temp;
	    }

	    counts[i][nj - 1] = 0.0;
	  }
	}
      }

      sumj[nj - 1] = 0.0;

      // do the both missing
      if (counts[ni - 1][nj - 1] > 0.0 && total_missing != sum) {
	for (i = 0; i < ni - 1; i++) {
	  for (j = 0; j < nj - 1; j++) {
	    temp = (counts_copy[i][j]/(sum - total_missing)) * 
	      counts_copy[ni - 1][nj - 1];
	    counts[i][j] += temp;
	    sumi[i] += temp;
	    sumj[j] += temp;
	  }
	}

	counts[ni - 1][nj - 1] = 0.0;
      }
    }

    return  ContingencyTables.symmetricalUncertainty(counts);
  }


  /**
   * Return a description of the evaluator
   * @return description as a string
   */
  public String toString () {
    StringBuffer text = new StringBuffer();

    if (m_trainInstances == null) {
      text.append("\tSymmetrical Uncertainty evaluator has not been built");
    }
    else {
      text.append("\tSymmetrical Uncertainty Ranking Filter");
      if (!m_missing_merge) {
	text.append("\n\tMissing values treated as seperate");
      }
    }

    text.append("\n");
    return  text.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5447 $");
  }

  // ============
  // Test method.
  // ============
  /**
   * Main method for testing this class.
   *
   * @param argv should contain the following arguments:
   * -t training file
   */
  public static void main (String[] argv) {
    runEvaluator(new SymmetricalUncertAttributeEval(), argv);
  }
}

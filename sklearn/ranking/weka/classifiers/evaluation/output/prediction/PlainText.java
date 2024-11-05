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
 * PlainText.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.evaluation.output.prediction;

import weka.classifiers.Classifier;
import weka.core.Instance;
import weka.core.Utils;
import weka.core.labelranking.PreferenceAttribute;

/**
 <!-- globalinfo-start -->
 * Outputs the predictions in plain text.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -p &lt;range&gt;
 *  The range of attributes to print in addition to the classification.
 *  (default: none)</pre>
 * 
 * <pre> -distribution
 *  Whether to turn on the output of the class distribution.
 *  Only for nominal class attributes.
 *  (default: off)</pre>
 * 
 * <pre> -decimals &lt;num&gt;
 *  The number of digits after the decimal point.
 *  (default: 3)</pre>
 * 
 * <pre> -file &lt;path&gt;
 *  The file to store the output in, instead of outputting it on stdout.
 *  Gets ignored if the supplied path is a directory.
 *  (default: .)</pre>
 * 
 * <pre> -suppress
 *  In case the data gets stored in a file, then this flag can be used
 *  to suppress the regular output.
 *  (default: not suppressed)</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6804 $
 */
public class PlainText
  extends AbstractOutput {
  
  /** for serialization. */
  private static final long serialVersionUID = 2033389864898242735L;
  
  /**
   * Returns a string describing the output generator.
   * 
   * @return 		a description suitable for
   * 			displaying in the GUI
   */
  public String globalInfo() {
    return "Outputs the predictions in plain text.";
  }
  
  /**
   * Returns a short display text, to be used in comboboxes.
   * 
   * @return 		a short display text
   */
  public String getDisplay() {
    return "Plain text";
  }

  /**
   * Performs the actual printing of the header.
   */
  protected void doPrintHeader() {
    if (m_Header.classAttribute().isNominal() || m_Header.classAttribute().isRanking())
      if (m_OutputDistribution)
	append(" inst#     actual  predicted error distribution");
      else
	append(" inst#     actual  predicted error prediction");
    else
      append(" inst#     actual  predicted      error");
    
    if (m_Attributes != null) {
      append(" (");
      boolean first = true;
      for (int i = 0; i < m_Header.numAttributes(); i++) {
        if (i == m_Header.classIndex())
          continue;

        if (m_Attributes.isInRange(i)) {
          if (!first)
            append(",");
          append(m_Header.attribute(i).name());
          first = false;
        }
      }
      append(")");
    }
    
    append("\n");
  }

  /**
   * Builds a string listing the attribute values in a specified range of indices,
   * separated by commas and enclosed in brackets.
   *
   * @param instance 	the instance to print the values from
   * @return 		a string listing values of the attributes in the range
   */
  protected String attributeValuesString(Instance instance) {
    StringBuffer text = new StringBuffer();
    if (m_Attributes != null) {
      boolean firstOutput = true;
      m_Attributes.setUpper(instance.numAttributes() - 1);
      for (int i=0; i<instance.numAttributes(); i++)
	if (m_Attributes.isInRange(i) && i != instance.classIndex()) {
	  if (firstOutput) text.append("(");
	  else text.append(",");
	  text.append(instance.toString(i));
	  firstOutput = false;
	}
      if (!firstOutput) text.append(")");
    }
    return text.toString();
  }

  /**
   * Store the prediction made by the classifier as a string.
   * 
   * @param classifier	the classifier to use
   * @param inst	the instance to generate text from
   * @param index	the index in the dataset
   * @throws Exception	if something goes wrong
   */
  protected void doPrintClassification(Classifier classifier, Instance inst, int index) throws Exception {
    int width = 7 + m_NumDecimals;
    int prec = m_NumDecimals;

    Instance withMissing = (Instance)inst.copy();
    withMissing.setDataset(inst.dataset());    
    inst = preProcessInstance(inst, withMissing, classifier);
    
    double predValue = classifier.classifyInstance(withMissing);

    // index
    append(Utils.padLeft("" + (index+1), 6));

    if (inst.dataset().classAttribute().isNumeric()) {
      // actual
      if (inst.classIsMissing())
	append(" " + Utils.padLeft("?", width));
      else
	append(" " + Utils.doubleToString(inst.classValue(), width, prec));
      // predicted
      if (Utils.isMissingValue(predValue))
	append(" " + Utils.padLeft("?", width));
      else
	append(" " + Utils.doubleToString(predValue, width, prec));
      // error
      if (Utils.isMissingValue(predValue) || inst.classIsMissing())
	append(" " + Utils.padLeft("?", width));
      else
	append(" " + Utils.doubleToString(predValue - inst.classValue(), width, prec));
    } 

    else {

//      else{
      // actual
      append(" " + Utils.padLeft(((int) inst.classValue()+1) + ":" + inst.toString(inst.classIndex()), width));
      // predicted
      if (Utils.isMissingValue(predValue))
	append(" " + Utils.padLeft("?", width));
      else
	append(" " + Utils.padLeft(((int) predValue+1) + ":" + inst.dataset().classAttribute().value((int)predValue), width));
      // error?
      if (!Utils.isMissingValue(predValue) && !inst.classIsMissing() && ((int) predValue+1 != (int) inst.classValue()+1))
	append(" " + "  +  ");
      else
	append(" " + "     ");
      // prediction/distribution
      if (m_OutputDistribution) {
	if (Utils.isMissingValue(predValue)) {
	  append(" " + "?");
	}
	else {
	  append(" ");
	  double[] dist = classifier.distributionForInstance(withMissing);
	  for (int n = 0; n < dist.length; n++) {
	    if (n > 0)
	      append(",");
	    if (n == (int) predValue)
	      append("*");
            append(Utils.doubleToString(dist[n], prec));
	  }
	}
      }
      else {
	if (Utils.isMissingValue(predValue))
	  append(" " + "?");
	else
	  append(" " + Utils.doubleToString(classifier.distributionForInstance(withMissing) [(int)predValue], prec));
      }
  }
    // attributes
    append(" " + attributeValuesString(withMissing) + "\n");
  }
  
  /**
   * Does nothing.
   */
  protected void doPrintFooter() {
  }
}

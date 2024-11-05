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
 * XML.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.evaluation.output.prediction;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Utils;
import weka.core.Version;
import weka.core.labelranking.PreferenceAttribute;
import weka.core.xml.XMLDocument;

/**
 <!-- globalinfo-start -->
 * Outputs the predictions in XML.<br/>
 * <br/>
 * The following DTD is used:<br/>
 * <br/>
 * &lt;!DOCTYPE predictions<br/>
 * [<br/>
 *   &lt;!ELEMENT predictions (prediction*)&gt;<br/>
 *   &lt;!ATTLIST predictions version CDATA "3.5.8"&gt;<br/>
 *   &lt;!ATTLIST predictions name CDATA #REQUIRED&gt;<br/>
 * <br/>
 *   &lt;!ELEMENT prediction ((actual_label,predicted_label,error,(prediction|distribution),attributes?)|(actual_value,predicted_value,error,attributes?))&gt;<br/>
 *   &lt;!ATTLIST prediction index CDATA #REQUIRED&gt;<br/>
 * <br/>
 *   &lt;!ELEMENT actual_label ANY&gt;<br/>
 *   &lt;!ATTLIST actual_label index CDATA #REQUIRED&gt;<br/>
 *   &lt;!ELEMENT predicted_label ANY&gt;<br/>
 *   &lt;!ATTLIST predicted_label index CDATA #REQUIRED&gt;<br/>
 *   &lt;!ELEMENT error ANY&gt;<br/>
 *   &lt;!ELEMENT prediction ANY&gt;<br/>
 *   &lt;!ELEMENT distribution (class_label+)&gt;<br/>
 *   &lt;!ELEMENT class_label ANY&gt;<br/>
 *   &lt;!ATTLIST class_label index CDATA #REQUIRED&gt;<br/>
 *   &lt;!ATTLIST class_label predicted (yes|no) "no"&gt;<br/>
 *   &lt;!ELEMENT actual_value ANY&gt;<br/>
 *   &lt;!ELEMENT predicted_value ANY&gt;<br/>
 *   &lt;!ELEMENT attributes (attribute+)&gt;<br/>
 *   &lt;!ELEMENT attribute ANY&gt;<br/>
 *   &lt;!ATTLIST attribute index CDATA #REQUIRED&gt;<br/>
 *   &lt;!ATTLIST attribute name CDATA #REQUIRED&gt;<br/>
 *   &lt;!ATTLIST attribute type (numeric|date|nominal|string|relational) #REQUIRED&gt;<br/>
 * ]<br/>
 * &gt;
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
public class XML
  extends AbstractOutput {
  
  /** for serialization. */
  private static final long serialVersionUID = -3165514277316824801L;

  /** the DocType definition. */
  public final static String DTD_DOCTYPE = XMLDocument.DTD_DOCTYPE;
  
  /** the Element definition. */
  public final static String DTD_ELEMENT = XMLDocument.DTD_ELEMENT;
  
  /** the AttList definition. */
  public final static String DTD_ATTLIST = XMLDocument.DTD_ATTLIST;
  
  /** the optional marker. */
  public final static String DTD_OPTIONAL = XMLDocument.DTD_OPTIONAL;
  
  /** the at least one marker. */
  public final static String DTD_AT_LEAST_ONE = XMLDocument.DTD_AT_LEAST_ONE;
  
  /** the zero or more marker. */
  public final static String DTD_ZERO_OR_MORE = XMLDocument.DTD_ZERO_OR_MORE;
  
  /** the option separator. */
  public final static String DTD_SEPARATOR = XMLDocument.DTD_SEPARATOR;
  
  /** the CDATA placeholder. */
  public final static String DTD_CDATA = XMLDocument.DTD_CDATA; 
  
  /** the ANY placeholder. */
  public final static String DTD_ANY = XMLDocument.DTD_ANY; 
  
  /** the #PCDATA placeholder. */
  public final static String DTD_PCDATA = XMLDocument.DTD_PCDATA; 
  
  /** the #IMPLIED placeholder. */
  public final static String DTD_IMPLIED = XMLDocument.DTD_IMPLIED; 
  
  /** the #REQUIRED placeholder. */
  public final static String DTD_REQUIRED = XMLDocument.DTD_REQUIRED; 

  /** the "version" attribute. */
  public final static String ATT_VERSION = XMLDocument.ATT_VERSION;
 
  /** the "name" attribute. */
  public final static String ATT_NAME = XMLDocument.ATT_NAME;
  
  /** the "type" attribute. */
  public final static String ATT_TYPE = "type";

  /** the value "yes". */
  public final static String VAL_YES = XMLDocument.VAL_YES;
  
  /** the value "no". */
  public final static String VAL_NO = XMLDocument.VAL_NO;
  
  /** the predictions tag. */
  public final static String TAG_PREDICTIONS = "predictions";
  
  /** the prediction tag. */
  public final static String TAG_PREDICTION = "prediction";

  /** the actual_nominal tag. */
  public final static String TAG_ACTUAL_LABEL = "actual_label";

  /** the predicted_nominal tag. */
  public final static String TAG_PREDICTED_LABEL = "predicted_label";

  /** the error tag. */
  public final static String TAG_ERROR = "error";

  /** the distribution tag. */
  public final static String TAG_DISTRIBUTION = "distribution";

  /** the class_label tag. */
  public final static String TAG_CLASS_LABEL = "class_label";

  /** the actual_numeric tag. */
  public final static String TAG_ACTUAL_VALUE = "actual_value";

  /** the predicted_numeric tag. */
  public final static String TAG_PREDICTED_VALUE = "predicted_value";

  /** the attributes tag. */
  public final static String TAG_ATTRIBUTES = "attributes";

  /** the attribute tag. */
  public final static String TAG_ATTRIBUTE = "attribute";

  /** the index attribute. */
  public final static String ATT_INDEX = "index";

  /** the predicted attribute. */
  public final static String ATT_PREDICTED = "predicted";
  
  /** the DTD. */
  public final static String DTD = 
    "<!" + DTD_DOCTYPE + " " + TAG_PREDICTIONS + "\n"
    + "[\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_PREDICTIONS + " (" + TAG_PREDICTION + DTD_ZERO_OR_MORE + ")" + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_PREDICTIONS + " " + ATT_VERSION + " " + DTD_CDATA + " \"" + Version.VERSION + "\"" + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_PREDICTIONS + " " + ATT_NAME + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_PREDICTION + " " 
             + "(" 
             + "(" + TAG_ACTUAL_LABEL + "," + TAG_PREDICTED_LABEL + "," + TAG_ERROR + "," + "(" + TAG_PREDICTION + DTD_SEPARATOR + TAG_DISTRIBUTION + ")" + "," + TAG_ATTRIBUTES + DTD_OPTIONAL + ")" 
             + DTD_SEPARATOR
             + "(" + TAG_ACTUAL_VALUE + "," + TAG_PREDICTED_VALUE + "," + TAG_ERROR + "," + TAG_ATTRIBUTES + DTD_OPTIONAL + ")"
             + ")" + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_PREDICTION + " " + ATT_INDEX + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_ACTUAL_LABEL + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_ACTUAL_LABEL + " " + ATT_INDEX + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_PREDICTED_LABEL + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_PREDICTED_LABEL + " " + ATT_INDEX + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_ERROR + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_PREDICTION + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_DISTRIBUTION + " (" + TAG_CLASS_LABEL + DTD_AT_LEAST_ONE + ")" + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_CLASS_LABEL + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_CLASS_LABEL + " " + ATT_INDEX + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_CLASS_LABEL + " " + ATT_PREDICTED + " (" + VAL_YES + DTD_SEPARATOR + VAL_NO + ") " + "\"" + VAL_NO + "\"" + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_ACTUAL_VALUE + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_PREDICTED_VALUE + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_ATTRIBUTES + " (" + TAG_ATTRIBUTE + DTD_AT_LEAST_ONE + ")" + ">\n"
    + "  <!" + DTD_ELEMENT + " " + TAG_ATTRIBUTE + " " + DTD_ANY + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_INDEX + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_NAME + " " + DTD_CDATA + " " + DTD_REQUIRED + ">\n"
    + "  <!" + DTD_ATTLIST + " " + TAG_ATTRIBUTE + " " + ATT_TYPE + " " + "(" + Attribute.typeToString(Attribute.NUMERIC) + DTD_SEPARATOR + Attribute.typeToString(Attribute.DATE) + DTD_SEPARATOR + Attribute.typeToString(Attribute.NOMINAL) + DTD_SEPARATOR + Attribute.typeToString(PreferenceAttribute.RANKING) + DTD_SEPARATOR + Attribute.typeToString(Attribute.STRING) + DTD_SEPARATOR + Attribute.typeToString(Attribute.RELATIONAL) + ")" + " " + DTD_REQUIRED + ">\n"
    + "]\n"
    + ">";
  
  /**
   * Returns a string describing the output generator.
   * 
   * @return 		a description suitable for
   * 			displaying in the GUI
   */
  public String globalInfo() {
    return 
        "Outputs the predictions in XML.\n\n"
      + "The following DTD is used:\n\n"
      + DTD;
  }
  
  /**
   * Returns a short display text, to be used in comboboxes.
   * 
   * @return 		a short display text
   */
  public String getDisplay() {
    return "XML";
  }

  /**
   * Replaces certain characters with their XML entities.
   * 
   * @param s		the string to process
   * @return		the processed string
   */
  protected String sanitize(String s) {
    String 	result;
    
    result = s;
    result = result.replaceAll("&", "&amp;");
    result = result.replaceAll("<", "&lt;");
    result = result.replaceAll(">", "&gt;");
    result = result.replaceAll("\"", "&quot;");
    
    return result;
  }
  
  /**
   * Performs the actual printing of the header.
   */
  protected void doPrintHeader() {
    append("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    append("\n");
    append(DTD + "\n\n");
    append("<" + TAG_PREDICTIONS + " " + ATT_VERSION + "=\"" + Version.VERSION + "\"" + " " + ATT_NAME + "=\"" + sanitize(m_Header.relationName()) + "\">\n");
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
      text.append("    <" + TAG_ATTRIBUTES + ">\n");
      m_Attributes.setUpper(instance.numAttributes() - 1);
      for (int i=0; i<instance.numAttributes(); i++) {
	if (m_Attributes.isInRange(i) && i != instance.classIndex()) {
	  text.append("      <" + TAG_ATTRIBUTE + " " + ATT_INDEX + "=\"" + (i+1) + "\"" + " " + ATT_NAME + "=\"" + sanitize(instance.attribute(i).name()) + "\"" + " " + ATT_TYPE + "=\"" + Attribute.typeToString(instance.attribute(i).type()) + "\"" + ">");
	  text.append(sanitize(instance.toString(i)));
	  text.append("</" + TAG_ATTRIBUTE + ">\n");
	}
      }
      text.append("    </" + TAG_ATTRIBUTES + ">\n");
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
    int prec = m_NumDecimals;

    Instance withMissing = (Instance)inst.copy();
    withMissing.setDataset(inst.dataset());
    inst = preProcessInstance(inst, withMissing, classifier);
    
    double predValue = classifier.classifyInstance(withMissing);

    // opening tag
    append("  <" + TAG_PREDICTION + " " + ATT_INDEX + "=\"" + (index+1) + "\">\n");

    if (inst.dataset().classAttribute().isNumeric()) {
      // actual
      append("    <" + TAG_ACTUAL_VALUE + ">");
      if (inst.classIsMissing())
	append("?");
      else
	append(Utils.doubleToString(inst.classValue(), prec));
      append("</" + TAG_ACTUAL_VALUE + ">\n");
      // predicted
      append("    <" + TAG_PREDICTED_VALUE + ">");
      if (inst.classIsMissing())
	append("?");
      else
	append(Utils.doubleToString(predValue, prec));
      append("</" + TAG_PREDICTED_VALUE + ">\n");
      // error
      append("    <" + TAG_ERROR + ">");
      if (Utils.isMissingValue(predValue) || inst.classIsMissing())
	append("?");
      else
	append(Utils.doubleToString(predValue - inst.classValue(), prec));
      append("</" + TAG_ERROR + ">\n");
    } else {
      // actual
      append("    <" + TAG_ACTUAL_LABEL + " " + ATT_INDEX + "=\"" + ((int) inst.classValue()+1) + "\"" + ">");
      append(sanitize(inst.toString(inst.classIndex())));
      append("</" + TAG_ACTUAL_LABEL + ">\n");
      // predicted
      append("    <" + TAG_PREDICTED_LABEL + " " + ATT_INDEX + "=\"" + ((int) predValue+1) + "\"" + ">");
      if (Utils.isMissingValue(predValue))
	append("?");
      else
	append(sanitize(inst.dataset().classAttribute().value((int)predValue)));
      append("</" + TAG_PREDICTED_LABEL + ">\n");
      // error?
      append("    <" + TAG_ERROR + ">");
      if (!Utils.isMissingValue(predValue) && !inst.classIsMissing() && ((int) predValue+1 != (int) inst.classValue()+1))
	append(VAL_YES);
      else
	append(VAL_NO);
      append("</" + TAG_ERROR + ">\n");
      // prediction/distribution
      if (m_OutputDistribution) {
	append("    <" + TAG_DISTRIBUTION + ">\n");
	double[] dist = classifier.distributionForInstance(withMissing);
	for (int n = 0; n < dist.length; n++) {
	  append("      <" + TAG_CLASS_LABEL + " " + ATT_INDEX + "=\"" + (n+1) + "\"");
	  if (!Utils.isMissingValue(predValue) && (n == (int) predValue))
	    append(" " + ATT_PREDICTED + "=\"" + VAL_YES + "\"");
	  append(">");
	  append(Utils.doubleToString(dist[n], prec));
	  append("</" + TAG_CLASS_LABEL + ">\n");
	}
	append("    </" + TAG_DISTRIBUTION + ">\n");
      }
      else {
	append("    <" + TAG_PREDICTION + ">");
	if (Utils.isMissingValue(predValue))
	  append("?");
	else
	  append(Utils.doubleToString(classifier.distributionForInstance(withMissing) [(int)predValue], prec));
	append("</" + TAG_PREDICTION + ">\n");
      }
    }

    // attributes
    if (m_Attributes != null)
      append(attributeValuesString(withMissing));
    
    // closing tag
    append("  </" + TAG_PREDICTION + ">\n");
  }
  
  /**
   * Does nothing.
   */
  protected void doPrintFooter() {
    append("</" + TAG_PREDICTIONS + ">\n");
  }
}

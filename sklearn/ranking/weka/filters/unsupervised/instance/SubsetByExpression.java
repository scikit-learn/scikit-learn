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
 * SubsetByExpression.java
 * Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.instance;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.SimpleBatchFilter;
import weka.filters.unsupervised.instance.subsetbyexpression.Parser;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Filters instances according to a user-specified expression.<br/>
 * <br/>
 * Grammar:<br/>
 * <br/>
 * boolexpr_list ::= boolexpr_list boolexpr_part | boolexpr_part;<br/>
 * <br/>
 * boolexpr_part ::= boolexpr:e {: parser.setResult(e); :} ;<br/>
 * <br/>
 * boolexpr ::=    BOOLEAN <br/>
 *               | true<br/>
 *               | false<br/>
 *               | expr &lt; expr<br/>
 *               | expr &lt;= expr<br/>
 *               | expr &gt; expr<br/>
 *               | expr &gt;= expr<br/>
 *               | expr = expr<br/>
 *               | ( boolexpr )<br/>
 *               | not boolexpr<br/>
 *               | boolexpr and boolexpr<br/>
 *               | boolexpr or boolexpr<br/>
 *               | ATTRIBUTE is STRING<br/>
 *               ;<br/>
 * <br/>
 * expr      ::=   NUMBER<br/>
 *               | ATTRIBUTE<br/>
 *               | ( expr )<br/>
 *               | opexpr<br/>
 *               | funcexpr<br/>
 *               ;<br/>
 * <br/>
 * opexpr    ::=   expr + expr<br/>
 *               | expr - expr<br/>
 *               | expr * expr<br/>
 *               | expr / expr<br/>
 *               ;<br/>
 * <br/>
 * funcexpr ::=    abs ( expr )<br/>
 *               | sqrt ( expr )<br/>
 *               | log ( expr )<br/>
 *               | exp ( expr )<br/>
 *               | sin ( expr )<br/>
 *               | cos ( expr )<br/>
 *               | tan ( expr )<br/>
 *               | rint ( expr )<br/>
 *               | floor ( expr )<br/>
 *               | pow ( expr for base , expr for exponent )<br/>
 *               | ceil ( expr )<br/>
 *               ;<br/>
 * <br/>
 * Notes:<br/>
 * - NUMBER<br/>
 *   any integer or floating point number <br/>
 *   (but not in scientific notation!)<br/>
 * - STRING<br/>
 *   any string surrounded by single quotes; <br/>
 *   the string may not contain a single quote though.<br/>
 * - ATTRIBUTE<br/>
 *   the following placeholders are recognized for <br/>
 *   attribute values:<br/>
 *   - CLASS for the class value in case a class attribute is set.<br/>
 *   - ATTxyz with xyz a number from 1 to # of attributes in the<br/>
 *     dataset, representing the value of indexed attribute.<br/>
 * <br/>
 * Examples:<br/>
 * - extracting only mammals and birds from the 'zoo' UCI dataset:<br/>
 *   (CLASS is 'mammal') or (CLASS is 'bird')<br/>
 * - extracting only animals with at least 2 legs from the 'zoo' UCI dataset:<br/>
 *   (ATT14 &gt;= 2)<br/>
 * - extracting only instances with non-missing 'wage-increase-second-year'<br/>
 *   from the 'labor' UCI dataset:<br/>
 *   not ismissing(ATT3)<br/>
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -E &lt;expr&gt;
 *  The expression to use for filtering
 *  (default: true).</pre>
 * 
 <!-- options-end -->
 *
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6113 $
 */
public class SubsetByExpression
  extends SimpleBatchFilter {

  /** for serialization. */
  private static final long serialVersionUID = 5628686110979589602L;
  
  /** the expresion to use for filtering. */
  protected String m_Expression = "true";
  
  /**
   * Returns a string describing this filter.
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Filters instances according to a user-specified expression.\n\n"
      + "Grammar:\n\n"
      + "boolexpr_list ::= boolexpr_list boolexpr_part | boolexpr_part;\n"
      + "\n"
      + "boolexpr_part ::= boolexpr:e {: parser.setResult(e); :} ;\n"
      + "\n"
      + "boolexpr ::=    BOOLEAN \n"
      + "              | true\n"
      + "              | false\n"
      + "              | expr < expr\n"
      + "              | expr <= expr\n"
      + "              | expr > expr\n"
      + "              | expr >= expr\n"
      + "              | expr = expr\n"
      + "              | ( boolexpr )\n"
      + "              | not boolexpr\n"
      + "              | boolexpr and boolexpr\n"
      + "              | boolexpr or boolexpr\n"
      + "              | ATTRIBUTE is STRING\n"
      + "              ;\n"
      + "\n"
      + "expr      ::=   NUMBER\n"
      + "              | ATTRIBUTE\n"
      + "              | ( expr )\n"
      + "              | opexpr\n"
      + "              | funcexpr\n"
      + "              ;\n"
      + "\n"
      + "opexpr    ::=   expr + expr\n"
      + "              | expr - expr\n"
      + "              | expr * expr\n"
      + "              | expr / expr\n"
      + "              ;\n"
      + "\n"
      + "funcexpr ::=    abs ( expr )\n"
      + "              | sqrt ( expr )\n"
      + "              | log ( expr )\n"
      + "              | exp ( expr )\n"
      + "              | sin ( expr )\n"
      + "              | cos ( expr )\n"
      + "              | tan ( expr )\n"
      + "              | rint ( expr )\n"
      + "              | floor ( expr )\n"
      + "              | pow ( expr for base , expr for exponent )\n"
      + "              | ceil ( expr )\n"
      + "              ;\n"
      + "\n"
      + "Notes:\n"
      + "- NUMBER\n"
      + "  any integer or floating point number \n"
      + "  (but not in scientific notation!)\n"
      + "- STRING\n"
      + "  any string surrounded by single quotes; \n"
      + "  the string may not contain a single quote though.\n"
      + "- ATTRIBUTE\n"
      + "  the following placeholders are recognized for \n"
      + "  attribute values:\n"
      + "  - CLASS for the class value in case a class attribute is set.\n"
      + "  - ATTxyz with xyz a number from 1 to # of attributes in the\n"
      + "    dataset, representing the value of indexed attribute.\n"
      + "\n"
      + "Examples:\n"
      + "- extracting only mammals and birds from the 'zoo' UCI dataset:\n"
      + "  (CLASS is 'mammal') or (CLASS is 'bird')\n"
      + "- extracting only animals with at least 2 legs from the 'zoo' UCI dataset:\n"
      + "  (ATT14 >= 2)\n"
      + "- extracting only instances with non-missing 'wage-increase-second-year'\n"
      + "  from the 'labor' UCI dataset:\n"
      + "  not ismissing(ATT3)\n"
      ;
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector	result;
    
    result = new Vector();

    result.addElement(new Option(
	"\tThe expression to use for filtering\n"
	+ "\t(default: true).",
	"E", 1, "-E <expr>"));

    return result.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -E &lt;expr&gt;
   *  The expression to use for filtering
   *  (default: true).</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('E', options);
    if (tmpStr.length() != 0)
      setExpression(tmpStr);
    else
      setExpression("true");
   
    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>	result;
    
    result = new Vector();

    result.add("-E");
    result.add("" + getExpression());

    return result.toArray(new String[result.size()]);
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
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
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
    return result;
  }

  /**
   * Sets the expression used for filtering.
   *
   * @param value	the expression
   */
  public void setExpression(String value) {
    m_Expression = value;
  }

  /**
   * Returns the expression used for filtering.
   *
   * @return 		the expression
   */
  public String getExpression() {
    return m_Expression;
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String expressionTipText() {
    return "The expression to used for filtering the dataset.";
  }

  /**
   * Determines the output format based on the input format and returns 
   * this.
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the determination goes wrong
   */
  protected Instances determineOutputFormat(Instances inputFormat)
      throws Exception {
    
    return new Instances(inputFormat, 0);
  }

  /**
   * Processes the given data (may change the provided dataset) and returns
   * the modified version. This method is called in batchFinished().
   *
   * @param instances   the data to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   * @see               #batchFinished()
   */
  protected Instances process(Instances instances) throws Exception {
    if (!isFirstBatchDone())
      return Parser.filter(m_Expression, instances);
    else
      return instances;
  }

  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6113 $");
  }

  /**
   * Main method for running this filter.
   *
   * @param args 	arguments for the filter: use -h for help
   */
  public static void main(String[] args) {
    runFilter(new SubsetByExpression(), args);
  }
}

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
 *    AddExpression.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.AttributeExpression;
import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Filter;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * An instance filter that creates a new attribute by applying a mathematical expression to existing attributes. The expression can contain attribute references and numeric constants. Supported operators are :<br/>
 * +, -, *, /, ^, log, abs, cos, exp, sqrt, floor, ceil, rint, tan, sin, (, )<br/>
 * Attributes are specified by prefixing with 'a', eg. a7 is attribute number 7 (starting from 1).<br/>
 * Example expression : a1^2*a5/log(a7*4.0).
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -E &lt;expression&gt;
 *  Specify the expression to apply. Eg a1^2*a5/log(a7*4.0).
 *  Supported opperators: ,+, -, *, /, ^, log, abs, cos, 
 *  exp, sqrt, floor, ceil, rint, tan, sin, (, )
 *  (default: a1^2)</pre>
 * 
 * <pre> -N &lt;name&gt;
 *  Specify the name for the new attribute. (default is the expression provided with -E)</pre>
 * 
 * <pre> -D
 *  Debug. Names attribute with the postfix parse of the expression.</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class AddExpression 
  extends Filter 
  implements UnsupervisedFilter, StreamableFilter, OptionHandler {

  /** for serialization */
  static final long serialVersionUID = 402130384261736245L;
  
  /** The infix expression */
  private String m_infixExpression = "a1^2";

  /** Name of the new attribute. "expression"  length string will use the 
      provided expression as the new attribute name */
  private String m_attributeName="expression";

  /** If true, makes the attribute name equal to the postfix parse of the
      expression */
  private boolean m_Debug = false;

  private AttributeExpression m_attributeExpression = null;

  /**
   * Returns a string describing this filter
   *
   * @return 		a description of the filter suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "An instance filter that creates a new attribute by applying a "
      + "mathematical expression to existing attributes. The expression "
      + "can contain attribute references and numeric constants. Supported "
      + "operators are :\n"
      + "+, -, *, /, ^, log, abs, cos, exp, sqrt, floor, ceil, rint, tan, "
      + "sin, (, )\n"
      + "Attributes are specified by prefixing with 'a', eg. a7 is "
      + "attribute number 7 (starting from 1).\n"
      + "Example expression : a1^2*a5/log(a7*4.0).";
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(3); 

    newVector.addElement(new Option(
	     "\tSpecify the expression to apply. Eg a1^2*a5/log(a7*4.0)."
	     +"\n\tSupported opperators: ,+, -, *, /, ^, log, abs, cos, "
	     +"\n\texp, sqrt, floor, ceil, rint, tan, sin, (, )"
	     +"\n\t(default: a1^2)",
	     "E",1,"-E <expression>"));

    newVector.addElement(new Option(
	     "\tSpecify the name for the new attribute. (default is the "
	     +"expression provided with -E)",
	     "N",1,"-N <name>"));

    newVector.addElement(new Option(
	     "\tDebug. Names attribute with the postfix parse of the "
	     +"expression.","D",0,"-D"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -E &lt;expression&gt;
   *  Specify the expression to apply. Eg a1^2*a5/log(a7*4.0).
   *  Supported opperators: ,+, -, *, /, ^, log, abs, cos, 
   *  exp, sqrt, floor, ceil, rint, tan, sin, (, )
   *  (default: a1^2)</pre>
   * 
   * <pre> -N &lt;name&gt;
   *  Specify the name for the new attribute. (default is the expression provided with -E)</pre>
   * 
   * <pre> -D
   *  Debug. Names attribute with the postfix parse of the expression.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String expString = Utils.getOption('E', options);
    if (expString.length() != 0) {
      setExpression(expString);
    } else {
      setExpression("a1^2");
    }

    String name = Utils.getOption('N',options);
    if (name.length() != 0) {
      setName(name);
    }

    setDebug(Utils.getFlag('D', options));
  }
  
  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [5];
    int current = 0;
    
    options[current++] = "-E"; options[current++] = getExpression();
    options[current++] = "-N"; options[current++] = getName();

    if (getDebug()) {
      options[current++] = "-D";
    }
    
    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String nameTipText() {
    return "Set the name of the new attribute.";
  }

  /**
   * Set the name for the new attribute. The string "expression" can
   * be used to make the name of the new attribute equal to the expression
   * provided.
   * @param name the name of the new attribute
   */
  public void setName(String name) {
    m_attributeName = name;
  }

  /**
   * Returns the name of the new attribute
   * @return the name of the new attribute
   */
  public String getName() {
    return m_attributeName;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String debugTipText() {
    return "Set debug mode. If true then the new attribute will be named with "
      +"the postfix parse of the supplied expression.";
  }
  
  /**
   * Set debug mode. Causes the new attribute to be named with the postfix
   * parse of the expression
   * @param d true if debug mode is to be used
   */
  public void setDebug(boolean d) {
    m_Debug = d;
  }

  /**
   * Gets whether debug is set
   * @return true if debug is set
   */
  public boolean getDebug() {
    return m_Debug;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String expressionTipText() {
    return "Set the math expression to apply. Eg. a1^2*a5/log(a7*4.0)";
  }

  /**
   * Set the expression to apply
   * @param expr a mathematical expression to apply
   */
  public void setExpression(String expr) {
    m_infixExpression = expr;
  }

  /**
   * Get the expression
   * @return the expression
   */
  public String getExpression() {
    return m_infixExpression;
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
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the format couldn't be set successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    m_attributeExpression = new AttributeExpression();
    m_attributeExpression.
      convertInfixToPostfix(new String(m_infixExpression));

    super.setInputFormat(instanceInfo);

    Instances outputFormat = new Instances(instanceInfo, 0);
    Attribute newAttribute;
    if (m_Debug) {
      newAttribute = 
        new Attribute(m_attributeExpression.getPostFixExpression());
    } else if (m_attributeName.compareTo("expression") != 0) {
      newAttribute = new Attribute(m_attributeName);
    } else {
      newAttribute = new Attribute(m_infixExpression);
    }
    outputFormat.insertAttributeAt(newAttribute, 
				   instanceInfo.numAttributes());
    setOutputFormat(outputFormat);
    return true;
  }

  /**
   * Input an instance for filtering. Ordinarily the instance is processed
   * and made available for output immediately. Some filters require all
   * instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been defined.
   * @throws Exception if there was a problem during the filtering.
   */
  public boolean input(Instance instance) throws Exception {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }

    double[] vals = new double[instance.numAttributes()+1];
    for(int i = 0; i < instance.numAttributes(); i++) {
      if (instance.isMissing(i)) {
	vals[i] = Utils.missingValue();
      } else {
	vals[i] = instance.value(i);
      }
    }

    m_attributeExpression.evaluateExpression(vals);

    Instance inst = null;
    if (instance instanceof SparseInstance) {
      inst = new SparseInstance(instance.weight(), vals);
    } 
    //RANKING BEGIN
    if(instance instanceof PreferenceDenseInstance){
    	PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	inst = new PreferenceDenseInstance(instance.weight(), vals, pdi.getHashMap());
    }
    //RANKING END
    else {
      inst = new DenseInstance(instance.weight(), vals);
    }

    inst.setDataset(getOutputFormat());
    copyValues(inst, false, instance.dataset(), getOutputFormat());
    inst.setDataset(getOutputFormat());
    push(inst);
    return true;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }
  
  /**
   * Main method for testing this class.
   *
   * @param args should contain arguments to the filter: use -h for help
   */
  public static void main(String [] args) {
    runFilter(new AddExpression(), args);
  }
}

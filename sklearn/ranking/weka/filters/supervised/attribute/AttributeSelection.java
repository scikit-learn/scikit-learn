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
 *    AttributeSelection.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.supervised.attribute;

import weka.attributeSelection.ASEvaluation;
import weka.attributeSelection.ASSearch;
import weka.attributeSelection.AttributeEvaluator;
import weka.attributeSelection.AttributeTransformer;
import weka.attributeSelection.BestFirst;
import weka.attributeSelection.CfsSubsetEval;
import weka.attributeSelection.Ranker;
import weka.attributeSelection.UnsupervisedAttributeEvaluator;
import weka.attributeSelection.UnsupervisedSubsetEvaluator;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * A supervised attribute filter that can be used to select attributes. It is very flexible and allows various search and evaluation methods to be combined.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;"Name of search class [search options]"&gt;
 *  Sets search method for subset evaluators.
 *  eg. -S "weka.attributeSelection.BestFirst -S 8"</pre>
 * 
 * <pre> -E &lt;"Name of attribute/subset evaluation class [evaluator options]"&gt;
 *  Sets attribute/subset evaluator.
 *  eg. -E "weka.attributeSelection.CfsSubsetEval -L"</pre>
 * 
 * <pre> 
 * Options specific to evaluator weka.attributeSelection.CfsSubsetEval:
 * </pre>
 * 
 * <pre> -M
 *  Treat missing values as a seperate value.</pre>
 * 
 * <pre> -L
 *  Don't include locally predictive attributes.</pre>
 * 
 * <pre> 
 * Options specific to search weka.attributeSelection.BestFirst:
 * </pre>
 * 
 * <pre> -P &lt;start set&gt;
 *  Specify a starting set of attributes.
 *  Eg. 1,3,5-7.</pre>
 * 
 * <pre> -D &lt;0 = backward | 1 = forward | 2 = bi-directional&gt;
 *  Direction of search. (default = 1).</pre>
 * 
 * <pre> -N &lt;num&gt;
 *  Number of non-improving nodes to
 *  consider before terminating search.</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Size of lookup cache for evaluated subsets.
 *  Expressed as a multiple of the number of
 *  attributes in the data set. (default = 1)</pre>
 * 
 <!-- options-end -->
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class AttributeSelection 
  extends Filter
  implements SupervisedFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -296211247688169716L;

  /** the attribute selection evaluation object */
  private weka.attributeSelection.AttributeSelection m_trainSelector;

  /** the attribute evaluator to use */
  private ASEvaluation m_ASEvaluator;

  /** the search method if any */
  private ASSearch m_ASSearch;

  /** holds a copy of the full set of valid  options passed to the filter */
  private String [] m_FilterOptions;

  /** holds the selected attributes  */
  private int [] m_SelectedAttributes;

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "A supervised attribute filter that can be used to select " 
      + "attributes. It is very flexible and allows various search " 
      + "and evaluation methods to be combined.";
  }

  /**
   * Constructor
   */
  public AttributeSelection () {
    
    resetOptions();
  }

  /**
   * Returns an enumeration describing the available options.
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    
    Vector newVector = new Vector(6);

    newVector.addElement(new Option(
	"\tSets search method for subset evaluators.\n"
	+ "\teg. -S \"weka.attributeSelection.BestFirst -S 8\"", 
	"S", 1,
	"-S <\"Name of search class [search options]\">"));

    newVector.addElement(new Option(
	"\tSets attribute/subset evaluator.\n"
	+ "\teg. -E \"weka.attributeSelection.CfsSubsetEval -L\"",
	"E", 1,
	"-E <\"Name of attribute/subset evaluation class [evaluator options]\">"));
    
    if ((m_ASEvaluator != null) && (m_ASEvaluator instanceof OptionHandler)) {
      Enumeration enu = ((OptionHandler)m_ASEvaluator).listOptions();
      
      newVector.addElement(new Option("", "", 0, "\nOptions specific to "
	   + "evaluator " + m_ASEvaluator.getClass().getName() + ":"));
      while (enu.hasMoreElements()) {
	newVector.addElement((Option)enu.nextElement());
      }
    }
  
    if ((m_ASSearch != null) && (m_ASSearch instanceof OptionHandler)) {
      Enumeration enu = ((OptionHandler)m_ASSearch).listOptions();
    
      newVector.addElement(new Option("", "", 0, "\nOptions specific to "
	      + "search " + m_ASSearch.getClass().getName() + ":"));
      while (enu.hasMoreElements()) {
	newVector.addElement((Option)enu.nextElement());
      }
    }
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;"Name of search class [search options]"&gt;
   *  Sets search method for subset evaluators.
   *  eg. -S "weka.attributeSelection.BestFirst -S 8"</pre>
   * 
   * <pre> -E &lt;"Name of attribute/subset evaluation class [evaluator options]"&gt;
   *  Sets attribute/subset evaluator.
   *  eg. -E "weka.attributeSelection.CfsSubsetEval -L"</pre>
   * 
   * <pre> 
   * Options specific to evaluator weka.attributeSelection.CfsSubsetEval:
   * </pre>
   * 
   * <pre> -M
   *  Treat missing values as a seperate value.</pre>
   * 
   * <pre> -L
   *  Don't include locally predictive attributes.</pre>
   * 
   * <pre> 
   * Options specific to search weka.attributeSelection.BestFirst:
   * </pre>
   * 
   * <pre> -P &lt;start set&gt;
   *  Specify a starting set of attributes.
   *  Eg. 1,3,5-7.</pre>
   * 
   * <pre> -D &lt;0 = backward | 1 = forward | 2 = bi-directional&gt;
   *  Direction of search. (default = 1).</pre>
   * 
   * <pre> -N &lt;num&gt;
   *  Number of non-improving nodes to
   *  consider before terminating search.</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Size of lookup cache for evaluated subsets.
   *  Expressed as a multiple of the number of
   *  attributes in the data set. (default = 1)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String optionString;
    resetOptions();

    if (Utils.getFlag('X',options)) {
	throw new Exception("Cross validation is not a valid option"
			    + " when using attribute selection as a Filter.");
    }

    optionString = Utils.getOption('E',options);
    if (optionString.length() != 0) {
      optionString = optionString.trim();
      // split a quoted evaluator name from its options (if any)
      int breakLoc = optionString.indexOf(' ');
      String evalClassName = optionString;
      String evalOptionsString = "";
      String [] evalOptions=null;
      if (breakLoc != -1) {
	evalClassName = optionString.substring(0, breakLoc);
	evalOptionsString = optionString.substring(breakLoc).trim();
	evalOptions = Utils.splitOptions(evalOptionsString);
      }
      setEvaluator(ASEvaluation.forName(evalClassName, evalOptions));
    }

    if (m_ASEvaluator instanceof AttributeEvaluator) {
      setSearch(new Ranker());
    }

    optionString = Utils.getOption('S',options);
    if (optionString.length() != 0) {
      optionString = optionString.trim();
      int breakLoc = optionString.indexOf(' ');
      String SearchClassName = optionString;
      String SearchOptionsString = "";
      String [] SearchOptions=null;
      if (breakLoc != -1) {
	SearchClassName = optionString.substring(0, breakLoc);
	SearchOptionsString = optionString.substring(breakLoc).trim();
	SearchOptions = Utils.splitOptions(SearchOptionsString);
      }
      setSearch(ASSearch.forName(SearchClassName, SearchOptions));
    }

    Utils.checkForRemainingOptions(options);
  }


  /**
   * Gets the current settings for the attribute selection (search, evaluator)
   * etc.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String [] getOptions() {
    String [] EvaluatorOptions = new String[0];
    String [] SearchOptions = new String[0];
    int current = 0;

    if (m_ASEvaluator instanceof OptionHandler) {
      EvaluatorOptions = ((OptionHandler)m_ASEvaluator).getOptions();
    }

    if (m_ASSearch instanceof OptionHandler) {
      SearchOptions = ((OptionHandler)m_ASSearch).getOptions();
    }

    String [] setOptions = new String [10];
    setOptions[current++]="-E";
    setOptions[current++]= getEvaluator().getClass().getName()
      +" "+Utils.joinOptions(EvaluatorOptions);

    setOptions[current++]="-S";
    setOptions[current++]=getSearch().getClass().getName() 
      + " "+Utils.joinOptions(SearchOptions);

    while (current < setOptions.length) {
      setOptions[current++] = "";
    }
    
    return setOptions;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String evaluatorTipText() {

    return "Determines how attributes/attribute subsets are evaluated.";
  }

  /**
   * set attribute/subset evaluator
   * 
   * @param evaluator the evaluator to use
   */
  public void setEvaluator(ASEvaluation evaluator) {
    m_ASEvaluator = evaluator;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String searchTipText() {

    return "Determines the search method.";
  }

  /**
   * Set search class
   * 
   * @param search the search class to use
   */
  public void setSearch(ASSearch search) {
    m_ASSearch = search;
  }

  /**
   * Get the name of the attribute/subset evaluator
   *
   * @return the name of the attribute/subset evaluator as a string
   */
  public ASEvaluation getEvaluator() {
    
      return m_ASEvaluator;
  }

  /**
   * Get the name of the search method
   *
   * @return the name of the search method as a string
   */
  public ASSearch getSearch() {
    
      return m_ASSearch;
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities	result;
    
    if (m_ASEvaluator == null) {
      result = super.getCapabilities();
      result.disableAll();
    } else {
      result = m_ASEvaluator.getCapabilities();
      // class index will be set if necessary, so we always allow the dataset
      // to have no class attribute set. see the following method:
      //   weka.attributeSelection.AttributeSelection.SelectAttributes(Instances)
      result.enable(Capability.NO_CLASS);
    }
    
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    result.disable(Capability.RANKING);
    
    result.setMinimumNumberInstances(0);
    
    return result;
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
   * @throws Exception if the input instance was not of the correct format 
   * or if there was a problem with the filtering.
   */
  public boolean input(Instance instance) throws Exception {
    
    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }

    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }

    if (isOutputFormatDefined()) {
      convertInstance(instance);
      return true;
    }

    bufferInput(instance);
    return false;
  }

  /**
   * Signify that this batch of input to the filter is finished. If the filter
   * requires all instances prior to filtering, output() may now be called
   * to retrieve the filtered instances.
   *
   * @return true if there are instances pending output.
   * @throws IllegalStateException if no input structure has been defined.
   * @throws Exception if there is a problem during the attribute selection.
   */
  public boolean batchFinished() throws Exception {
    
    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }

    if (!isOutputFormatDefined()) {
      m_trainSelector.setEvaluator(m_ASEvaluator);
      m_trainSelector.setSearch(m_ASSearch);
      m_trainSelector.SelectAttributes(getInputFormat());
      //      System.out.println(m_trainSelector.toResultsString());

      m_SelectedAttributes = m_trainSelector.selectedAttributes();
      if (m_SelectedAttributes == null) {
	throw new Exception("No selected attributes\n");
      }
     
      setOutputFormat();
      
      // Convert pending input instances
      for (int i = 0; i < getInputFormat().numInstances(); i++) {
	convertInstance(getInputFormat().instance(i));
      }
      flushInput();
    }
    
    m_NewBatch = true;
    return (numPendingOutput() != 0);
  }

  /**
   * Set the output format. Takes the currently defined attribute set 
   * m_InputFormat and calls setOutputFormat(Instances) appropriately.
   * 
   * @throws Exception if something goes wrong
   */
  protected void setOutputFormat() throws Exception {
    Instances informat;

    if (m_SelectedAttributes == null) {
      setOutputFormat(null);
      return;
    }

    FastVector attributes = new FastVector(m_SelectedAttributes.length);

    int i;
    if (m_ASEvaluator instanceof AttributeTransformer) {
      informat = ((AttributeTransformer)m_ASEvaluator).transformedHeader();
    } else {
      informat = getInputFormat();
    }

    for (i=0;i < m_SelectedAttributes.length;i++) {
      attributes.
	addElement(informat.attribute(m_SelectedAttributes[i]).copy());
    }

    Instances outputFormat = 
      new Instances(getInputFormat().relationName(), attributes, 0);


    if (!(m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) &&
	!(m_ASEvaluator instanceof UnsupervisedAttributeEvaluator)) {
      outputFormat.setClassIndex(m_SelectedAttributes.length - 1);
    }
    
    setOutputFormat(outputFormat);  
  }

  /**
   * Convert a single instance over. Selected attributes only are transfered.
   * The converted instance is added to the end of
   * the output queue.
   *
   * @param instance the instance to convert
   * @throws Exception if something goes wrong
   */
  protected void convertInstance(Instance instance) throws Exception {
    double[] newVals = new double[getOutputFormat().numAttributes()];

    if (m_ASEvaluator instanceof AttributeTransformer) {
      Instance tempInstance = ((AttributeTransformer)m_ASEvaluator).
	convertInstance(instance);
      for (int i = 0; i < m_SelectedAttributes.length; i++) {
	int current = m_SelectedAttributes[i];
	newVals[i] = tempInstance.value(current);
      }
    } else {
      for (int i = 0; i < m_SelectedAttributes.length; i++) {
	int current = m_SelectedAttributes[i];
	newVals[i] = instance.value(current);
      }
    }
    if (instance instanceof SparseInstance) {
      push(new SparseInstance(instance.weight(), newVals));
    } else {
      push(new DenseInstance(instance.weight(), newVals));
    }
  }

  /**
   * set options to their default values
   */
  protected void resetOptions() {

    m_trainSelector = new weka.attributeSelection.AttributeSelection();
    setEvaluator(new CfsSubsetEval());
    setSearch(new BestFirst());
    m_SelectedAttributes = null;
    m_FilterOptions = null;
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
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new AttributeSelection(), argv);
  }
}


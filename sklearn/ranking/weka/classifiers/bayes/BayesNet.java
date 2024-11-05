/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * BayesNet.java
 * Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 * 
 */
package weka.classifiers.bayes;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.bayes.net.ADNode;
import weka.classifiers.bayes.net.BIFReader;
import weka.classifiers.bayes.net.ParentSet;
import weka.classifiers.bayes.net.estimate.BayesNetEstimator;
import weka.classifiers.bayes.net.estimate.DiscreteEstimatorBayes;
import weka.classifiers.bayes.net.estimate.SimpleEstimator;
import weka.classifiers.bayes.net.search.SearchAlgorithm;
import weka.classifiers.bayes.net.search.local.K2;
import weka.classifiers.bayes.net.search.local.LocalScoreSearchAlgorithm;
import weka.classifiers.bayes.net.search.local.Scoreable;
import weka.core.AdditionalMeasureProducer;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Drawable;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceAttribute;
import weka.estimators.Estimator;
import weka.filters.Filter;
import weka.filters.supervised.attribute.Discretize;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Bayes Network learning using various search algorithms and quality measures.<br/>
 * Base class for a Bayes Network classifier. Provides datastructures (network structure, conditional probability distributions, etc.) and facilities common to Bayes Network learning algorithms like K2 and B.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * http://sourceforge.net/projects/weka/files/documentation/WekaManual-3-7-0.pdf/download
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Do not use ADTree data structure
 * </pre>
 * 
 * <pre> -B &lt;BIF file&gt;
 *  BIF file to compare with
 * </pre>
 * 
 * <pre> -Q weka.classifiers.bayes.net.search.SearchAlgorithm
 *  Search algorithm
 * </pre>
 * 
 * <pre> -E weka.classifiers.bayes.net.estimate.SimpleEstimator
 *  Estimator algorithm
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 5928 $
 */
public class BayesNet
  extends AbstractClassifier
  implements OptionHandler, WeightedInstancesHandler, Drawable, 
             AdditionalMeasureProducer {

  /** for serialization */
  static final long serialVersionUID = 746037443258775954L;


  /**
   * The parent sets.
   */
  protected ParentSet[] m_ParentSets;

  /**
   * The attribute estimators containing CPTs.
   */
  public Estimator[][] m_Distributions;


  /** filter used to quantize continuous variables, if any **/
  protected Discretize m_DiscretizeFilter = null;

  /** attribute index of a non-nominal attribute */
  int m_nNonDiscreteAttribute = -1;

  /** filter used to fill in missing values, if any **/
  protected ReplaceMissingValues m_MissingValuesFilter = null;	

  /**
   * The number of classes
   */
  protected int m_NumClasses;

  /**
   * The dataset header for the purposes of printing out a semi-intelligible
   * model
   */
  public Instances m_Instances;

  /**
   * Datastructure containing ADTree representation of the database.
   * This may result in more efficient access to the data.
   */
  ADNode m_ADTree;

  /**
   * Bayes network to compare the structure with.
   */
  protected BIFReader m_otherBayesNet = null;

  /**
   * Use the experimental ADTree datastructure for calculating contingency tables
   */
  boolean m_bUseADTree = false;

  /**
   * Search algorithm used for learning the structure of a network.
   */
  SearchAlgorithm m_SearchAlgorithm = new K2();

  /**
   * Search algorithm used for learning the structure of a network.
   */
  BayesNetEstimator m_BayesNetEstimator = new SimpleEstimator();

  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);

    // instances
    result.setMinimumNumberInstances(0);

    return result;
  }

  /**
   * Generates the classifier.
   * 
   * @param instances set of instances serving as training data
   * @throws Exception if the classifier has not been generated
   * successfully
   */
  public void buildClassifier(Instances instances) throws Exception {

    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();

    // ensure we have a data set with discrete variables only and with no missing values
    instances = normalizeDataSet(instances);

    // Copy the instances
    m_Instances = new Instances(instances);

    // sanity check: need more than 1 variable in datat set
    m_NumClasses = instances.numClasses();

    // initialize ADTree
    if (m_bUseADTree) {
      m_ADTree = ADNode.makeADTree(instances);
      //      System.out.println("Oef, done!");
    }

    // build the network structure
    initStructure();

    // build the network structure
    buildStructure();

    // build the set of CPTs
    estimateCPTs();

    // Save space
    // m_Instances = new Instances(m_Instances, 0);
    m_ADTree = null;
  } // buildClassifier

  /** ensure that all variables are nominal and that there are no missing values
   * @param instances data set to check and quantize and/or fill in missing values
   * @return filtered instances
   * @throws Exception if a filter (Discretize, ReplaceMissingValues) fails
   */
  protected Instances normalizeDataSet(Instances instances) throws Exception {
    m_DiscretizeFilter = null;
    m_MissingValuesFilter = null;

    boolean bHasNonNominal = false;
    boolean bHasMissingValues = false;

    Enumeration enu = instances.enumerateAttributes();		
    while (enu.hasMoreElements()) {
      Attribute attribute = (Attribute) enu.nextElement();
      if (attribute.type() != Attribute.NOMINAL && attribute.type() != PreferenceAttribute.RANKING) {
	m_nNonDiscreteAttribute = attribute.index();
	bHasNonNominal = true;
	//throw new UnsupportedAttributeTypeException("BayesNet handles nominal variables only. Non-nominal variable in dataset detected.");
      }
      Enumeration enum2 = instances.enumerateInstances();
      while (enum2.hasMoreElements()) {
	if (((Instance) enum2.nextElement()).isMissing(attribute)) {
	  bHasMissingValues = true;
	  // throw new NoSupportForMissingValuesException("BayesNet: no missing values, please.");
	}
      }
    }

    if (bHasNonNominal) {
      System.err.println("Warning: discretizing data set");
      m_DiscretizeFilter = new Discretize();
      m_DiscretizeFilter.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_DiscretizeFilter);
    }

    if (bHasMissingValues) {
      System.err.println("Warning: filling in missing values in data set");
      m_MissingValuesFilter = new ReplaceMissingValues();
      m_MissingValuesFilter.setInputFormat(instances);
      instances = Filter.useFilter(instances, m_MissingValuesFilter);
    }
    return instances;
  } // normalizeDataSet

  /** ensure that all variables are nominal and that there are no missing values
   * @param instance instance to check and quantize and/or fill in missing values
   * @return filtered instance
   * @throws Exception if a filter (Discretize, ReplaceMissingValues) fails
   */
  protected Instance normalizeInstance(Instance instance) throws Exception {
    if ((m_DiscretizeFilter != null) &&
	(instance.attribute(m_nNonDiscreteAttribute).type() != Attribute.NOMINAL && instance.attribute(m_nNonDiscreteAttribute).type() != PreferenceAttribute.RANKING)) {
      m_DiscretizeFilter.input(instance);
      instance = m_DiscretizeFilter.output();
    }
    if (m_MissingValuesFilter != null) {
      m_MissingValuesFilter.input(instance);
      instance = m_MissingValuesFilter.output();
    } else {
      // is there a missing value in this instance?
      // this can happen when there is no missing value in the training set
      for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
	if (iAttribute != instance.classIndex() && instance.isMissing(iAttribute)) {
	  System.err.println("Warning: Found missing value in test set, filling in values.");
	  m_MissingValuesFilter = new ReplaceMissingValues();
	  m_MissingValuesFilter.setInputFormat(m_Instances);
	  Filter.useFilter(m_Instances, m_MissingValuesFilter);
	  m_MissingValuesFilter.input(instance);
	  instance = m_MissingValuesFilter.output();
	  iAttribute = m_Instances.numAttributes();
	}
      }
    }
    return instance;
  } // normalizeInstance

  /**
   * Init structure initializes the structure to an empty graph or a Naive Bayes
   * graph (depending on the -N flag).
   * 
   * @throws Exception in case of an error
   */
  public void initStructure() throws Exception {

    // initialize topological ordering
    //    m_nOrder = new int[m_Instances.numAttributes()];
    //    m_nOrder[0] = m_Instances.classIndex();

    int nAttribute = 0;

    for (int iOrder = 1; iOrder < m_Instances.numAttributes(); iOrder++) {
      if (nAttribute == m_Instances.classIndex()) {
	nAttribute++;
      }

      //      m_nOrder[iOrder] = nAttribute++;
    }

    // reserve memory
    m_ParentSets = new ParentSet[m_Instances.numAttributes()];

    for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
      m_ParentSets[iAttribute] = new ParentSet(m_Instances.numAttributes());
    }
  } // initStructure

  /**
   * buildStructure determines the network structure/graph of the network.
   * The default behavior is creating a network where all nodes have the first
   * node as its parent (i.e., a BayesNet that behaves like a naive Bayes classifier).
   * This method can be overridden by derived classes to restrict the class
   * of network structures that are acceptable.
   * 
   * @throws Exception in case of an error
   */
  public void buildStructure() throws Exception {
    m_SearchAlgorithm.buildStructure(this, m_Instances);
  } // buildStructure

  /**
   * estimateCPTs estimates the conditional probability tables for the Bayes
   * Net using the network structure.
   * 
   * @throws Exception in case of an error
   */
  public void estimateCPTs() throws Exception {
    m_BayesNetEstimator.estimateCPTs(this);
  } // estimateCPTs

  /**
   * initializes the conditional probabilities
   * 
   * @throws Exception in case of an error
   */
  public void initCPTs() throws Exception {
    m_BayesNetEstimator.initCPTs(this);
  } // estimateCPTs

  /**
   * Updates the classifier with the given instance.
   * 
   * @param instance the new training instance to include in the model
   * @throws Exception if the instance could not be incorporated in
   * the model.
   */
  public void updateClassifier(Instance instance) throws Exception {
    instance = normalizeInstance(instance);
    m_BayesNetEstimator.updateClassifier(this, instance);
  } // updateClassifier

  /**
   * Calculates the class membership probabilities for the given test
   * instance.
   * 
   * @param instance the instance to be classified
   * @return predicted class probability distribution
   * @throws Exception if there is a problem generating the prediction
   */
  public double[] distributionForInstance(Instance instance) throws Exception {
    instance = normalizeInstance(instance);
    return m_BayesNetEstimator.distributionForInstance(this, instance);
  } // distributionForInstance

  /**
   * Calculates the counts for Dirichlet distribution for the 
   * class membership probabilities for the given test instance.
   * 
   * @param instance the instance to be classified
   * @return counts for Dirichlet distribution for class probability 
   * @throws Exception if there is a problem generating the prediction
   */
  public double[] countsForInstance(Instance instance) throws Exception {
    double[] fCounts = new double[m_NumClasses];

    for (int iClass = 0; iClass < m_NumClasses; iClass++) {
      fCounts[iClass] = 0.0;
    }

    for (int iClass = 0; iClass < m_NumClasses; iClass++) {
      double fCount = 0;

      for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
	double iCPT = 0;

	for (int iParent = 0; iParent < m_ParentSets[iAttribute].getNrOfParents(); iParent++) {
	  int nParent = m_ParentSets[iAttribute].getParent(iParent);

	  if (nParent == m_Instances.classIndex()) {
	    iCPT = iCPT * m_NumClasses + iClass;
	  } else {
	    iCPT = iCPT * m_Instances.attribute(nParent).numValues() + instance.value(nParent);
	  }
	}

	if (iAttribute == m_Instances.classIndex()) {
	  fCount += ((DiscreteEstimatorBayes) m_Distributions[iAttribute][(int) iCPT]).getCount(iClass);
	} else {
	  fCount
	  += ((DiscreteEstimatorBayes) m_Distributions[iAttribute][(int) iCPT]).getCount(
	      instance.value(iAttribute));
	}
      }

      fCounts[iClass] += fCount;
    }
    return fCounts;
  } // countsForInstance

  /**
   * Returns an enumeration describing the available options
   * 
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector newVector = new Vector(4);

    newVector.addElement(new Option("\tDo not use ADTree data structure\n", "D", 0, "-D"));
    newVector.addElement(new Option("\tBIF file to compare with\n", "B", 1, "-B <BIF file>"));
    newVector.addElement(new Option("\tSearch algorithm\n", "Q", 1, "-Q weka.classifiers.bayes.net.search.SearchAlgorithm"));
    newVector.addElement(new Option("\tEstimator algorithm\n", "E", 1, "-E weka.classifiers.bayes.net.estimate.SimpleEstimator"));

    return newVector.elements();
  } // listOptions

  /**
   * Parses a given list of options. <p>
   * 
     <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Do not use ADTree data structure
   * </pre>
   * 
   * <pre> -B &lt;BIF file&gt;
   *  BIF file to compare with
   * </pre>
   * 
   * <pre> -Q weka.classifiers.bayes.net.search.SearchAlgorithm
   *  Search algorithm
   * </pre>
   * 
   * <pre> -E weka.classifiers.bayes.net.estimate.SimpleEstimator
   *  Estimator algorithm
   * </pre>
   * 
     <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    m_bUseADTree = !(Utils.getFlag('D', options));

    String sBIFFile = Utils.getOption('B', options);
    if (sBIFFile != null && !sBIFFile.equals("")) {
      setBIFFile(sBIFFile);
    }

    String searchAlgorithmName = Utils.getOption('Q', options);
    if (searchAlgorithmName.length() != 0) {
      setSearchAlgorithm(
	  (SearchAlgorithm) Utils.forName(
	      SearchAlgorithm.class,
	      searchAlgorithmName,
	      partitionOptions(options)));
    }
    else {
      setSearchAlgorithm(new K2());
    }


    String estimatorName = Utils.getOption('E', options);
    if (estimatorName.length() != 0) {
      setEstimator(
	  (BayesNetEstimator) Utils.forName(
	      BayesNetEstimator.class,
	      estimatorName,
	      Utils.partitionOptions(options)));
    }
    else {
      setEstimator(new SimpleEstimator());
    }

    Utils.checkForRemainingOptions(options);
  } // setOptions

  /**
   * Returns the secondary set of options (if any) contained in
   * the supplied options array. The secondary set is defined to
   * be any options after the first "--" but before the "-E". These 
   * options are removed from the original options array.
   *
   * @param options the input array of options
   * @return the array of secondary options
   */
  public static String [] partitionOptions(String [] options) {

    for (int i = 0; i < options.length; i++) {
      if (options[i].equals("--")) {
	// ensure it follows by a -E option
	int j = i;
	while ((j < options.length) && !(options[j].equals("-E"))) {
	  j++;
	}
        /*	if (j >= options.length) {
	  return new String[0];
          } */
	options[i++] = "";
	String [] result = new String [options.length - i];
	j = i;
	while ((j < options.length) && !(options[j].equals("-E"))) {
	  result[j - i] = options[j];
	  options[j] = "";
	  j++;
	}
	while(j < options.length) {
	  result[j - i] = "";
	  j++;
	}		 
	return result;
      }
    }
    return new String [0];
  }


  /**
   * Gets the current settings of the classifier.
   * 
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    String[] searchOptions = m_SearchAlgorithm.getOptions();
    String[] estimatorOptions = m_BayesNetEstimator.getOptions();
    String[] options = new String[11 + searchOptions.length + estimatorOptions.length];
    int current = 0;

    if (!m_bUseADTree) {
      options[current++] = "-D";
    }

    if (m_otherBayesNet != null) {
      options[current++] = "-B";
      options[current++] = ((BIFReader) m_otherBayesNet).getFileName();
    }

    options[current++] = "-Q";
    options[current++] = "" + getSearchAlgorithm().getClass().getName();
    options[current++] = "--";
    for (int iOption = 0; iOption < searchOptions.length; iOption++) {
      options[current++] = searchOptions[iOption];
    }

    options[current++] = "-E";
    options[current++] = "" + getEstimator().getClass().getName();
    options[current++] = "--";
    for (int iOption = 0; iOption < estimatorOptions.length; iOption++) {
      options[current++] = estimatorOptions[iOption];
    }

    // Fill up rest with empty strings, not nulls!
    while (current < options.length) {
      options[current++] = "";
    }

    return options;
  } // getOptions

  /**
   * Set the SearchAlgorithm used in searching for network structures. 
   * @param newSearchAlgorithm the SearchAlgorithm to use.
   */
  public void setSearchAlgorithm(SearchAlgorithm newSearchAlgorithm) {
    m_SearchAlgorithm = newSearchAlgorithm;
  }

  /**
   * Get the SearchAlgorithm used as the search algorithm
   * @return the SearchAlgorithm used as the search algorithm
   */
  public SearchAlgorithm getSearchAlgorithm() {
    return m_SearchAlgorithm;
  }

  /**
   * Set the Estimator Algorithm used in calculating the CPTs 
   * @param newBayesNetEstimator the Estimator to use.
   */
  public void setEstimator(BayesNetEstimator newBayesNetEstimator) {
    m_BayesNetEstimator = newBayesNetEstimator;
  }

  /**
   * Get the BayesNetEstimator used for calculating the CPTs
   * @return the BayesNetEstimator used.
   */
  public BayesNetEstimator getEstimator() {
    return m_BayesNetEstimator;
  }

  /**
   * Set whether ADTree structure is used or not
   * @param bUseADTree true if an ADTree structure is used
   */
  public void setUseADTree(boolean bUseADTree) {
    m_bUseADTree = bUseADTree;
  }

  /**
   * Method declaration
   * @return whether ADTree structure is used or not
   */
  public boolean getUseADTree() {
    return m_bUseADTree;
  }

  /**
   * Set name of network in BIF file to compare with
   * @param sBIFFile the name of the BIF file
   */
  public void setBIFFile(String sBIFFile) {
    try {
      m_otherBayesNet = new BIFReader().processFile(sBIFFile);
    } catch (Throwable t) {
      m_otherBayesNet = null;
    }
  }

  /**
   * Get name of network in BIF file to compare with
   * @return BIF file name
   */
  public String getBIFFile() {
    if (m_otherBayesNet != null) {
      return m_otherBayesNet.getFileName();
    }
    return "";
  }


  /**
   * Returns a description of the classifier.
   * 
   * @return a description of the classifier as a string.
   */
  public String toString() {
    StringBuffer text = new StringBuffer();

    text.append("Bayes Network Classifier");
    text.append("\n" + (m_bUseADTree ? "Using " : "not using ") + "ADTree");

    if (m_Instances == null) {
      text.append(": No model built yet.");
    } else {

      // flatten BayesNet down to text
      text.append("\n#attributes=");
      text.append(m_Instances.numAttributes());
      text.append(" #classindex=");
      text.append(m_Instances.classIndex());
      text.append("\nNetwork structure (nodes followed by parents)\n");

      for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
	text.append(
	    m_Instances.attribute(iAttribute).name()
	    + "("
		+ m_Instances.attribute(iAttribute).numValues()
		+ "): ");

	for (int iParent = 0; iParent < m_ParentSets[iAttribute].getNrOfParents(); iParent++) {
	  text.append(m_Instances.attribute(m_ParentSets[iAttribute].getParent(iParent)).name() + " ");
	}

	text.append("\n");

	// Description of distributions tends to be too much detail, so it is commented out here
	// for (int iParent = 0; iParent < m_ParentSets[iAttribute].GetCardinalityOfParents(); iParent++) {
	// text.append('(' + m_Distributions[iAttribute][iParent].toString() + ')');
	// }
	// text.append("\n");
      }

      text.append("LogScore Bayes: " + measureBayesScore() + "\n");
      text.append("LogScore BDeu: " + measureBDeuScore() + "\n");
      text.append("LogScore MDL: " + measureMDLScore() + "\n");
      text.append("LogScore ENTROPY: " + measureEntropyScore() + "\n");
      text.append("LogScore AIC: " + measureAICScore() + "\n");

      if (m_otherBayesNet != null) {
	text.append(
	    "Missing: "
	    + m_otherBayesNet.missingArcs(this)
	    + " Extra: "
	    + m_otherBayesNet.extraArcs(this)
	    + " Reversed: "
	    + m_otherBayesNet.reversedArcs(this)
	    + "\n");
	text.append("Divergence: " + m_otherBayesNet.divergence(this) + "\n");
      }
    }

    return text.toString();
  } // toString


  /**
   *  Returns the type of graph this classifier
   *  represents.
   *  @return Drawable.TREE
   */   
  public int graphType() {
    return Drawable.BayesNet;
  }

  /**
   * Returns a BayesNet graph in XMLBIF ver 0.3 format.
   * @return String representing this BayesNet in XMLBIF ver  0.3
   * @throws Exception in case BIF generation fails
   */
  public String graph() throws Exception {
    return toXMLBIF03();
  }

  public String getBIFHeader() {
    StringBuffer text = new StringBuffer();
    text.append("<?xml version=\"1.0\"?>\n");
    text.append("<!-- DTD for the XMLBIF 0.3 format -->\n");
    text.append("<!DOCTYPE BIF [\n");
    text.append("	<!ELEMENT BIF ( NETWORK )*>\n");
    text.append("	      <!ATTLIST BIF VERSION CDATA #REQUIRED>\n");
    text.append("	<!ELEMENT NETWORK ( NAME, ( PROPERTY | VARIABLE | DEFINITION )* )>\n");
    text.append("	<!ELEMENT NAME (#PCDATA)>\n");
    text.append("	<!ELEMENT VARIABLE ( NAME, ( OUTCOME |  PROPERTY )* ) >\n");
    text.append("	      <!ATTLIST VARIABLE TYPE (nature|decision|utility) \"nature\">\n");
    text.append("	<!ELEMENT OUTCOME (#PCDATA)>\n");
    text.append("	<!ELEMENT DEFINITION ( FOR | GIVEN | TABLE | PROPERTY )* >\n");
    text.append("	<!ELEMENT FOR (#PCDATA)>\n");
    text.append("	<!ELEMENT GIVEN (#PCDATA)>\n");
    text.append("	<!ELEMENT TABLE (#PCDATA)>\n");
    text.append("	<!ELEMENT PROPERTY (#PCDATA)>\n");
    text.append("]>\n");
    return text.toString();
  } // getBIFHeader

  /**
   * Returns a description of the classifier in XML BIF 0.3 format.
   * See http://www-2.cs.cmu.edu/~fgcozman/Research/InterchangeFormat/
   * for details on XML BIF.
   * @return an XML BIF 0.3 description of the classifier as a string.
   */
  public String toXMLBIF03() {
    if (m_Instances == null) {
      return("<!--No model built yet-->");
    }

    StringBuffer text = new StringBuffer();
    text.append(getBIFHeader());
    text.append("\n");
    text.append("\n");
    text.append("<BIF VERSION=\"0.3\">\n");
    text.append("<NETWORK>\n");
    text.append("<NAME>" + XMLNormalize(m_Instances.relationName()) + "</NAME>\n");
    for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
      text.append("<VARIABLE TYPE=\"nature\">\n");
      text.append("<NAME>" + XMLNormalize(m_Instances.attribute(iAttribute).name()) + "</NAME>\n");
      for (int iValue = 0; iValue < m_Instances.attribute(iAttribute).numValues(); iValue++) {
	text.append("<OUTCOME>" + XMLNormalize(m_Instances.attribute(iAttribute).value(iValue)) + "</OUTCOME>\n");
      }
      text.append("</VARIABLE>\n");
    }

    for (int iAttribute = 0; iAttribute < m_Instances.numAttributes(); iAttribute++) {
      text.append("<DEFINITION>\n");
      text.append("<FOR>" + XMLNormalize(m_Instances.attribute(iAttribute).name()) + "</FOR>\n");
      for (int iParent = 0; iParent < m_ParentSets[iAttribute].getNrOfParents(); iParent++) {
	text.append("<GIVEN>"
	    + XMLNormalize(m_Instances.attribute(m_ParentSets[iAttribute].getParent(iParent)).name()) +
	"</GIVEN>\n");
      }
      text.append("<TABLE>\n");
      for (int iParent = 0; iParent < m_ParentSets[iAttribute].getCardinalityOfParents(); iParent++) {
	for (int iValue = 0; iValue < m_Instances.attribute(iAttribute).numValues(); iValue++) {
	  text.append(m_Distributions[iAttribute][iParent].getProbability(iValue));
	  text.append(' ');
	}
	text.append('\n');
      }
      text.append("</TABLE>\n");
      text.append("</DEFINITION>\n");
    }
    text.append("</NETWORK>\n");
    text.append("</BIF>\n");
    return text.toString();
  } // toXMLBIF03


  /** XMLNormalize converts the five standard XML entities in a string
   * g.e. the string V&D's is returned as V&amp;D&apos;s
   * @param sStr string to normalize
   * @return normalized string
   */
  protected String XMLNormalize(String sStr) {
    StringBuffer sStr2 = new StringBuffer();
    for (int iStr = 0; iStr < sStr.length(); iStr++) {
      char c = sStr.charAt(iStr);
      switch (c) {
	case '&': sStr2.append("&amp;"); break;
	case '\'': sStr2.append("&apos;"); break;
	case '\"': sStr2.append("&quot;"); break;
	case '<': sStr2.append("&lt;"); break;
	case '>': sStr2.append("&gt;"); break;
	default:
	  sStr2.append(c);
      }
    }
    return sStr2.toString();
  } // XMLNormalize


  /**
   * @return a string to describe the UseADTreeoption.
   */
  public String useADTreeTipText() {
    return "When ADTree (the data structure for increasing speed on counts,"
    + " not to be confused with the classifier under the same name) is used"
    + " learning time goes down typically. However, because ADTrees are memory"
    + " intensive, memory problems may occur. Switching this option off makes"
    + " the structure learning algorithms slower, and run with less memory."
    + " By default, ADTrees are used.";
  }

  /**
   * @return a string to describe the SearchAlgorithm.
   */
  public String searchAlgorithmTipText() {
    return "Select method used for searching network structures.";
  }

  /**
   * This will return a string describing the BayesNetEstimator.
   * @return The string.
   */
  public String estimatorTipText() {
    return "Select Estimator algorithm for finding the conditional probability tables"
    + " of the Bayes Network.";
  }

  /**
   * @return a string to describe the BIFFile.
   */
  public String BIFFileTipText() {
    return "Set the name of a file in BIF XML format. A Bayes network learned"
    + " from data can be compared with the Bayes network represented by the BIF file."
    + " Statistics calculated are o.a. the number of missing and extra arcs.";
  }

  /**
   * This will return a string describing the classifier.
   * @return The string.
   */
  public String globalInfo() {
    return 
    "Bayes Network learning using various search algorithms and "
    + "quality measures.\n"
    + "Base class for a Bayes Network classifier. Provides "
    + "datastructures (network structure, conditional probability "
    + "distributions, etc.) and facilities common to Bayes Network "
    + "learning algorithms like K2 and B.\n\n"
    + "For more information see:\n\n"
    + "http://www.cs.waikato.ac.nz/~remco/weka.pdf";
  }

  /**
   * Main method for testing this class.
   * 
   * @param argv the options
   */
  public static void main(String[] argv) {
    runClassifier(new BayesNet(), argv);
  } // main

  /** get name of the Bayes network
   * @return name of the Bayes net
   */
  public String getName() {
    return m_Instances.relationName();
  }

  /** get number of nodes in the Bayes network
   * @return number of nodes
   */
  public int getNrOfNodes() {
    return m_Instances.numAttributes();
  }

  /** get name of a node in the Bayes network
   * @param iNode index of the node
   * @return name of the specified node
   */
  public String getNodeName(int iNode) {
    return m_Instances.attribute(iNode).name();
  }

  /** get number of values a node can take
   * @param iNode index of the node
   * @return cardinality of the specified node
   */
  public int getCardinality(int iNode) {
    return m_Instances.attribute(iNode).numValues();
  }

  /** get name of a particular value of a node
   * @param iNode index of the node
   * @param iValue index of the value
   * @return cardinality of the specified node
   */
  public String getNodeValue(int iNode, int iValue) {
    return m_Instances.attribute(iNode).value(iValue);
  }

  /** get number of parents of a node in the network structure
   * @param iNode index of the node
   * @return number of parents of the specified node
   */
  public int getNrOfParents(int iNode) {
    return m_ParentSets[iNode].getNrOfParents();
  }

  /** get node index of a parent of a node in the network structure
   * @param iNode index of the node
   * @param iParent index of the parents, e.g., 0 is the first parent, 1 the second parent, etc.
   * @return node index of the iParent's parent of the specified node
   */
  public int getParent(int iNode, int iParent) {
    return m_ParentSets[iNode].getParent(iParent);
  }

  /** Get full set of parent sets.
   * @return parent sets;
   */
  public ParentSet[] getParentSets() { 
    return m_ParentSets;
  }

  /** Get full set of estimators.
   * @return estimators;
   */
  public Estimator[][] getDistributions() {
    return m_Distributions;
  }

  /** get number of values the collection of parents of a node can take
   * @param iNode index of the node
   * @return cardinality of the parent set of the specified node
   */
  public int getParentCardinality(int iNode) {
    return m_ParentSets[iNode].getCardinalityOfParents();
  }

  /** get particular probability of the conditional probability distribtion
   * of a node given its parents.
   * @param iNode index of the node
   * @param iParent index of the parent set, 0 <= iParent <= getParentCardinality(iNode)
   * @param iValue index of the value, 0 <= iValue <= getCardinality(iNode)
   * @return probability
   */
  public double getProbability(int iNode, int iParent, int iValue) {
    return m_Distributions[iNode][iParent].getProbability(iValue);
  }

  /** get the parent set of a node 
   * @param iNode index of the node
   * @return Parent set of the specified node.
   */
  public ParentSet getParentSet(int iNode) {
    return m_ParentSets[iNode];
  }

  /** get ADTree strucrture containing efficient representation of counts.
   * @return ADTree strucrture
   */
  public ADNode getADTree() { return m_ADTree;}

  // implementation of AdditionalMeasureProducer interface
  /**
   * Returns an enumeration of the measure names. Additional measures
   * must follow the naming convention of starting with "measure", eg.
   * double measureBlah()
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector(4);
    newVector.addElement("measureExtraArcs");
    newVector.addElement("measureMissingArcs");
    newVector.addElement("measureReversedArcs");
    newVector.addElement("measureDivergence");
    newVector.addElement("measureBayesScore");
    newVector.addElement("measureBDeuScore");
    newVector.addElement("measureMDLScore");
    newVector.addElement("measureAICScore");
    newVector.addElement("measureEntropyScore");
    return newVector.elements();
  } // enumerateMeasures

  public double measureExtraArcs() {
    if (m_otherBayesNet != null) {
      return m_otherBayesNet.extraArcs(this); 
    }
    return 0;
  } // measureExtraArcs

  public double measureMissingArcs() {
    if (m_otherBayesNet != null) {
      return m_otherBayesNet.missingArcs(this); 
    }
    return 0;
  } // measureMissingArcs

  public double measureReversedArcs() {
    if (m_otherBayesNet != null) {
      return m_otherBayesNet.reversedArcs(this); 
    }
    return 0;
  } // measureReversedArcs

  public double measureDivergence() {
    if (m_otherBayesNet != null) {
      return m_otherBayesNet.divergence(this); 
    }
    return 0;
  } // measureDivergence

  public double measureBayesScore() {
    LocalScoreSearchAlgorithm s = new LocalScoreSearchAlgorithm(this, m_Instances);
    return s.logScore(Scoreable.BAYES);
  } // measureBayesScore

  public double measureBDeuScore() {
    LocalScoreSearchAlgorithm s = new LocalScoreSearchAlgorithm(this, m_Instances);
    return s.logScore(Scoreable.BDeu);
  } // measureBDeuScore

  public double measureMDLScore() {
    LocalScoreSearchAlgorithm s = new LocalScoreSearchAlgorithm(this, m_Instances);
    return s.logScore(Scoreable.MDL);
  } // measureMDLScore

  public double measureAICScore() {
    LocalScoreSearchAlgorithm s = new LocalScoreSearchAlgorithm(this, m_Instances);
    return s.logScore(Scoreable.AIC);
  } // measureAICScore

  public double measureEntropyScore() {
    LocalScoreSearchAlgorithm s = new LocalScoreSearchAlgorithm(this, m_Instances);
    return s.logScore(Scoreable.ENTROPY);
  } // measureEntropyScore

  /**
   * Returns the value of the named measure
   * @param measureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws IllegalArgumentException if the named measure is not supported
   */
  public double getMeasure(String measureName) {
    if (measureName.equals("measureExtraArcs")) {
      return measureExtraArcs();
    }
    if (measureName.equals("measureMissingArcs")) {
      return measureMissingArcs();
    }
    if (measureName.equals("measureReversedArcs")) {
      return measureReversedArcs();
    }
    if (measureName.equals("measureDivergence")) {
      return measureDivergence();
    }
    if (measureName.equals("measureBayesScore")) {
      return measureBayesScore();
    }
    if (measureName.equals("measureBDeuScore")) {
      return measureBDeuScore();
    }
    if (measureName.equals("measureMDLScore")) {
      return measureMDLScore();
    }
    if (measureName.equals("measureAICScore")) {
      return measureAICScore();
    }
    if (measureName.equals("measureEntropyScore")) {
      return measureEntropyScore();
    }
    return 0;
  } // getMeasure

  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5928 $");
  }
} // class BayesNet

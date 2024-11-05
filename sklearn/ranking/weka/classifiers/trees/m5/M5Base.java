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
 *    M5Base.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.m5;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.functions.LinearRegression;
import weka.core.AdditionalMeasureProducer;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.supervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.RemoveUseless;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * M5Base. Implements base routines
 * for generating M5 Model trees and rules. <p>
 * 
 * The original algorithm M5 was invented by Quinlan: <br/>
 * 
 * Quinlan J. R. (1992). Learning with continuous classes. Proceedings of
 * the Australian Joint Conference on Artificial Intelligence. 343--348.
 * World Scientific, Singapore. <p/>
 * 
 * Yong Wang made improvements and created M5': <br/>
 * 
 * Wang, Y and Witten, I. H. (1997). Induction of model trees for
 * predicting continuous classes. Proceedings of the poster papers of the
 * European Conference on Machine Learning. University of Economics,
 * Faculty of Informatics and Statistics, Prague. <p/>
 *
 * Valid options are:<p>
 * 
 * -U <br>
 * Use unsmoothed predictions. <p>
 *
 * -R <br>
 * Build regression tree/rule rather than model tree/rule
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 6237 $
 */
public abstract class M5Base 
  extends AbstractClassifier 
  implements AdditionalMeasureProducer,
	     TechnicalInformationHandler {

  /** for serialization */
  private static final long serialVersionUID = -4022221950191647679L;

  /**
   * the instances covered by the tree/rules
   */
  private Instances m_instances;

  /**
   * the rule set
   */
  protected FastVector m_ruleSet;

  /**
   * generate a decision list instead of a single tree.
   */
  private boolean m_generateRules;

  /**
   * use unsmoothed predictions
   */
  private boolean m_unsmoothedPredictions;

  /**
   * filter to fill in missing values
   */
  private ReplaceMissingValues m_replaceMissing;

  /**
   * filter to convert nominal attributes to binary
   */
  private NominalToBinary m_nominalToBinary;
  
  /**
   * for removing useless attributes 
   */
  private RemoveUseless m_removeUseless;

  /**
   * Save instances at each node in an M5 tree for visualization purposes.
   */
  protected boolean m_saveInstances = false;

  /**
   * Make a regression tree/rule instead of a model tree/rule
   */
  protected boolean m_regressionTree;

  /**
   * Do not prune tree/rules
   */
  protected boolean m_useUnpruned = false;

  /**
   * The minimum number of instances to allow at a leaf node
   */
  protected double m_minNumInstances = 4;

  /**
   * Constructor
   */
  public M5Base() {
    m_generateRules = false;
    m_unsmoothedPredictions = false;
    m_useUnpruned = false;
    m_minNumInstances = 4;
  }

  /**
   * returns information about the classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "M5Base. Implements base routines for generating M5 Model trees and " 
      + "rules\n"
      + "The original algorithm M5 was invented by R. Quinlan and Yong Wang "
      + "made improvements.\n\n"
      + "For more information see:\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    TechnicalInformation 	additional;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Ross J. Quinlan");
    result.setValue(Field.TITLE, "Learning with Continuous Classes");
    result.setValue(Field.BOOKTITLE, "5th Australian Joint Conference on Artificial Intelligence");
    result.setValue(Field.YEAR, "1992");
    result.setValue(Field.PAGES, "343-348");
    result.setValue(Field.PUBLISHER, "World Scientific");
    result.setValue(Field.ADDRESS, "Singapore");
    
    additional = result.add(Type.INPROCEEDINGS);
    additional.setValue(Field.AUTHOR, "Y. Wang and I. H. Witten");
    additional.setValue(Field.TITLE, "Induction of model trees for predicting continuous classes");
    additional.setValue(Field.BOOKTITLE, "Poster papers of the 9th European Conference on Machine Learning");
    additional.setValue(Field.YEAR, "1997");
    additional.setValue(Field.PUBLISHER, "Springer");
    
    return result;
  }

  /**
   * Returns an enumeration describing the available options
   * 
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector newVector = new Vector(4);

    newVector.addElement(new Option("\tUse unpruned tree/rules", 
				    "N", 0, "-N"));

    newVector.addElement(new Option("\tUse unsmoothed predictions", 
				    "U", 0, "-U"));

    newVector.addElement(new Option("\tBuild regression tree/rule rather "
				    +"than a model tree/rule", 
				    "R", 0, "-R"));

    newVector.addElement(new Option("\tSet minimum number of instances "
				    +"per leaf\n\t(default 4)",
				    "M",1,"-M <minimum number of instances>"));
    return newVector.elements();
  } 

  /**
   * Parses a given list of options. <p/>
   * 
   * Valid options are:<p>
   * 
   * -U <br>
   * Use unsmoothed predictions. <p>
   *
   * -R <br>
   * Build a regression tree rather than a model tree. <p>
   * 
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    setUnpruned(Utils.getFlag('N', options));
    setUseUnsmoothed(Utils.getFlag('U', options));
    setBuildRegressionTree(Utils.getFlag('R', options));
    String optionString = Utils.getOption('M', options);
    if (optionString.length() != 0) {
      setMinNumInstances((new Double(optionString)).doubleValue());
    }
    Utils.checkForRemainingOptions(options);
  } 

  /**
   * Gets the current settings of the classifier.
   * 
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    String[] options = new String[5];
    int current = 0;

    if (getUnpruned()) {
      options[current++] = "-N";
    }

    if (getUseUnsmoothed()) {
      options[current++] = "-U";
    } 

    if (getBuildRegressionTree()) {
      options[current++] = "-R";
    }

    options[current++] = "-M"; 
    options[current++] = ""+getMinNumInstances();

    while (current < options.length) {
      options[current++] = "";
    } 
    return options;
  } 

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String unprunedTipText() {
    return "Whether unpruned tree/rules are to be generated.";
  }

  /**
   * Use unpruned tree/rules
   *
   * @param unpruned true if unpruned tree/rules are to be generated
   */
  public void setUnpruned(boolean unpruned) {
    m_useUnpruned = unpruned;
  }

  /**
   * Get whether unpruned tree/rules are being generated
   *
   * @return true if unpruned tree/rules are to be generated
   */
  public boolean getUnpruned() {
    return m_useUnpruned;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String generateRulesTipText() {
    return "Whether to generate rules (decision list) rather than a tree.";
  }

  /**
   * Generate rules (decision list) rather than a tree
   * 
   * @param u true if rules are to be generated
   */
  protected void setGenerateRules(boolean u) {
    m_generateRules = u;
  } 

  /**
   * get whether rules are being generated rather than a tree
   * 
   * @return true if rules are to be generated
   */
  protected boolean getGenerateRules() {
    return m_generateRules;
  } 

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String useUnsmoothedTipText() {
    return "Whether to use unsmoothed predictions.";
  }

  /**
   * Use unsmoothed predictions
   * 
   * @param s true if unsmoothed predictions are to be used
   */
  public void setUseUnsmoothed(boolean s) {
    m_unsmoothedPredictions = s;
  } 

  /**
   * Get whether or not smoothing is being used
   * 
   * @return true if unsmoothed predictions are to be used
   */
  public boolean getUseUnsmoothed() {
    return m_unsmoothedPredictions;
  } 

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String buildRegressionTreeTipText() {
    return "Whether to generate a regression tree/rule instead of a model tree/rule.";
  }

  /**
   * Get the value of regressionTree.
   *
   * @return Value of regressionTree.
   */
  public boolean getBuildRegressionTree() {
    
    return m_regressionTree;
  }
  
  /**
   * Set the value of regressionTree.
   *
   * @param newregressionTree Value to assign to regressionTree.
   */
  public void setBuildRegressionTree(boolean newregressionTree) {
    
    m_regressionTree = newregressionTree;
  }

  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String minNumInstancesTipText() {
    return "The minimum number of instances to allow at a leaf node.";
  }

  /**
   * Set the minimum number of instances to allow at a leaf node
   *
   * @param minNum the minimum number of instances
   */
  public void setMinNumInstances(double minNum) {
    m_minNumInstances = minNum;
  }

  /**
   * Get the minimum number of instances to allow at a leaf node
   *
   * @return a <code>double</code> value
   */
  public double getMinNumInstances() {
    return m_minNumInstances;
  }

  /**
   * Returns default capabilities of the classifier, i.e., of LinearRegression.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    return new LinearRegression().getCapabilities();
  }

  /**
   * Generates the classifier.
   * 
   * @param data set of instances serving as training data
   * @throws Exception if the classifier has not been generated
   * successfully
   */
  public void buildClassifier(Instances data) throws Exception {
    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    data = new Instances(data);
    data.deleteWithMissingClass();
    
    m_instances = new Instances(data);

    m_replaceMissing = new ReplaceMissingValues();
    m_replaceMissing.setInputFormat(m_instances);
    m_instances = Filter.useFilter(m_instances, m_replaceMissing);

    m_nominalToBinary = new NominalToBinary();
    m_nominalToBinary.setInputFormat(m_instances);
    m_instances = Filter.useFilter(m_instances, m_nominalToBinary);

    m_removeUseless = new RemoveUseless();
    m_removeUseless.setInputFormat(m_instances);
    m_instances = Filter.useFilter(m_instances, m_removeUseless);
    
    m_instances.randomize(new Random(1));

    m_ruleSet = new FastVector();

    Rule tempRule;

    if (m_generateRules) {
      Instances tempInst = m_instances;
     
      do {
	tempRule = new Rule();
	tempRule.setSmoothing(!m_unsmoothedPredictions);
	tempRule.setRegressionTree(m_regressionTree);
	tempRule.setUnpruned(m_useUnpruned);
	tempRule.setSaveInstances(false);
	tempRule.setMinNumInstances(m_minNumInstances);
	tempRule.buildClassifier(tempInst);
	m_ruleSet.addElement(tempRule);
	//	System.err.println("Built rule : "+tempRule.toString());
	tempInst = tempRule.notCoveredInstances();
	tempRule.freeNotCoveredInstances();
      } while (tempInst.numInstances() > 0);
    } else {
      // just build a single tree
      tempRule = new Rule();

      tempRule.setUseTree(true);
      //      tempRule.setGrowFullTree(true);
      tempRule.setSmoothing(!m_unsmoothedPredictions);
      tempRule.setSaveInstances(m_saveInstances);
      tempRule.setRegressionTree(m_regressionTree);
      tempRule.setUnpruned(m_useUnpruned);
      tempRule.setMinNumInstances(m_minNumInstances);

      Instances temp_train;

      temp_train = m_instances;

      tempRule.buildClassifier(temp_train);

      m_ruleSet.addElement(tempRule);      

      //      System.err.print(tempRule.m_topOfTree.treeToString(0));
    }
    
    // save space
    m_instances = new Instances(m_instances, 0);
  } 

  /**
   * Calculates a prediction for an instance using a set of rules
   * or an M5 model tree
   * 
   * @param inst the instance whos class value is to be predicted
   * @return the prediction
   * @throws Exception if a prediction can't be made.
   */
  public double classifyInstance(Instance inst) throws Exception {
    Rule   temp;
    double prediction = 0;
    boolean success = false;

    m_replaceMissing.input(inst);
    inst = m_replaceMissing.output();
    m_nominalToBinary.input(inst);
    inst = m_nominalToBinary.output();
    m_removeUseless.input(inst);
    inst = m_removeUseless.output();

    if (m_ruleSet == null) {
      throw new Exception("Classifier has not been built yet!");
    } 

    if (!m_generateRules) {
      temp = (Rule) m_ruleSet.elementAt(0);
      return temp.classifyInstance(inst);
    } 

    boolean cont;
    int     i;

    for (i = 0; i < m_ruleSet.size(); i++) {
      cont = false;
      temp = (Rule) m_ruleSet.elementAt(i);

      try {
	prediction = temp.classifyInstance(inst);
	success = true;
      } catch (Exception e) {
	cont = true;
      } 

      if (!cont) {
	break;
      } 
    } 

    if (!success) {
      System.out.println("Error in predicting (DecList)");
    } 
    return prediction;
  } 

  /**
   * Returns a description of the classifier
   * 
   * @return a description of the classifier as a String
   */
  public String toString() {
    StringBuffer text = new StringBuffer();
    Rule	 temp;

    if (m_ruleSet == null) {
      return "Classifier hasn't been built yet!";
    } 

    if (m_generateRules) {
      text.append("M5 "
		  + ((m_useUnpruned == true)
		     ? "unpruned "
		     : "pruned ")
		  + ((m_regressionTree == true) 
		     ?  "regression "
		     : "model ")
		  + "rules ");

      if (!m_unsmoothedPredictions) {
	text.append("\n(using smoothed linear models) ");
      }

      text.append(":\n");

      text.append("Number of Rules : " + m_ruleSet.size() + "\n\n");

      for (int j = 0; j < m_ruleSet.size(); j++) {
	temp = (Rule) m_ruleSet.elementAt(j);

	text.append("Rule: " + (j + 1) + "\n");
	text.append(temp.toString());
      } 
    } else {
      temp = (Rule) m_ruleSet.elementAt(0);
      text.append(temp.toString());
    } 
    return text.toString();
  } 

  /**
   * Returns an enumeration of the additional measure names
   * @return an enumeration of the measure names
   */
  public Enumeration enumerateMeasures() {
    Vector newVector = new Vector(1);
    newVector.addElement("measureNumRules");
    return newVector.elements();
  }

  /**
   * Returns the value of the named measure
   * @param additionalMeasureName the name of the measure to query for its value
   * @return the value of the named measure
   * @throws Exception if the named measure is not supported
   */
  public double getMeasure(String additionalMeasureName) 
    {
    if (additionalMeasureName.compareToIgnoreCase("measureNumRules") == 0) {
      return measureNumRules();
    } else {
      throw new IllegalArgumentException(additionalMeasureName 
					 + " not supported (M5)");
    }
  }

  /**
   * return the number of rules
   * @return the number of rules (same as # linear models &
   * # leaves in the tree)
   */
  public double measureNumRules() {
    if (m_generateRules) {
      return m_ruleSet.size();
    }
    return ((Rule)m_ruleSet.elementAt(0)).m_topOfTree.numberOfLinearModels();
  }

  public RuleNode getM5RootNode() {
    Rule temp = (Rule) m_ruleSet.elementAt(0);
    return temp.getM5RootNode();
  }
}

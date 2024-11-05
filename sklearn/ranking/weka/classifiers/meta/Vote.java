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
 *    Vote.java
 *    Copyright (C) 2000 University of Waikato
 *    Copyright (C) 2006 Roberto Perdisci
 *
 */

package weka.classifiers.meta;

import weka.classifiers.RandomizableMultipleClassifiersCombiner;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Class for combining classifiers. Different combinations of probability estimates for classification are available.<br/>
 * <br/>
 * For more information see:<br/>
 * <br/>
 * Ludmila I. Kuncheva (2004). Combining Pattern Classifiers: Methods and Algorithms. John Wiley and Sons, Inc..<br/>
 * <br/>
 * J. Kittler, M. Hatef, Robert P.W. Duin, J. Matas (1998). On combining classifiers. IEEE Transactions on Pattern Analysis and Machine Intelligence. 20(3):226-239.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;num&gt;
 *  Random number seed.
 *  (default 1)</pre>
 * 
 * <pre> -B &lt;classifier specification&gt;
 *  Full class name of classifier to include, followed
 *  by scheme options. May be specified multiple times.
 *  (default: "weka.classifiers.rules.ZeroR")</pre>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 * <pre> -R &lt;AVG|PROD|MAJ|MIN|MAX|MED&gt;
 *  The combination rule to use
 *  (default: AVG)</pre>
 * 
 <!-- options-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;book{Kuncheva2004,
 *    author = {Ludmila I. Kuncheva},
 *    publisher = {John Wiley and Sons, Inc.},
 *    title = {Combining Pattern Classifiers: Methods and Algorithms},
 *    year = {2004}
 * }
 * 
 * &#64;article{Kittler1998,
 *    author = {J. Kittler and M. Hatef and Robert P.W. Duin and J. Matas},
 *    journal = {IEEE Transactions on Pattern Analysis and Machine Intelligence},
 *    number = {3},
 *    pages = {226-239},
 *    title = {On combining classifiers},
 *    volume = {20},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 * @author Alexander K. Seewald (alex@seewald.at)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Roberto Perdisci (roberto.perdisci@gmail.com)
 * @version $Revision: 5987 $
 */
public class Vote
  extends RandomizableMultipleClassifiersCombiner
  implements TechnicalInformationHandler {
    
  /** for serialization */
  static final long serialVersionUID = -637891196294399624L;
  
  /** combination rule: Average of Probabilities */
  public static final int AVERAGE_RULE = 1;
  /** combination rule: Product of Probabilities (only nominal classes) */
  public static final int PRODUCT_RULE = 2;
  /** combination rule: Majority Voting (only nominal classes) */
  public static final int MAJORITY_VOTING_RULE = 3;
  /** combination rule: Minimum Probability */
  public static final int MIN_RULE = 4;
  /** combination rule: Maximum Probability */
  public static final int MAX_RULE = 5;
  /** combination rule: Median Probability (only numeric class) */
  public static final int MEDIAN_RULE = 6;
  /** combination rules */
  public static final Tag[] TAGS_RULES = {
    new Tag(AVERAGE_RULE, "AVG", "Average of Probabilities"),
    new Tag(PRODUCT_RULE, "PROD", "Product of Probabilities"),
    new Tag(MAJORITY_VOTING_RULE, "MAJ", "Majority Voting"),
    new Tag(MIN_RULE, "MIN", "Minimum Probability"),
    new Tag(MAX_RULE, "MAX", "Maximum Probability"),
    new Tag(MEDIAN_RULE, "MED", "Median")
  };
  
  /** Combination Rule variable */
  protected int m_CombinationRule = AVERAGE_RULE;
  
  /** the random number generator used for breaking ties in majority voting
   * @see #distributionForInstanceMajorityVoting(Instance) */
  protected Random m_Random;
  
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Class for combining classifiers. Different combinations of "
      + "probability estimates for classification are available.\n\n"
      + "For more information see:\n\n"
      + getTechnicalInformation().toString();
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Enumeration 	enm;
    Vector		result;
    
    result = new Vector();
    
    enm = super.listOptions();
    while (enm.hasMoreElements())
      result.addElement(enm.nextElement());

    result.addElement(new Option(
	"\tThe combination rule to use\n"
	+ "\t(default: AVG)",
	"R", 1, "-R " + Tag.toOptionList(TAGS_RULES)));
    
    return result.elements();
  }
  
  /**
   * Gets the current settings of Vote.
   *
   * @return an array of strings suitable for passing to setOptions()
   */
  public String [] getOptions() {
    int       	i;
    Vector    	result;
    String[]  	options;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    result.add("-R");
    result.add("" + getCombinationRule());

    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;num&gt;
   *  Random number seed.
   *  (default 1)</pre>
   * 
   * <pre> -B &lt;classifier specification&gt;
   *  Full class name of classifier to include, followed
   *  by scheme options. May be specified multiple times.
   *  (default: "weka.classifiers.rules.ZeroR")</pre>
   * 
   * <pre> -D
   *  If set, classifier is run in debug mode and
   *  may output additional info to the console</pre>
   * 
   * <pre> -R &lt;AVG|PROD|MAJ|MIN|MAX|MED&gt;
   *  The combination rule to use
   *  (default: AVG)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String 	tmpStr;
    
    tmpStr = Utils.getOption('R', options);
    if (tmpStr.length() != 0) 
      setCombinationRule(new SelectedTag(tmpStr, TAGS_RULES));
    else
      setCombinationRule(new SelectedTag(AVERAGE_RULE, TAGS_RULES));

    super.setOptions(options);
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
    
    result = new TechnicalInformation(Type.BOOK);
    result.setValue(Field.AUTHOR, "Ludmila I. Kuncheva");
    result.setValue(Field.TITLE, "Combining Pattern Classifiers: Methods and Algorithms");
    result.setValue(Field.YEAR, "2004");
    result.setValue(Field.PUBLISHER, "John Wiley and Sons, Inc.");

    additional = result.add(Type.ARTICLE);
    additional.setValue(Field.AUTHOR, "J. Kittler and M. Hatef and Robert P.W. Duin and J. Matas");
    additional.setValue(Field.YEAR, "1998");
    additional.setValue(Field.TITLE, "On combining classifiers");
    additional.setValue(Field.JOURNAL, "IEEE Transactions on Pattern Analysis and Machine Intelligence");
    additional.setValue(Field.VOLUME, "20");
    additional.setValue(Field.NUMBER, "3");
    additional.setValue(Field.PAGES, "226-239");
    
    return result;
  }
  
  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();

    // class
    if (    (m_CombinationRule == PRODUCT_RULE) 
	 || (m_CombinationRule == MAJORITY_VOTING_RULE) ) {
      result.disableAllClasses();
      result.disableAllClassDependencies();
      result.enable(Capability.NOMINAL_CLASS);
      result.enableDependency(Capability.NOMINAL_CLASS);
    }
    else if (m_CombinationRule == MEDIAN_RULE) {
      result.disableAllClasses();
      result.disableAllClassDependencies();
      result.enable(Capability.NUMERIC_CLASS);
      result.enableDependency(Capability.NUMERIC_CLASS);
    }
    
    //RANKING BEGIN
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    result.disable(Capability.RANKING);
    //RANKING END
    
    return result;
  }

  /**
   * Buildclassifier selects a classifier from the set of classifiers
   * by minimising error on the training data.
   *
   * @param data the training data to be used for generating the
   * boosted classifier.
   * @throws Exception if the classifier could not be built successfully
   */
  public void buildClassifier(Instances data) throws Exception {

    // can classifier handle the data?
    getCapabilities().testWithFail(data);

    // remove instances with missing class
    Instances newData = new Instances(data);
    newData.deleteWithMissingClass();

    m_Random = new Random(getSeed());
    
    for (int i = 0; i < m_Classifiers.length; i++) {
      getClassifier(i).buildClassifier(newData);
    }
  }

  /**
   * Classifies the given test instance.
   *
   * @param instance the instance to be classified
   * @return the predicted most likely class for the instance or 
   * Utils.missingValue() if no prediction is made
   * @throws Exception if an error occurred during the prediction
   */
  public double classifyInstance(Instance instance) throws Exception {
    double result;
    double[] dist;
    int index;
    
    switch (m_CombinationRule) {
      case AVERAGE_RULE:
      case PRODUCT_RULE:
      case MAJORITY_VOTING_RULE:
      case MIN_RULE:
      case MAX_RULE:
	dist = distributionForInstance(instance);
	if (instance.classAttribute().isNominal() || instance.classAttribute().isRanking()) {
	  index = Utils.maxIndex(dist);
	  if (dist[index] == 0)
	    result = Utils.missingValue();
	  else
	    result = index;
	}
	else if (instance.classAttribute().isNumeric()){
	  result = dist[0];
	}
	else {
	  result = Utils.missingValue();
	}
	break;
      case MEDIAN_RULE:
	result = classifyInstanceMedian(instance);
	break;
      default:
	throw new IllegalStateException("Unknown combination rule '" + m_CombinationRule + "'!");
    }
    
    return result;
  }

  /**
   * Classifies the given test instance, returning the median from all
   * classifiers.
   *
   * @param instance the instance to be classified
   * @return the predicted most likely class for the instance or 
   * Utils.missingValue() if no prediction is made
   * @throws Exception if an error occurred during the prediction
   */
  protected double classifyInstanceMedian(Instance instance) throws Exception {
    double[] results = new double[m_Classifiers.length];
    double result;

    for (int i = 0; i < results.length; i++)
      results[i] = m_Classifiers[i].classifyInstance(instance);
    
    if (results.length == 0)
      result = 0;
    else if (results.length == 1)
      result = results[0];
    else
      result = Utils.kthSmallestValue(results, results.length / 2);
    
    return result;
  }

  /**
   * Classifies a given instance using the selected combination rule.
   *
   * @param instance the instance to be classified
   * @return the distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  public double[] distributionForInstance(Instance instance) throws Exception {
    double[] result = new double[instance.numClasses()];
    
    switch (m_CombinationRule) {
      case AVERAGE_RULE:
	result = distributionForInstanceAverage(instance);
	break;
      case PRODUCT_RULE:
	result = distributionForInstanceProduct(instance);
	break;
      case MAJORITY_VOTING_RULE:
	result = distributionForInstanceMajorityVoting(instance);
	break;
      case MIN_RULE:
	result = distributionForInstanceMin(instance);
	break;
      case MAX_RULE:
	result = distributionForInstanceMax(instance);
	break;
      case MEDIAN_RULE:
	result[0] = classifyInstance(instance);
	break;
      default:
	throw new IllegalStateException("Unknown combination rule '" + m_CombinationRule + "'!");
    }
    
    if (!instance.classAttribute().isNumeric() && (Utils.sum(result) > 0))
      Utils.normalize(result);
    
    return result;
  }
  
  /**
   * Classifies a given instance using the Average of Probabilities 
   * combination rule.
   *
   * @param instance the instance to be classified
   * @return the distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  protected double[] distributionForInstanceAverage(Instance instance) throws Exception {

    double[] probs = getClassifier(0).distributionForInstance(instance);
    for (int i = 1; i < m_Classifiers.length; i++) {
      double[] dist = getClassifier(i).distributionForInstance(instance);
      for (int j = 0; j < dist.length; j++) {
    	  probs[j] += dist[j];
      }
    }
    for (int j = 0; j < probs.length; j++) {
      probs[j] /= (double)m_Classifiers.length;
    }
    return probs;
  }
  
  /**
   * Classifies a given instance using the Product of Probabilities 
   * combination rule.
   *
   * @param instance the instance to be classified
   * @return the distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  protected double[] distributionForInstanceProduct(Instance instance) throws Exception {

    double[] probs = getClassifier(0).distributionForInstance(instance);
    for (int i = 1; i < m_Classifiers.length; i++) {
      double[] dist = getClassifier(i).distributionForInstance(instance);
      for (int j = 0; j < dist.length; j++) {
    	  probs[j] *= dist[j];
      }
    }
    
    return probs;
  }
  
  /**
   * Classifies a given instance using the Majority Voting combination rule.
   *
   * @param instance the instance to be classified
   * @return the distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  protected double[] distributionForInstanceMajorityVoting(Instance instance) throws Exception {

    double[] probs = new double[instance.classAttribute().numValues()];
    double[] votes = new double[probs.length];
    
    for (int i = 0; i < m_Classifiers.length; i++) {
      probs = getClassifier(i).distributionForInstance(instance);
      int maxIndex = 0;
      for(int j = 0; j<probs.length; j++) {
          if(probs[j] > probs[maxIndex])
        	  maxIndex = j;
      }
      
      // Consider the cases when multiple classes happen to have the same probability
      for (int j=0; j<probs.length; j++) {
	if (probs[j] == probs[maxIndex])
	  votes[j]++;
      }
    }
    
    int tmpMajorityIndex = 0;
    for (int k = 1; k < votes.length; k++) {
      if (votes[k] > votes[tmpMajorityIndex])
	tmpMajorityIndex = k;
    }
    
    // Consider the cases when multiple classes receive the same amount of votes
    Vector<Integer> majorityIndexes = new Vector<Integer>();
    for (int k = 0; k < votes.length; k++) {
      if (votes[k] == votes[tmpMajorityIndex])
	majorityIndexes.add(k);
     }
    // Resolve the ties according to a uniform random distribution
    int majorityIndex = majorityIndexes.get(m_Random.nextInt(majorityIndexes.size()));
    
    //set probs to 0
    for (int k = 0; k<probs.length; k++)
      probs[k] = 0;
    probs[majorityIndex] = 1; //the class that have been voted the most receives 1
    
    return probs;
  }
  
  /**
   * Classifies a given instance using the Maximum Probability combination rule.
   *
   * @param instance the instance to be classified
   * @return the distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  protected double[] distributionForInstanceMax(Instance instance) throws Exception {

    double[] max = getClassifier(0).distributionForInstance(instance);
    for (int i = 1; i < m_Classifiers.length; i++) {
      double[] dist = getClassifier(i).distributionForInstance(instance);
      for (int j = 0; j < dist.length; j++) {
    	  if(max[j]<dist[j])
    		  max[j]=dist[j];
      }
    }
    
    return max;
  }
  
  /**
   * Classifies a given instance using the Minimum Probability combination rule.
   *
   * @param instance the instance to be classified
   * @return the distribution
   * @throws Exception if instance could not be classified
   * successfully
   */
  protected double[] distributionForInstanceMin(Instance instance) throws Exception {

    double[] min = getClassifier(0).distributionForInstance(instance);
    for (int i = 1; i < m_Classifiers.length; i++) {
      double[] dist = getClassifier(i).distributionForInstance(instance);
      for (int j = 0; j < dist.length; j++) {
    	  if(dist[j]<min[j])
    		  min[j]=dist[j];
      }
    }
    
    return min;
  } 
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String combinationRuleTipText() {
    return "The combination rule used.";
  }
  
  /**
   * Gets the combination rule used
   *
   * @return 		the combination rule used
   */
  public SelectedTag getCombinationRule() {
    return new SelectedTag(m_CombinationRule, TAGS_RULES);
  }

  /**
   * Sets the combination rule to use. Values other than
   *
   * @param newRule 	the combination rule method to use
   */
  public void setCombinationRule(SelectedTag newRule) {
    if (newRule.getTags() == TAGS_RULES)
      m_CombinationRule = newRule.getSelectedTag().getID();
  }
  
  /**
   * Output a representation of this classifier
   * 
   * @return a string representation of the classifier
   */
  public String toString() {

    if (m_Classifiers == null) {
      return "Vote: No model built yet.";
    }

    String result = "Vote combines";
    result += " the probability distributions of these base learners:\n";
    for (int i = 0; i < m_Classifiers.length; i++) {
      result += '\t' + getClassifierSpec(i) + '\n';
    }
    result += "using the '";
    
    switch (m_CombinationRule) {
      case AVERAGE_RULE:
	result += "Average of Probabilities";
	break;
	
      case PRODUCT_RULE:
	result += "Product of Probabilities";
	break;
	
      case MAJORITY_VOTING_RULE:
	result += "Majority Voting";
	break;
	
      case MIN_RULE:
	result += "Minimum Probability";
	break;
	
      case MAX_RULE:
	result += "Maximum Probability";
	break;
	
      case MEDIAN_RULE:
	result += "Median Probability";
	break;
	
      default:
	throw new IllegalStateException("Unknown combination rule '" + m_CombinationRule + "'!");
    }
    
    result += "' combination rule \n";

    return result;
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
   * @param argv should contain the following arguments:
   * -t training file [-T test file] [-c class index]
   */
  public static void main(String [] argv) {
    runClassifier(new Vote(), argv);
  }
}

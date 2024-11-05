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
 *    DecisionStump.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Sourcable;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.ContingencyTables;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;

/**
 <!-- globalinfo-start -->
 * Class for building and using a decision stump. Usually used in conjunction with a boosting algorithm. Does regression (based on mean-squared error) or classification (based on entropy). Missing is treated as a separate value.
 * <p/>
 <!-- globalinfo-end -->
 *
 * Typical usage: <p>
 * <code>java weka.classifiers.meta.LogitBoost -I 100 -W weka.classifiers.trees.DecisionStump 
 * -t training_data </code><p>
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  If set, classifier is run in debug mode and
 *  may output additional info to the console</pre>
 * 
 <!-- options-end -->
 * 
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5928 $
 */
public class DecisionStump 
  extends AbstractClassifier 
  implements WeightedInstancesHandler, Sourcable {

  /** for serialization */
  static final long serialVersionUID = 1618384535950391L;
  
  /** The attribute used for classification. */
  private int m_AttIndex;

  /** The split point (index respectively). */
  private double m_SplitPoint;

  /** The distribution of class values or the means in each subset. */
  private double[][] m_Distribution;

  /** The instances used for training. */
  private Instances m_Instances;

  /** a ZeroR model in case no model can be built from the data */
  private Classifier m_ZeroR;
    
  /**
   * Returns a string describing classifier
   * @return a description suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return  "Class for building and using a decision stump. Usually used in "
      + "conjunction with a boosting algorithm. Does regression (based on "
      + "mean-squared error) or classification (based on entropy). Missing "
      + "is treated as a separate value.";
  }

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
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);

    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    
    return result;
  }

  /**
   * Generates the classifier.
   *
   * @param instances set of instances serving as training data 
   * @throws Exception if the classifier has not been generated successfully
   */
  public void buildClassifier(Instances instances) throws Exception {
    
    double bestVal = Double.MAX_VALUE, currVal;
    double bestPoint = -Double.MAX_VALUE;
    int bestAtt = -1, numClasses;

    // can classifier handle the data?
    getCapabilities().testWithFail(instances);

    // remove instances with missing class
    instances = new Instances(instances);
    instances.deleteWithMissingClass();
    
    // only class? -> build ZeroR model
    if (instances.numAttributes() == 1) {
      System.err.println(
	  "Cannot build model (only class attribute present in data!), "
	  + "using ZeroR model instead!");
      m_ZeroR = new weka.classifiers.rules.ZeroR();
      m_ZeroR.buildClassifier(instances);
      return;
    }
    else {
      m_ZeroR = null;
    }
    
    double[][] bestDist = new double[3][instances.numClasses()];

    m_Instances = new Instances(instances);

    if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
      numClasses = m_Instances.numClasses();
    } else {
      numClasses = 1;
    }

    // For each attribute
    boolean first = true;
    for (int i = 0; i < m_Instances.numAttributes(); i++) {
      if (i != m_Instances.classIndex()) {

	// Reserve space for distribution.
	m_Distribution = new double[3][numClasses];

	// Compute value of criterion for best split on attribute
	if (m_Instances.attribute(i).isNominal() || m_Instances.attribute(i).isRanking()) {
	  currVal = findSplitNominal(i);
	} else {
	  currVal = findSplitNumeric(i);
	}
	if ((first) || (currVal < bestVal)) {
	  bestVal = currVal;
	  bestAtt = i;
	  bestPoint = m_SplitPoint;
	  for (int j = 0; j < 3; j++) {
	    System.arraycopy(m_Distribution[j], 0, bestDist[j], 0, 
			     numClasses);
	  }
	}
	
	// First attribute has been investigated
	first = false;
      }
    }
    
    // Set attribute, split point and distribution.
    m_AttIndex = bestAtt;
    m_SplitPoint = bestPoint;
    m_Distribution = bestDist;
    if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
      for (int i = 0; i < m_Distribution.length; i++) {
	double sumCounts = Utils.sum(m_Distribution[i]);
	if (sumCounts == 0) { // This means there were only missing attribute values
	  System.arraycopy(m_Distribution[2], 0, m_Distribution[i], 0, 
			   m_Distribution[2].length);
	  Utils.normalize(m_Distribution[i]);
	} else {
	  Utils.normalize(m_Distribution[i], sumCounts); 
	}
      }
    }
    
    // Save memory
    m_Instances = new Instances(m_Instances, 0);
  }

  /**
   * Calculates the class membership probabilities for the given test instance.
   *
   * @param instance the instance to be classified
   * @return predicted class probability distribution
   * @throws Exception if distribution can't be computed
   */
  public double[] distributionForInstance(Instance instance) throws Exception {

    // default model?
    if (m_ZeroR != null) {
      return m_ZeroR.distributionForInstance(instance);
    }
    
    return m_Distribution[whichSubset(instance)];
  }

  /**
   * Returns the decision tree as Java source code.
   *
   * @param className the classname of the generated code
   * @return the tree as Java source code
   * @throws Exception if something goes wrong
   */
  public String toSource(String className) throws Exception {

    StringBuffer text = new StringBuffer("class ");
    Attribute c = m_Instances.classAttribute();
    text.append(className)
      .append(" {\n"
	      +"  public static double classify(Object[] i) {\n");
    text.append("    /* " + m_Instances.attribute(m_AttIndex).name() + " */\n");
    text.append("    if (i[").append(m_AttIndex);
    text.append("] == null) { return ");
    text.append(sourceClass(c, m_Distribution[2])).append(";");
    if (m_Instances.attribute(m_AttIndex).isNominal() || m_Instances.attribute(m_AttIndex).isRanking()) {
      text.append(" } else if (((String)i[").append(m_AttIndex);
      text.append("]).equals(\"");
      text.append(m_Instances.attribute(m_AttIndex).value((int)m_SplitPoint));
      text.append("\")");
    } else {
      text.append(" } else if (((Double)i[").append(m_AttIndex);
      text.append("]).doubleValue() <= ").append(m_SplitPoint);
    }
    text.append(") { return ");
    text.append(sourceClass(c, m_Distribution[0])).append(";");
    text.append(" } else { return ");
    text.append(sourceClass(c, m_Distribution[1])).append(";");
    text.append(" }\n  }\n}\n");
    return text.toString();
  }

  /**
   * Returns the value as string out of the given distribution
   * 
   * @param c the attribute to get the value for
   * @param dist the distribution to extract the value
   * @return the value
   */
  private String sourceClass(Attribute c, double []dist) {

    if (c.isNominal() || c.isRanking()) {
      return Integer.toString(Utils.maxIndex(dist));
    } else {
      return Double.toString(dist[0]);
    }
  }

  /**
   * Returns a description of the classifier.
   *
   * @return a description of the classifier as a string.
   */
  public String toString(){

    // only ZeroR model?
    if (m_ZeroR != null) {
      StringBuffer buf = new StringBuffer();
      buf.append(this.getClass().getName().replaceAll(".*\\.", "") + "\n");
      buf.append(this.getClass().getName().replaceAll(".*\\.", "").replaceAll(".", "=") + "\n\n");
      buf.append("Warning: No model could be built, hence ZeroR model is used:\n\n");
      buf.append(m_ZeroR.toString());
      return buf.toString();
    }
    
    if (m_Instances == null) {
      return "Decision Stump: No model built yet.";
    }
    try {
      StringBuffer text = new StringBuffer();
      
      text.append("Decision Stump\n\n");
      text.append("Classifications\n\n");
      Attribute att = m_Instances.attribute(m_AttIndex);
      if (att.isNominal() || att.isRanking()) {
	text.append(att.name() + " = " + att.value((int)m_SplitPoint) + 
		    " : ");
	text.append(printClass(m_Distribution[0]));
	text.append(att.name() + " != " + att.value((int)m_SplitPoint) + 
		    " : ");
	text.append(printClass(m_Distribution[1]));
      } else {
	text.append(att.name() + " <= " + m_SplitPoint + " : ");
	text.append(printClass(m_Distribution[0]));
	text.append(att.name() + " > " + m_SplitPoint + " : ");
	text.append(printClass(m_Distribution[1]));
      }
      text.append(att.name() + " is missing : ");
      text.append(printClass(m_Distribution[2]));

      if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
	text.append("\nClass distributions\n\n");
	if (att.isNominal() || att.isRanking()) {
	  text.append(att.name() + " = " + att.value((int)m_SplitPoint) + 
		      "\n");
	  text.append(printDist(m_Distribution[0]));
	  text.append(att.name() + " != " + att.value((int)m_SplitPoint) + 
		      "\n");
	  text.append(printDist(m_Distribution[1]));
	} else {
	  text.append(att.name() + " <= " + m_SplitPoint + "\n");
	  text.append(printDist(m_Distribution[0]));
	  text.append(att.name() + " > " + m_SplitPoint + "\n");
	  text.append(printDist(m_Distribution[1]));
	}
	text.append(att.name() + " is missing\n");
	text.append(printDist(m_Distribution[2]));
      }

      return text.toString();
    } catch (Exception e) {
      return "Can't print decision stump classifier!";
    }
  }

  /** 
   * Prints a class distribution.
   *
   * @param dist the class distribution to print
   * @return the distribution as a string
   * @throws Exception if distribution can't be printed
   */
  private String printDist(double[] dist) throws Exception {

    StringBuffer text = new StringBuffer();
    
    if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
      for (int i = 0; i < m_Instances.numClasses(); i++) {
	text.append(m_Instances.classAttribute().value(i) + "\t");
      }
      text.append("\n");
      for (int i = 0; i < m_Instances.numClasses(); i++) {
	text.append(dist[i] + "\t");
      }
      text.append("\n");
    }
    
    return text.toString();
  }

  /** 
   * Prints a classification.
   *
   * @param dist the class distribution
   * @return the classificationn as a string
   * @throws Exception if the classification can't be printed
   */
  private String printClass(double[] dist) throws Exception {

    StringBuffer text = new StringBuffer();
    
    if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
      text.append(m_Instances.classAttribute().value(Utils.maxIndex(dist)));
    } else {
      text.append(dist[0]);
    }
    
    return text.toString() + "\n";
  }

  /**
   * Finds best split for nominal attribute and returns value.
   *
   * @param index attribute index
   * @return value of criterion for the best split
   * @throws Exception if something goes wrong
   */
  private double findSplitNominal(int index) throws Exception {

    if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
      return findSplitNominalNominal(index);
    } else {
      return findSplitNominalNumeric(index);
    }
  }

  /**
   * Finds best split for nominal attribute and nominal class
   * and returns value.
   *
   * @param index attribute index
   * @return value of criterion for the best split
   * @throws Exception if something goes wrong
   */
  private double findSplitNominalNominal(int index) throws Exception {

    double bestVal = Double.MAX_VALUE, currVal;
    double[][] counts = new double[m_Instances.attribute(index).numValues() 
				  + 1][m_Instances.numClasses()];
    double[] sumCounts = new double[m_Instances.numClasses()];
    double[][] bestDist = new double[3][m_Instances.numClasses()];
    int numMissing = 0;

    // Compute counts for all the values
    for (int i = 0; i < m_Instances.numInstances(); i++) {
      Instance inst = m_Instances.instance(i);
      if (inst.isMissing(index)) {
	numMissing++;
	counts[m_Instances.attribute(index).numValues()]
	  [(int)inst.classValue()] += inst.weight();
      } else {
	counts[(int)inst.value(index)][(int)inst.classValue()] += inst
	  .weight();
      }
    }

    // Compute sum of counts
    for (int i = 0; i < m_Instances.attribute(index).numValues(); i++) {
      for (int j = 0; j < m_Instances.numClasses(); j++) {
	sumCounts[j] += counts[i][j];
      }
    }
    
    // Make split counts for each possible split and evaluate
    System.arraycopy(counts[m_Instances.attribute(index).numValues()], 0,
		     m_Distribution[2], 0, m_Instances.numClasses());
    for (int i = 0; i < m_Instances.attribute(index).numValues(); i++) {
      for (int j = 0; j < m_Instances.numClasses(); j++) {
	m_Distribution[0][j] = counts[i][j];
	m_Distribution[1][j] = sumCounts[j] - counts[i][j];
      }
      currVal = ContingencyTables.entropyConditionedOnRows(m_Distribution);
      if (currVal < bestVal) {
	bestVal = currVal;
	m_SplitPoint = (double)i;
	for (int j = 0; j < 3; j++) {
	  System.arraycopy(m_Distribution[j], 0, bestDist[j], 0, 
			   m_Instances.numClasses());
	}
      }
    }

    // No missing values in training data.
    if (numMissing == 0) {
      System.arraycopy(sumCounts, 0, bestDist[2], 0, 
		       m_Instances.numClasses());
    }
   
    m_Distribution = bestDist;
    return bestVal;
  }

  /**
   * Finds best split for nominal attribute and numeric class
   * and returns value.
   *
   * @param index attribute index
   * @return value of criterion for the best split
   * @throws Exception if something goes wrong
   */
  private double findSplitNominalNumeric(int index) throws Exception {

    double bestVal = Double.MAX_VALUE, currVal;
    double[] sumsSquaresPerValue = 
      new double[m_Instances.attribute(index).numValues()], 
      sumsPerValue = new double[m_Instances.attribute(index).numValues()], 
      weightsPerValue = new double[m_Instances.attribute(index).numValues()];
    double totalSumSquaresW = 0, totalSumW = 0, totalSumOfWeightsW = 0,
      totalSumOfWeights = 0, totalSum = 0;
    double[] sumsSquares = new double[3], sumOfWeights = new double[3];
    double[][] bestDist = new double[3][1];

    // Compute counts for all the values
    for (int i = 0; i < m_Instances.numInstances(); i++) {
      Instance inst = m_Instances.instance(i);
      if (inst.isMissing(index)) {
	m_Distribution[2][0] += inst.classValue() * inst.weight();
	sumsSquares[2] += inst.classValue() * inst.classValue() 
	  * inst.weight();
	sumOfWeights[2] += inst.weight();
      } else {
	weightsPerValue[(int)inst.value(index)] += inst.weight();
	sumsPerValue[(int)inst.value(index)] += inst.classValue() 
	  * inst.weight();
	sumsSquaresPerValue[(int)inst.value(index)] += 
	  inst.classValue() * inst.classValue() * inst.weight();
      }
      totalSumOfWeights += inst.weight();
      totalSum += inst.classValue() * inst.weight();
    }

    // Check if the total weight is zero
    if (totalSumOfWeights <= 0) {
      return bestVal;
    }

    // Compute sum of counts without missing ones
    for (int i = 0; i < m_Instances.attribute(index).numValues(); i++) {
      totalSumOfWeightsW += weightsPerValue[i];
      totalSumSquaresW += sumsSquaresPerValue[i];
      totalSumW += sumsPerValue[i];
    }
    
    // Make split counts for each possible split and evaluate
    for (int i = 0; i < m_Instances.attribute(index).numValues(); i++) {
      
      m_Distribution[0][0] = sumsPerValue[i];
      sumsSquares[0] = sumsSquaresPerValue[i];
      sumOfWeights[0] = weightsPerValue[i];
      m_Distribution[1][0] = totalSumW - sumsPerValue[i];
      sumsSquares[1] = totalSumSquaresW - sumsSquaresPerValue[i];
      sumOfWeights[1] = totalSumOfWeightsW - weightsPerValue[i];

      currVal = variance(m_Distribution, sumsSquares, sumOfWeights);
      
      if (currVal < bestVal) {
	bestVal = currVal;
	m_SplitPoint = (double)i;
	for (int j = 0; j < 3; j++) {
	  if (sumOfWeights[j] > 0) {
	    bestDist[j][0] = m_Distribution[j][0] / sumOfWeights[j];
	  } else {
	    bestDist[j][0] = totalSum / totalSumOfWeights;
	  }
	}
      }
    }

    m_Distribution = bestDist;
    return bestVal;
  }

  /**
   * Finds best split for numeric attribute and returns value.
   *
   * @param index attribute index
   * @return value of criterion for the best split
   * @throws Exception if something goes wrong
   */
  private double findSplitNumeric(int index) throws Exception {

    if (m_Instances.classAttribute().isNominal() || m_Instances.classAttribute().isRanking()) {
      return findSplitNumericNominal(index);
    } else {
      return findSplitNumericNumeric(index);
    }
  }

  /**
   * Finds best split for numeric attribute and nominal class
   * and returns value.
   *
   * @param index attribute index
   * @return value of criterion for the best split
   * @throws Exception if something goes wrong
   */
  private double findSplitNumericNominal(int index) throws Exception {

    double bestVal = Double.MAX_VALUE, currVal, currCutPoint;
    int numMissing = 0;
    double[] sum = new double[m_Instances.numClasses()];
    double[][] bestDist = new double[3][m_Instances.numClasses()];

    // Compute counts for all the values
    for (int i = 0; i < m_Instances.numInstances(); i++) {
      Instance inst = m_Instances.instance(i);
      if (!inst.isMissing(index)) {
	m_Distribution[1][(int)inst.classValue()] += inst.weight();
      } else {
	m_Distribution[2][(int)inst.classValue()] += inst.weight();
	numMissing++;
      }
    }
    System.arraycopy(m_Distribution[1], 0, sum, 0, m_Instances.numClasses());

    // Save current distribution as best distribution
    for (int j = 0; j < 3; j++) {
      System.arraycopy(m_Distribution[j], 0, bestDist[j], 0, 
		       m_Instances.numClasses());
    }

    // Sort instances
    m_Instances.sort(index);
    
    // Make split counts for each possible split and evaluate
    for (int i = 0; i < m_Instances.numInstances() - (numMissing + 1); i++) {
      Instance inst = m_Instances.instance(i);
      Instance instPlusOne = m_Instances.instance(i + 1);
      m_Distribution[0][(int)inst.classValue()] += inst.weight();
      m_Distribution[1][(int)inst.classValue()] -= inst.weight();
      if (inst.value(index) < instPlusOne.value(index)) {
	currCutPoint = (inst.value(index) + instPlusOne.value(index)) / 2.0;
	currVal = ContingencyTables.entropyConditionedOnRows(m_Distribution);
	if (currVal < bestVal) {
	  m_SplitPoint = currCutPoint;
	  bestVal = currVal;
	  for (int j = 0; j < 3; j++) {
	    System.arraycopy(m_Distribution[j], 0, bestDist[j], 0, 
			     m_Instances.numClasses());
	  }
	}
      }
    }

    // No missing values in training data.
    if (numMissing == 0) {
      System.arraycopy(sum, 0, bestDist[2], 0, m_Instances.numClasses());
    }
 
    m_Distribution = bestDist;
    return bestVal;
  }

  /**
   * Finds best split for numeric attribute and numeric class
   * and returns value.
   *
   * @param index attribute index
   * @return value of criterion for the best split
   * @throws Exception if something goes wrong
   */
  private double findSplitNumericNumeric(int index) throws Exception {

    double bestVal = Double.MAX_VALUE, currVal, currCutPoint;
    int numMissing = 0;
    double[] sumsSquares = new double[3], sumOfWeights = new double[3];
    double[][] bestDist = new double[3][1];
    double totalSum = 0, totalSumOfWeights = 0;

    // Compute counts for all the values
    for (int i = 0; i < m_Instances.numInstances(); i++) {
      Instance inst = m_Instances.instance(i);
      if (!inst.isMissing(index)) {
	m_Distribution[1][0] += inst.classValue() * inst.weight();
	sumsSquares[1] += inst.classValue() * inst.classValue() 
	  * inst.weight();
	sumOfWeights[1] += inst.weight();
      } else {
	m_Distribution[2][0] += inst.classValue() * inst.weight();
	sumsSquares[2] += inst.classValue() * inst.classValue() 
	  * inst.weight();
	sumOfWeights[2] += inst.weight();
	numMissing++;
      }
      totalSumOfWeights += inst.weight();
      totalSum += inst.classValue() * inst.weight();
    }

    // Check if the total weight is zero
    if (totalSumOfWeights <= 0) {
      return bestVal;
    }

    // Sort instances
    m_Instances.sort(index);
    
    // Make split counts for each possible split and evaluate
    for (int i = 0; i < m_Instances.numInstances() - (numMissing + 1); i++) {
      Instance inst = m_Instances.instance(i);
      Instance instPlusOne = m_Instances.instance(i + 1);
      m_Distribution[0][0] += inst.classValue() * inst.weight();
      sumsSquares[0] += inst.classValue() * inst.classValue() * inst.weight();
      sumOfWeights[0] += inst.weight();
      m_Distribution[1][0] -= inst.classValue() * inst.weight();
      sumsSquares[1] -= inst.classValue() * inst.classValue() * inst.weight();
      sumOfWeights[1] -= inst.weight();
      if (inst.value(index) < instPlusOne.value(index)) {
	currCutPoint = (inst.value(index) + instPlusOne.value(index)) / 2.0;
	currVal = variance(m_Distribution, sumsSquares, sumOfWeights);
	if (currVal < bestVal) {
	  m_SplitPoint = currCutPoint;
	  bestVal = currVal;
	  for (int j = 0; j < 3; j++) {
	    if (sumOfWeights[j] > 0) {
	      bestDist[j][0] = m_Distribution[j][0] / sumOfWeights[j];
	    } else {
	      bestDist[j][0] = totalSum / totalSumOfWeights;
	    }
	  }
	}
      }
    }

    m_Distribution = bestDist;
    return bestVal;
  }

  /**
   * Computes variance for subsets.
   * 
   * @param s
   * @param sS
   * @param sumOfWeights
   * @return the variance
   */
  private double variance(double[][] s,double[] sS,double[] sumOfWeights) {

    double var = 0;

    for (int i = 0; i < s.length; i++) {
      if (sumOfWeights[i] > 0) {
	var += sS[i] - ((s[i][0] * s[i][0]) / (double) sumOfWeights[i]);
      }
    }
    
    return var;
  }

  /**
   * Returns the subset an instance falls into.
   * 
   * @param instance the instance to check
   * @return the subset the instance falls into
   * @throws Exception if something goes wrong
   */
  private int whichSubset(Instance instance) throws Exception {

    if (instance.isMissing(m_AttIndex)) {
      return 2;
    } else if (instance.attribute(m_AttIndex).isNominal() || instance.attribute(m_AttIndex).isRanking()) {
      if ((int)instance.value(m_AttIndex) == m_SplitPoint) {
	return 0;
      } else {
	return 1;
      }
    } else {
      if (instance.value(m_AttIndex) <= m_SplitPoint) {
	return 0;
      } else {
	return 1;
      }
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5928 $");
  }
 
  /**
   * Main method for testing this class.
   *
   * @param argv the options
   */
  public static void main(String [] argv) {
    runClassifier(new DecisionStump(), argv);
  }
}

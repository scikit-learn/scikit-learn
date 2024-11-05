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
 *    NBTreeSplit.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.classifiers.bayes.NaiveBayesUpdateable;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.filters.Filter;
import weka.filters.supervised.attribute.Discretize;

import java.util.Random;

/**
 * Class implementing a NBTree split on an attribute.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 6088 $
 */
public class NBTreeSplit
  extends ClassifierSplitModel{

  /** for serialization */
  private static final long serialVersionUID = 8922627123884975070L;

  /** Desired number of branches. */
  private int m_complexityIndex;  

  /** Attribute to split on. */
  private int m_attIndex;         

  /** Minimum number of objects in a split.   */
  private int m_minNoObj;         

  /** Value of split point. */
  private double m_splitPoint;   

  /** The sum of the weights of the instances. */
  private double m_sumOfWeights;  

  /** The weight of the instances incorrectly classified by the 
      naive bayes models arising from this split*/
  private double m_errors;

  private C45Split m_c45S;

  /** The global naive bayes model for this node */
  NBTreeNoSplit m_globalNB;

  /**
   * Initializes the split model.
   */
  public NBTreeSplit(int attIndex, int minNoObj, double sumOfWeights) {
    
    // Get index of attribute to split on.
    m_attIndex = attIndex;
        
    // Set minimum number of objects.
    m_minNoObj = minNoObj;

    // Set the sum of the weights
    m_sumOfWeights = sumOfWeights;
    
  }

  /**
   * Creates a NBTree-type split on the given data. Assumes that none of
   * the class values is missing.
   *
   * @exception Exception if something goes wrong
   */
  public void buildClassifier(Instances trainInstances) 
       throws Exception {

    // Initialize the remaining instance variables.
    m_numSubsets = 0;
    m_splitPoint = Double.MAX_VALUE;
    m_errors = 0;
    if (m_globalNB != null) {
      m_errors = m_globalNB.getErrors();
    }

    // Different treatment for enumerated and numeric
    // attributes.
    if (trainInstances.attribute(m_attIndex).isNominal() || trainInstances.attribute(m_attIndex).isRanking()) {
      m_complexityIndex = trainInstances.attribute(m_attIndex).numValues();
      handleEnumeratedAttribute(trainInstances);
    }else{
      m_complexityIndex = 2;
      trainInstances.sort(trainInstances.attribute(m_attIndex));
      handleNumericAttribute(trainInstances);
    }
  }

  /**
   * Returns index of attribute for which split was generated.
   */
  public final int attIndex() {
    
    return m_attIndex;
  }

  /**
   * Creates split on enumerated attribute.
   *
   * @exception Exception if something goes wrong
   */
  private void handleEnumeratedAttribute(Instances trainInstances)
       throws Exception {

    m_c45S = new C45Split(m_attIndex, 2, m_sumOfWeights, true);
    m_c45S.buildClassifier(trainInstances);
    if (m_c45S.numSubsets() == 0) {
      return;
    }
    m_errors = 0;
    Instance instance;

    Instances [] trainingSets = new Instances [m_complexityIndex];
    for (int i = 0; i < m_complexityIndex; i++) {
      trainingSets[i] = new Instances(trainInstances, 0);
    }
    /*    m_distribution = new Distribution(m_complexityIndex,
	  trainInstances.numClasses()); */
    int subset;
    for (int i = 0; i < trainInstances.numInstances(); i++) {
      instance = trainInstances.instance(i);
      subset = m_c45S.whichSubset(instance);
      if (subset > -1) {
	trainingSets[subset].add((Instance)instance.copy());
      } else {
	double [] weights = m_c45S.weights(instance);
	for (int j = 0; j < m_complexityIndex; j++) {
	  try {
	    Instance temp = (Instance) instance.copy();
	    if (weights.length == m_complexityIndex) {
	      temp.setWeight(temp.weight() * weights[j]);
	    } else {
	      temp.setWeight(temp.weight() / m_complexityIndex);
	    }
	    trainingSets[j].add(temp);
	  } catch (Exception ex) {
	    ex.printStackTrace();
	    System.err.println("*** "+m_complexityIndex);
	    System.err.println(weights.length);
	    System.exit(1);
	  }
	}
      }
    }

    /*    // compute weights (weights of instances per subset
    m_weights = new double [m_complexityIndex];
    for (int i = 0; i < m_complexityIndex; i++) {
      m_weights[i] = trainingSets[i].sumOfWeights();
    }
    Utils.normalize(m_weights); */

    /*
    // Only Instances with known values are relevant.
    Enumeration enu = trainInstances.enumerateInstances();
    while (enu.hasMoreElements()) {
      instance = (Instance) enu.nextElement();
      if (!instance.isMissing(m_attIndex)) {
	//	m_distribution.add((int)instance.value(m_attIndex),instance);
	trainingSets[(int)instances.value(m_attIndex)].add(instance);
      } else {
	// add these to the error count
	m_errors += instance.weight();
      }
      } */

    Random r = new Random(1);
    int minNumCount = 0;
    for (int i = 0; i < m_complexityIndex; i++) {
      if (trainingSets[i].numInstances() >= 5) {
	minNumCount++;
	// Discretize the sets
	Discretize disc = new Discretize();
	disc.setInputFormat(trainingSets[i]);
	trainingSets[i] = Filter.useFilter(trainingSets[i], disc);

	trainingSets[i].randomize(r);
	trainingSets[i].stratify(5);
	NaiveBayesUpdateable fullModel = new NaiveBayesUpdateable();
	fullModel.buildClassifier(trainingSets[i]);

	// add the errors for this branch of the split
	m_errors += NBTreeNoSplit.crossValidate(fullModel, trainingSets[i], r);
      } else {
	// if fewer than min obj then just count them as errors
	for (int j = 0; j < trainingSets[i].numInstances(); j++) {
	  m_errors += trainingSets[i].instance(j).weight();
	}
      }
    }
    
    // Check if there are at least five instances in at least two of the subsets
    // subsets.
    if (minNumCount > 1) {
      m_numSubsets = m_complexityIndex;
    }
  }

  /**
   * Creates split on numeric attribute.
   *
   * @exception Exception if something goes wrong
   */
  private void handleNumericAttribute(Instances trainInstances)
       throws Exception {

    m_c45S = new C45Split(m_attIndex, 2, m_sumOfWeights, true);
    m_c45S.buildClassifier(trainInstances);
    if (m_c45S.numSubsets() == 0) {
      return;
    }
    m_errors = 0;

    Instances [] trainingSets = new Instances [m_complexityIndex];
    trainingSets[0] = new Instances(trainInstances, 0);
    trainingSets[1] = new Instances(trainInstances, 0);
    int subset = -1;
    
    // populate the subsets
    for (int i = 0; i < trainInstances.numInstances(); i++) {
      Instance instance = trainInstances.instance(i);
      subset = m_c45S.whichSubset(instance);
      if (subset != -1) {
	trainingSets[subset].add((Instance)instance.copy());
      } else {
	double [] weights = m_c45S.weights(instance);
	for (int j = 0; j < m_complexityIndex; j++) {
	  Instance temp = (Instance)instance.copy();
	  if (weights.length == m_complexityIndex) {
	    temp.setWeight(temp.weight() * weights[j]);
	  } else {
	    temp.setWeight(temp.weight() / m_complexityIndex);
	  }
	  trainingSets[j].add(temp); 
	}
      }
    }
    
    /*    // compute weights (weights of instances per subset
    m_weights = new double [m_complexityIndex];
    for (int i = 0; i < m_complexityIndex; i++) {
      m_weights[i] = trainingSets[i].sumOfWeights();
    }
    Utils.normalize(m_weights); */

    Random r = new Random(1);
    int minNumCount = 0;
    for (int i = 0; i < m_complexityIndex; i++) {
      if (trainingSets[i].numInstances() > 5) {
	minNumCount++;
	// Discretize the sets
		Discretize disc = new Discretize();
	disc.setInputFormat(trainingSets[i]);
	trainingSets[i] = Filter.useFilter(trainingSets[i], disc);

	trainingSets[i].randomize(r);
	trainingSets[i].stratify(5);
	NaiveBayesUpdateable fullModel = new NaiveBayesUpdateable();
	fullModel.buildClassifier(trainingSets[i]);

	// add the errors for this branch of the split
	m_errors += NBTreeNoSplit.crossValidate(fullModel, trainingSets[i], r);
      } else {
	for (int j = 0; j < trainingSets[i].numInstances(); j++) {
	  m_errors += trainingSets[i].instance(j).weight();
	}
      }
    }
    
    // Check if minimum number of Instances in at least two
    // subsets.
    if (minNumCount > 1) {
      m_numSubsets = m_complexityIndex;
    }
  }

  /**
   * Returns index of subset instance is assigned to.
   * Returns -1 if instance is assigned to more than one subset.
   *
   * @exception Exception if something goes wrong
   */
  public final int whichSubset(Instance instance) 
    throws Exception {
    
    return m_c45S.whichSubset(instance);
  }

  /**
   * Returns weights if instance is assigned to more than one subset.
   * Returns null if instance is only assigned to one subset.
   */
  public final double [] weights(Instance instance) {
    return m_c45S.weights(instance);
    //     return m_weights;
  }

  /**
   * Returns a string containing java source code equivalent to the test
   * made at this node. The instance being tested is called "i".
   *
   * @param index index of the nominal value tested
   * @param data the data containing instance structure info
   * @return a value of type 'String'
   */
  public final String sourceExpression(int index, Instances data) {
    return m_c45S.sourceExpression(index, data);
  }

  /**
   * Prints the condition satisfied by instances in a subset.
   *
   * @param index of subset 
   * @param data training set.
   */
  public final String rightSide(int index,Instances data) {
    return m_c45S.rightSide(index, data);
  }

  /**
   * Prints left side of condition..
   *
   * @param data training set.
   */
  public final String leftSide(Instances data) {

    return m_c45S.leftSide(data);
  }

  /**
   * Return the probability for a class value
   *
   * @param classIndex the index of the class value
   * @param instance the instance to generate a probability for
   * @param theSubset the subset to consider
   * @return a probability
   * @exception Exception if an error occurs
   */
  public double classProb(int classIndex, Instance instance, int theSubset) 
    throws Exception {

    // use the global naive bayes model
    if (theSubset > -1) {
      return m_globalNB.classProb(classIndex, instance, theSubset);
    } else {
      throw new Exception("This shouldn't happen!!!");
    }
  }

  /**
   * Return the global naive bayes model for this node
   *
   * @return a <code>NBTreeNoSplit</code> value
   */
  public NBTreeNoSplit getGlobalModel() {
    return m_globalNB;
  }

  /**
   * Set the global naive bayes model for this node
   *
   * @param global a <code>NBTreeNoSplit</code> value
   */
  public void setGlobalModel(NBTreeNoSplit global) {
    m_globalNB = global;
  }

  /**
   * Return the errors made by the naive bayes models arising
   * from this split.
   *
   * @return a <code>double</code> value
   */
  public double getErrors() {
    return m_errors;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6088 $");
  }
}

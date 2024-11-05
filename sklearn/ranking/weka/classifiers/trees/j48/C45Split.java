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
 *    C45Split.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;

/**
 * Class implementing a C4.5-type split on an attribute.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6073 $
 */
public class C45Split
  extends ClassifierSplitModel{

  /** for serialization */
  private static final long serialVersionUID = 3064079330067903161L;

  /** Desired number of branches. */
  private int m_complexityIndex;  

  /** Attribute to split on. */
  private int m_attIndex;         

  /** Minimum number of objects in a split.   */
  private int m_minNoObj;         

  /** Use MDL correction? */
  private boolean m_useMDLcorrection;         

  /** Value of split point. */
  private double m_splitPoint;   

  /** InfoGain of split. */ 
  private double m_infoGain; 

  /** GainRatio of split.  */
  private double m_gainRatio;  

  /** The sum of the weights of the instances. */
  private double m_sumOfWeights;  

  /** Number of split points. */
  private int m_index;            

  /** Static reference to splitting criterion. */
  private static InfoGainSplitCrit infoGainCrit = new InfoGainSplitCrit();

  /** Static reference to splitting criterion. */
  private static GainRatioSplitCrit gainRatioCrit = new GainRatioSplitCrit();

  /**
   * Initializes the split model.
   */
  public C45Split(int attIndex,int minNoObj, double sumOfWeights,
                  boolean useMDLcorrection) {

    // Get index of attribute to split on.
    m_attIndex = attIndex;
        
    // Set minimum number of objects.
    m_minNoObj = minNoObj;

    // Set the sum of the weights
    m_sumOfWeights = sumOfWeights;

    // Whether to use the MDL correction for numeric attributes
    m_useMDLcorrection = useMDLcorrection;
  }

  /**
   * Creates a C4.5-type split on the given data. Assumes that none of
   * the class values is missing.
   *
   * @exception Exception if something goes wrong
   */
  public void buildClassifier(Instances trainInstances) 
       throws Exception {

    // Initialize the remaining instance variables.
    m_numSubsets = 0;
    m_splitPoint = Double.MAX_VALUE;
    m_infoGain = 0;
    m_gainRatio = 0;

    // Different treatment for enumerated and numeric
    // attributes.
    if (trainInstances.attribute(m_attIndex).isNominal() || trainInstances.attribute(m_attIndex).isRanking()) {
      m_complexityIndex = trainInstances.attribute(m_attIndex).numValues();
      m_index = m_complexityIndex;
      handleEnumeratedAttribute(trainInstances);
    }else{
      m_complexityIndex = 2;
      m_index = 0;
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
   * Returns the split point (numeric attribute only).
   * 
   * @return the split point used for a test on a numeric attribute
   */
  public double splitPoint() {
    return m_splitPoint;
  }

  /**
   * Gets class probability for instance.
   *
   * @exception Exception if something goes wrong
   */
  public final double classProb(int classIndex,Instance instance,
				int theSubset) throws Exception {

    if (theSubset <= -1) {
      double [] weights = weights(instance);
      if (weights == null) {
	return m_distribution.prob(classIndex);
      } else {
	double prob = 0;
	for (int i = 0; i < weights.length; i++) {
	  prob += weights[i] * m_distribution.prob(classIndex, i);
	}
	return prob;
      }
    } else {
      if (Utils.gr(m_distribution.perBag(theSubset), 0)) {
	return m_distribution.prob(classIndex, theSubset);
      } else {
	return m_distribution.prob(classIndex);
      }
    }
  }
 
  /**
   * Returns coding cost for split (used in rule learner).
   */
  public final double codingCost() {

    return Utils.log2(m_index);
  }
 
  /**
   * Returns (C4.5-type) gain ratio for the generated split.
   */
  public final double gainRatio() {
    return m_gainRatio;
  }

  /**
   * Creates split on enumerated attribute.
   *
   * @exception Exception if something goes wrong
   */
  private void handleEnumeratedAttribute(Instances trainInstances)
       throws Exception {
    
    Instance instance;

    m_distribution = new Distribution(m_complexityIndex,
			      trainInstances.numClasses());
    
    // Only Instances with known values are relevant.
    Enumeration enu = trainInstances.enumerateInstances();
    while (enu.hasMoreElements()) {
      instance = (Instance) enu.nextElement();
      if (!instance.isMissing(m_attIndex))
	m_distribution.add((int)instance.value(m_attIndex),instance);
    }
    
    // Check if minimum number of Instances in at least two
    // subsets.
    if (m_distribution.check(m_minNoObj)) {
      m_numSubsets = m_complexityIndex;
      m_infoGain = infoGainCrit.
	splitCritValue(m_distribution,m_sumOfWeights);
      m_gainRatio = 
	gainRatioCrit.splitCritValue(m_distribution,m_sumOfWeights,
				     m_infoGain);
    }
  }
  
  /**
   * Creates split on numeric attribute.
   *
   * @exception Exception if something goes wrong
   */
  private void handleNumericAttribute(Instances trainInstances)
       throws Exception {
  
    int firstMiss;
    int next = 1;
    int last = 0;
    int splitIndex = -1;
    double currentInfoGain;
    double defaultEnt;
    double minSplit;
    Instance instance;
    int i;

    // Current attribute is a numeric attribute.
    m_distribution = new Distribution(2,trainInstances.numClasses());
    
    // Only Instances with known values are relevant.
    Enumeration enu = trainInstances.enumerateInstances();
    i = 0;
    while (enu.hasMoreElements()) {
      instance = (Instance) enu.nextElement();
      if (instance.isMissing(m_attIndex))
	break;
      m_distribution.add(1,instance);
      i++;
    }
    firstMiss = i;
	
    // Compute minimum number of Instances required in each
    // subset.
    minSplit =  0.1*(m_distribution.total())/
      ((double)trainInstances.numClasses());
    if (Utils.smOrEq(minSplit,m_minNoObj)) 
      minSplit = m_minNoObj;
    else
      if (Utils.gr(minSplit,25)) 
	minSplit = 25;
	
    // Enough Instances with known values?
    if (Utils.sm((double)firstMiss,2*minSplit))
      return;
    
    // Compute values of criteria for all possible split
    // indices.
    defaultEnt = infoGainCrit.oldEnt(m_distribution);
    while (next < firstMiss) {
	  
      if (trainInstances.instance(next-1).value(m_attIndex)+1e-5 < 
	  trainInstances.instance(next).value(m_attIndex)) { 
	
	// Move class values for all Instances up to next 
	// possible split point.
	m_distribution.shiftRange(1,0,trainInstances,last,next);
	
	// Check if enough Instances in each subset and compute
	// values for criteria.
	if (Utils.grOrEq(m_distribution.perBag(0),minSplit) &&
	    Utils.grOrEq(m_distribution.perBag(1),minSplit)) {
	  currentInfoGain = infoGainCrit.
	    splitCritValue(m_distribution,m_sumOfWeights,
			   defaultEnt);
	  if (Utils.gr(currentInfoGain,m_infoGain)) {
	    m_infoGain = currentInfoGain;
	    splitIndex = next-1;
	  }
	  m_index++;
	}
	last = next;
      }
      next++;
    }
    
    // Was there any useful split?
    if (m_index == 0)
      return;
    
    // Compute modified information gain for best split.
    if (m_useMDLcorrection) {
      m_infoGain = m_infoGain-(Utils.log2(m_index)/m_sumOfWeights);
    }
    if (Utils.smOrEq(m_infoGain,0))
      return;
    
    // Set instance variables' values to values for
    // best split.
    m_numSubsets = 2;
    m_splitPoint = 
      (trainInstances.instance(splitIndex+1).value(m_attIndex)+
       trainInstances.instance(splitIndex).value(m_attIndex))/2;

    // In case we have a numerical precision problem we need to choose the
    // smaller value
    if (m_splitPoint == trainInstances.instance(splitIndex + 1).value(m_attIndex)) {
      m_splitPoint = trainInstances.instance(splitIndex).value(m_attIndex);
    }

    // Restore distributioN for best split.
    m_distribution = new Distribution(2,trainInstances.numClasses());
    m_distribution.addRange(0,trainInstances,0,splitIndex+1);
    m_distribution.addRange(1,trainInstances,splitIndex+1,firstMiss);

    // Compute modified gain ratio for best split.
    m_gainRatio = gainRatioCrit.
      splitCritValue(m_distribution,m_sumOfWeights,
		     m_infoGain);
  }

  /**
   * Returns (C4.5-type) information gain for the generated split.
   */
  public final double infoGain() {

    return m_infoGain;
  }

  /**
   * Prints left side of condition..
   *
   * @param data training set.
   */
  public final String leftSide(Instances data) {

    return data.attribute(m_attIndex).name();
  }

  /**
   * Prints the condition satisfied by instances in a subset.
   *
   * @param index of subset 
   * @param data training set.
   */
  public final String rightSide(int index,Instances data) {

    StringBuffer text;

    text = new StringBuffer();
    if (data.attribute(m_attIndex).isNominal() || data.attribute(m_attIndex).isRanking())
      text.append(" = "+
		  data.attribute(m_attIndex).value(index));
    else
      if (index == 0)
	text.append(" <= "+
		    Utils.doubleToString(m_splitPoint,6));
      else
	text.append(" > "+
		    Utils.doubleToString(m_splitPoint,6));
    return text.toString();
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

    StringBuffer expr = null;
    if (index < 0) {
      return "i[" + m_attIndex + "] == null";
    }
    if (data.attribute(m_attIndex).isNominal() || data.attribute(m_attIndex).isRanking()) {
      expr = new StringBuffer("i[");
      expr.append(m_attIndex).append("]");
      expr.append(".equals(\"").append(data.attribute(m_attIndex)
				     .value(index)).append("\")");
    } else {
      expr = new StringBuffer("((Double) i[");
      expr.append(m_attIndex).append("])");
      if (index == 0) {
	expr.append(".doubleValue() <= ").append(m_splitPoint);
      } else {
	expr.append(".doubleValue() > ").append(m_splitPoint);
      }
    }
    return expr.toString();
  }  

  /**
   * Sets split point to greatest value in given data smaller or equal to
   * old split point.
   * (C4.5 does this for some strange reason).
   */
  public final void setSplitPoint(Instances allInstances) {
    
    double newSplitPoint = -Double.MAX_VALUE;
    double tempValue;
    Instance instance;
    
    if ((allInstances.attribute(m_attIndex).isNumeric()) &&
	(m_numSubsets > 1)) {
      Enumeration enu = allInstances.enumerateInstances();
      while (enu.hasMoreElements()) {
	instance = (Instance) enu.nextElement();
	if (!instance.isMissing(m_attIndex)) {
	  tempValue = instance.value(m_attIndex);
	  if (Utils.gr(tempValue,newSplitPoint) && 
	      Utils.smOrEq(tempValue,m_splitPoint))
	    newSplitPoint = tempValue;
	}
      }
      m_splitPoint = newSplitPoint;
    }
  }
  
  /**
   * Returns the minsAndMaxs of the index.th subset.
   */
  public final double [][] minsAndMaxs(Instances data, double [][] minsAndMaxs,
				       int index) {

    double [][] newMinsAndMaxs = new double[data.numAttributes()][2];
    
    for (int i = 0; i < data.numAttributes(); i++) {
      newMinsAndMaxs[i][0] = minsAndMaxs[i][0];
      newMinsAndMaxs[i][1] = minsAndMaxs[i][1];
      if (i == m_attIndex)
	if (data.attribute(m_attIndex).isNominal() || data.attribute(m_attIndex).isRanking())
	  newMinsAndMaxs[m_attIndex][1] = 1;
	else
	  newMinsAndMaxs[m_attIndex][1-index] = m_splitPoint;
    }

    return newMinsAndMaxs;
  }
  
  /**
   * Sets distribution associated with model.
   */
  public void resetDistribution(Instances data) throws Exception {
    
    Instances insts = new Instances(data, data.numInstances());
    for (int i = 0; i < data.numInstances(); i++) {
      if (whichSubset(data.instance(i)) > -1) {
	insts.add(data.instance(i));
      }
    }
    Distribution newD = new Distribution(insts, this);
    newD.addInstWithUnknown(data, m_attIndex);
    m_distribution = newD;
  }

  /**
   * Returns weights if instance is assigned to more than one subset.
   * Returns null if instance is only assigned to one subset.
   */
  public final double [] weights(Instance instance) {
    
    double [] weights;
    int i;
    
    if (instance.isMissing(m_attIndex)) {
      weights = new double [m_numSubsets];
      for (i=0;i<m_numSubsets;i++)
	weights [i] = m_distribution.perBag(i)/m_distribution.total();
      return weights;
    }else{
      return null;
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
    
    if (instance.isMissing(m_attIndex))
      return -1;
    else{
      if (instance.attribute(m_attIndex).isNominal() || instance.attribute(m_attIndex).isRanking())
	return (int)instance.value(m_attIndex);
      else
	if (Utils.smOrEq(instance.value(m_attIndex),m_splitPoint))
	  return 0;
	else
	  return 1;
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6073 $");
  }
}

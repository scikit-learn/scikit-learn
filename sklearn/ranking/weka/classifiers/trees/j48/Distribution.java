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
 *    Distribution.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.trees.j48;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.Serializable;
import java.util.Enumeration;

/**
 * Class for handling a distribution of class values.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.12 $
 */
public class Distribution
  implements Cloneable, Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = 8526859638230806576L;

  /** Weight of instances per class per bag. */
  private double m_perClassPerBag[][]; 

  /** Weight of instances per bag. */
  private double m_perBag[];           

  /** Weight of instances per class. */
  private double m_perClass[];         

  /** Total weight of instances. */
  private double totaL;            

  /**
   * Creates and initializes a new distribution.
   */
  public Distribution(int numBags,int numClasses) {

    int i;

    m_perClassPerBag = new double [numBags][0];
    m_perBag = new double [numBags];
    m_perClass = new double [numClasses];
    for (i=0;i<numBags;i++)
      m_perClassPerBag[i] = new double [numClasses];
    totaL = 0;
  }

  /**
   * Creates and initializes a new distribution using the given
   * array. WARNING: it just copies a reference to this array.
   */
  public Distribution(double [][] table) {

    int i, j;

    m_perClassPerBag = table;
    m_perBag = new double [table.length];
    m_perClass = new double [table[0].length];
    for (i = 0; i < table.length; i++) 
      for (j  = 0; j < table[i].length; j++) {
	m_perBag[i] += table[i][j];
	m_perClass[j] += table[i][j];
	totaL += table[i][j];
      }
  }

  /**
   * Creates a distribution with only one bag according
   * to instances in source.
   *
   * @exception Exception if something goes wrong
   */
  public Distribution(Instances source) throws Exception {
    
    m_perClassPerBag = new double [1][0];
    m_perBag = new double [1];
    totaL = 0;
    m_perClass = new double [source.numClasses()];
    m_perClassPerBag[0] = new double [source.numClasses()];
    Enumeration enu = source.enumerateInstances();
    while (enu.hasMoreElements())
      add(0,(Instance) enu.nextElement());
  }

  /**
   * Creates a distribution according to given instances and
   * split model.
   *
   * @exception Exception if something goes wrong
   */

  public Distribution(Instances source, 
		      ClassifierSplitModel modelToUse)
       throws Exception {

    int index;
    Instance instance;
    double[] weights;

    m_perClassPerBag = new double [modelToUse.numSubsets()][0];
    m_perBag = new double [modelToUse.numSubsets()];
    totaL = 0;
    m_perClass = new double [source.numClasses()];
    for (int i = 0; i < modelToUse.numSubsets(); i++)
      m_perClassPerBag[i] = new double [source.numClasses()];
    Enumeration enu = source.enumerateInstances();
    while (enu.hasMoreElements()) {
      instance = (Instance) enu.nextElement();
      index = modelToUse.whichSubset(instance);
      if (index != -1)
	add(index, instance);
      else {
	weights = modelToUse.weights(instance);
	addWeights(instance, weights);
      }
    }
  }

  /**
   * Creates distribution with only one bag by merging all
   * bags of given distribution.
   */
  public Distribution(Distribution toMerge) {

    totaL = toMerge.totaL;
    m_perClass = new double [toMerge.numClasses()];
    System.arraycopy(toMerge.m_perClass,0,m_perClass,0,toMerge.numClasses());
    m_perClassPerBag = new double [1] [0];
    m_perClassPerBag[0] = new double [toMerge.numClasses()];
    System.arraycopy(toMerge.m_perClass,0,m_perClassPerBag[0],0,
		     toMerge.numClasses());
    m_perBag = new double [1];
    m_perBag[0] = totaL;
  }

  /**
   * Creates distribution with two bags by merging all bags apart of
   * the indicated one.
   */
  public Distribution(Distribution toMerge, int index) {

    int i;

    totaL = toMerge.totaL;
    m_perClass = new double [toMerge.numClasses()];
    System.arraycopy(toMerge.m_perClass,0,m_perClass,0,toMerge.numClasses());
    m_perClassPerBag = new double [2] [0];
    m_perClassPerBag[0] = new double [toMerge.numClasses()];
    System.arraycopy(toMerge.m_perClassPerBag[index],0,m_perClassPerBag[0],0,
		     toMerge.numClasses());
    m_perClassPerBag[1] = new double [toMerge.numClasses()];
    for (i=0;i<toMerge.numClasses();i++)
      m_perClassPerBag[1][i] = toMerge.m_perClass[i]-m_perClassPerBag[0][i];
    m_perBag = new double [2];
    m_perBag[0] = toMerge.m_perBag[index];
    m_perBag[1] = totaL-m_perBag[0];
  }
  
  /**
   * Returns number of non-empty bags of distribution.
   */
  public final int actualNumBags() {
    
    int returnValue = 0;
    int i;

    for (i=0;i<m_perBag.length;i++)
      if (Utils.gr(m_perBag[i],0))
	returnValue++;
    
    return returnValue;
  }

  /**
   * Returns number of classes actually occuring in distribution.
   */
  public final int actualNumClasses() {

    int returnValue = 0;
    int i;

    for (i=0;i<m_perClass.length;i++)
      if (Utils.gr(m_perClass[i],0))
	returnValue++;
    
    return returnValue;
  }

  /**
   * Returns number of classes actually occuring in given bag.
   */
  public final int actualNumClasses(int bagIndex) {

    int returnValue = 0;
    int i;

    for (i=0;i<m_perClass.length;i++)
      if (Utils.gr(m_perClassPerBag[bagIndex][i],0))
	returnValue++;
    
    return returnValue;
  }

  /**
   * Adds given instance to given bag.
   *
   * @exception Exception if something goes wrong
   */
  public final void add(int bagIndex,Instance instance) 
       throws Exception {
    
    int classIndex;
    double weight;

    classIndex = (int)instance.classValue();
    weight = instance.weight();
    m_perClassPerBag[bagIndex][classIndex] = 
      m_perClassPerBag[bagIndex][classIndex]+weight;
    m_perBag[bagIndex] = m_perBag[bagIndex]+weight;
    m_perClass[classIndex] = m_perClass[classIndex]+weight;
    totaL = totaL+weight;
  }

  /**
   * Subtracts given instance from given bag.
   *
   * @exception Exception if something goes wrong
   */
  public final void sub(int bagIndex,Instance instance) 
       throws Exception {
    
    int classIndex;
    double weight;

    classIndex = (int)instance.classValue();
    weight = instance.weight();
    m_perClassPerBag[bagIndex][classIndex] = 
      m_perClassPerBag[bagIndex][classIndex]-weight;
    m_perBag[bagIndex] = m_perBag[bagIndex]-weight;
    m_perClass[classIndex] = m_perClass[classIndex]-weight;
    totaL = totaL-weight;
  }

  /**
   * Adds counts to given bag.
   */
  public final void add(int bagIndex, double[] counts) {
    
    double sum = Utils.sum(counts);

    for (int i = 0; i < counts.length; i++)
      m_perClassPerBag[bagIndex][i] += counts[i];
    m_perBag[bagIndex] = m_perBag[bagIndex]+sum;
    for (int i = 0; i < counts.length; i++)
      m_perClass[i] = m_perClass[i]+counts[i];
    totaL = totaL+sum;
  }

  /**
   * Adds all instances with unknown values for given attribute, weighted
   * according to frequency of instances in each bag.
   *
   * @exception Exception if something goes wrong
   */
  public final void addInstWithUnknown(Instances source,
				       int attIndex)
       throws Exception {

    double [] probs;
    double weight,newWeight;
    int classIndex;
    Instance instance;
    int j;

    probs = new double [m_perBag.length];
    for (j=0;j<m_perBag.length;j++) {
      if (Utils.eq(totaL, 0)) {
	probs[j] = 1.0 / probs.length;
      } else {
	probs[j] = m_perBag[j]/totaL;
      }
    }
    Enumeration enu = source.enumerateInstances();
    while (enu.hasMoreElements()) {
      instance = (Instance) enu.nextElement();
      if (instance.isMissing(attIndex)) {
	classIndex = (int)instance.classValue();
	weight = instance.weight();
	m_perClass[classIndex] = m_perClass[classIndex]+weight;
	totaL = totaL+weight;
	for (j = 0; j < m_perBag.length; j++) {
	  newWeight = probs[j]*weight;
	  m_perClassPerBag[j][classIndex] = m_perClassPerBag[j][classIndex]+
	    newWeight;
	  m_perBag[j] = m_perBag[j]+newWeight;
	}
      }
    }
  }

  /**
   * Adds all instances in given range to given bag.
   *
   * @exception Exception if something goes wrong
   */
  public final void addRange(int bagIndex,Instances source,
			     int startIndex, int lastPlusOne)
       throws Exception {

    double sumOfWeights = 0;
    int classIndex;
    Instance instance;
    int i;

    for (i = startIndex; i < lastPlusOne; i++) {
      instance = (Instance) source.instance(i);
      classIndex = (int)instance.classValue();
      sumOfWeights = sumOfWeights+instance.weight();
      m_perClassPerBag[bagIndex][classIndex] += instance.weight();
      m_perClass[classIndex] += instance.weight();
    }
    m_perBag[bagIndex] += sumOfWeights;
    totaL += sumOfWeights;
  }

  /**
   * Adds given instance to all bags weighting it according to given weights.
   *
   * @exception Exception if something goes wrong
   */
  public final void addWeights(Instance instance, 
			       double [] weights)
       throws Exception {

    int classIndex;
    int i;

    classIndex = (int)instance.classValue();
    for (i=0;i<m_perBag.length;i++) {
      double weight = instance.weight() * weights[i];
      m_perClassPerBag[i][classIndex] = m_perClassPerBag[i][classIndex] + weight;
      m_perBag[i] = m_perBag[i] + weight;
      m_perClass[classIndex] = m_perClass[classIndex] + weight;
      totaL = totaL + weight;
    }
  }

  /**
   * Checks if at least two bags contain a minimum number of instances.
   */
  public final boolean check(double minNoObj) {

    int counter = 0;
    int i;

    for (i=0;i<m_perBag.length;i++)
      if (Utils.grOrEq(m_perBag[i],minNoObj))
	counter++;
    if (counter > 1)
      return true;
    else
      return false;
  }

  /**
   * Clones distribution (Deep copy of distribution).
   */
  public final Object clone() {

    int i,j;

    Distribution newDistribution = new Distribution (m_perBag.length,
						     m_perClass.length);
    for (i=0;i<m_perBag.length;i++) {
      newDistribution.m_perBag[i] = m_perBag[i];
      for (j=0;j<m_perClass.length;j++)
	newDistribution.m_perClassPerBag[i][j] = m_perClassPerBag[i][j];
    }
    for (j=0;j<m_perClass.length;j++)
      newDistribution.m_perClass[j] = m_perClass[j];
    newDistribution.totaL = totaL;
  
    return newDistribution;
  }

  /**
   * Deletes given instance from given bag.
   *
   * @exception Exception if something goes wrong
   */
  public final void del(int bagIndex,Instance instance) 
       throws Exception {

    int classIndex;
    double weight;

    classIndex = (int)instance.classValue();
    weight = instance.weight();
    m_perClassPerBag[bagIndex][classIndex] = 
      m_perClassPerBag[bagIndex][classIndex]-weight;
    m_perBag[bagIndex] = m_perBag[bagIndex]-weight;
    m_perClass[classIndex] = m_perClass[classIndex]-weight;
    totaL = totaL-weight;
  }

  /**
   * Deletes all instances in given range from given bag.
   *
   * @exception Exception if something goes wrong
   */
  public final void delRange(int bagIndex,Instances source,
			     int startIndex, int lastPlusOne)
       throws Exception {

    double sumOfWeights = 0;
    int classIndex;
    Instance instance;
    int i;

    for (i = startIndex; i < lastPlusOne; i++) {
      instance = (Instance) source.instance(i);
      classIndex = (int)instance.classValue();
      sumOfWeights = sumOfWeights+instance.weight();
      m_perClassPerBag[bagIndex][classIndex] -= instance.weight();
      m_perClass[classIndex] -= instance.weight();
    }
    m_perBag[bagIndex] -= sumOfWeights;
    totaL -= sumOfWeights;
  }

  /**
   * Prints distribution.
   */
  
  public final String dumpDistribution() {

    StringBuffer text;
    int i,j;

    text = new StringBuffer();
    for (i=0;i<m_perBag.length;i++) {
      text.append("Bag num "+i+"\n");
      for (j=0;j<m_perClass.length;j++)
	text.append("Class num "+j+" "+m_perClassPerBag[i][j]+"\n");
    }
    return text.toString();
  }

  /**
   * Sets all counts to zero.
   */
  public final void initialize() {

    for (int i = 0; i < m_perClass.length; i++) 
      m_perClass[i] = 0;
    for (int i = 0; i < m_perBag.length; i++)
      m_perBag[i] = 0;
    for (int i = 0; i < m_perBag.length; i++)
      for (int j = 0; j < m_perClass.length; j++)
	m_perClassPerBag[i][j] = 0;
    totaL = 0;
  }

  /**
   * Returns matrix with distribution of class values.
   */
  public final double[][] matrix() {

    return m_perClassPerBag;
  }
  
  /**
   * Returns index of bag containing maximum number of instances.
   */
  public final int maxBag() {

    double max;
    int maxIndex;
    int i;
    
    max = 0;
    maxIndex = -1;
    for (i=0;i<m_perBag.length;i++)
      if (Utils.grOrEq(m_perBag[i],max)) {
	max = m_perBag[i];
	maxIndex = i;
      }
    return maxIndex;
  }

  /**
   * Returns class with highest frequency over all bags.
   */
  public final int maxClass() {

    double maxCount = 0;
    int maxIndex = 0;
    int i;

    for (i=0;i<m_perClass.length;i++)
      if (Utils.gr(m_perClass[i],maxCount)) {
	maxCount = m_perClass[i];
	maxIndex = i;
      }

    return maxIndex;
  }

  /**
   * Returns class with highest frequency for given bag.
   */
  public final int maxClass(int index) {

    double maxCount = 0;
    int maxIndex = 0;
    int i;

    if (Utils.gr(m_perBag[index],0)) {
      for (i=0;i<m_perClass.length;i++)
	if (Utils.gr(m_perClassPerBag[index][i],maxCount)) {
	  maxCount = m_perClassPerBag[index][i];
	  maxIndex = i;
	}
      return maxIndex;
    }else
      return maxClass();
  }

  /**
   * Returns number of bags.
   */
  public final int numBags() {
    
    return m_perBag.length;
  }

  /**
   * Returns number of classes.
   */
  public final int numClasses() {

    return m_perClass.length;
  }

  /**
   * Returns perClass(maxClass()).
   */
  public final double numCorrect() {

    return m_perClass[maxClass()];
  }

  /**
   * Returns perClassPerBag(index,maxClass(index)).
   */
  public final double numCorrect(int index) {

    return m_perClassPerBag[index][maxClass(index)];
  }

  /**
   * Returns total-numCorrect().
   */
  public final double numIncorrect() {

    return totaL-numCorrect();
  }

  /**
   * Returns perBag(index)-numCorrect(index).
   */
  public final double numIncorrect(int index) {

    return m_perBag[index]-numCorrect(index);
  }

  /**
   * Returns number of (possibly fractional) instances of given class in 
   * given bag.
   */
  public final double perClassPerBag(int bagIndex, int classIndex) {

    return m_perClassPerBag[bagIndex][classIndex];
  }

  /**
   * Returns number of (possibly fractional) instances in given bag.
   */
  public final double perBag(int bagIndex) {

    return m_perBag[bagIndex];
  }

  /**
   * Returns number of (possibly fractional) instances of given class.
   */
  public final double perClass(int classIndex) {

    return m_perClass[classIndex];
  }

  /**
   * Returns relative frequency of class over all bags with
   * Laplace correction.
   */
  public final double laplaceProb(int classIndex) {

    return (m_perClass[classIndex] + 1) / 
      (totaL + (double) m_perClass.length);
  }

  /**
   * Returns relative frequency of class for given bag.
   */
  public final double laplaceProb(int classIndex, int intIndex) {

	  if (Utils.gr(m_perBag[intIndex],0))
		return (m_perClassPerBag[intIndex][classIndex] + 1.0) /
	           (m_perBag[intIndex] + (double) m_perClass.length);
	  else
	    return laplaceProb(classIndex);
	  
  }

  /**
   * Returns relative frequency of class over all bags.
   */
  public final double prob(int classIndex) {

    if (!Utils.eq(totaL, 0)) {
      return m_perClass[classIndex]/totaL;
    } else {
      return 0;
    }
  }

  /**
   * Returns relative frequency of class for given bag.
   */
  public final double prob(int classIndex,int intIndex) {

    if (Utils.gr(m_perBag[intIndex],0))
      return m_perClassPerBag[intIndex][classIndex]/m_perBag[intIndex];
    else
      return prob(classIndex);
  }

  /** 
   * Subtracts the given distribution from this one. The results
   * has only one bag.
   */
  public final Distribution subtract(Distribution toSubstract) {

    Distribution newDist = new Distribution(1,m_perClass.length);

    newDist.m_perBag[0] = totaL-toSubstract.totaL;
    newDist.totaL = newDist.m_perBag[0];
    for (int i = 0; i < m_perClass.length; i++) {
      newDist.m_perClassPerBag[0][i] = m_perClass[i] - toSubstract.m_perClass[i];
      newDist.m_perClass[i] = newDist.m_perClassPerBag[0][i];
    }
    return newDist;
  }

  /**
   * Returns total number of (possibly fractional) instances.
   */
  public final double total() {

    return totaL;
  }

  /**
   * Shifts given instance from one bag to another one.
   *
   * @exception Exception if something goes wrong
   */
  public final void shift(int from,int to,Instance instance) 
       throws Exception {
    
    int classIndex;
    double weight;

    classIndex = (int)instance.classValue();
    weight = instance.weight();
    m_perClassPerBag[from][classIndex] -= weight;
    m_perClassPerBag[to][classIndex] += weight;
    m_perBag[from] -= weight;
    m_perBag[to] += weight;
  }

  /**
   * Shifts all instances in given range from one bag to another one.
   *
   * @exception Exception if something goes wrong
   */
  public final void shiftRange(int from,int to,Instances source,
			       int startIndex,int lastPlusOne) 
       throws Exception {
    
    int classIndex;
    double weight;
    Instance instance;
    int i;

    for (i = startIndex; i < lastPlusOne; i++) {
      instance = (Instance) source.instance(i);
      classIndex = (int)instance.classValue();
      weight = instance.weight();
      m_perClassPerBag[from][classIndex] -= weight;
      m_perClassPerBag[to][classIndex] += weight;
      m_perBag[from] -= weight;
      m_perBag[to] += weight;
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.12 $");
  }
}

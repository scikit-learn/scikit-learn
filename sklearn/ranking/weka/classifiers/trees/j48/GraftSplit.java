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
 *  GraftSplit.java
 *  Copyright (C) 2007 Geoff Webb & Janice Boughton
 *  a split object for nodes added to a tree during grafting.
 *  (used in classifier J48g).
 */

package weka.classifiers.trees.j48;

import weka.core.*;

/**
 * Class implementing a split for nodes added to a tree during grafting.
 *
 * @author Janice Boughton (jrbought@infotech.monash.edu.au)
 * @version $Revision 1.0 $
 */
public class GraftSplit
  extends ClassifierSplitModel
  implements Comparable {

  /** for serialzation. */
  private static final long serialVersionUID = 722773260393182051L;

  /** the distribution for graft values, from cases in atbop */
  private Distribution m_graftdistro;
	
  /** the attribute we are splitting on */
  private int m_attIndex;

  /** value of split point (if numeric attribute) */
  private double m_splitPoint;

  /** dominant class of the subset specified by m_testType */
  private int m_maxClass;

  /** dominant class of the subset not specified by m_testType */
  private int m_otherLeafMaxClass;

  /** laplace value of the subset specified by m_testType for m_maxClass */
  private double m_laplace;

  /** leaf for the subset specified by m_testType */
  private Distribution m_leafdistro;

  /** 
   * type of test:
   * 0: <= test
   * 1: > test
   * 2: = test
   * 3: != test
   */
  private int m_testType;


  /**
   * constructor
   *
   * @param a the attribute to split on
   * @param v the value of a where split occurs
   * @param t the test type (0 is <=, 1 is >, 2 is =, 3 is !)
   * @param c the class to label the leaf node pointed to by test as.
   * @param l the laplace value (needed when sorting GraftSplits)
   */
  public GraftSplit(int a, double v, int t, double c, double l) {

    m_attIndex = a;
    m_splitPoint = v;
    m_testType = t;
    m_maxClass = (int)c;
    m_laplace = l;
  }


  /**
   * constructor
   *
   * @param a the attribute to split on
   * @param v the value of a where split occurs
   * @param t the test type (0 is <=, 1 is >, 2 is =, 3 is !=)
   * @param oC the class to label the leaf node not pointed to by test as.
   * @param counts the distribution for this split
   */
  public GraftSplit(int a, double v, int t, double oC, double [][] counts)
                                                           throws Exception {
    m_attIndex = a;
    m_splitPoint = v;
    m_testType = t;
    m_otherLeafMaxClass = (int)oC;

    // only deal with binary cuts (<= and >; = and !=)
    m_numSubsets = 2;

    // which subset are we looking at for the graft?
    int subset = subsetOfInterest();  // this is the subset for m_leaf

    // create graft distribution, based on counts
    m_distribution = new Distribution(counts);

    // create a distribution object for m_leaf
    double [][] lcounts = new double[1][m_distribution.numClasses()];
    for(int c = 0; c < lcounts[0].length; c++) {
       lcounts[0][c] = counts[subset][c];
    }
    m_leafdistro = new Distribution(lcounts);

    // set the max class
    m_maxClass = m_distribution.maxClass(subset);
 
    // set the laplace value (assumes binary class) for subset of interest
    m_laplace = (m_distribution.perClassPerBag(subset, m_maxClass) + 1.0) 
               / (m_distribution.perBag(subset) + 2.0);
  }


  /**
   * deletes the cases in data that belong to leaf pointed to by
   * the test (i.e. the subset of interest).  this is useful so
   * the instances belonging to that leaf aren't passed down the
   * other branch.
   *
   * @param data the instances to delete from
   */
  public void deleteGraftedCases(Instances data) {

    int subOfInterest = subsetOfInterest();
    for(int x = 0; x < data.numInstances(); x++) {
       if(whichSubset(data.instance(x)) == subOfInterest) {
          data.delete(x--);
       }
    }
  }


  /**
   * builds m_graftdistro using the passed data
   *
   * @param data the instances to use when creating the distribution
   */
  public void buildClassifier(Instances data) throws Exception {

    // distribution for the graft, not counting cases in atbop, only orig leaf
    m_graftdistro = new Distribution(2, data.numClasses());
 
    // which subset are we looking at for the graft?
    int subset = subsetOfInterest();  // this is the subset for m_leaf

    double thisNodeCount = 0;
    double knownCases = 0;
    boolean allKnown = true;
    // populate distribution
    for(int x = 0; x < data.numInstances(); x++) {
       Instance instance = data.instance(x);
       if(instance.isMissing(m_attIndex)) {
          allKnown = false;
          continue;
       }
       knownCases += instance.weight();
       int subst = whichSubset(instance);
       if(subst == -1)
          continue;
       m_graftdistro.add(subst, instance);
       if(subst == subset) {  // instance belongs at m_leaf
          thisNodeCount += instance.weight();
       }
    }
    double factor = (knownCases == 0) ? (1.0 / (double)2.0)
                                      : (thisNodeCount / knownCases);
    if(!allKnown) {
       for(int x = 0; x < data.numInstances(); x++) {
          if(data.instance(x).isMissing(m_attIndex)) {
             Instance instance = data.instance(x);
             int subst = whichSubset(instance);
             if(subst == -1)
                continue;
             instance.setWeight(instance.weight() * factor);
             m_graftdistro.add(subst, instance);
          }
       }
    }

    // if there are no cases at the leaf, make sure the desired
    // class is chosen, by setting counts to 0.01
    if(m_graftdistro.perBag(subset) == 0) {
       double [] counts = new double[data.numClasses()];
       counts[m_maxClass] = 0.01;
       m_graftdistro.add(subset, counts);
    }
    if(m_graftdistro.perBag((subset == 0) ? 1 : 0) == 0) {
       double [] counts = new double[data.numClasses()];
       counts[(int)m_otherLeafMaxClass] = 0.01;
       m_graftdistro.add((subset == 0) ? 1 : 0, counts);
    }
  }


  /**
   * @return the NoSplit object for the leaf pointed to by m_testType branch
   */
  public NoSplit getLeaf() {
    return new NoSplit(m_leafdistro);
  }


  /**
   * @return the NoSplit object for the leaf not pointed to by m_testType branch
   */
  public NoSplit getOtherLeaf() {

    // the bag (subset) that isn't pointed to by m_testType branch
    int bag = (subsetOfInterest() == 0) ? 1 : 0;

    double [][] counts = new double[1][m_graftdistro.numClasses()];
    double totals = 0;
    for(int c = 0; c < counts[0].length; c++) {
       counts[0][c] = m_graftdistro.perClassPerBag(bag, c);
       totals += counts[0][c];
    }
    // if empty, make sure proper class gets chosen
    if(totals == 0) {
       counts[0][m_otherLeafMaxClass] += 0.01;
    }
    return new NoSplit(new Distribution(counts));
  }


  /**
   * Prints label for subset index of instances (eg class).
   *
   * @param index the bag to dump label for
   * @param data to get attribute names and such
   * @return the label as a string
   * @exception Exception if something goes wrong
   */
  public final String dumpLabelG(int index, Instances data) throws Exception {

    StringBuffer text;

    text = new StringBuffer();
    text.append(((Instances)data).classAttribute().
       value((index==subsetOfInterest()) ? m_maxClass : m_otherLeafMaxClass));
    text.append(" ("+Utils.roundDouble(m_graftdistro.perBag(index),1));
    if(Utils.gr(m_graftdistro.numIncorrect(index),0))
       text.append("/"
        +Utils.roundDouble(m_graftdistro.numIncorrect(index),2));

    // show the graft values, only if this is subsetOfInterest()
    if(index == subsetOfInterest()) {
       text.append("|"+Utils.roundDouble(m_distribution.perBag(index),2));
       if(Utils.gr(m_distribution.numIncorrect(index),0))
          text.append("/"
             +Utils.roundDouble(m_distribution.numIncorrect(index),2));
    }
    text.append(")");
    return text.toString();
  }


  /**
   * @return the subset that is specified by the test type
   */
  public int subsetOfInterest() {
    if(m_testType == 2)
       return 0;
    if(m_testType == 3)
       return 1;
    return m_testType;
  }


  /**
   * @return the number of positive cases in the subset of interest
   */
  public double positivesForSubsetOfInterest() {
    return (m_distribution.perClassPerBag(subsetOfInterest(), m_maxClass));
  }


  /**
   * @param subset the subset to get the positives for
   * @return the number of positive cases in the specified subset
   */
  public double positives(int subset) {
    return (m_distribution.perClassPerBag(subset, 
                                    m_distribution.maxClass(subset)));
  }


  /**
   * @return the number of instances in the subset of interest
   */
  public double totalForSubsetOfInterest() {
    return (m_distribution.perBag(subsetOfInterest()));
  }

  
  /**
   * @param subset the index of the bag to get the total for
   * @return the number of instances in the subset
   */
  public double totalForSubset(int subset) {
    return (m_distribution.perBag(subset));
  }


  /**
   * Prints left side of condition satisfied by instances.
   *
   * @param data the data.
   */
  public String leftSide(Instances data) {
    return data.attribute(m_attIndex).name();
  }


  /**
   * @return the index of the attribute to split on
   */ 
  public int attribute() {
    return m_attIndex;
  }


  /**
   * Prints condition satisfied by instances in subset index.
   */
  public final String rightSide(int index, Instances data) {

    StringBuffer text;

    text = new StringBuffer();
    if(data.attribute(m_attIndex).isNominal() || data.attribute(m_attIndex).isRanking())
       if(index == 0)
          text.append(" = "+
                      data.attribute(m_attIndex).value((int)m_splitPoint));
       else
          text.append(" != "+
                      data.attribute(m_attIndex).value((int)m_splitPoint));
    else
       if(index == 0)
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
    if(index < 0) {
       return "i[" + m_attIndex + "] == null";
    }
    if(data.attribute(m_attIndex).isNominal() || data.attribute(m_attIndex).isRanking()) {
       if(index == 0)
          expr = new StringBuffer("i[");
       else
          expr = new StringBuffer("!i[");
       expr.append(m_attIndex).append("]");
       expr.append(".equals(\"").append(data.attribute(m_attIndex)
                                      .value((int)m_splitPoint)).append("\")");
    } else {
       expr = new StringBuffer("((Double) i[");
       expr.append(m_attIndex).append("])");
       if(index == 0) {
          expr.append(".doubleValue() <= ").append(m_splitPoint);
       } else {
          expr.append(".doubleValue() > ").append(m_splitPoint);
       }
    }
    return expr.toString();
  }


  /**
   * @param instance the instance to produce the weights for
   * @return a double array of weights, null if only belongs to one subset
   */
  public double [] weights(Instance instance) {

    double [] weights;
    int i;

    if(instance.isMissing(m_attIndex)) {
       weights = new double [m_numSubsets];
       for(i=0;i<m_numSubsets;i++) {
          weights [i] = m_graftdistro.perBag(i)/m_graftdistro.total();
       }
       return weights;
    } else {
       return null;
    }
  }


  /**
   * @param instance the instance for which to determine the subset
   * @return an int indicating the subset this instance belongs to
   */
  public int whichSubset(Instance instance) {

    if(instance.isMissing(m_attIndex))
       return -1;

    if(instance.attribute(m_attIndex).isNominal() || instance.attribute(m_attIndex).isRanking()) {
       // in the case of nominal, m_splitPoint is the = value, all else is !=
       if(instance.value(m_attIndex) == m_splitPoint)
          return 0;
       else
          return 1;
    } else {
       if(Utils.smOrEq(instance.value(m_attIndex), m_splitPoint))
          return 0;
       else
          return 1;
    }
  }


  /**
   * @return the value of the split point
   */
  public double splitPoint() {
    return m_splitPoint;
  }

  /**
   * @return the dominate class for the subset of interest
   */
  public int maxClassForSubsetOfInterest() {
    return m_maxClass;
  }

  /**
   * @return the laplace value for maxClass of subset of interest
   */
  public double laplaceForSubsetOfInterest() {
    return m_laplace;
  }

  /**
   * returns the test type
   * @return value of testtype
   */
  public int testType() {
    return m_testType;
  }

  /**
   * method needed for sorting a collection of GraftSplits by laplace value
   * @param g the graft split to compare to this one
   * @return -1, 0, or 1 if this GraftSplit laplace is <, = or > than that of g
   */
  public int compareTo(Object g) {

    if(m_laplace > ((GraftSplit)g).laplaceForSubsetOfInterest())
       return 1;
    if(m_laplace < ((GraftSplit)g).laplaceForSubsetOfInterest())
       return -1;
    return 0;
  }

  /**
   * returns the probability for instance for the specified class
   * @param classIndex the index of the class
   * @param instance the instance to get the probability for
   * @param theSubset the subset
   */
  public final double classProb(int classIndex, Instance instance, 
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
   * method for returning information about this GraftSplit
   * @param data instances for determining names of attributes and values
   * @return a string showing this GraftSplit's information
   */
  public String toString(Instances data) {

    String theTest;
    if(m_testType == 0)
       theTest = " <= ";
    else if(m_testType == 1)
       theTest = " > ";
    else if(m_testType == 2)
       theTest = " = ";
    else
       theTest = " != ";

    if(data.attribute(m_attIndex).isNominal() || data.attribute(m_attIndex).isRanking())
       theTest += data.attribute(m_attIndex).value((int)m_splitPoint);
    else
       theTest += Double.toString(m_splitPoint);

    return data.attribute(m_attIndex).name() + theTest
           + " (" + Double.toString(m_laplace) + ") --> " 
           + data.attribute(data.classIndex()).value(m_maxClass);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.2 $");
  }
}

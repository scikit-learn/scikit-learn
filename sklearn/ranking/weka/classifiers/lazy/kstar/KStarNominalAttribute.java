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

/**
 *    KStarNominalAttribute.java
 *    Copyright (C) 1995 Univeristy of Waikato
 *    Java port to Weka by Abdelaziz Mahoui (am14@cs.waikato.ac.nz).
 *
 */


package weka.classifiers.lazy.kstar;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.labelranking.PreferenceAttribute;

/**
 * A custom class which provides the environment for computing the
 * transformation probability of a specified test instance nominal
 * attribute to a specified train instance nominal attribute.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @author Abdelaziz Mahoui (am14@cs.waikato.ac.nz)
 * @version $Revision 1.0 $
 */
public class KStarNominalAttribute
  implements KStarConstants, RevisionHandler {
  
  /** The training instances used for classification. */
  protected Instances m_TrainSet;

  /** The test instance */
  protected Instance m_Test;

  /** The train instance */
  protected Instance m_Train;

  /** The index of the nominal attribute in the test and train instances */
  protected int m_AttrIndex;

  /** The stop parameter */
  protected double m_Stop = 1.0;

  /** Probability of test attribute transforming into train attribute 
      with missing value */
  protected double m_MissingProb = 1.0;

  /** Average probability of test attribute transforming into train 
      attribute */
  protected double m_AverageProb = 1.0;

  /** Smallest probability of test attribute transforming into 
      train attribute */
  protected double m_SmallestProb = 1.0;

  /** Number of trai instances with no missing attribute values */
  protected int m_TotalCount;

  /** Distribution of the attribute value in the train dataset */
  protected int [] m_Distribution;

  /** Set of colomns: each colomn representing a randomised version 
      of the train dataset class colomn */
  protected int [][] m_RandClassCols;

  /** A cache for storing attribute values and their corresponding 
      stop parameters */
  protected KStarCache m_Cache;

  // KStar Global settings

  /** The number of instances in the dataset */
  protected int m_NumInstances;

  /** The number of class values */
  protected int m_NumClasses;

  /** The number of attributes */
  protected int m_NumAttributes;

  /** The class attribute type */
  protected int m_ClassType;

  /** missing value treatment */
  protected int m_MissingMode = M_AVERAGE;

  /** B_SPHERE = use specified blend, B_ENTROPY = entropic blend setting */
  protected int m_BlendMethod = B_SPHERE ;

  /** default sphere of influence blend setting */
  protected int m_BlendFactor = 20;
  
  /**
   * Constructor
   */
  public KStarNominalAttribute(Instance test, Instance train, int attrIndex,
			       Instances trainSet, int [][] randClassCol, 
			       KStarCache cache)
  {
    m_Test = test;
    m_Train = train;
    m_AttrIndex = attrIndex;
    m_TrainSet = trainSet;
    m_RandClassCols = randClassCol;
    m_Cache = cache;
    init();
  }

  /**
   * Initializes the m_Attributes of the class.
   */
  private void init() {
    try {
      m_NumInstances  = m_TrainSet.numInstances();
      m_NumClasses    = m_TrainSet.numClasses();
      m_NumAttributes = m_TrainSet.numAttributes();
      m_ClassType     = m_TrainSet.classAttribute().type();
    } catch(Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Calculates the probability of the indexed nominal attribute of the test
   * instance transforming into the indexed nominal attribute of the training 
   * instance.
   *
   * @return the value of the transformation probability.
   */
  public double transProb() {
    String debug = "(KStarNominalAttribute.transProb) ";
    double transProb = 0.0;
    // check if the attribute value has been encountred before
    // in which case it should be in the nominal cache
    if (m_Cache.containsKey(m_Test.value(m_AttrIndex))) {
      KStarCache.TableEntry te = 
	m_Cache.getCacheValues(m_Test.value(m_AttrIndex));
      m_Stop = te.value;
      m_MissingProb = te.pmiss;
    }
    else {
      generateAttrDistribution();
      // we have to compute the parameters
      if (m_BlendMethod == B_ENTROPY) {
	m_Stop = stopProbUsingEntropy();
      }
      else { // default is B_SPHERE
	m_Stop = stopProbUsingBlend();
      }
      // store the values in cache
      m_Cache.store( m_Test.value(m_AttrIndex), m_Stop, m_MissingProb );
    }
    // we've got our m_Stop, then what?
    if (m_Train.isMissing(m_AttrIndex)) {
      transProb = m_MissingProb;
    }
    else {
      try {
	transProb = (1.0 - m_Stop) / m_Test.attribute(m_AttrIndex).numValues();
	if ( (int)m_Test.value(m_AttrIndex) == 
	     (int)m_Train.value(m_AttrIndex) )
	  {
	    transProb += m_Stop;
	  }
      } catch (Exception e) {
	e.printStackTrace();
      }
    }
    return transProb;
  }
  
  /**
   * Calculates the "stop parameter" for this attribute using
   * the entropy method: the value is computed using a root finder
   * algorithm. The method takes advantage of the calculation to
   * compute the smallest and average transformation probabilities
   * once the stop factor is obtained. It also sets the transformation
   * probability to an attribute with a missing value.
   *
   * @return the value of the stop parameter.
   *
   */
  private double stopProbUsingEntropy() {
    String debug = "(KStarNominalAttribute.stopProbUsingEntropy)";
    if ( m_ClassType != Attribute.NOMINAL && m_ClassType != PreferenceAttribute.RANKING ) {
      System.err.println("Error: "+debug+" attribute class must be nominal!");
      System.exit(1);
    }
    int itcount = 0;
    double stopProb;
    double lower, upper, pstop;
    double bestminprob = 0.0, bestpsum = 0.0;
    double bestdiff = 0.0, bestpstop = 0.0;
    double currentdiff, lastdiff, stepsize, delta;
    
    KStarWrapper botvals = new KStarWrapper();
    KStarWrapper upvals = new KStarWrapper();
    KStarWrapper vals = new KStarWrapper();

    // Initial values for root finder
    lower = 0.0 + ROOT_FINDER_ACCURACY/2.0;
    upper = 1.0 - ROOT_FINDER_ACCURACY/2.0;
    
    // Find (approx) entropy ranges
    calculateEntropy(upper, upvals);
    calculateEntropy(lower, botvals);
    
    if (upvals.avgProb == 0) {
      // When there are no training instances with the test value:
      // doesn't matter what exact value we use for pstop, just acts as
      // a constant scale factor in this case.
      calculateEntropy(lower, vals);
    }
    else
      {
	// Optimise the scale factor
	if ( (upvals.randEntropy - upvals.actEntropy < 
	      botvals.randEntropy - botvals.actEntropy) &&
	     (botvals.randEntropy - botvals.actEntropy > FLOOR) )
	  {
	    bestpstop = pstop = lower;
	    stepsize = INITIAL_STEP;
	    bestminprob = botvals.minProb;
	    bestpsum = botvals.avgProb;
	  }
	else {
	  bestpstop = pstop = upper;
	  stepsize = -INITIAL_STEP;
	  bestminprob = upvals.minProb;
	  bestpsum = upvals.avgProb;
	}
	bestdiff = currentdiff = FLOOR;
	itcount = 0;
	/* Enter the root finder */
	while (true)
	  {
	    itcount++;	
	    lastdiff = currentdiff;
	    pstop += stepsize;
	    if (pstop <= lower) {
	      pstop = lower;
	      currentdiff = 0.0;
	      delta = -1.0;
	    }
	    else if (pstop >= upper) {
	      pstop = upper;
	      currentdiff = 0.0;
	      delta = -1.0;
	    }
	    else {
	      calculateEntropy(pstop, vals);
	      currentdiff = vals.randEntropy - vals.actEntropy;

	      if (currentdiff < FLOOR) {
		currentdiff = FLOOR;
		if ((Math.abs(stepsize) < INITIAL_STEP) && 
		    (bestdiff == FLOOR)) {
		  bestpstop = lower;
		  bestminprob = botvals.minProb;
		  bestpsum = botvals.avgProb;
		  break;
		}
	      }
	      delta = currentdiff - lastdiff;
	    }
	    if (currentdiff > bestdiff) {
	      bestdiff = currentdiff;
	      bestpstop = pstop;
	      bestminprob = vals.minProb;
	      bestpsum = vals.avgProb;
	    }
	    if (delta < 0) {
	      if (Math.abs(stepsize) < ROOT_FINDER_ACCURACY) {
		break;
	      }
	      else {
		stepsize /= -2.0;
	      }
	    }
	    if (itcount > ROOT_FINDER_MAX_ITER) {
	      break;
	    }
	  }
      }
    
    m_SmallestProb = bestminprob;
    m_AverageProb = bestpsum;
    // Set the probability of transforming to a missing value
    switch ( m_MissingMode )
      {
      case M_DELETE:
	m_MissingProb = 0.0;
	break;
      case M_NORMAL:
	m_MissingProb = 1.0;
	break;
      case M_MAXDIFF:
	m_MissingProb = m_SmallestProb;
	break;
      case M_AVERAGE:
	m_MissingProb = m_AverageProb;
	break;
      }

    if ( Math.abs(bestpsum - (double)m_TotalCount) < EPSILON) { 
      // No difference in the values
      stopProb = 1.0;
    }
    else {
      stopProb = bestpstop;
    }
    return stopProb;
  }

  /**
   * Calculates the entropy of the actual class prediction
   * and the entropy for random class prediction. It also
   * calculates the smallest and average transformation probabilities.
   *
   * @param stop the stop parameter
   * @param params the object wrapper for the parameters:
   * actual entropy, random entropy, average probability and smallest 
   * probability.
   * @return the values are returned in the object "params".
   *
   */
  private void calculateEntropy( double stop, KStarWrapper params) {
    String debug = "(KStarNominalAttribute.calculateEntropy)";
    int i,j,k;
    Instance train;
    double actent = 0.0, randent=0.0;
    double pstar, tprob, psum=0.0, minprob=1.0;
    double actClassProb, randClassProb;
    double [][] pseudoClassProb = new double[NUM_RAND_COLS+1][m_NumClasses];
    // init ...
    for(j = 0; j <= NUM_RAND_COLS; j++) {
      for(i = 0; i < m_NumClasses; i++) {
	pseudoClassProb[j][i] = 0.0;
      }
    }
    for (i=0; i < m_NumInstances; i++) {
      train = m_TrainSet.instance(i);
      if (!train.isMissing(m_AttrIndex)) {
	pstar = PStar(m_Test, train, m_AttrIndex, stop);
	tprob = pstar / m_TotalCount;
	if (pstar < minprob) {
	  minprob = pstar;
	}
	psum += tprob;
	// filter instances with same class value
	for (k=0 ; k <= NUM_RAND_COLS ; k++) {
	  // instance i is assigned a random class value in colomn k;
	  // colomn k = NUM_RAND_COLS contains the original mapping: 
	  // instance -> class vlaue
	  pseudoClassProb[k][ m_RandClassCols[k][i] ] += tprob;
	}
      }
    }
    // compute the actual entropy using the class probs
    // with the original class value mapping (colomn NUM_RAND_COLS)
    for (j=m_NumClasses-1; j>=0; j--) {
      actClassProb = pseudoClassProb[NUM_RAND_COLS][j] / psum;
      if (actClassProb > 0) {
    	actent -= actClassProb * Math.log(actClassProb) / LOG2;
      }
    }
    // compute a random entropy using the pseudo class probs
    // excluding the colomn NUM_RAND_COLS
    for (k=0; k < NUM_RAND_COLS;k++) {
      for (i = m_NumClasses-1; i >= 0; i--) {
  	randClassProb = pseudoClassProb[k][i] / psum;
  	if (randClassProb > 0) {
  	  randent -= randClassProb * Math.log(randClassProb) / LOG2;
	}
      }
    }
    randent /= NUM_RAND_COLS;
    // return the results ... Yuk !!!
    params.actEntropy = actent;
    params.randEntropy = randent;
    params.avgProb = psum;
    params.minProb = minprob;
  }
  
  /**
   * Calculates the "stop parameter" for this attribute using
   * the blend method: the value is computed using a root finder
   * algorithm. The method takes advantage of this calculation to
   * compute the smallest and average transformation probabilities
   * once the stop factor is obtained. It also sets the transformation
   * probability to an attribute with a missing value.
   *
   * @return the value of the stop parameter.
   *
   */
  private double stopProbUsingBlend() {
    String debug = "(KStarNominalAttribute.stopProbUsingBlend) ";
    int itcount = 0;
    double stopProb, aimfor;
    double lower, upper, tstop;

    KStarWrapper botvals = new KStarWrapper();
    KStarWrapper upvals = new KStarWrapper();
    KStarWrapper vals = new KStarWrapper();

    int testvalue = (int)m_Test.value(m_AttrIndex);
    aimfor = (m_TotalCount - m_Distribution[testvalue]) * 
      (double)m_BlendFactor / 100.0 + m_Distribution[testvalue];

    // Initial values for root finder
    tstop = 1.0 - (double)m_BlendFactor / 100.0;
    lower = 0.0 + ROOT_FINDER_ACCURACY/2.0;
    upper = 1.0 - ROOT_FINDER_ACCURACY/2.0;

    // Find out function border values
    calculateSphereSize(testvalue, lower, botvals);
    botvals.sphere -= aimfor;
    calculateSphereSize(testvalue, upper, upvals);
    upvals.sphere -= aimfor;
    
    if (upvals.avgProb == 0) {
      // When there are no training instances with the test value:
      // doesn't matter what exact value we use for tstop, just acts as
      // a constant scale factor in this case.
      calculateSphereSize(testvalue, tstop, vals);
    }
    else if (upvals.sphere > 0) {
      // Can't include aimfor instances, going for min possible
      tstop = upper;
      vals.avgProb = upvals.avgProb;
    }
    else {
      // Enter the root finder
      for (;;) {
	itcount++;
	calculateSphereSize(testvalue, tstop, vals);
	vals.sphere -= aimfor;
	if ( Math.abs(vals.sphere) <= ROOT_FINDER_ACCURACY ||
	     itcount >= ROOT_FINDER_MAX_ITER )
	  {
	    break;
	  }
	if (vals.sphere > 0.0) {
	  lower = tstop;
	  tstop = (upper + lower) / 2.0;
	}
	else {
	  upper = tstop;
	  tstop = (upper + lower) / 2.0;
	}
      }
    }

    m_SmallestProb = vals.minProb;
    m_AverageProb = vals.avgProb;
    // Set the probability of transforming to a missing value
    switch ( m_MissingMode )
      {
      case M_DELETE:
	m_MissingProb = 0.0;
	break;
      case M_NORMAL:
	m_MissingProb = 1.0;
	break;
      case M_MAXDIFF:
	m_MissingProb = m_SmallestProb;
	break;
      case M_AVERAGE:
	m_MissingProb = m_AverageProb;
	break;
      }
    
    if ( Math.abs(vals.avgProb - m_TotalCount) < EPSILON) { 
      // No difference in the values
      stopProb = 1.0;
    }
    else {
      stopProb = tstop;
    }
    return stopProb;
  }
  
  /**
   * Calculates the size of the "sphere of influence" defined as:
   * sphere = sum(P^2)/sum(P)^2
   * P(i|j) = (1-tstop)*P(i) + ((i==j)?tstop:0).
   * This method takes advantage of the calculation to compute the values of
   * the "smallest" and "average" transformation probabilities when using
   * the specified stop parameter.
   *
   * @param testValue the value of the test instance
   * @param stop the stop parameter
   * @param params a wrapper of the parameters to be computed:
   * "sphere" the sphere size
   * "avgprob" the average transformation probability
   * "minProb" the smallest transformation probability
   * @return the values are returned in "params" object.
   *
   */
  private void calculateSphereSize(int testvalue, double stop, 
				   KStarWrapper params) {
    String debug = "(KStarNominalAttribute.calculateSphereSize) ";
    int i, thiscount;
    double tprob, tval = 0.0, t1 = 0.0;
    double sphere, minprob = 1.0, transprob = 0.0;

    for(i = 0; i < m_Distribution.length; i++) {
      thiscount = m_Distribution[i];
      if ( thiscount != 0 ) {
	if ( testvalue == i ) {
	  tprob = (stop + (1 - stop) / m_Distribution.length) / m_TotalCount;
	  tval += tprob * thiscount;
	  t1 += tprob * tprob * thiscount;
	}
	else {
	  tprob = ((1 - stop) / m_Distribution.length) / m_TotalCount;
	  tval += tprob * thiscount;
	  t1 += tprob * tprob * thiscount;
	}
	if ( minprob > tprob * m_TotalCount ) {
	  minprob = tprob * m_TotalCount;
	}
      }
    }
    transprob = tval;
    sphere = (t1 == 0) ? 0 : ((tval * tval) / t1);
    // return values ... Yck!!!
    params.sphere = sphere;
    params.avgProb = transprob;
    params.minProb = minprob;
  }
  
  /**
   * Calculates the nominal probability function defined as:
   * P(i|j) = (1-stop) * P(i) + ((i==j) ? stop : 0)
   * In this case, it calculates the transformation probability of the
   * indexed test attribute to the indexed train attribute.
   *
   * @param test the test instance
   * @param train the train instance
   * @param col the attribute index
   * @return the value of the tranformation probability.
   *
   */
  private double PStar(Instance test, Instance train, int col, double stop) {
    String debug = "(KStarNominalAttribute.PStar) ";
    double pstar;
    int numvalues = 0;
    try {
      numvalues = test.attribute(col).numValues();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    if ( (int)test.value(col) == (int)train.value(col) ) {
      pstar = stop + (1 - stop) / numvalues;
    }
    else {
      pstar = (1 - stop) / numvalues;
    }
    return pstar;
  }
  
  /**
   * Calculates the distribution, in the dataset, of the indexed nominal
   * attribute values. It also counts the actual number of training instances
   * that contributed (those with non-missing values) to calculate the 
   * distribution.
   */
  private void generateAttrDistribution() {
    String debug = "(KStarNominalAttribute.generateAttrDistribution)";
    m_Distribution = new int[ m_TrainSet.attribute(m_AttrIndex).numValues() ];
    int i;
    Instance train;
    for (i=0; i < m_NumInstances; i++) {
      train = m_TrainSet.instance(i);
      if ( !train.isMissing(m_AttrIndex) ) {
	m_TotalCount++;
	m_Distribution[(int)train.value(m_AttrIndex)]++;
      }
    }
  }

  /**
   * Sets the options.
   *
   */
  public void setOptions(int missingmode, int blendmethod, int blendfactor) {
    m_MissingMode = missingmode;
    m_BlendMethod = blendmethod;
    m_BlendFactor = blendfactor;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.7 $");
  }
} // class

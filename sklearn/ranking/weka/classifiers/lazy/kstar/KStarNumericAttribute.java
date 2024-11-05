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
 *    KStarNumericAttribute.java
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
 * transformation probability of a specified test instance numeric
 * attribute to a specified train instance numeric attribute.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @author Abdelaziz Mahoui (am14@cs.waikato.ac.nz)
 * @version $Revision 1.0 $
 */
public class KStarNumericAttribute
  implements KStarConstants, RevisionHandler {

  /** The training instances used for classification. */
  protected Instances m_TrainSet;

  /** The test instance */
  protected Instance m_Test;

  /** The train instance */
  protected Instance m_Train;

  /** The index of the attribute in the test and train instances */
  protected int m_AttrIndex;

  /** The scale parameter */
  protected double m_Scale = 1.0;

  /** Probability of test attribute transforming into train attribute 
      with missing value */
  protected double m_MissingProb = 1.0;

  /** Average probability of test attribute transforming into train 
      attribute */
  protected double m_AverageProb = 1.0;

  /** Smallest probability of test attribute transforming into train 
      attribute */
  protected double m_SmallestProb = 1.0;

  /** The set of disctances from the test attribute to the set of train 
      attributes */
  protected double [] m_Distances;

  /** Set of colomns: each colomn representing a randomised version of 
      the train dataset class colomn */
  protected int [][] m_RandClassCols;

  /** The number of train instances with no missing attribute values */
  protected int m_ActualCount = 0;

  /** A cache for storing attribute values and their corresponding scale 
      parameters */
  protected KStarCache m_Cache;

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

  /** 0 = use specified blend, 1 = entropic blend setting */
  protected int m_BlendMethod = B_SPHERE ;

  /** default sphere of influence blend setting */
  protected int m_BlendFactor = 20;
  
  /**
   * Constructor
   */
  public KStarNumericAttribute(Instance test, Instance train, int attrIndex,
			       Instances trainSet, 
			       int [][] randClassCols, 
			       KStarCache cache)
  {
    m_Test      = test;
    m_Train     = train;
    m_AttrIndex = attrIndex;
    m_TrainSet  = trainSet;
    m_RandClassCols = randClassCols;
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
   * Calculates the transformation probability of the attribute indexed
   * "m_AttrIndex" in test instance "m_Test" to the same attribute in
   * the train instance "m_Train".
   *
   * @return the probability value
   */
  public double transProb() {
    String debug = "(KStarNumericAttribute.transProb) ";
    double transProb, distance, scale;
    // check if the attribute value has been encountred before
    // in which case it should be in the numeric cache
    if ( m_Cache.containsKey(m_Test.value(m_AttrIndex))) {
      KStarCache.TableEntry te = 
	m_Cache.getCacheValues( m_Test.value(m_AttrIndex) );
      m_Scale = te.value;
      m_MissingProb = te.pmiss;
    }
    else {
      if (m_BlendMethod == B_ENTROPY) {
	m_Scale = scaleFactorUsingEntropy();
      }
      else { // default is B_SPHERE
	m_Scale = scaleFactorUsingBlend();
      }
      m_Cache.store( m_Test.value(m_AttrIndex), m_Scale, m_MissingProb );
    }
    // now what???
    if (m_Train.isMissing(m_AttrIndex)) {
      transProb = m_MissingProb;
    }
    else {
      distance = 
	Math.abs( m_Test.value(m_AttrIndex) - m_Train.value(m_AttrIndex) );
      transProb = PStar( distance, m_Scale );
    }
    return transProb;
  }
  
  /**
   * Calculates the scale factor for the attribute indexed
   * "m_AttrIndex" in test instance "m_Test" using a global
   * blending factor (default value is 20%).
   *
   * @return the scale factor value
   */
  private double scaleFactorUsingBlend() {
    String debug = "(KStarNumericAttribute.scaleFactorUsingBlend)";
    int i, j, lowestcount = 0, count = 0;
    double lowest = -1.0, nextlowest = -1.0;
    double root, broot, up, bot;
    double aimfor, min_val = 9e300, scale = 1.0;
    double avgprob = 0.0, minprob = 0.0, min_pos = 0.0;

    KStarWrapper botvals = new KStarWrapper();
    KStarWrapper upvals = new KStarWrapper();
    KStarWrapper vals = new KStarWrapper();

    m_Distances = new double [m_NumInstances];

    for (j=0; j<m_NumInstances; j++) {
      if ( m_TrainSet.instance(j).isMissing(m_AttrIndex) ) {
	// mark the train instance with a missing value by setting 
	// the distance to -1.0
	m_Distances[j] = -1.0;
      }
      else {
	m_Distances[j] = Math.abs(m_TrainSet.instance(j).value(m_AttrIndex) - 
				  m_Test.value(m_AttrIndex));
	if ( (m_Distances[j]+1e-5) < nextlowest || nextlowest == -1.0 ) {
	  if ( (m_Distances[j]+1e-5) < lowest || lowest == -1.0 ) {
	    nextlowest = lowest;
	    lowest = m_Distances[j];
	    lowestcount = 1;
	  }
	  else if ( Math.abs(m_Distances[j]-lowest) < 1e-5 ) {
	    // record the number training instances (number n0) at
	    // the smallest distance from test instance
	    lowestcount++;
	  }
	  else {
	    nextlowest = m_Distances[j];
	  }
	}
	// records the actual number of instances with no missing value
	m_ActualCount++;
      }
    }
    
    if (nextlowest == -1 || lowest == -1) { // Data values are all the same
      scale = 1.0;
      m_SmallestProb = m_AverageProb = 1.0;
      return scale;
    }
    else {
      // starting point for root
      root = 1.0 / (nextlowest - lowest);
      i = 0;
      // given the expression: n0 <= E(scale) <= N
      // E(scale) =  (N - n0) * b + n0  with blending factor: 0 <= b <= 1
      // aimfor = (N - n0) * b + n0
      aimfor = (m_ActualCount - lowestcount) * 
	(double)m_BlendFactor / 100.0 + lowestcount;
      if (m_BlendFactor == 0) {
	aimfor += 1.0;
      }
      // root is bracketed in interval [bot,up]
      bot = 0.0 + ROOT_FINDER_ACCURACY / 2.0;
      up = root * 16;     // This is bodgy
      // E(bot)
      calculateSphereSize(bot, botvals);
      botvals.sphere -= aimfor;
      // E(up)
      calculateSphereSize(up, upvals);
      upvals.sphere -= aimfor;
      
      if (botvals.sphere < 0) {    // Couldn't include that many 
	                           // instances - going for max possible
	min_pos = bot;
	avgprob = botvals.avgProb;
	minprob = botvals.minProb;
      }
      else if (upvals.sphere > 0) { // Couldn't include that few, 
	                            // going for min possible
	min_pos = up;
	avgprob = upvals.avgProb;
	minprob = upvals.minProb;
      }
      else {
	// Root finding Algorithm starts here !
	for (;;) {
	  calculateSphereSize(root, vals);
	  vals.sphere -= aimfor;
	  if ( Math.abs(vals.sphere) < min_val ) {
	    min_val = Math.abs(vals.sphere);
	    min_pos = root;
	    avgprob = vals.avgProb;
	    minprob = vals.minProb;
	  }
	  if ( Math.abs(vals.sphere) <= ROOT_FINDER_ACCURACY ) {
	    break;        // converged to a solution, done!
	  }
	  if (vals.sphere > 0.0) {
	    broot = (root + up) / 2.0;
	    bot = root;
	    root = broot;
	  }
	  else {
	    broot = (root + bot) / 2.0;
	    up = root;
	    root = broot;
	  }
	  i++;
	  if (i > ROOT_FINDER_MAX_ITER) {
	    //	    System.err.println("Warning: "+debug+" 
	    // ROOT_FINDER_MAX_ITER exceeded");
	    root = min_pos;
	    break;
	  }
	}
      }

      m_SmallestProb = minprob;
      m_AverageProb = avgprob;
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
      // set the scale factor value
      scale = min_pos;
      return scale;
    }
  }
  
  /**
   * Calculates the size of the "sphere of influence" defined as:
   * sphere = sum(P)^2/sum(P^2) where
   * P(i) = root*exp(-2*i*root).
   * Since there are n different training instances we multiply P(i) by 1/n.
   */
  private void calculateSphereSize(double scale, KStarWrapper params) {
    String debug = "(KStarNumericAttribute.calculateSphereSize)";
    int i;
    double sphereSize, minprob = 1.0;
    double pstar;                // P*(b|a)
    double pstarSum = 0.0;       // sum(P*)
    double pstarSquareSum = 0.0; // sum(P*^2)
    double inc;
    for (i = 0; i < m_NumInstances; i++) {
      if (m_Distances[i] < 0) {
	// instance with missing value
	continue;
      }
      else {
	pstar = PStar( m_Distances[i], scale );
	if (minprob > pstar) {
	  minprob = pstar;
	}
	inc = pstar / m_ActualCount;
	pstarSum += inc;
	pstarSquareSum += inc * inc;
      }
    }
    sphereSize = (pstarSquareSum == 0 ? 0 
		  : pstarSum * pstarSum / pstarSquareSum);
    // return the values
    params.sphere = sphereSize;
    params.avgProb = pstarSum;
    params.minProb = minprob;
  }
  
  /**
   * Calculates the scale factor using entropy.
   *
   * @return the scale factor value
   */
  private double scaleFactorUsingEntropy() {
    String debug = "(KStarNumericAttribute.scaleFactorUsingEntropy)";
    if ( m_ClassType != Attribute.NOMINAL && m_ClassType != PreferenceAttribute.RANKING) {
      System.err.println("Error: "+debug+" attribute class must be nominal!");
      System.exit(1);
    }
    int i,j, lowestcount = 0, count, itcount;
    double lowest = -1.0, nextlowest = -1.0;
    double root, up, bot, stepsize, delta;
    double actentropy = 0.0, randentropy = 0.0, actscale, randscale;
    double minrand = 0.0, minact = 0.0, maxrand = 0.0, maxact = 0.0;
    double bestdiff, bestroot, currentdiff, lastdiff;
    double bestpsum, bestminprob, scale = 1.0;

    KStarWrapper botvals = new KStarWrapper();
    KStarWrapper upvals = new KStarWrapper();
    KStarWrapper vals = new KStarWrapper();

    m_Distances = new double [m_NumInstances];

    for (j=0; j<m_NumInstances; j++) {
      if ( m_TrainSet.instance(j).isMissing(m_AttrIndex) ) {
	// mark the train instance with a missing value by setting 
	// the distance to -1.0
	m_Distances[j] = -1.0;
      }
      else {
	m_Distances[j] = Math.abs(m_TrainSet.instance(j).value(m_AttrIndex) - 
				  m_Test.value(m_AttrIndex));
	
	if ( (m_Distances[j]+1e-5) < nextlowest || nextlowest == -1.0 ) {
	  if ( (m_Distances[j]+1e-5) < lowest || lowest == -1.0 ) {
	    nextlowest = lowest;
	    lowest = m_Distances[j];
	    lowestcount = 1;
	  }
	  else if ( Math.abs(m_Distances[j]-lowest) < 1e-5 ) {
	    // record the number training instances (number n0) at
	    // the smallest distance from test instance
	    lowestcount++;
	  }
	  else {
	    nextlowest = m_Distances[j];
	  }
	}
	// records the actual number of instances with no missing value
	m_ActualCount++;
      }
    } // for
    
    if (nextlowest == -1 || lowest == -1) { // Data values are all the same
      scale = 1.0;
      m_SmallestProb = m_AverageProb = 1.0;
      return scale;
    }
    else {
      // starting point for root
      root = 1.0 / (nextlowest - lowest);
      // root is bracketed in interval [bot,up]
      bot = 0.0 + ROOT_FINDER_ACCURACY / 2;  
      up = root * 8; // This is bodgy
      // Find (approx) entropy ranges
      calculateEntropy(up, upvals);
      calculateEntropy(bot, botvals);
      actscale = botvals.actEntropy - upvals.actEntropy;
      randscale = botvals.randEntropy - upvals.randEntropy;
      // Optimise the scale factor
      bestroot = root = bot;
      bestdiff = currentdiff = FLOOR1;
      bestpsum = botvals.avgProb;
      bestminprob = botvals.minProb;
      stepsize = (up - bot) / 20.0;
      itcount = 0;
      // Root finding algorithm starts here!
      while (true)
	{
	  itcount++;
	  lastdiff = currentdiff;
	  root += Math.log(root + 1.0) * stepsize;
	  if (root <= bot) {
	    root = bot;
	    currentdiff = 0.0;
	    delta = -1.0;
	  }
	  else if (root >= up) {
	    root = up;
	    currentdiff = 0.0;
	    delta = -1.0;
	  }
	  else {
	    calculateEntropy(root, vals);
	    // Normalise entropies
	    vals.randEntropy = (vals.randEntropy - upvals.randEntropy) / 
	      randscale;
	    vals.actEntropy = (vals.actEntropy - upvals.actEntropy) / 
	      randscale;
	    currentdiff = vals.randEntropy - vals.actEntropy;

	    if (currentdiff < FLOOR1) {
	      currentdiff = FLOOR1;
	      if (stepsize < 0) { 
		// If we've hit the end and turned around we can't 
		// have found any peaks
		bestdiff = currentdiff;
		bestroot = bot;
		bestpsum = botvals.avgProb;
		bestminprob = botvals.minProb;
		break;
	      }
	    }
	    delta = currentdiff - lastdiff;
	  }
	  if (currentdiff > bestdiff) {
	    bestdiff = currentdiff;
	    bestroot = root;
	    bestminprob = vals.minProb;
	    bestpsum = vals.avgProb;
	  }
	  if (delta < 0) {
	    if (Math.abs(stepsize) < ROOT_FINDER_ACCURACY) {
	      break;
	    }
	    else {
	      stepsize /= -4.0;
	    }
	  }
	  if (itcount > ROOT_FINDER_MAX_ITER) {
	    //  System.err.println("Warning: "+debug+" ROOT_FINDER_MAX_ITER 
	    // exceeded");
	    break;
	  }
	} // while

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
      // set scale factor
      scale = bestroot;
    } // else
    return scale;
  }

  /**
   * Calculates several parameters aside from the entropy: for a specified
   * scale factor, calculates the actual entropy, a random entropy using a
   * randomized set of class value colomns, and records the average and
   * smallest probabilities (for use in missing value case).
   */
  private void calculateEntropy(double scale, KStarWrapper params) {    
    String debug = "(KStarNumericAttribute.calculateEntropy)";
    int i,j,k;
    double actent = 0.0, randent = 0.0;
    double pstar, tprob, avgprob = 0.0, minprob = 1.0;
    double actClassProb, randClassProb;
    double [][] pseudoClassProbs = new double[NUM_RAND_COLS+1][m_NumClasses];
    // init
    for(j = 0; j <= NUM_RAND_COLS; j++) {
      for(i = 0; i < m_NumClasses; i++) {
	pseudoClassProbs[j][i] = 0.0;
      }
    }
    for (i=0; i < m_NumInstances; i++) {
      if (m_Distances[i] < 0) {
	// train instance has mising value
	continue;
      }
      else {
	pstar = PStar(m_Distances[i], scale);
	tprob = pstar / m_ActualCount;
	avgprob += tprob;
	if (pstar < minprob) {
	  minprob = pstar;
	}
	// filter instances with same class value
	for (k=0; k <= NUM_RAND_COLS; k++) {
	  // instance i is assigned a random class value in colomn k;
	  // colomn k = NUM_RAND_COLS contains the original mapping: 
	  // instance -> class vlaue
	  pseudoClassProbs[k][ m_RandClassCols[k][i] ] += tprob;
	}
      }
    }
    // compute the actual entropy using the class probabilities
    // with the original class value mapping (colomn NUM_RAND_COLS)
    for (j = m_NumClasses-1; j >= 0; j--) {
      actClassProb = pseudoClassProbs[NUM_RAND_COLS][j] / avgprob;
      if (actClassProb > 0) {
    	actent -= actClassProb * Math.log(actClassProb) / LOG2;
      }
    }
    // compute a random entropy using the pseudo class probs
    // excluding the colomn NUM_RAND_COLS
    for (k=0; k < NUM_RAND_COLS; k++) {
      for (i = m_NumClasses-1; i >= 0; i--) {
  	randClassProb = pseudoClassProbs[k][i] / avgprob;
  	if (randClassProb > 0) {
  	  randent -= randClassProb * Math.log(randClassProb) / LOG2;
	}
      }
    }
    randent /= NUM_RAND_COLS;
    // return the values
    params.actEntropy = actent;
    params.randEntropy = randent;
    params.avgProb = avgprob;
    params.minProb = minprob;
  }

  /**
   * Calculates the value of P for a given value x using the expression:
   * P(x) = scale * exp( -2.0 * x * scale )
   *
   * @param x input value
   * @param scale the scale factor
   * @return output of the function P(x)
   */          
  private double PStar(double x, double scale) {
    return scale * Math.exp( -2.0 * x * scale );
  }

  /**
   * Set options.
   * @param missingmode the missing value treatment to use
   * @param blendmethod the blending method to use
   * @param blendfactor the level of blending to use
   */
  public void setOptions(int missingmode, int blendmethod, int blendfactor) {
    m_MissingMode = missingmode;
    m_BlendMethod = blendmethod;
    m_BlendFactor = blendfactor;
  }

  /**
   * Set the missing value mode.
   * @param mode the type of missing value treatment to use
   */
  public void setMissingMode(int mode) {
    m_MissingMode = mode;
  }

  /**
   * Set the blending method
   * @param method the blending method to use
   */
  public void setBlendMethod(int method) {
    m_BlendMethod = method;
  }

  /**
   * Set the blending factor
   * @param factor the level of blending to use
   */
  public void setBlendFactor(int factor) {
    m_BlendFactor = factor;
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

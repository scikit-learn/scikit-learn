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
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */
 
package weka.classifiers.bayes.net.estimate;

import weka.classifiers.bayes.BayesNet;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.estimators.Estimator;

import java.util.Enumeration;

/** 
 <!-- globalinfo-start -->
 * SimpleEstimator is used for estimating the conditional probability tables of a Bayes network once the structure has been learned. Estimates probabilities directly from data.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -A &lt;alpha&gt;
 *  Initial count (alpha)
 * </pre>
 * 
 <!-- options-end -->
 * 
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.6 $
 */
public class SimpleEstimator 
    extends BayesNetEstimator {

    /** for serialization */
    static final long serialVersionUID = 5874941612331806172L;
    
    /**
     * Returns a string describing this object
     * @return a description of the classifier suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
      return 
          "SimpleEstimator is used for estimating the conditional probability "
      	+ "tables of a Bayes network once the structure has been learned. "
      	+ "Estimates probabilities directly from data.";
    }
  
    /**
     * estimateCPTs estimates the conditional probability tables for the Bayes
     * Net using the network structure.
     * 
     * @param bayesNet the bayes net to use
     * @throws Exception if something goes wrong
     */
    public void estimateCPTs(BayesNet bayesNet) throws Exception {
            initCPTs(bayesNet);

            // Compute counts
            Enumeration enumInsts = bayesNet.m_Instances.enumerateInstances();
            while (enumInsts.hasMoreElements()) {
                Instance instance = (Instance) enumInsts.nextElement();

                updateClassifier(bayesNet, instance);
            }
    } // estimateCPTs

    /**
     * Updates the classifier with the given instance.
     * 
     * @param bayesNet the bayes net to use
     * @param instance the new training instance to include in the model
     * @throws Exception if the instance could not be incorporated in
     * the model.
     */
    public void updateClassifier(BayesNet bayesNet, Instance instance) throws Exception {
        for (int iAttribute = 0; iAttribute < bayesNet.m_Instances.numAttributes(); iAttribute++) {
            double iCPT = 0;

            for (int iParent = 0; iParent < bayesNet.getParentSet(iAttribute).getNrOfParents(); iParent++) {
                int nParent = bayesNet.getParentSet(iAttribute).getParent(iParent);

                iCPT = iCPT * bayesNet.m_Instances.attribute(nParent).numValues() + instance.value(nParent);
            }

            bayesNet.m_Distributions[iAttribute][(int) iCPT].addValue(instance.value(iAttribute), instance.weight());
        }
    } // updateClassifier


    /** 
     * initCPTs reserves space for CPTs and set all counts to zero
     * 
     * @param bayesNet the bayes net to use
     * @throws Exception if something goes wrong
     */
    public void initCPTs(BayesNet bayesNet) throws Exception {
        Instances instances = bayesNet.m_Instances;
        
        // Reserve space for CPTs
        int nMaxParentCardinality = 1;
        for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
            if (bayesNet.getParentSet(iAttribute).getCardinalityOfParents() > nMaxParentCardinality) {
                nMaxParentCardinality = bayesNet.getParentSet(iAttribute).getCardinalityOfParents();
            }
        }
	
        // Reserve plenty of memory
        bayesNet.m_Distributions = new Estimator[instances.numAttributes()][nMaxParentCardinality];
	
        // estimate CPTs
        for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
            for (int iParent = 0; iParent < bayesNet.getParentSet(iAttribute).getCardinalityOfParents(); iParent++) {
                bayesNet.m_Distributions[iAttribute][iParent] =
                    new DiscreteEstimatorBayes(instances.attribute(iAttribute).numValues(), m_fAlpha);
            }
        }
    } // initCPTs

    /**
     * Calculates the class membership probabilities for the given test
     * instance.
     * 
     * @param bayesNet the bayes net to use
     * @param instance the instance to be classified
     * @return predicted class probability distribution
     * @throws Exception if there is a problem generating the prediction
     */
    public double[] distributionForInstance(BayesNet bayesNet, Instance instance) throws Exception {
        Instances instances = bayesNet.m_Instances;
        int nNumClasses = instances.numClasses();
        double[] fProbs = new double[nNumClasses];

        for (int iClass = 0; iClass < nNumClasses; iClass++) {
            fProbs[iClass] = 1.0;
        }

        for (int iClass = 0; iClass < nNumClasses; iClass++) {
            double logfP = 0;

            for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
                double iCPT = 0;

                for (int iParent = 0; iParent < bayesNet.getParentSet(iAttribute).getNrOfParents(); iParent++) {
                    int nParent = bayesNet.getParentSet(iAttribute).getParent(iParent);

                    if (nParent == instances.classIndex()) {
                        iCPT = iCPT * nNumClasses + iClass;
                    } else {
                        iCPT = iCPT * instances.attribute(nParent).numValues() + instance.value(nParent);
                    }
                }

                if (iAttribute == instances.classIndex()) {
                    //	  fP *= 
                    //	    m_Distributions[iAttribute][(int) iCPT].getProbability(iClass);
                    logfP += Math.log(bayesNet.m_Distributions[iAttribute][(int) iCPT].getProbability(iClass));
                } else {
                    //	  fP *= 
                    //	    m_Distributions[iAttribute][(int) iCPT]
                    //	      .getProbability(instance.value(iAttribute));
                    logfP
                        += Math.log(bayesNet.m_Distributions[iAttribute][(int) iCPT].getProbability(instance.value(iAttribute)));
                }
            }

            //      fProbs[iClass] *= fP;
            fProbs[iClass] += logfP;
        }

        // Find maximum
        double fMax = fProbs[0];
        for (int iClass = 0; iClass < nNumClasses; iClass++) {
            if (fProbs[iClass] > fMax) {
                fMax = fProbs[iClass];
            }
        }
        // transform from log-space to normal-space
        for (int iClass = 0; iClass < nNumClasses; iClass++) {
            fProbs[iClass] = Math.exp(fProbs[iClass] - fMax);
        }

        // Display probabilities
        Utils.normalize(fProbs);

        return fProbs;
    } // distributionForInstance
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 1.6 $");
    }

} // SimpleEstimator

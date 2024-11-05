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

package weka.classifiers.bayes.net.estimate;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.search.local.K2;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Statistics;
import weka.core.Utils;
import weka.estimators.Estimator;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Multinomial BMA Estimator.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -k2
 *  Whether to use K2 prior.
 * </pre>
 * 
 * <pre> -A &lt;alpha&gt;
 *  Initial count (alpha)
 * </pre>
 * 
 <!-- options-end -->
 * 
 * @version $Revision: 5987 $
 * @author Remco Bouckaert (rrb@xm.co.nz)
 */
public class MultiNomialBMAEstimator 
    extends BayesNetEstimator {

    /** for serialization */
    static final long serialVersionUID = 8330705772601586313L;
    
    /** whether to use K2 prior */
    protected boolean m_bUseK2Prior = true;
    
    /**
     * Returns a string describing this object
     * @return a description of the classifier suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
      return 
      	  "Multinomial BMA Estimator.";
    }

    /**
     * estimateCPTs estimates the conditional probability tables for the Bayes
     * Net using the network structure.
     * 
     * @param bayesNet the bayes net to use
     * @throws Exception if number of parents doesn't fit (more than 1)
     */
    public void estimateCPTs(BayesNet bayesNet) throws Exception {
        initCPTs(bayesNet);
        
        // sanity check to see if nodes have not more than one parent
        for (int iAttribute = 0; iAttribute < bayesNet.m_Instances.numAttributes(); iAttribute++) {
            if (bayesNet.getParentSet(iAttribute).getNrOfParents() > 1) {
                throw new Exception("Cannot handle networks with nodes with more than 1 parent (yet).");
            }
        }

		// filter data to binary
        Instances instances = new Instances(bayesNet.m_Instances);
        while (instances.numInstances() > 0) {
            instances.delete(0);
        }
        for (int iAttribute = instances.numAttributes() - 1; iAttribute >= 0; iAttribute--) {
            if (iAttribute != instances.classIndex()) {
                FastVector values = new FastVector();
                values.addElement("0");
                values.addElement("1");
                Attribute a = new Attribute(instances.attribute(iAttribute).name(), (FastVector) values);
                instances.deleteAttributeAt(iAttribute);
                instances.insertAttributeAt(a,iAttribute);
            }
        }
        
        for (int iInstance = 0; iInstance < bayesNet.m_Instances.numInstances(); iInstance++) {
            Instance instanceOrig = bayesNet.m_Instances.instance(iInstance);
            Instance instance = new DenseInstance(instances.numAttributes());
            for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
                if (iAttribute != instances.classIndex()) {
                    if (instanceOrig.value(iAttribute) > 0) {
                        instance.setValue(iAttribute, 1);
                    }
                } else {
                    instance.setValue(iAttribute, instanceOrig.value(iAttribute));
                }
            }
        }
        // ok, now all data is binary, except the class attribute
        // now learn the empty and tree network

        BayesNet EmptyNet = new BayesNet();
        K2 oSearchAlgorithm = new K2();
        oSearchAlgorithm.setInitAsNaiveBayes(false);
        oSearchAlgorithm.setMaxNrOfParents(0);
        EmptyNet.setSearchAlgorithm(oSearchAlgorithm);
        EmptyNet.buildClassifier(instances);

        BayesNet NBNet = new BayesNet();
        oSearchAlgorithm.setInitAsNaiveBayes(true);
        oSearchAlgorithm.setMaxNrOfParents(1);
        NBNet.setSearchAlgorithm(oSearchAlgorithm);
        NBNet.buildClassifier(instances);

        // estimate CPTs
        for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
            if (iAttribute != instances.classIndex()) {
                  double w1 = 0.0, w2 = 0.0;
                  int nAttValues = instances.attribute(iAttribute).numValues();
                  if (m_bUseK2Prior == true) {
                      // use Cooper and Herskovitz's metric
                      for (int iAttValue = 0; iAttValue < nAttValues; iAttValue++) {
                        w1 += Statistics.lnGamma(1 + ((DiscreteEstimatorBayes)EmptyNet.m_Distributions[iAttribute][0]).getCount(iAttValue))
                              - Statistics.lnGamma(1);
                      }
                      w1 += Statistics.lnGamma(nAttValues) - Statistics.lnGamma(nAttValues + instances.numInstances());

                      for (int iParent = 0; iParent < bayesNet.getParentSet(iAttribute).getCardinalityOfParents(); iParent++) {
                        int nTotal = 0;
                          for (int iAttValue = 0; iAttValue < nAttValues; iAttValue++) {
                            double nCount = ((DiscreteEstimatorBayes)NBNet.m_Distributions[iAttribute][iParent]).getCount(iAttValue);
                            w2 += Statistics.lnGamma(1 + nCount)
                                  - Statistics.lnGamma(1);
                            nTotal += nCount;
                          }
                        w2 += Statistics.lnGamma(nAttValues) - Statistics.lnGamma(nAttValues + nTotal);
                      }
                  } else {
                      // use BDe metric
                      for (int iAttValue = 0; iAttValue < nAttValues; iAttValue++) {
                        w1 += Statistics.lnGamma(1.0/nAttValues + ((DiscreteEstimatorBayes)EmptyNet.m_Distributions[iAttribute][0]).getCount(iAttValue))
                              - Statistics.lnGamma(1.0/nAttValues);
                      }
                      w1 += Statistics.lnGamma(1) - Statistics.lnGamma(1 + instances.numInstances());

					  int nParentValues = bayesNet.getParentSet(iAttribute).getCardinalityOfParents();
                      for (int iParent = 0; iParent < nParentValues; iParent++) {
                        int nTotal = 0;
                          for (int iAttValue = 0; iAttValue < nAttValues; iAttValue++) {
                            double nCount = ((DiscreteEstimatorBayes)NBNet.m_Distributions[iAttribute][iParent]).getCount(iAttValue);
                            w2 += Statistics.lnGamma(1.0/(nAttValues * nParentValues) + nCount)
                                  - Statistics.lnGamma(1.0/(nAttValues * nParentValues));
                            nTotal += nCount;
                          }
                        w2 += Statistics.lnGamma(1) - Statistics.lnGamma(1 + nTotal);
                      }
                  }
		
//    System.out.println(w1 + " " + w2 + " " + (w2 - w1));
                  // normalize weigths
                  if (w1 < w2) {
                    w2 = w2 - w1;
                    w1 = 0;
                    w1 = 1 / (1 + Math.exp(w2));
                    w2 = Math.exp(w2) / (1 + Math.exp(w2));
                  } else {
                    w1 = w1 - w2;
                    w2 = 0;
                    w2 = 1 / (1 + Math.exp(w1));
                    w1 = Math.exp(w1) / (1 + Math.exp(w1));
                  }
		
                  for (int iParent = 0; iParent < bayesNet.getParentSet(iAttribute).getCardinalityOfParents(); iParent++) {
                      bayesNet.m_Distributions[iAttribute][iParent] = 
                      new DiscreteEstimatorFullBayes(
                        instances.attribute(iAttribute).numValues(), 
                        w1, w2,
                        (DiscreteEstimatorBayes) EmptyNet.m_Distributions[iAttribute][0],
                        (DiscreteEstimatorBayes) NBNet.m_Distributions[iAttribute][iParent],
                        m_fAlpha
                       );
                  } 
            }
        }
        int iAttribute = instances.classIndex();
        bayesNet.m_Distributions[iAttribute][0] = EmptyNet.m_Distributions[iAttribute][0];
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
        throw new Exception("updateClassifier does not apply to BMA estimator");
    } // updateClassifier

    /** 
     * initCPTs reserves space for CPTs and set all counts to zero
     * 
     * @param bayesNet the bayes net to use
     * @throws Exception doesn't apply
     */
    public void initCPTs(BayesNet bayesNet) throws Exception {
        // Reserve sufficient memory
        bayesNet.m_Distributions = new Estimator[bayesNet.m_Instances.numAttributes()][2];
    } // initCPTs


    /**
     * @return boolean
     */
    public boolean isUseK2Prior() {
        return m_bUseK2Prior;
    }

    /**
     * Sets the UseK2Prior.
     * 
     * @param bUseK2Prior The bUseK2Prior to set
     */
    public void setUseK2Prior(boolean bUseK2Prior) {
        m_bUseK2Prior = bUseK2Prior;
    }

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
                    logfP += Math.log(bayesNet.m_Distributions[iAttribute][(int) iCPT].getProbability(iClass));
                } else {
                    logfP += instance.value(iAttribute) * Math.log(
                      bayesNet.m_Distributions[iAttribute][(int) iCPT].getProbability(instance.value(1)));
                }
            }

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
     * Returns an enumeration describing the available options
     * 
     * @return an enumeration of all the available options
     */
    public Enumeration listOptions() {
        Vector newVector = new Vector(1);

        newVector.addElement(new Option(
            "\tWhether to use K2 prior.\n", 
            "k2", 0, "-k2"));

        Enumeration enu = super.listOptions();
        while (enu.hasMoreElements()) {
                newVector.addElement(enu.nextElement());
        }

        return newVector.elements();
    } // listOptions

    /**
     * Parses a given list of options. <p/>
     *
     <!-- options-start -->
     * Valid options are: <p/>
     * 
     * <pre> -k2
     *  Whether to use K2 prior.
     * </pre>
     * 
     * <pre> -A &lt;alpha&gt;
     *  Initial count (alpha)
     * </pre>
     * 
     <!-- options-end -->
     * 
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    public void setOptions(String[] options) throws Exception {
        setUseK2Prior(Utils.getFlag("k2", options));

        super.setOptions(options);
    } // setOptions

    /**
     * Gets the current settings of the classifier.
     * 
     * @return an array of strings suitable for passing to setOptions
     */
    public String[] getOptions() {
        String[] superOptions = super.getOptions();
        String[] options = new String[1 + superOptions.length];
        int current = 0;

        if (isUseK2Prior())
          options[current++] = "-k2";

        // insert options from parent class
        for (int iOption = 0; iOption < superOptions.length; iOption++) {
                options[current++] = superOptions[iOption];
        }

        // Fill up rest with empty strings, not nulls!
        while (current < options.length) {
                options[current++] = "";
        }

        return options;
    } // getOptions
    
    /**
     * Returns the revision string.
     * 
     * @return		the revision
     */
    public String getRevision() {
      return RevisionUtils.extract("$Revision: 5987 $");
    }
} // class MultiNomialBMAEstimator

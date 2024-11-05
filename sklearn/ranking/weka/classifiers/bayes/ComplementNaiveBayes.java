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
 *    ComplementNaiveBayes.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.bayes;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;


/**
 <!-- globalinfo-start -->
 * Class for building and using a Complement class Naive Bayes classifier.<br/>
 * <br/>
 * For more information see, <br/>
 * <br/>
 * Jason D. Rennie, Lawrence Shih, Jaime Teevan, David R. Karger: Tackling the Poor Assumptions of Naive Bayes Text Classifiers. In: ICML, 616-623, 2003.<br/>
 * <br/>
 * P.S.: TF, IDF and length normalization transforms, as described in the paper, can be performed through weka.filters.unsupervised.StringToWordVector.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Rennie2003,
 *    author = {Jason D. Rennie and Lawrence Shih and Jaime Teevan and David R. Karger},
 *    booktitle = {ICML},
 *    pages = {616-623},
 *    publisher = {AAAI Press},
 *    title = {Tackling the Poor Assumptions of Naive Bayes Text Classifiers},
 *    year = {2003}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N
 *  Normalize the word weights for each class
 * </pre>
 * 
 * <pre> -S
 *  Smoothing value to avoid zero WordGivenClass probabilities (default=1.0).
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz)
 * @version $Revision: 5928 $ 
 */
public class ComplementNaiveBayes extends AbstractClassifier
    implements OptionHandler, WeightedInstancesHandler, TechnicalInformationHandler {
    
    /** for serialization */
    static final long serialVersionUID = 7246302925903086397L;
  
    /**
      Weight of words for each class. The weight is actually the
      log of the probability of a word (w) given a class (c) 
      (i.e. log(Pr[w|c])). The format of the matrix is: 
      wordWeights[class][wordAttribute]
    */
    private double[][] wordWeights;
    
    /** Holds the smoothing value to avoid word probabilities of zero.<br>
        P.S.: According to the paper this is the Alpha i parameter 
     */
    private double smoothingParameter = 1.0;
    
    /** True if the words weights are to be normalized */
    private boolean m_normalizeWordWeights = false;
    
    /** Holds the number of Class values present in the set of specified 
        instances */
    private int numClasses;
    
    /** The instances header that'll be used in toString */
    private Instances header;

    
    /**
     * Returns an enumeration describing the available options.
     *
     * @return an enumeration of all the available options.
     */
    public java.util.Enumeration listOptions() {
        FastVector newVector = new FastVector(2);
        newVector.addElement(
        new Option("\tNormalize the word weights for each class\n",
                   "N", 0,"-N"));
        newVector.addElement(
        new Option("\tSmoothing value to avoid zero WordGivenClass"+
                   " probabilities (default=1.0).\n",
                   "S", 1,"-S"));
        
        return newVector.elements();
    }
    
    /**
     * Gets the current settings of the classifier.
     *
     * @return an array of strings suitable for passing to setOptions
     */
    public String[] getOptions() {
        String options[] = new String[4];
        int current=0;
        
        if(getNormalizeWordWeights())
            options[current++] = "-N";
        
        options[current++] = "-S";
        options[current++] = Double.toString(smoothingParameter);
        
        while (current < options.length) {
            options[current++] = "";
        }
        
        return options;
    }        

    /**
     * Parses a given list of options. <p/>
     *
     <!-- options-start -->
     * Valid options are: <p/>
     * 
     * <pre> -N
     *  Normalize the word weights for each class
     * </pre>
     * 
     * <pre> -S
     *  Smoothing value to avoid zero WordGivenClass probabilities (default=1.0).
     * </pre>
     * 
     <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    public void setOptions(String[] options) throws Exception {
        
        setNormalizeWordWeights(Utils.getFlag('N', options));
        
        String val = Utils.getOption('S', options);
        if(val.length()!=0)
          setSmoothingParameter(Double.parseDouble(val));
        else
          setSmoothingParameter(1.0);
    }
    
    /**
     * Returns true if the word weights for each class are to be normalized
     * 
     * @return true if the word weights are normalized
     */
    public boolean getNormalizeWordWeights() {
        return m_normalizeWordWeights;
    }
    
    /**
     * Sets whether if the word weights for each class should be normalized
     * 
     * @param doNormalize whether the word weights are to be normalized
     */
    public void setNormalizeWordWeights(boolean doNormalize) {
        m_normalizeWordWeights = doNormalize;
    }
    
    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String normalizeWordWeightsTipText() {
        return "Normalizes the word weights for each class.";
    }
    
    /**
     * Gets the smoothing value to be used to avoid zero WordGivenClass
     * probabilities.
     * 
     * @return the smoothing value
     */
    public double getSmoothingParameter() {
        return smoothingParameter;
    }

    /**
     * Sets the smoothing value used to avoid zero WordGivenClass probabilities
     * 
     * @param val the new smooting value
     */
    public void setSmoothingParameter(double val) {
        smoothingParameter = val;
    }
        
    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String smoothingParameterTipText() {
        return "Sets the smoothing parameter to avoid zero WordGivenClass "+
               "probabilities (default=1.0).";
    }

    /**
     * Returns a string describing this classifier
     * @return a description of the classifier suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
        
        return  "Class for building and using a Complement class Naive Bayes "+
                "classifier.\n\nFor more information see, \n\n"+
                getTechnicalInformation().toString() + "\n\n" +
                "P.S.: TF, IDF and length normalization transforms, as "+
                "described in the paper, can be performed through "+
                "weka.filters.unsupervised.StringToWordVector.";
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
      
      result = new TechnicalInformation(Type.INPROCEEDINGS);
      result.setValue(Field.AUTHOR, "Jason D. Rennie and Lawrence Shih and Jaime Teevan and David R. Karger");
      result.setValue(Field.TITLE, "Tackling the Poor Assumptions of Naive Bayes Text Classifiers");
      result.setValue(Field.BOOKTITLE, "ICML");
      result.setValue(Field.YEAR, "2003");
      result.setValue(Field.PAGES, "616-623");
      result.setValue(Field.PUBLISHER, "AAAI Press");
      
      return result;
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
      result.enable(Capability.NUMERIC_ATTRIBUTES);
      result.enable(Capability.MISSING_VALUES);

      // class
      result.enable(Capability.NOMINAL_CLASS);
      result.enable(Capability.MISSING_CLASS_VALUES);
      
      return result;
    }
    
    /**
     * Generates the classifier.
     *
     * @param instances set of instances serving as training data 
     * @throws Exception if the classifier has not been built successfully
     */
    public void buildClassifier(Instances instances) throws Exception {

      // can classifier handle the data?
      getCapabilities().testWithFail(instances);

      // remove instances with missing class
      instances = new Instances(instances);
      instances.deleteWithMissingClass();
      
        numClasses = instances.numClasses();
	int numAttributes = instances.numAttributes();
        
        header = new Instances(instances, 0);
	double [][] ocrnceOfWordInClass = new double[numClasses][numAttributes];        
        wordWeights = new double[numClasses][numAttributes];
        //double [] docsPerClass = new double[numClasses];
	double[] wordsPerClass = new double[numClasses];
        double totalWordOccurrences = 0;
        double sumOfSmoothingParams = (numAttributes-1)*smoothingParameter;
        int classIndex = instances.instance(0).classIndex();        
	Instance instance;
	int docClass;
	double numOccurrences;        
        
        java.util.Enumeration enumInsts = instances.enumerateInstances();
	while (enumInsts.hasMoreElements()) {
		instance = (Instance) enumInsts.nextElement();
		docClass = (int)instance.value(classIndex);
		//docsPerClass[docClass] += instance.weight();
		
		for(int a = 0; a<instance.numValues(); a++)
		    if(instance.index(a) != instance.classIndex()) {
			    if(!instance.isMissing(a)) {
				    numOccurrences = instance.valueSparse(a) * instance.weight();
				    if(numOccurrences < 0)
					throw new Exception("Numeric attribute"+
                                                  " values must all be greater"+
                                                  " or equal to zero.");
                                    totalWordOccurrences += numOccurrences;
				    wordsPerClass[docClass] += numOccurrences;
				    ocrnceOfWordInClass[docClass]
                                          [instance.index(a)] += numOccurrences;
                                    //For the time being wordweights[0][i] 
                                    //will hold the total occurrence of word
                                    // i over all classes
                                    wordWeights[0]
                                          [instance.index(a)] += numOccurrences;
                            }
                    }
        }

	//Calculating the complement class probability for all classes except 0        
	for(int c=1; c<numClasses; c++) {
            //total occurrence of words in classes other than c
            double totalWordOcrnces = totalWordOccurrences - wordsPerClass[c];

            for(int w=0; w<numAttributes; w++) {
                if(w != classIndex ) {
                     //occurrence of w in classes other that c
                    double ocrncesOfWord = 
                                wordWeights[0][w] - ocrnceOfWordInClass[c][w];

                    wordWeights[c][w] = 
                        Math.log((ocrncesOfWord+smoothingParameter) / 
                                (totalWordOcrnces+sumOfSmoothingParams));
                }
            }
        }
        
	//Now calculating the complement class probability for class 0
        for(int w=0; w<numAttributes; w++) {
            if(w != classIndex) {
                //occurrence of w in classes other that c
                double ocrncesOfWord = wordWeights[0][w] - ocrnceOfWordInClass[0][w];
                //total occurrence of words in classes other than c
                double totalWordOcrnces = totalWordOccurrences - wordsPerClass[0];
                
                wordWeights[0][w] =
                Math.log((ocrncesOfWord+smoothingParameter) /
                (totalWordOcrnces+sumOfSmoothingParams));
            }            
        }
        
       	//Normalizing weights
        if(m_normalizeWordWeights==true)
            for(int c=0; c<numClasses; c++) {
                double sum=0;
                for(int w=0; w<numAttributes; w++) {
                    if(w!=classIndex)
                        sum += Math.abs(wordWeights[c][w]);
                }
                for(int w=0; w<numAttributes; w++) {
                    if(w!=classIndex) {
                        wordWeights[c][w] = wordWeights[c][w]/sum;
                    }
                }
            }

    }

    
    /**
     * Classifies a given instance. <p>
     *
     * The classification rule is: <br>
     *     MinC(forAllWords(ti*Wci)) <br>
     *      where <br>
     *         ti is the frequency of word i in the given instance <br>
     *         Wci is the weight of word i in Class c. <p>
     *
     * For more information see section 4.4 of the paper mentioned above
     * in the classifiers description.
     *
     * @param instance the instance to classify
     * @return the index of the class the instance is most likely to belong.
     * @throws Exception if the classifier has not been built yet.
     */
    public double classifyInstance(Instance instance) throws Exception {

        if(wordWeights==null)
            throw new Exception("Error. The classifier has not been built "+
                                "properly.");
        
        double [] valueForClass = new double[numClasses];
	double sumOfClassValues=0;
	
	for(int c=0; c<numClasses; c++) {
	    double sumOfWordValues=0;
	    for(int w=0; w<instance.numValues(); w++) {
                if(instance.index(w)!=instance.classIndex()) {
                    double freqOfWordInDoc = instance.valueSparse(w);
                    sumOfWordValues += freqOfWordInDoc * 
                                  wordWeights[c][instance.index(w)];
                }
	    }
	    //valueForClass[c] = Math.log(probOfClass[c]) - sumOfWordValues;
	    valueForClass[c] = sumOfWordValues;
	    sumOfClassValues += valueForClass[c];
	}

        int minidx=0;
	for(int i=0; i<numClasses; i++)
	    if(valueForClass[i]<valueForClass[minidx])
		minidx = i;
	
	return minidx;
    }


    /**
     * Prints out the internal model built by the classifier. In this case
     * it prints out the word weights calculated when building the classifier.
     */
    public String toString() {
        if(wordWeights==null) {            
            return "The classifier hasn't been built yet.";
        }
        
        int numAttributes = header.numAttributes();
        StringBuffer result = new StringBuffer("The word weights for each class are: \n"+
                                               "------------------------------------\n\t");
        
        for(int c = 0; c<numClasses; c++)
            result.append(header.classAttribute().value(c)).append("\t");
        
        result.append("\n");
        
        for(int w = 0; w<numAttributes; w++) {
            result.append(header.attribute(w).name()).append("\t");
            for(int c = 0; c<numClasses; c++)
                result.append(Double.toString(wordWeights[c][w])).append("\t");
            result.append("\n");
        }
        
        return result.toString();
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
      runClassifier(new ComplementNaiveBayes(), argv);
    }        
}


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
 *    SimpleLogistic.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions;

import weka.classifiers.Classifier;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.trees.lmt.LogisticBase;
import weka.core.AdditionalMeasureProducer;
import weka.core.Capabilities;
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
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Classifier for building linear logistic regression models. LogitBoost with simple regression functions as base learners is used for fitting the logistic models. The optimal number of LogitBoost iterations to perform is cross-validated, which leads to automatic attribute selection. For more information see:<br/>
 * Niels Landwehr, Mark Hall, Eibe Frank (2005). Logistic Model Trees.<br/>
 * <br/>
 * Marc Sumner, Eibe Frank, Mark Hall: Speeding up Logistic Model Tree Induction. In: 9th European Conference on Principles and Practice of Knowledge Discovery in Databases, 675-683, 2005.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Landwehr2005,
 *    author = {Niels Landwehr and Mark Hall and Eibe Frank},
 *    booktitle = {Machine Learning},
 *    number = {1-2},
 *    pages = {161-205},
 *    title = {Logistic Model Trees},
 *    volume = {95},
 *    year = {2005}
 * }
 * 
 * &#64;inproceedings{Sumner2005,
 *    author = {Marc Sumner and Eibe Frank and Mark Hall},
 *    booktitle = {9th European Conference on Principles and Practice of Knowledge Discovery in Databases},
 *    pages = {675-683},
 *    publisher = {Springer},
 *    title = {Speeding up Logistic Model Tree Induction},
 *    year = {2005}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -I &lt;iterations&gt;
 *  Set fixed number of iterations for LogitBoost</pre>
 * 
 * <pre> -S
 *  Use stopping criterion on training set (instead of
 *  cross-validation)</pre>
 * 
 * <pre> -P
 *  Use error on probabilities (rmse) instead of
 *  misclassification error for stopping criterion</pre>
 * 
 * <pre> -M &lt;iterations&gt;
 *  Set maximum number of boosting iterations</pre>
 * 
 * <pre> -H &lt;iterations&gt;
 *  Set parameter for heuristic for early stopping of
 *  LogitBoost.
 *  If enabled, the minimum is selected greedily, stopping
 *  if the current minimum has not changed for iter iterations.
 *  By default, heuristic is enabled with value 50. Set to
 *  zero to disable heuristic.</pre>
 * 
 * <pre> -W &lt;beta&gt;
 *  Set beta for weight trimming for LogitBoost. Set to 0 for no weight trimming.
 * </pre>
 * 
 * <pre> -A
 *  The AIC is used to choose the best iteration (instead of CV or training error).
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Niels Landwehr 
 * @author Marc Sumner 
 * @version $Revision: 5928 $
 */
public class SimpleLogistic 
  extends AbstractClassifier 
  implements OptionHandler, AdditionalMeasureProducer, WeightedInstancesHandler,
             TechnicalInformationHandler {

    /** for serialization */
    static final long serialVersionUID = 7397710626304705059L;
  
    /**The actual logistic regression model */
    protected LogisticBase m_boostedModel;
    
    /**Filter for converting nominal attributes to binary ones*/
    protected NominalToBinary m_NominalToBinary = null;

    /**Filter for replacing missing values*/
    protected ReplaceMissingValues m_ReplaceMissingValues = null;
    
    /**If non-negative, use this as fixed number of LogitBoost iterations*/ 
    protected int m_numBoostingIterations;
    
    /**Maximum number of iterations for LogitBoost*/
    protected int m_maxBoostingIterations = 500;
    
    /**Parameter for the heuristic for early stopping of LogitBoost*/
    protected int m_heuristicStop = 50;

    /**If true, cross-validate number of LogitBoost iterations*/
    protected boolean m_useCrossValidation;

    /**If true, use minimize error on probabilities instead of misclassification error*/
    protected boolean m_errorOnProbabilities;
    
    /**Threshold for trimming weights. Instances with a weight lower than this (as a percentage
     * of total weights) are not included in the regression fit.
     */
    protected double m_weightTrimBeta = 0;
    
    /** If true, the AIC is used to choose the best iteration*/
    private boolean m_useAIC = false;

    /**
     * Constructor for creating SimpleLogistic object with standard options.
     */
    public SimpleLogistic() {
	m_numBoostingIterations = 0;
	m_useCrossValidation = true;
	m_errorOnProbabilities = false;
        m_weightTrimBeta = 0;
        m_useAIC = false;
    }

    /**
     * Constructor for creating SimpleLogistic object.
     * @param numBoostingIterations if non-negative, use this as fixed number of iterations for LogitBoost
     * @param useCrossValidation cross-validate number of LogitBoost iterations.
     * @param errorOnProbabilities minimize error on probabilities instead of misclassification error
     */
    public SimpleLogistic(int numBoostingIterations, boolean useCrossValidation, 
			      boolean errorOnProbabilities) { 
  	m_numBoostingIterations = numBoostingIterations;
	m_useCrossValidation = useCrossValidation;
	m_errorOnProbabilities = errorOnProbabilities;
        m_weightTrimBeta = 0;
        m_useAIC = false;
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
      result.enable(Capability.MISSING_CLASS_VALUES);
      
      return result;
    }

    /**
     * Builds the logistic regression using LogitBoost.
     * @param data the training data
     * @throws Exception if something goes wrong 
     */
    public void buildClassifier(Instances data) throws Exception {

      // can classifier handle the data?
      getCapabilities().testWithFail(data);

      // remove instances with missing class
      data = new Instances(data);
      data.deleteWithMissingClass();

	//replace missing values
	m_ReplaceMissingValues = new ReplaceMissingValues();
	m_ReplaceMissingValues.setInputFormat(data);
	data = Filter.useFilter(data, m_ReplaceMissingValues);
	
	//convert nominal attributes
	m_NominalToBinary = new NominalToBinary();
	m_NominalToBinary.setInputFormat(data);
	data = Filter.useFilter(data, m_NominalToBinary);
	
	//create actual logistic model
	m_boostedModel = new LogisticBase(m_numBoostingIterations, m_useCrossValidation, m_errorOnProbabilities);
	m_boostedModel.setMaxIterations(m_maxBoostingIterations);
	m_boostedModel.setHeuristicStop(m_heuristicStop);
        m_boostedModel.setWeightTrimBeta(m_weightTrimBeta);
        m_boostedModel.setUseAIC(m_useAIC);
	
	//build logistic model
	m_boostedModel.buildClassifier(data);
    }
    
    /** 
     * Returns class probabilities for an instance.
     *
     * @param inst the instance to compute the probabilities for
     * @return the probabilities
     * @throws Exception if distribution can't be computed successfully
     */
    public double[] distributionForInstance(Instance inst) 
	throws Exception {
	
	//replace missing values / convert nominal atts
	m_ReplaceMissingValues.input(inst);
	inst = m_ReplaceMissingValues.output();
	m_NominalToBinary.input(inst);
	inst = m_NominalToBinary.output();	
	
	//obtain probs from logistic model
	return m_boostedModel.distributionForInstance(inst);	
    }

    /**
     * Returns an enumeration describing the available options.
     *
     * @return an enumeration of all the available options.
     */
    public Enumeration listOptions() {
	Vector newVector = new Vector();
	
	newVector.addElement(new Option(
	    "\tSet fixed number of iterations for LogitBoost",
	    "I",1,"-I <iterations>"));
	
	newVector.addElement(new Option(
	    "\tUse stopping criterion on training set (instead of\n"
	    + "\tcross-validation)",
	    "S",0,"-S"));
	
	newVector.addElement(new Option(
	    "\tUse error on probabilities (rmse) instead of\n"
	    + "\tmisclassification error for stopping criterion",
	    "P",0,"-P"));

	newVector.addElement(new Option(
	    "\tSet maximum number of boosting iterations",
	    "M",1,"-M <iterations>"));

	newVector.addElement(new Option(
	    "\tSet parameter for heuristic for early stopping of\n"
	    + "\tLogitBoost.\n"
	    + "\tIf enabled, the minimum is selected greedily, stopping\n"
	    + "\tif the current minimum has not changed for iter iterations.\n"
	    + "\tBy default, heuristic is enabled with value 50. Set to\n"
	    + "\tzero to disable heuristic.",
	    "H",1,"-H <iterations>"));
        
        newVector.addElement(new Option("\tSet beta for weight trimming for LogitBoost. Set to 0 for no weight trimming.\n",
                                        "W",1,"-W <beta>"));
        
        newVector.addElement(new Option("\tThe AIC is used to choose the best iteration (instead of CV or training error).\n",
                                        "A", 0, "-A"));
	
	return newVector.elements();
    } 
    

    /**
     * Parses a given list of options. <p/>
     *
     <!-- options-start -->
     * Valid options are: <p/>
     * 
     * <pre> -I &lt;iterations&gt;
     *  Set fixed number of iterations for LogitBoost</pre>
     * 
     * <pre> -S
     *  Use stopping criterion on training set (instead of
     *  cross-validation)</pre>
     * 
     * <pre> -P
     *  Use error on probabilities (rmse) instead of
     *  misclassification error for stopping criterion</pre>
     * 
     * <pre> -M &lt;iterations&gt;
     *  Set maximum number of boosting iterations</pre>
     * 
     * <pre> -H &lt;iterations&gt;
     *  Set parameter for heuristic for early stopping of
     *  LogitBoost.
     *  If enabled, the minimum is selected greedily, stopping
     *  if the current minimum has not changed for iter iterations.
     *  By default, heuristic is enabled with value 50. Set to
     *  zero to disable heuristic.</pre>
     * 
     * <pre> -W &lt;beta&gt;
     *  Set beta for weight trimming for LogitBoost. Set to 0 for no weight trimming.
     * </pre>
     * 
     * <pre> -A
     *  The AIC is used to choose the best iteration (instead of CV or training error).
     * </pre>
     * 
     <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    public void setOptions(String[] options) throws Exception {

	String optionString = Utils.getOption('I', options);
	if (optionString.length() != 0) {
	    setNumBoostingIterations((new Integer(optionString)).intValue());
	}
		
	setUseCrossValidation(!Utils.getFlag('S', options));
	setErrorOnProbabilities(Utils.getFlag('P', options));
	
	optionString = Utils.getOption('M', options);
	if (optionString.length() != 0) {
	    setMaxBoostingIterations((new Integer(optionString)).intValue());
	}

	optionString = Utils.getOption('H', options);
	if (optionString.length() != 0) {
	    setHeuristicStop((new Integer(optionString)).intValue());
	}
        
        optionString = Utils.getOption('W', options);
        if (optionString.length() != 0) {
            setWeightTrimBeta((new Double(optionString)).doubleValue());
        }
        
        setUseAIC(Utils.getFlag('A', options));        

	Utils.checkForRemainingOptions(options);
    } 

    /**
     * Gets the current settings of the Classifier.
     *
     * @return an array of strings suitable for passing to setOptions
     */
    public String[] getOptions() {
	String[] options = new String[11];
	int current = 0;
		
	options[current++] = "-I"; 
	options[current++] = ""+getNumBoostingIterations();
	
	if (!getUseCrossValidation()) {
	    options[current++] = "-S";
	} 

	if (getErrorOnProbabilities()) {
	    options[current++] = "-P";
	} 

	options[current++] = "-M"; 
	options[current++] = ""+getMaxBoostingIterations();
	
	options[current++] = "-H"; 
	options[current++] = ""+getHeuristicStop();
        
        options[current++] = "-W";
        options[current++] = ""+getWeightTrimBeta();
        
        if (getUseAIC()) {
            options[current++] = "-A";
        }

	while (current < options.length) {
	    options[current++] = "";
	} 
	return options;
    } 

    /**
     * Get the value of numBoostingIterations.
     * 
     * @return the number of boosting iterations
     */
    public int getNumBoostingIterations(){
	return m_numBoostingIterations;
    }
    /**
     * Get the value of useCrossValidation.
     * 
     * @return true if cross-validation is used
     */
    public boolean getUseCrossValidation(){
	return m_useCrossValidation;
    }

    /**
     * Get the value of errorOnProbabilities.
     * 
     * @return 	If true, use minimize error on probabilities instead of 
     * 		misclassification error
     */
    public boolean getErrorOnProbabilities(){
	return m_errorOnProbabilities;
    }
    
    /**
     * Get the value of maxBoostingIterations.
     * 
     * @return the maximum number of boosting iterations
     */
    public int getMaxBoostingIterations(){
	return m_maxBoostingIterations;
    }

    /**
     * Get the value of heuristicStop.
     * 
     * @return the value of heuristicStop
     */
    public int getHeuristicStop(){
	return m_heuristicStop;
    }
    
    /**
     * Get the value of weightTrimBeta.
     */
    public double getWeightTrimBeta(){
        return m_weightTrimBeta;
    }
    
    /**
     * Get the value of useAIC.
     *
     * @return Value of useAIC.
     */
    public boolean getUseAIC(){
        return m_useAIC;
    }
    
    /**
     * Set the value of numBoostingIterations.
     * 
     * @param n the number of boosting iterations
     */
    public void setNumBoostingIterations(int n){
	m_numBoostingIterations = n;
    }

    /**
     * Set the value of useCrossValidation.
     * 
     * @param l whether to use cross-validation
     */
    public void setUseCrossValidation(boolean l){
	m_useCrossValidation = l;
    }

    /**
     * Set the value of errorOnProbabilities.
     * 
     * @param l If true, use minimize error on probabilities instead of 
     * 		misclassification error
     */
    public void setErrorOnProbabilities(boolean l){
	m_errorOnProbabilities = l;
    }

    /**
     * Set the value of maxBoostingIterations.
     * 
     * @param n the maximum number of boosting iterations
     */
    public void setMaxBoostingIterations(int n){
	m_maxBoostingIterations = n;
    } 

    /**
     * Set the value of heuristicStop.
     * 
     * @param n the value of heuristicStop
     */
    public void setHeuristicStop(int n){
	if (n == 0) 
	  m_heuristicStop = m_maxBoostingIterations; 
	else 
	  m_heuristicStop = n;
    }
    
    /**
     * Set the value of weightTrimBeta.
     */
    public void setWeightTrimBeta(double n){
        m_weightTrimBeta = n;
    }
    
    /**
     * Set the value of useAIC.
     *
     * @param c Value to assign to useAIC.
     */
    public void setUseAIC(boolean c){
        m_useAIC = c;
    }

    /**
     * Get the number of LogitBoost iterations performed (= the number of 
     * regression functions fit by LogitBoost).
     * 
     * @return the number of LogitBoost iterations performed
     */
    public int getNumRegressions(){
	return m_boostedModel.getNumRegressions();
    }

    /**
     * Returns a description of the logistic model (attributes/coefficients).
     * 
     * @return the model as string
     */
    public String toString(){
	if (m_boostedModel == null) return "No model built";
	return "SimpleLogistic:\n" + m_boostedModel.toString();
    }

    /**
     * Returns the fraction of all attributes in the data that are used in the 
     * logistic model (in percent). An attribute is used in the model if it is 
     * used in any of the models for the different classes.
     * 
     * @return percentage of attributes used in the model
     */
    public double measureAttributesUsed(){
	return m_boostedModel.percentAttributesUsed();
    }
       
     /**
     * Returns an enumeration of the additional measure names
     * @return an enumeration of the measure names
     */
    public Enumeration enumerateMeasures() {
	Vector newVector = new Vector(3);
	newVector.addElement("measureAttributesUsed");
	newVector.addElement("measureNumIterations");
	return newVector.elements();
    }
    
    /**
     * Returns the value of the named measure
     * @param additionalMeasureName the name of the measure to query for its value
     * @return the value of the named measure
     * @throws IllegalArgumentException if the named measure is not supported
     */
    public double getMeasure(String additionalMeasureName) {
	if (additionalMeasureName.compareToIgnoreCase("measureAttributesUsed") == 0) {
	    return measureAttributesUsed();
      	} else if(additionalMeasureName.compareToIgnoreCase("measureNumIterations") == 0){
	    return getNumRegressions();
	} else {
	    throw new IllegalArgumentException(additionalMeasureName 
					       + " not supported (SimpleLogistic)");
	}
    }    


    /**
     * Returns a string describing classifier
     * @return a description suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "Classifier for building linear logistic regression models. LogitBoost with simple regression "
	    +"functions as base learners is used for fitting the logistic models. The optimal number of LogitBoost "
	    +"iterations to perform is cross-validated, which leads to automatic attribute selection. "
	    +"For more information see:\n"
	    + getTechnicalInformation().toString();
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
      TechnicalInformation 	additional;
      
      result = new TechnicalInformation(Type.ARTICLE);
      result.setValue(Field.AUTHOR, "Niels Landwehr and Mark Hall and Eibe Frank");
      result.setValue(Field.TITLE, "Logistic Model Trees");
      result.setValue(Field.BOOKTITLE, "Machine Learning");
      result.setValue(Field.YEAR, "2005");
      result.setValue(Field.VOLUME, "95");
      result.setValue(Field.PAGES, "161-205");
      result.setValue(Field.NUMBER, "1-2");
      
      additional = result.add(Type.INPROCEEDINGS);
      additional.setValue(Field.AUTHOR, "Marc Sumner and Eibe Frank and Mark Hall");
      additional.setValue(Field.TITLE, "Speeding up Logistic Model Tree Induction");
      additional.setValue(Field.BOOKTITLE, "9th European Conference on Principles and Practice of Knowledge Discovery in Databases");
      additional.setValue(Field.YEAR, "2005");
      additional.setValue(Field.PAGES, "675-683");
      additional.setValue(Field.PUBLISHER, "Springer");
      
      return result;
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String numBoostingIterationsTipText() {
	return "Set fixed number of iterations for LogitBoost. If >= 0, this sets the number of LogitBoost iterations "
	    +"to perform. If < 0, the number is cross-validated or a stopping criterion on the training set is used "
	    +"(depending on the value of useCrossValidation).";
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String useCrossValidationTipText() {
	return "Sets whether the number of LogitBoost iterations is to be cross-validated or the stopping criterion "
	    +"on the training set should be used. If not set (and no fixed number of iterations was given), "
	    +"the number of LogitBoost iterations is used that minimizes the error on the training set "
	    +"(misclassification error or error on probabilities depending on errorOnProbabilities).";
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String errorOnProbabilitiesTipText() {
	return "Use error on the probabilties as error measure when determining the best number of LogitBoost iterations. "
	    +"If set, the number of LogitBoost iterations is chosen that minimizes the root mean squared error "
	    +"(either on the training set or in the cross-validation, depending on useCrossValidation).";
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String maxBoostingIterationsTipText() {
	return "Sets the maximum number of iterations for LogitBoost. Default value is 500, for very small/large "
	    +"datasets a lower/higher value might be preferable.";
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String heuristicStopTipText() {
	return "If heuristicStop > 0, the heuristic for greedy stopping while cross-validating the number of "
	    +"LogitBoost iterations is enabled. This means LogitBoost is stopped if no new error minimum "
	    +"has been reached in the last heuristicStop iterations. It is recommended to use this heuristic, "
	    +"it gives a large speed-up especially on small datasets. The default value is 50.";
    }    
    
    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String weightTrimBetaTipText() {
        return "Set the beta value used for weight trimming in LogitBoost. "
        +"Only instances carrying (1 - beta)% of the weight from previous iteration "
        +"are used in the next iteration. Set to 0 for no weight trimming. "
        +"The default value is 0.";
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String useAICTipText() {
        return "The AIC is used to determine when to stop LogitBoost iterations "
        +"(instead of cross-validation or training error).";
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
     * Main method for testing this class
     *
     * @param argv commandline options 
     */
    public static void main(String[] argv) {	
        runClassifier(new SimpleLogistic(), argv);
    }
}

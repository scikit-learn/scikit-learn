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
 *    NominalToBinary.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.supervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.UnassignedClassException;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Converts all nominal attributes into binary numeric attributes. An attribute with k values is transformed into k binary attributes if the class is nominal (using the one-attribute-per-value approach). Binary attributes are left binary, if option '-A' is not given.If the class is numeric, k - 1 new binary attributes are generated in the manner described in "Classification and Regression Trees" by Breiman et al. (i.e. taking the average class value associated with each attribute value into account)<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * L. Breiman, J.H. Friedman, R.A. Olshen, C.J. Stone (1984). Classification and Regression Trees. Wadsworth Inc.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;book{Breiman1984,
 *    author = {L. Breiman and J.H. Friedman and R.A. Olshen and C.J. Stone},
 *    publisher = {Wadsworth Inc},
 *    title = {Classification and Regression Trees},
 *    year = {1984},
 *    ISBN = {0412048418}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N
 *  Sets if binary attributes are to be coded as nominal ones.</pre>
 * 
 * <pre> -A
 *  For each nominal value a new attribute is created, 
 *  not only if there are more than 2 values.</pre>
 * 
 <!-- options-end -->
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz) 
 * @version $Revision: 5987 $ 
 */
public class NominalToBinary 
  extends Filter 
  implements SupervisedFilter, OptionHandler, TechnicalInformationHandler {
  
  /** for serialization */
  static final long serialVersionUID = -5004607029857673950L;

  /** The sorted indices of the attribute values. */
  private int[][] m_Indices = null;

  /** Are the new attributes going to be nominal or numeric ones? */
  private boolean m_Numeric = true;

  /** Are all values transformed into new attributes? */
  private boolean m_TransformAll = false;

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Converts all nominal attributes into binary numeric attributes. An "
      + "attribute with k values is transformed into k binary attributes if "
      + "the class is nominal (using the one-attribute-per-value approach). "
      + "Binary attributes are left binary, if option '-A' is not given."
      + "If the class is numeric, k - 1 new binary attributes are generated "
      + "in the manner described in \"Classification and Regression "
      + "Trees\" by Breiman et al. (i.e. taking the average class value associated "
      + "with each attribute value into account)\n\n"
      + "For more information, see:\n\n"
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
    
    result = new TechnicalInformation(Type.BOOK);
    result.setValue(Field.AUTHOR, "L. Breiman and J.H. Friedman and R.A. Olshen and C.J. Stone");
    result.setValue(Field.TITLE, "Classification and Regression Trees");
    result.setValue(Field.YEAR, "1984");
    result.setValue(Field.PUBLISHER, "Wadsworth Inc");
    result.setValue(Field.ISBN, "0412048418");
    
    return result;
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.NOMINAL_CLASS);
    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input 
   * instance structure (any instances contained in the object are 
   * ignored - only the structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the input format can't be set 
   * successfully
   */
  public boolean setInputFormat(Instances instanceInfo) 
       throws Exception {

    super.setInputFormat(instanceInfo);
    if (instanceInfo.classIndex() < 0) {
      throw new UnassignedClassException("No class has been assigned to the instances");
    }
    setOutputFormat();
    m_Indices = null;
    if (instanceInfo.classAttribute().isNominal()) {
      return true;
    } else {
      return false;
    }
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been set
   */
  public boolean input(Instance instance) {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    if ((m_Indices != null) || 
	(getInputFormat().classAttribute().isNominal())) {
      convertInstance(instance);
      return true;
    }
    bufferInput(instance);
    return false;
  }

  /**
   * Signify that this batch of input to the filter is finished. 
   * If the filter requires all instances prior to filtering,
   * output() may now be called to retrieve the filtered instances.
   *
   * @return true if there are instances pending output
   * @throws IllegalStateException if no input structure has been defined
   */
  public boolean batchFinished() {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if ((m_Indices == null) && 
	(getInputFormat().classAttribute().isNumeric())) {
      computeAverageClassValues();
      setOutputFormat();

      // Convert pending input instances

      for(int i = 0; i < getInputFormat().numInstances(); i++) {
	convertInstance(getInputFormat().instance(i));
      }
    } 
    flushInput();

    m_NewBatch = true;
    return (numPendingOutput() != 0);
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(1);

    newVector.addElement(new Option(
	"\tSets if binary attributes are to be coded as nominal ones.",
	"N", 0, "-N"));
    
    newVector.addElement(new Option(
	"\tFor each nominal value a new attribute is created, \n"
	+ "\tnot only if there are more than 2 values.",
	"A", 0, "-A"));

    return newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -N
   *  Sets if binary attributes are to be coded as nominal ones.</pre>
   * 
   * <pre> -A
   *  For each nominal value a new attribute is created, 
   *  not only if there are more than 2 values.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    setBinaryAttributesNominal(Utils.getFlag('N', options));

    setTransformAllValues(Utils.getFlag('A', options));

    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [1];
    int current = 0;

    if (getBinaryAttributesNominal()) {
      options[current++] = "-N";
    }

    if (getTransformAllValues()) {
      options[current++] = "-A";
    }

    while (current < options.length) {
      options[current++] = "";
    }
    return options;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String binaryAttributesNominalTipText() {
    return "Whether resulting binary attributes will be nominal.";
  }

  /**
   * Gets if binary attributes are to be treated as nominal ones.
   *
   * @return true if binary attributes are to be treated as nominal ones
   */
  public boolean getBinaryAttributesNominal() {

    return !m_Numeric;
  }

  /**
   * Sets if binary attributes are to be treates as nominal ones.
   *
   * @param bool true if binary attributes are to be treated as nominal ones
   */
  public void setBinaryAttributesNominal(boolean bool) {

    m_Numeric = !bool;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String transformAllValuesTipText() {
    return "Whether all nominal values are turned into new attributes, not only if there are more than 2.";
  }

  /**
   * Gets if all nominal values are turned into new attributes, not only if
   * there are more than 2.
   *
   * @return true all nominal values are transformed into new attributes
   */
  public boolean getTransformAllValues() {

    return m_TransformAll;
  }

  /**
   * Sets whether all nominal values are transformed into new attributes, not
   * just if there are more than 2.
   *
   * @param bool true if all nominal value are transformed into new attributes
   */
  public void setTransformAllValues(boolean bool) {

    m_TransformAll = bool;
  }

  /** Computes average class values for each attribute and value */
  private void computeAverageClassValues() {

    double totalCounts, sum;
    Instance instance;
    double [] counts;

    double [][] avgClassValues = new double[getInputFormat().numAttributes()][0];
    m_Indices = new int[getInputFormat().numAttributes()][0];
    for (int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if (att.isNominal()) {
	avgClassValues[j] = new double [att.numValues()];
	counts = new double [att.numValues()];
	for (int i = 0; i < getInputFormat().numInstances(); i++) {
	  instance = getInputFormat().instance(i);
	  if (!instance.classIsMissing() && 
	      (!instance.isMissing(j))) {
	    counts[(int)instance.value(j)] += instance.weight();
	    avgClassValues[j][(int)instance.value(j)] += 
	      instance.weight() * instance.classValue();
	  }
	}
	sum = Utils.sum(avgClassValues[j]);
	totalCounts = Utils.sum(counts);
	if (Utils.gr(totalCounts, 0)) {
	  for (int k = 0; k < att.numValues(); k++) {
	    if (Utils.gr(counts[k], 0)) {
	      avgClassValues[j][k] /= (double)counts[k];
	    } else {
	      avgClassValues[j][k] = sum / (double)totalCounts;
	    }
	  }
	}
	m_Indices[j] = Utils.sort(avgClassValues[j]);
      }
    }
  }

  /** Set the output format. */
  private void setOutputFormat() {

    if (getInputFormat().classAttribute().isNominal()) {
      setOutputFormatNominal();
    } else {
      setOutputFormatNumeric();
    }
  }

  /**
   * Convert a single instance over. The converted instance is 
   * added to the end of the output queue.
   *
   * @param instance the instance to convert
   */
  private void convertInstance(Instance inst) {

    if (getInputFormat().classAttribute().isNominal()) {
      convertInstanceNominal(inst);
    } else {
      convertInstanceNumeric(inst);
    }
  }

  /**
   * Set the output format if the class is nominal.
   */
  private void setOutputFormatNominal() {

    FastVector newAtts;
    int newClassIndex;
    StringBuffer attributeName;
    Instances outputFormat;
    FastVector vals;

    // Compute new attributes

    newClassIndex = getInputFormat().classIndex();
    newAtts = new FastVector();
    for (int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if ((!att.isNominal()) || 
	  (j == getInputFormat().classIndex())) {
	newAtts.addElement(att.copy());
      } else {
	if ( (att.numValues() <= 2) && (!m_TransformAll) ) {
	  if (m_Numeric) {
	    newAtts.addElement(new Attribute(att.name()));
	  } else {
	    newAtts.addElement(att.copy());
	  }
	} else {

	  if (j < getInputFormat().classIndex()) {
	    newClassIndex += att.numValues() - 1;
	  }

	  // Compute values for new attributes
	  for (int k = 0; k < att.numValues(); k++) {
	    attributeName = 
	      new StringBuffer(att.name() + "=");
	    attributeName.append(att.value(k));
	    if (m_Numeric) {
	      newAtts.
		addElement(new Attribute(attributeName.toString()));
	    } else {
	      vals = new FastVector(2);
	      vals.addElement("f"); vals.addElement("t");
	      newAtts.
		addElement(new Attribute(attributeName.toString(), vals));
	    }
	  }
	}
      }
    }
    outputFormat = new Instances(getInputFormat().relationName(),
				 newAtts, 0);
    outputFormat.setClassIndex(newClassIndex);
    setOutputFormat(outputFormat);
  }

  /**
   * Set the output format if the class is numeric.
   */
  private void setOutputFormatNumeric() {

    if (m_Indices == null) {
      setOutputFormat(null);
      return;
    }
    FastVector newAtts;
    int newClassIndex;
    StringBuffer attributeName;
    Instances outputFormat;
    FastVector vals;

    // Compute new attributes

    newClassIndex = getInputFormat().classIndex();
    newAtts = new FastVector();
    for (int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if ((!att.isNominal()) || 
	  (j == getInputFormat().classIndex())) {
	newAtts.addElement(att.copy());
      } else {
	if (j < getInputFormat().classIndex())
	  newClassIndex += att.numValues() - 2;
	  
	// Compute values for new attributes
	  
	for (int k = 1; k < att.numValues(); k++) {
	  attributeName = 
	    new StringBuffer(att.name() + "=");
	  for (int l = k; l < att.numValues(); l++) {
	    if (l > k) {
	      attributeName.append(',');
	    }
	    attributeName.append(att.value(m_Indices[j][l]));
	  }
	  if (m_Numeric) {
	    newAtts.
	      addElement(new Attribute(attributeName.toString()));
	  } else {
	    vals = new FastVector(2);
	    vals.addElement("f"); vals.addElement("t");
	    newAtts.
	      addElement(new Attribute(attributeName.toString(), vals));
	  }
	}
      }
    }
    outputFormat = new Instances(getInputFormat().relationName(),
				 newAtts, 0);
    outputFormat.setClassIndex(newClassIndex);
    setOutputFormat(outputFormat);
  }

  /**
   * Convert a single instance over if the class is nominal. The converted 
   * instance is added to the end of the output queue.
   *
   * @param instance the instance to convert
   */
  private void convertInstanceNominal(Instance instance) {

    double [] vals = new double [outputFormatPeek().numAttributes()];
    int attSoFar = 0;

    for(int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if ((!att.isNominal()) || (j == getInputFormat().classIndex())) {
	vals[attSoFar] = instance.value(j);
	attSoFar++;
      } else {
	if ( (att.numValues() <= 2) && (!m_TransformAll) ) {
	  vals[attSoFar] = instance.value(j);
	  attSoFar++;
	} else {
	  if (instance.isMissing(j)) {
	    for (int k = 0; k < att.numValues(); k++) {
              vals[attSoFar + k] = instance.value(j);
	    }
	  } else {
	    for (int k = 0; k < att.numValues(); k++) {
	      if (k == (int)instance.value(j)) {
                vals[attSoFar + k] = 1;
	      } else {
                vals[attSoFar + k] = 0;
	      }
	    }
	  }
	  attSoFar += att.numValues();
	}
      }
    }
    Instance inst = null;
    if (instance instanceof SparseInstance) {
      inst = new SparseInstance(instance.weight(), vals);
    } else {
      inst = new DenseInstance(instance.weight(), vals);
    }
    inst.setDataset(getOutputFormat());
    copyValues(inst, false, instance.dataset(), getOutputFormat());
    inst.setDataset(getOutputFormat());
    push(inst);
  }

  /**
   * Convert a single instance over if the class is numeric. The converted 
   * instance is added to the end of the output queue.
   *
   * @param instance the instance to convert
   */
  private void convertInstanceNumeric(Instance instance) {

    double [] vals = new double [outputFormatPeek().numAttributes()];
    int attSoFar = 0;

    for(int j = 0; j < getInputFormat().numAttributes(); j++) {
      Attribute att = getInputFormat().attribute(j);
      if ((!att.isNominal()) || (j == getInputFormat().classIndex())) {
	vals[attSoFar] = instance.value(j);
	attSoFar++;
      } else {
	if (instance.isMissing(j)) {
	  for (int k = 0; k < att.numValues() - 1; k++) {
            vals[attSoFar + k] = instance.value(j);
	  }
	} else {
	  int k = 0;
	  while ((int)instance.value(j) != m_Indices[j][k]) {
            vals[attSoFar + k] = 1;
	    k++;
	  }
	  while (k < att.numValues() - 1) {
            vals[attSoFar + k] = 0;
	    k++;
	  }
	}
	attSoFar += att.numValues() - 1;
      }
    }
    Instance inst = null;
    if (instance instanceof SparseInstance) {
      inst = new SparseInstance(instance.weight(), vals);
    } else {
      inst = new DenseInstance(instance.weight(), vals);
    }
    inst.setDataset(getOutputFormat());
    copyValues(inst, false, instance.dataset(), getOutputFormat());
    inst.setDataset(getOutputFormat());
    push(inst);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: 
   * use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new NominalToBinary(), argv);
  }
}

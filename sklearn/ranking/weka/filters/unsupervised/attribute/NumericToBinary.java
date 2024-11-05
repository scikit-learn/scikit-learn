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
 *    NumericToBinary.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

/** 
 <!-- globalinfo-start -->
 * Converts all numeric attributes into binary attributes (apart from the class attribute, if set): if the value of the numeric attribute is exactly zero, the value of the new attribute will be zero. If the value of the numeric attribute is missing, the value of the new attribute will be missing. Otherwise, the value of the new attribute will be one. The new attributes will be nominal.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -unset-class-temporarily
 *  Unsets the class index temporarily before the filter is
 *  applied to the data.
 *  (default: no)</pre>
 * 
 <!-- options-end -->
 * 
 * @author Eibe Frank (eibe@cs.waikato.ac.nz) 
 * @version $Revision: 5987 $ 
 */
public class NumericToBinary 
  extends PotentialClassIgnorer
  implements UnsupervisedFilter, StreamableFilter {

  /** for serialization */
  static final long serialVersionUID = 2616879323359470802L;
  
  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Converts all numeric attributes into binary attributes (apart from "
      + "the class attribute, if set): if the value of the numeric attribute is "
      + "exactly zero, the value of the new attribute will be zero. If the "
      + "value of the numeric attribute is missing, the value of the new "
      + "attribute will be missing. Otherwise, the value of the new "
      + "attribute will be one. The new attributes will be nominal.";
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
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
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
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    super.setInputFormat(instanceInfo);
    setOutputFormat();
    return true;
  }

  /**
   * Input an instance for filtering.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been defined.
   */
  public boolean input(Instance instance) {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    convertInstance(instance);
    return true;
  }

  /** 
   * Set the output format. 
   */
  private void setOutputFormat() {

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
      if ((j == newClassIndex) || (!att.isNumeric())) {
	newAtts.addElement(att.copy());
      } else {
	attributeName = new StringBuffer(att.name() + "_binarized");
	vals = new FastVector(2);
	vals.addElement("0"); vals.addElement("1");
	newAtts.addElement(new Attribute(attributeName.toString(), vals));
      }
    }
    outputFormat = new Instances(getInputFormat().relationName(), newAtts, 0);
    outputFormat.setClassIndex(newClassIndex);
    setOutputFormat(outputFormat);
  }

  /**
   * Convert a single instance over. The converted instance is 
   * added to the end of the output queue.
   *
   * @param instance the instance to convert
   */
  private void convertInstance(Instance instance) {
  
    Instance inst = null;
    if (instance instanceof SparseInstance) {
      double[] vals = new double[instance.numValues()];
      int[] newIndices = new int[instance.numValues()];
      for (int j = 0; j < instance.numValues(); j++) {
	Attribute att = getInputFormat().attribute(instance.index(j));
	if ((!att.isNumeric()) || (instance.index(j) == getInputFormat().classIndex())) {
	  vals[j] = instance.valueSparse(j);
	} else {
	  if (instance.isMissingSparse(j)) {
	    vals[j] = instance.valueSparse(j);
	  } else {
	    vals[j] = 1;
	  }
	} 
	newIndices[j] = instance.index(j);
      }
      inst = new SparseInstance(instance.weight(), vals, newIndices, 
                                outputFormatPeek().numAttributes());
    } else {
      double[] vals = new double[outputFormatPeek().numAttributes()];
      for (int j = 0; j < getInputFormat().numAttributes(); j++) {
	Attribute att = getInputFormat().attribute(j);
	if ((!att.isNumeric()) || (j == getInputFormat().classIndex())) {
	  vals[j] = instance.value(j);
	} else {
	  if (instance.isMissing(j) || (instance.value(j) == 0)) {
	    vals[j] = instance.value(j);
	  } else {
	    vals[j] = 1;
	  }
	} 
      }
      if(instance instanceof PreferenceDenseInstance){
    	  PreferenceDenseInstance pdi = (PreferenceDenseInstance) instance;
    	  inst = new PreferenceDenseInstance(instance.weight(), vals, pdi.getHashMap());
      }
      else
    	  inst = new DenseInstance(instance.weight(), vals);
    }
    inst.setDataset(instance.dataset());
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
    runFilter(new NumericToBinary(), argv);
  }
}

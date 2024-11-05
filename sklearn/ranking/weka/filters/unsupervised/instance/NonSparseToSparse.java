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
 *    NonSparseToSparse.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.unsupervised.instance;

import java.util.Enumeration;
import java.util.Vector;

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
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

/** 
 <!-- globalinfo-start -->
 * An instance filter that converts all incoming instances into sparse format.
 * <p/>
 <!-- globalinfo-end -->
 * 
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5987 $
 */
public class NonSparseToSparse 
  extends Filter
  implements UnsupervisedFilter, StreamableFilter, OptionHandler {

  /** for serialization */
  static final long serialVersionUID = 4694489111366063852L;
  
  protected boolean m_encodeMissingAsZero = false;
  
  protected boolean m_insertDummyNominalFirstValue = false;
  
  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "An instance filter that converts all incoming instances"
      + " into sparse format.";
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
    
    //RANKING BEGIN
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    //RANKING END
    
    return result;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result;
    
    result = new Vector();
    result.add(new Option("\tTreat missing values as zero.",
        "M", 0, "-M"));
    result.add(new Option("\tAdd a dummy first value for nominal attributes.",
        "F", 0, "-F"));
    
    return result.elements();
  }
  
  public void setOptions(String[] options) throws Exception {
    m_encodeMissingAsZero = Utils.getFlag('M', options);
    m_insertDummyNominalFirstValue = Utils.getFlag('F', options);
  }
  
  public String[] getOptions() {
    Vector result = new Vector();
    
    if (m_encodeMissingAsZero) {
      result.add("-M");      
    }
    
    if (m_insertDummyNominalFirstValue) {
      result.add("-F");
    }
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * Set whether missing values should be treated in the same
   * way as zeros
   * 
   * @param m true if missing values are to be treated the same
   * as zeros
   */
  public void setTreatMissingValuesAsZero(boolean m) {
    m_encodeMissingAsZero = m;
  }
  
  /**
   * Get whether missing values are to be treated in the same
   * way as zeros
   * 
   * @return true if missing values are to be treated in the
   * same way as zeros
   */
  public boolean getTreatMissingValuesAsZero() {
    return m_encodeMissingAsZero;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return            tip text for this property suitable for
   *                    displaying in the explorer/experimenter gui
   */
  public String treatMissingValuesAsZeroTipText() {
    return "Treat missing values in the same way as zeros.";
  }
  
  /**
   * Set whether to insert a dummy first value in the definition
   * for each nominal attribute or not.
   * 
   * @param d true if a dummy value is to be inserted for
   * each nominal attribute.
   */
  public void setInsertDummyNominalFirstValue(boolean d) {
    m_insertDummyNominalFirstValue = d;
  }
  
  /**
   * Get whether a dummy first value will be inserted in the definition
   * of each nominal attribute.
   * 
   * @return true if a dummy first value will be inserted for each nominal
   * attribute.
   */
  public boolean getInsertDummyNominalFirstValue() {
    return m_insertDummyNominalFirstValue;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return            tip text for this property suitable for
   *                    displaying in the explorer/experimenter gui
   */
  public String insertDummyNominalFirstValueTipText() {
    return "Insert a dummy value before the first declared value "
    + "for all nominal attributes. Useful when converting market "
    + "basket data that has been encoded for Apriori to sparse format. "
    + "Typically used in conjuction with treat missing values as zero.";
    
    		
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if format cannot be processed
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    super.setInputFormat(instanceInfo);
    Instances instNew = instanceInfo;
    
    if (m_insertDummyNominalFirstValue) {
      FastVector atts = new FastVector();
      for (int i = 0; i < instanceInfo.numAttributes(); i++) {
        if (instanceInfo.attribute(i).isNominal() || instanceInfo.attribute(i).isRanking()) {
          FastVector labels = new FastVector();
          labels.addElement("_d");
          for (int j = 0; j < instanceInfo.attribute(j).numValues(); j++) {
            labels.addElement(instanceInfo.attribute(i).value(j));
          }
          Attribute newAtt = new Attribute(instanceInfo.attribute(i).name(), 
              labels);
          atts.addElement(newAtt);
        } else {
          atts.addElement(instanceInfo.attribute(i));
        }
      }
      instNew = new Instances(instanceInfo.relationName(), atts, 0);
    }
    
    setOutputFormat(instNew);
    return true;
  }


  /**
   * Input an instance for filtering. Ordinarily the instance is processed
   * and made available for output immediately. Some filters require all
   * instances be read before producing output.
   *
   * @param instance the input instance.
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been set.
   */
  public boolean input(Instance instance) {

    Instance newInstance = null;
    
    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    
    if (m_encodeMissingAsZero && !m_insertDummyNominalFirstValue) {
      Instance tempInst = (Instance)instance.copy();
      tempInst.setDataset(getInputFormat());
      
      for (int i = 0; i < tempInst.numAttributes(); i++) {
        if (tempInst.isMissing(i)) {
          tempInst.setValue(i, 0);
        }
      }
      instance = tempInst;
    }
    
    if (m_insertDummyNominalFirstValue) {
      double[] values = instance.toDoubleArray();      
      for (int i = 0; i < instance.numAttributes(); i++) {
        if (instance.attribute(i).isNominal() || instance.attribute(i).isRanking()) {
          if (!Utils.isMissingValue(values[i])) {
            values[i]++;
          }
        }
        if (m_encodeMissingAsZero && Utils.isMissingValue(values[i])) {
          values[i] = 0;
        }
      }
      newInstance = new SparseInstance(instance.weight(), values);
      newInstance.setDataset(getOutputFormat());
      push(newInstance);
    } else {
      newInstance = new SparseInstance(instance);
      newInstance.setDataset(instance.dataset());
      push(newInstance);
    }
    
    /*Instance inst = new SparseInstance(instance);
    inst.setDataset(instance.dataset());
    push(inst); */
    return true;
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
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new NonSparseToSparse(), argv);
  }
}

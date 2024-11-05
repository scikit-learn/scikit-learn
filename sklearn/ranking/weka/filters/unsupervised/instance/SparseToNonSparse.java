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
 *    SparseToNonSparse.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.unsupervised.instance;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.StreamableFilter;
import weka.filters.UnsupervisedFilter;

/** 
 <!-- globalinfo-start -->
 * An instance filter that converts all incoming sparse instances into non-sparse format.
 * <p/>
 <!-- globalinfo-end -->
 * 
 * @author Len Trigg (len@reeltwo.com)
 * @version $Revision: 5987 $ 
 */
public class SparseToNonSparse 
  extends Filter 
  implements UnsupervisedFilter, StreamableFilter {

  /** for serialization */
  static final long serialVersionUID = 2481634184210236074L;
  
  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "An instance filter that converts all incoming sparse instances"
      + " into non-sparse format.";
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
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if format cannot be processed
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {

    super.setInputFormat(instanceInfo);
    setOutputFormat(instanceInfo);
    return true;
  }


  /**
   * Input an instance for filtering. Ordinarily the instance is processed
   * and made available for output immediately. Some filters require all
   * instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input structure has been defined
   */
  public boolean input(Instance instance) {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    Instance inst = null;
    if (instance instanceof SparseInstance) {
      inst = new DenseInstance(instance.weight(), instance.toDoubleArray());
      inst.setDataset(instance.dataset());
    } else {
      inst = instance;
    }
    push(inst);
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
    runFilter(new SparseToNonSparse(), argv);
  }
}

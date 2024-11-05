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
 *    ReplaceMissingValues.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.unsupervised.attribute;

import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.filters.Sourcable;
import weka.filters.UnsupervisedFilter;

/** 
 <!-- globalinfo-start -->
 * Replaces all missing values for nominal and numeric attributes in a dataset with the modes and means from the training data.
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
public class ReplaceMissingValues 
  extends PotentialClassIgnorer
  implements UnsupervisedFilter, Sourcable {

  /** for serialization */
  static final long serialVersionUID = 8349568310991609867L;
  
  /** The modes and means */
  private double[] m_ModesAndMeans = null;

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Replaces all missing values for nominal and numeric attributes in a "
      + "dataset with the modes and means from the training data.";
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
  public boolean setInputFormat(Instances instanceInfo) 
       throws Exception {

    super.setInputFormat(instanceInfo);
    setOutputFormat(instanceInfo);
    m_ModesAndMeans = null;
    return true;
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been set.
   */
  public boolean input(Instance instance) {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    if (m_ModesAndMeans == null) {
      bufferInput(instance);
      return false;
    } else {
      convertInstance(instance);
      return true;
    }
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

    if (m_ModesAndMeans == null) {
      // Compute modes and means
      double sumOfWeights =  getInputFormat().sumOfWeights();
      double[][] counts = new double[getInputFormat().numAttributes()][];
      for (int i = 0; i < getInputFormat().numAttributes(); i++) {
	if (getInputFormat().attribute(i).isNominal() || getInputFormat().attribute(i).isRanking()) {
	  counts[i] = new double[getInputFormat().attribute(i).numValues()];
          if (counts[i].length > 0)
            counts[i][0] = sumOfWeights;
	}
      }
      double[] sums = new double[getInputFormat().numAttributes()];
      for (int i = 0; i < sums.length; i++) {
	sums[i] = sumOfWeights;
      }
      double[] results = new double[getInputFormat().numAttributes()];
      for (int j = 0; j < getInputFormat().numInstances(); j++) {
	Instance inst = getInputFormat().instance(j);
	for (int i = 0; i < inst.numValues(); i++) {
	  if (!inst.isMissingSparse(i)) {
	    double value = inst.valueSparse(i);
	    if (inst.attributeSparse(i).isNominal() || inst.attributeSparse(i).isRanking()) {
              if (counts[inst.index(i)].length > 0) {
                counts[inst.index(i)][(int)value] += inst.weight();
                counts[inst.index(i)][0] -= inst.weight();
              }
	    } else if (inst.attributeSparse(i).isNumeric()) {
	      results[inst.index(i)] += inst.weight() * inst.valueSparse(i);
	    }
	  } else {
	    if (inst.attributeSparse(i).isNominal() || inst.attributeSparse(i).isRanking()) {
              if (counts[inst.index(i)].length > 0) {
	        counts[inst.index(i)][0] -= inst.weight();
              }
	    } else if (inst.attributeSparse(i).isNumeric()) {
	      sums[inst.index(i)] -= inst.weight();
	    }
	  }
	}
      }
      m_ModesAndMeans = new double[getInputFormat().numAttributes()];
      for (int i = 0; i < getInputFormat().numAttributes(); i++) {
	if (getInputFormat().attribute(i).isNominal() || getInputFormat().attribute(i).isRanking()) {
          if (counts[i].length == 0)
            m_ModesAndMeans[i] = Utils.missingValue();
          else
	    m_ModesAndMeans[i] = (double)Utils.maxIndex(counts[i]);
	} else if (getInputFormat().attribute(i).isNumeric()) {
	  if (Utils.gr(sums[i], 0)) {
	    m_ModesAndMeans[i] = results[i] / sums[i];
	  }
	}
      }

      // Convert pending input instances
      for(int i = 0; i < getInputFormat().numInstances(); i++) {
	convertInstance(getInputFormat().instance(i));
      }
    } 
    // Free memory
    flushInput();

    m_NewBatch = true;
    return (numPendingOutput() != 0);
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
      double []vals = new double[instance.numValues()];
      int []indices = new int[instance.numValues()];
      int num = 0;
      for (int j = 0; j < instance.numValues(); j++) {
	if (instance.isMissingSparse(j) &&
	    (getInputFormat().classIndex() != instance.index(j)) &&
	    (instance.attributeSparse(j).isNominal() ||
	     instance.attributeSparse(j).isRanking() ||
	     instance.attributeSparse(j).isNumeric())) {
	  if (m_ModesAndMeans[instance.index(j)] != 0.0) {
	    vals[num] = m_ModesAndMeans[instance.index(j)];
	    indices[num] = instance.index(j);
	    num++;
	  } 
	} else {
	  vals[num] = instance.valueSparse(j);
	  indices[num] = instance.index(j);
	  num++;
	}
      } 
      if (num == instance.numValues()) {
	inst = new SparseInstance(instance.weight(), vals, indices,
                                  instance.numAttributes());
      } else {
	double []tempVals = new double[num];
	int []tempInd = new int[num];
	System.arraycopy(vals, 0, tempVals, 0, num);
	System.arraycopy(indices, 0, tempInd, 0, num);
	inst = new SparseInstance(instance.weight(), tempVals, tempInd,
                                  instance.numAttributes());
      }
    } else {
      double []vals = new double[getInputFormat().numAttributes()];
      for (int j = 0; j < instance.numAttributes(); j++) {
	if (instance.isMissing(j) &&
	    (getInputFormat().classIndex() != j) &&
	    (getInputFormat().attribute(j).isNominal() ||
	     getInputFormat().attribute(j).isRanking() ||
	     getInputFormat().attribute(j).isNumeric())) {
	  vals[j] = m_ModesAndMeans[j]; 
	} else {
	  vals[j] = instance.value(j);
	}
      } 
      
      //RANKING BEGIN
      PreferenceDenseInstance prefCopy = null;
      if(instance instanceof PreferenceDenseInstance){
    	  prefCopy = (PreferenceDenseInstance)instance;
    	  inst = new PreferenceDenseInstance(prefCopy.weight(),vals,prefCopy.getHashMap());
      }
      else
        inst = new DenseInstance(instance.weight(), vals);     
      //RANKING END
    } 
    inst.setDataset(instance.dataset());
    push(inst);
  }
  
  /**
   * Returns a string that describes the filter as source. The
   * filter will be contained in a class with the given name (there may
   * be auxiliary classes),
   * and will contain two methods with these signatures:
   * <pre><code>
   * // converts one row
   * public static Object[] filter(Object[] i);
   * // converts a full dataset (first dimension is row index)
   * public static Object[][] filter(Object[][] i);
   * </code></pre>
   * where the array <code>i</code> contains elements that are either
   * Double, String, with missing values represented as null. The generated
   * code is public domain and comes with no warranty.
   *
   * @param className   the name that should be given to the source class.
   * @param data	the dataset used for initializing the filter
   * @return            the object source described by a string
   * @throws Exception  if the source can't be computed
   */
  public String toSource(String className, Instances data) throws Exception {
    StringBuffer        result;
    boolean[]		numeric;
    boolean[]		nominal;
    String[]		modes;
    double[]		means;
    int			i;
    
    result = new StringBuffer();
    
    // determine what attributes were processed
    numeric = new boolean[data.numAttributes()];
    nominal = new boolean[data.numAttributes()];
    modes   = new String[data.numAttributes()];
    means   = new double[data.numAttributes()];
    for (i = 0; i < data.numAttributes(); i++) {
      numeric[i] = (data.attribute(i).isNumeric() && (i != data.classIndex()));
      nominal[i] = ((data.attribute(i).isNominal() || data.attribute(i).isRanking()) && (i != data.classIndex()));
      
      if (numeric[i])
	means[i] = m_ModesAndMeans[i];
      else
	means[i] = Double.NaN;

      if (nominal[i])
	modes[i] = data.attribute(i).value((int) m_ModesAndMeans[i]);
      else
	modes[i] = null;
    }
    
    result.append("class " + className + " {\n");
    result.append("\n");
    result.append("  /** lists which numeric attributes will be processed */\n");
    result.append("  protected final static boolean[] NUMERIC = new boolean[]{" + Utils.arrayToString(numeric) + "};\n");
    result.append("\n");
    result.append("  /** lists which nominal attributes will be processed */\n");
    result.append("  protected final static boolean[] NOMINAL = new boolean[]{" + Utils.arrayToString(nominal) + "};\n");
    result.append("\n");
    result.append("  /** the means */\n");
    result.append("  protected final static double[] MEANS = new double[]{" + Utils.arrayToString(means).replaceAll("NaN", "Double.NaN") + "};\n");
    result.append("\n");
    result.append("  /** the modes */\n");
    result.append("  protected final static String[] MODES = new String[]{");
    for (i = 0; i < modes.length; i++) {
      if (i > 0)
	result.append(",");
      if (nominal[i])
	result.append("\"" + Utils.quote(modes[i]) + "\"");
      else
	result.append(modes[i]);
    }
    result.append("};\n");
    result.append("\n");
    result.append("  /**\n");
    result.append("   * filters a single row\n");
    result.append("   * \n");
    result.append("   * @param i the row to process\n");
    result.append("   * @return the processed row\n");
    result.append("   */\n");
    result.append("  public static Object[] filter(Object[] i) {\n");
    result.append("    Object[] result;\n");
    result.append("\n");
    result.append("    result = new Object[i.length];\n");
    result.append("    for (int n = 0; n < i.length; n++) {\n");
    result.append("      if (i[n] == null) {\n");
    result.append("        if (NUMERIC[n])\n");
    result.append("          result[n] = MEANS[n];\n");
    result.append("        else if (NOMINAL[n])\n");
    result.append("          result[n] = MODES[n];\n");
    result.append("        else\n");
    result.append("          result[n] = i[n];\n");
    result.append("      }\n");
    result.append("      else {\n");
    result.append("        result[n] = i[n];\n");
    result.append("      }\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    return result;\n");
    result.append("  }\n");
    result.append("\n");
    result.append("  /**\n");
    result.append("   * filters multiple rows\n");
    result.append("   * \n");
    result.append("   * @param i the rows to process\n");
    result.append("   * @return the processed rows\n");
    result.append("   */\n");
    result.append("  public static Object[][] filter(Object[][] i) {\n");
    result.append("    Object[][] result;\n");
    result.append("\n");
    result.append("    result = new Object[i.length][];\n");
    result.append("    for (int n = 0; n < i.length; n++) {\n");
    result.append("      result[n] = filter(i[n]);\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    return result;\n");
    result.append("  }\n");
    result.append("}\n");
    
    return result.toString();
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
    runFilter(new ReplaceMissingValues(), argv);
  }
}

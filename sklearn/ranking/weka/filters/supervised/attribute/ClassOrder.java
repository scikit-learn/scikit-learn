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
 *    ClassOrder.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */
package weka.filters.supervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Changes the order of the classes so that the class values are no longer of in the order specified in the header. The values will be in the order specified by the user -- it could be either in ascending/descending order by the class frequency or in random order. Note that this filter currently does not change the header, only the class values of the instances, so there is not much point in using it in conjunction with the FilteredClassifier. The value can also be converted back using 'originalValue(double value)' procedure.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -R &lt;seed&gt;
 *  Specify the seed of randomization
 *  used to randomize the class
 *  order (default: 1)</pre>
 * 
 * <pre> -C &lt;order&gt;
 *  Specify the class order to be
 *  sorted, could be 0: ascending
 *  1: descending and 2: random.(default: 0)</pre>
 * 
 <!-- options-end -->
 *
 * @author Xin Xu (xx5@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5491 $
 */
public class ClassOrder 
  extends Filter 
  implements SupervisedFilter, OptionHandler {
    
  /** for serialization */
  static final long serialVersionUID = -2116226838887628411L;
  
  /** The seed of randomization */
  private long m_Seed = 1;
    
  /** The random object */
  private Random m_Random = null;
    
  /** 
   * The 1-1 converting table from the original class values
   * to the new values
   */
  private int[] m_Converter = null;
    
  /** Class attribute of the data */
  private Attribute m_ClassAttribute = null;

  /** The class order to be sorted */
  private int m_ClassOrder = 0;
    
  /** The class values are sorted in ascending order based on their frequencies */
  public static final int FREQ_ASCEND = 0;

  /** The class values are sorted in descending order based on their frequencies */
  public static final int FREQ_DESCEND = 1;
   
  /** The class values are sorted in random order*/
  public static final int RANDOM =2;

  /** This class can provide the class distribution in the sorted order 
   *  as side effect */
  private double[] m_ClassCounts = null;

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Changes the order of the classes so that the class values are " 
      + "no longer of in the order specified in the header. "
      + "The values will be in the order specified by the user "
      + "-- it could be either in ascending/descending order by the class "
      + "frequency or in random order. Note that this filter currently does not "
      + "change the header, only the class values of the instances, "
      + "so there is not much point in using it in conjunction with the "
      + "FilteredClassifier. The value can also be converted back using "
      + "'originalValue(double value)' procedure.";
  }
    
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
	
    Vector newVector = new Vector(1);
	
    newVector.addElement(new Option("\tSpecify the seed of randomization\n"
				    + "\tused to randomize the class\n"
				    + "\torder (default: 1)",
				    "R", 1, "-R <seed>"));
	
    newVector.addElement(new Option("\tSpecify the class order to be\n"
				    + "\tsorted, could be 0: ascending\n"
				    + "\t1: descending and 2: random.(default: 0)",
				    "C", 1, "-C <order>"));
	
    return newVector.elements();
  }
    
    
  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -R &lt;seed&gt;
   *  Specify the seed of randomization
   *  used to randomize the class
   *  order (default: 1)</pre>
   * 
   * <pre> -C &lt;order&gt;
   *  Specify the class order to be
   *  sorted, could be 0: ascending
   *  1: descending and 2: random.(default: 0)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
	
    String seedString = Utils.getOption('R', options);
    if (seedString.length() != 0)
      m_Seed = Long.parseLong(seedString);
    else 
      m_Seed = 1;  
	
    String orderString = Utils.getOption('C', options);
    if (orderString.length() != 0)
      m_ClassOrder = Integer.parseInt(orderString);
    else 
      m_ClassOrder = FREQ_ASCEND;   	
	
    if (getInputFormat() != null)
      setInputFormat(getInputFormat()); 	
	
    m_Random = null;
  }
        
  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {
	
    String [] options = new String [4];
    int current = 0;
	
    options[current++] = "-R"; 
    options[current++] = "" + m_Seed;
    options[current++] = "-C"; 
    options[current++] = "" + m_ClassOrder;
	
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
  public String seedTipText() {
    return "Specify the seed of randomization of the class order";
  }
    
  /**
   * Get the current randomization seed
   *
   * @return a seed
   */
  public long getSeed() {	
    return m_Seed;
  }

  /**
   * Set randomization seed
   *
   * @param seed the set seed
   */
  public void setSeed(long seed){
    m_Seed = seed;
    m_Random = null;
  }
    
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String classOrderTipText() {
    return "Specify the class order after the filtering";
  }
    
  /**
   * Get the wanted class order
   *
   * @return class order
   */
  public int getClassOrder() {	
    return m_ClassOrder;
  }
    
  /**
   * Set the wanted class order
   *
   * @param order the class order
   */
  public void setClassOrder(int order){
    m_ClassOrder = order;
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
    result.enable(Capability.NOMINAL_CLASS);
    
    return result;
  }
    
  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if no class index set or class not nominal
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {     

    super.setInputFormat(new Instances(instanceInfo, 0));	

    m_ClassAttribute = instanceInfo.classAttribute();	
    m_Random = new Random(m_Seed);
    m_Converter = null;
    
    int numClasses = instanceInfo.numClasses();
    m_ClassCounts = new double[numClasses];	
    return false;
  }    
    
  /**
   * Input an instance for filtering. Ordinarily the instance is processed
   * and made available for output immediately. Some filters require all
   * instances be read before producing output.
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
    
    // In case some one use this routine in testing, 
    // although he/she should not do so
    if(m_Converter != null){
      Instance datum = (Instance)instance.copy();
      if (!datum.isMissing(m_ClassAttribute)){
	datum.setClassValue((double)m_Converter[(int)datum.classValue()]);
      }
      push(datum);
      return true;
    }
    
    if (!instance.isMissing(m_ClassAttribute)) {
      m_ClassCounts[(int)instance.classValue()] += instance.weight();
    }

    bufferInput(instance);
    return false;
  }
  
  /**
   * Signify that this batch of input to the filter is finished. If
   * the filter requires all instances prior to filtering, output()
   * may now be called to retrieve the filtered instances. Any
   * subsequent instances filtered should be filtered based on setting
   * obtained from the first batch (unless the inputFormat has been
   * re-assigned or new options have been set). This implementation 
   * sorts the class values and provide class counts in the output format
   *
   * @return true if there are instances pending output
   * @throws IllegalStateException if no input structure has been defined,
   * @throws Exception if there was a problem finishing the batch.
   */
  public boolean batchFinished() throws Exception {

    Instances data = getInputFormat();
    if (data == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_Converter == null) {

      // Get randomized indices and class counts 
      int[] randomIndices = new int[m_ClassCounts.length];
      for (int i = 0; i < randomIndices.length; i++) {
	randomIndices[i] = i;
      }
      for (int j = randomIndices.length - 1; j > 0; j--) {
	int toSwap = m_Random.nextInt(j + 1);
	int tmpIndex = randomIndices[j];
	randomIndices[j] = randomIndices[toSwap];
	randomIndices[toSwap] = tmpIndex;
      }
      
      double[] randomizedCounts = new double[m_ClassCounts.length];
      for (int i = 0; i < randomizedCounts.length; i++) {
	randomizedCounts[i] = m_ClassCounts[randomIndices[i]];
      } 

      // Create new order. For the moment m_Converter converts new indices
      // into old ones.
      if (m_ClassOrder == RANDOM) {
	m_Converter = randomIndices;
	m_ClassCounts = randomizedCounts;
      } else {
	int[] sorted = Utils.sort(randomizedCounts);
	m_Converter = new int[sorted.length];
	if (m_ClassOrder == FREQ_ASCEND) {
	  for (int i = 0; i < sorted.length; i++) {
	    m_Converter[i] = randomIndices[sorted[i]];
	  }
	} else if (m_ClassOrder == FREQ_DESCEND) {
	  for (int i = 0; i < sorted.length; i++) {
	    m_Converter[i] = randomIndices[sorted[sorted.length - i - 1]];
	  }
	} else {
	  throw new IllegalArgumentException("Class order not defined!");
	}
	
	// Change class counts
	double[] tmp2 = new double[m_ClassCounts.length];
	for (int i = 0; i < m_Converter.length; i++) {
	  tmp2[i] = m_ClassCounts[m_Converter[i]];
	}
	m_ClassCounts = tmp2;
      }
      
      // Change the class values
      FastVector values = new FastVector(data.classAttribute().numValues());
      for (int i = 0; i < data.numClasses(); i++) {
	values.addElement(data.classAttribute().value(m_Converter[i]));
      }
      FastVector newVec = new FastVector(data.numAttributes());
      for (int i = 0; i < data.numAttributes(); i++) {
	if (i == data.classIndex()) {
	  newVec.addElement(new Attribute(data.classAttribute().name(), values, 
					  data.classAttribute().getMetadata()));
	} else {
	  newVec.addElement(data.attribute(i));
	}
      }
      Instances newInsts = new Instances(data.relationName(), newVec, 0);
      newInsts.setClassIndex(data.classIndex());
      setOutputFormat(newInsts);

      // From now on we need m_Converter to convert old indices into new ones
      int[] temp = new int[m_Converter.length];
      for (int i = 0; i < temp.length; i++) {
	temp[m_Converter[i]] = i;
      }
      m_Converter = temp;

      // Process all instances
      for(int xyz=0; xyz<data.numInstances(); xyz++){
	Instance datum = data.instance(xyz);
	if (!datum.isMissing(datum.classIndex())) {
	  datum.setClassValue((double)m_Converter[(int)datum.classValue()]);
	}
	push(datum);
      }
    }
    flushInput();
    m_NewBatch = true;
    return (numPendingOutput() != 0);
  }
    
  /**
   * Get the class distribution of the sorted class values.  If class is numeric
   * it returns null 
   *
   * @return the class counts
   */
  public double[] getClassCounts(){ 

    if(m_ClassAttribute.isNominal())
      return m_ClassCounts; 
    else
      return null;
  }

  /**
   * Convert the given class distribution back to the distributions
   * with the original internal class index
   * 
   * @param before the given class distribution
   * @return the distribution converted back
   */
  public double[] distributionsByOriginalIndex (double[] before){

    double[] after = new double[m_Converter.length];
    for(int i=0; i < m_Converter.length; i++) 
      after[i] = before[m_Converter[i]];
    
    return after;
  }

  /**
   * Return the original internal class value given the randomized 
   * class value, i.e. the string presentations of the two indices
   * are the same.  It's useful when the filter is used within a classifier  
   * so that the filtering procedure should be transparent to the 
   * evaluation
   *
   * @param value the given value
   * @return the original internal value, -1 if not found
   * @throws Exception if the coverter table is not set yet
   */
  public double originalValue(double value)throws Exception{

    if(m_Converter == null)
      throw new IllegalStateException("Coverter table not defined yet!");
	
    for(int i=0; i < m_Converter.length; i++)
      if((int)value == m_Converter[i])
	return (double)i;

    return -1;
  }   
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5491 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new ClassOrder(), argv);
  }
}

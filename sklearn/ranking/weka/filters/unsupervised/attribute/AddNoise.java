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
 *    AddNoise.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 */

package weka.filters.unsupervised.attribute;

import weka.core.Capabilities;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SingleIndex;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * An instance filter that changes a percentage of a given attributes values. The attribute must be nominal. Missing value can be treated as value itself.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;col&gt;
 *  Index of the attribute to be changed 
 *  (default last attribute)</pre>
 * 
 * <pre> -M
 *  Treat missing values as an extra value 
 * </pre>
 * 
 * <pre> -P &lt;num&gt;
 *  Specify the percentage of noise introduced 
 *  to the data (default 10)</pre>
 * 
 * <pre> -S &lt;num&gt;
 *  Specify the random number seed (default 1)</pre>
 * 
 <!-- options-end -->
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @version $Revision: 5987 $ 
 */
public class AddNoise 
  extends Filter 
  implements UnsupervisedFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -8499673222857299082L;

  /** The attribute's index setting. */
  private SingleIndex m_AttIndex = new SingleIndex("last"); 

  /** Flag if missing values are taken as value. */
  private boolean m_UseMissing = false;

  /** The subsample size, percent of original set, default 10% */
  private int m_Percent = 10;
  
  /** The random number generator seed */
  private int m_RandomSeed = 1;

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "An instance filter that changes a percentage of a given"
           + " attributes values. The attribute must be nominal."
           + " Missing value can be treated as value itself.";
  }

  /**
   * Returns an enumeration describing the available options
   *
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(4);

    newVector.addElement(new Option(
              "\tIndex of the attribute to be changed \n"
              +"\t(default last attribute)",
              "C", 1, "-C <col>"));
    newVector.addElement(new Option(
              "\tTreat missing values as an extra value \n",
              "M", 1, "-M"));
    newVector.addElement(new Option(
              "\tSpecify the percentage of noise introduced \n"
              +"\tto the data (default 10)",
              "P", 1, "-P <num>"));
    newVector.addElement(new Option(
              "\tSpecify the random number seed (default 1)",
              "S", 1, "-S <num>"));

    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -C &lt;col&gt;
   *  Index of the attribute to be changed 
   *  (default last attribute)</pre>
   * 
   * <pre> -M
   *  Treat missing values as an extra value 
   * </pre>
   * 
   * <pre> -P &lt;num&gt;
   *  Specify the percentage of noise introduced 
   *  to the data (default 10)</pre>
   * 
   * <pre> -S &lt;num&gt;
   *  Specify the random number seed (default 1)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {

    String indexString = Utils.getOption('C', options);
    if (indexString.length() != 0) {
      setAttributeIndex(indexString);
    } else {
      setAttributeIndex("last");
    }

    if (Utils.getFlag('M', options)) {
      setUseMissing(true);
    }

    String percentString = Utils.getOption('P', options);
    if (percentString.length() != 0) {
      setPercent((int) Double.valueOf(percentString).doubleValue());
    } else {
      setPercent(10);
    }

    String seedString = Utils.getOption('S', options);
    if (seedString.length() != 0) {
      setRandomSeed(Integer.parseInt(seedString));
    } else {
      setRandomSeed(1);
    }

  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [7];
    int current = 0;

    options[current++] = "-C"; options[current++] = "" + getAttributeIndex();

    if (getUseMissing()) {
      options[current++] = "-M";
    }

    options[current++] = "-P"; options[current++] = "" + getPercent();

    options[current++] = "-S"; options[current++] = "" + getRandomSeed();

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
  public String useMissingTipText() {

    return "Flag to set if missing values are used.";
  }

  /**
   * Gets the flag if missing values are treated as extra values.
   *
   * @return the flag missing values.
   */
  public boolean getUseMissing() {

    return m_UseMissing;
  }

  /**
   * Sets the flag if missing values are treated as extra values.
   *
   * @param newUseMissing the new flag value.
   */
  public void setUseMissing(boolean newUseMissing) {

    m_UseMissing = newUseMissing;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String randomSeedTipText() {

    return "Random number seed.";
  }

  /**
   * Gets the random number seed.
   *
   * @return the random number seed.
   */
  public int getRandomSeed() {

    return m_RandomSeed;
  }
  
  /**
   * Sets the random number seed.
   *
   * @param newSeed the new random number seed.
   */
  public void setRandomSeed(int newSeed) {

    m_RandomSeed = newSeed;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String percentTipText() {

    return "Percentage of introduced noise to data.";
  }

  /**
   * Gets the size of noise data as a percentage of the original set.
   *
   * @return the noise data size
   */
  public int getPercent() {

    return m_Percent;
  }
  
  /**
   * Sets the size of noise data, as a percentage of the original set.
   *
   * @param newPercent the subsample set size, between 0 and 100.
   */
  public void setPercent(int newPercent) {

    m_Percent = newPercent;
  }
  
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String attributeIndexTipText() {

    return "Index of the attribute that is to changed.";
  }

  /**
   * Get the index of the attribute used.
   *
   * @return the index of the attribute
   */
  public String getAttributeIndex() {

    return m_AttIndex.getSingleIndex();
  }

  /**
   * Sets index of the attribute used.
   *
   * @param attIndex the index of the attribute
   */
  public void setAttributeIndex(String attIndex) {
    
    m_AttIndex.setSingleIndex(attIndex);
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
    // set input format
    //m_InputFormat = new Instances(instanceInfo, 0);
    m_AttIndex.setUpper(getInputFormat().numAttributes() - 1);
    // set index of attribute to be changed

    // test if nominal 
    if (!getInputFormat().attribute(m_AttIndex.getIndex()).isNominal() && !getInputFormat().attribute(m_AttIndex.getIndex()).isRanking()) {
      throw new Exception("Adding noise is not possible:"
                          + "Chosen attribute is numeric.");
      }

    // test if two values are given
    if ((getInputFormat().attribute(m_AttIndex.getIndex()).numValues() < 2)
        && (!m_UseMissing)) {
      throw new Exception("Adding noise is not possible:"
                          + "Chosen attribute has less than two values.");
    }
 
    setOutputFormat(getInputFormat());
    m_NewBatch = true; 
    return false;
  }

  /**
   * Input an instance for filtering. 
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws Exception if the input format was not set
   */
  public boolean input(Instance instance) throws Exception {

    // check if input format is defined
    if (getInputFormat() == null) {
      throw new Exception("No input instance format defined");
    }
    
    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }

    if (isFirstBatchDone()) {
      push(instance);
      return true;
    } else {
      bufferInput(instance);
      return false;
    }
  }

  /**
   * Signify that this batch of input to the filter is finished. 
   * If the filter requires all instances prior to filtering,
   * output() may now be called to retrieve the filtered instances.
   *
   * @return true if there are instances pending output
   * @throws Exception if no input structure has been defined
   */
  public boolean batchFinished() throws Exception {

    if (getInputFormat() == null) {
      throw new Exception("No input instance format defined");
    }

    // Do the subsample, and clear the input instances.
    addNoise (getInputFormat(), m_RandomSeed, m_Percent, m_AttIndex.getIndex(), 
              m_UseMissing);

    for(int i=0; i<getInputFormat().numInstances(); i++) {
      push ((Instance)getInputFormat().instance(i).copy());
    }

    flushInput();

    m_NewBatch = true;
    m_FirstBatchDone = true;
    return (numPendingOutput() != 0);
  }

  /**
   * add noise to the dataset
   * 
   * a given percentage of the instances are changed in the  way, that 
   * a set of instances are randomly selected using seed. The attribute 
   * given by its index is changed from its current value to one of the
   * other possibly ones, also randomly. This is done with leaving the
   * apportion the same.  
   * if m_UseMissing is true, missing value is  used as a value of its own
   * @param instances is the dataset
   * @param seed used for random function
   * @param percent percentage of instances that are changed
   * @param attIndex index of the attribute changed
   * @param useMissing if true missing values are treated as extra value
   */
  public void addNoise (Instances instances, 
                         int seed, 
                         int percent,
                         int attIndex,
                         boolean useMissing) {
    int indexList [];
    int partition_count [];
    int partition_max [];
    double splitPercent = (double) percent; // percentage used for splits

    // fill array with the indexes
    indexList = new int [instances.numInstances()];
    for (int i=0; i<instances.numInstances(); i++) {
      indexList[i] = i;
      }

    // randomize list of indexes
    Random random = new Random(seed);
    for (int i=instances.numInstances()-1; i>=0; i--) {
      int hValue = indexList[i];
      int hIndex = (int)(random.nextDouble()*(double) i);
      indexList[i] = indexList[hIndex];
      indexList[hIndex] = hValue;
      }
 
    // initialize arrays that are used to count instances
    // of each value and to keep the amount of instances of that value 
    // that has to be changed
    // this is done for the missing values in the two variables
    // missing_count and missing_max
    int numValues = instances.attribute(attIndex).numValues();

    partition_count = new int[numValues];
    partition_max = new int[numValues];
    int missing_count = 0;;
    int missing_max = 0;;

    for (int i = 0; i < numValues; i++) {
      partition_count[i] = 0;
      partition_max[i] = 0;
      }

    // go through the dataset and count all occurrences of values 
    // and all missing values using temporarily .._max arrays and
    // variable missing_max
    for (Enumeration e = instances.enumerateInstances();
         e.hasMoreElements();) {
      Instance instance = (Instance) e.nextElement(); 
      if (instance.isMissing(attIndex)) {
        missing_max++;
      }
      else {
        int j = (int) instance.value(attIndex);
        partition_max[(int) instance.value(attIndex)]++; 
      }
    }
      
    // use given percentage to calculate 
    // how many have to be changed per split and
    // how many of the missing values
    if (!useMissing) {
      missing_max = missing_count;
    } else {
      missing_max = (int) (((double)missing_max/100) * splitPercent + 0.5);
    }
    int sum_max = missing_max;
    for (int i=0; i<numValues; i++) {
      partition_max[i]=(int) (((double)partition_max[i]/100) * splitPercent 
                              + 0.5);
      sum_max = sum_max + partition_max[i];
      }

    // initialize sum_count to zero, use this variable to see if 
    // everything is done already
    int sum_count = 0;
  
    // add noise
    // using the randomized index-array
    // 
    Random randomValue = new Random (seed);
    int numOfValues = instances.attribute(attIndex).numValues();
    for(int i=0; i<instances.numInstances(); i++) {
       if (sum_count >= sum_max) { break; } // finished
       Instance currInstance = instances.instance(indexList[i]);
       // if value is missing then...
       if (currInstance.isMissing(attIndex)) {
         if (missing_count < missing_max) {
           changeValueRandomly (randomValue, 
                                numOfValues,
                                attIndex, 
                                currInstance,
                                useMissing); 
           missing_count++;
           sum_count++;
         }
         
       } else {
         int vIndex = (int) currInstance.value(attIndex);
         if (partition_count[vIndex] < partition_max[vIndex]) {
           changeValueRandomly (randomValue,
                                numOfValues,
                                attIndex,     
                                currInstance, 
                                useMissing);           
           partition_count[vIndex]++;
           sum_count++;
         }
       }
    }

  }

  /**
   * method to set a new value
   *
   * @param r random function
   * @param numOfValues 
   * @param instance
   * @param useMissing
   */
  private void changeValueRandomly(Random r, int numOfValues,
                                   int indexOfAtt, 
                                   Instance instance, 
                                   boolean useMissing) {
    int currValue;

    // get current value 
    // if value is missing set current value to number of values
    // whiche is the highest possible value plus one 
    if (instance.isMissing(indexOfAtt)) {
      currValue = numOfValues;
    } else {
      currValue = (int) instance.value(indexOfAtt);
    }

    // with only two possible values it is easier
    if ((numOfValues == 2) && (!instance.isMissing(indexOfAtt))) {
	instance.setValue(indexOfAtt, (double) ((currValue+1)% 2));
    } else {
      // get randomly a new value not equal to the current value
      // if missing values are used as values they must be treated
      // in a special way
      while (true) {
	  int newValue;
        if (useMissing) {
          newValue = (int) (r.nextDouble() * (double) (numOfValues + 1));
        } else {
          newValue = (int) (r.nextDouble() * (double) numOfValues);
        }
        // have we found a new value?
        if (newValue != currValue) { 
          // the value 1 above the highest possible value (=numOfValues)
          // is used as missing value
          if (newValue == numOfValues) { instance.setMissing(indexOfAtt); }
          else { instance.setValue(indexOfAtt, (double) newValue); }
          break;
        }
      }
    }
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
    runFilter(new AddNoise(), argv);
  }
}

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
 *    SpreadSubsample.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.filters.supervised.instance;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.UnassignedClassException;
import weka.core.UnsupportedClassTypeException;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Produces a random subsample of a dataset. The original dataset must fit entirely in memory. This filter allows you to specify the maximum "spread" between the rarest and most common class. For example, you may specify that there be at most a 2:1 difference in class frequencies. When used in batch mode, subsequent batches are NOT resampled.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;num&gt;
 *  Specify the random number seed (default 1)</pre>
 * 
 * <pre> -M &lt;num&gt;
 *  The maximum class distribution spread.
 *  0 = no maximum spread, 1 = uniform distribution, 10 = allow at most
 *  a 10:1 ratio between the classes (default 0)</pre>
 * 
 * <pre> -W
 *  Adjust weights so that total weight per class is maintained.
 *  Individual instance weighting is not preserved. (default no
 *  weights adjustment</pre>
 * 
 * <pre> -X &lt;num&gt;
 *  The maximum count for any class value (default 0 = unlimited).
 * </pre>
 * 
 <!-- options-end -->
 *
 * @author Stuart Inglis (stuart@reeltwo.com)
 * @version $Revision: 5492 $ 
 **/
public class SpreadSubsample 
  extends Filter 
  implements SupervisedFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = -3947033795243930016L;

  /** The random number generator seed */
  private int m_RandomSeed = 1;

  /** The maximum count of any class */
  private int m_MaxCount;

  /** True if the first batch has been done */
  private double m_DistributionSpread = 0;

  /**
   * True if instance weights will be adjusted to maintain
   * total weight per class.
   */
  private boolean m_AdjustWeights = false;

  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Produces a random subsample of a dataset. The original dataset must "
      + "fit entirely in memory. This filter allows you to specify the maximum "
      + "\"spread\" between the rarest and most common class. For example, you may "
      + "specify that there be at most a 2:1 difference in class frequencies. "
      + "When used in batch mode, subsequent batches are NOT resampled.";

  }
    
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String adjustWeightsTipText() {
    return "Wether instance weights will be adjusted to maintain total weight per "
      + "class.";
  }
  
  /**
   * Returns true if instance  weights will be adjusted to maintain
   * total weight per class.
   *
   * @return true if instance weights will be adjusted to maintain
   * total weight per class.
   */
  public boolean getAdjustWeights() {

    return m_AdjustWeights;
  }
  
  /**
   * Sets whether the instance weights will be adjusted to maintain
   * total weight per class.
   *
   * @param newAdjustWeights whether to adjust weights
   */
  public void setAdjustWeights(boolean newAdjustWeights) {

    m_AdjustWeights = newAdjustWeights;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(4);

    newVector.addElement(new Option(
              "\tSpecify the random number seed (default 1)",
              "S", 1, "-S <num>"));
    newVector.addElement(new Option(
              "\tThe maximum class distribution spread.\n"
              +"\t0 = no maximum spread, 1 = uniform distribution, 10 = allow at most\n"
	      +"\ta 10:1 ratio between the classes (default 0)",
              "M", 1, "-M <num>"));
    newVector.addElement(new Option(
              "\tAdjust weights so that total weight per class is maintained.\n"
              +"\tIndividual instance weighting is not preserved. (default no\n"
              +"\tweights adjustment",
              "W", 0, "-W"));
    newVector.addElement(new Option(
	      "\tThe maximum count for any class value (default 0 = unlimited).\n",
              "X", 0, "-X <num>"));

    return newVector.elements();
  }


  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -S &lt;num&gt;
   *  Specify the random number seed (default 1)</pre>
   * 
   * <pre> -M &lt;num&gt;
   *  The maximum class distribution spread.
   *  0 = no maximum spread, 1 = uniform distribution, 10 = allow at most
   *  a 10:1 ratio between the classes (default 0)</pre>
   * 
   * <pre> -W
   *  Adjust weights so that total weight per class is maintained.
   *  Individual instance weighting is not preserved. (default no
   *  weights adjustment</pre>
   * 
   * <pre> -X &lt;num&gt;
   *  The maximum count for any class value (default 0 = unlimited).
   * </pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    
    String seedString = Utils.getOption('S', options);
    if (seedString.length() != 0) {
      setRandomSeed(Integer.parseInt(seedString));
    } else {
      setRandomSeed(1);
    }

    String maxString = Utils.getOption('M', options);
    if (maxString.length() != 0) {
      setDistributionSpread(Double.valueOf(maxString).doubleValue());
    } else {
      setDistributionSpread(0);
    }

    String maxCount = Utils.getOption('X', options);
    if (maxCount.length() != 0) {
      setMaxCount(Double.valueOf(maxCount).doubleValue());
    } else {
      setMaxCount(0);
    }

    setAdjustWeights(Utils.getFlag('W', options));

    if (getInputFormat() != null) {
      setInputFormat(getInputFormat());
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

    options[current++] = "-M"; 
    options[current++] = "" + getDistributionSpread();

    options[current++] = "-X"; 
    options[current++] = "" + getMaxCount();

    options[current++] = "-S"; 
    options[current++] = "" + getRandomSeed();

    if (getAdjustWeights()) {
      options[current++] = "-W";
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
  public String distributionSpreadTipText() {
    return "The maximum class distribution spread. "
      + "(0 = no maximum spread, 1 = uniform distribution, 10 = allow at most a "
      + "10:1 ratio between the classes).";
  }
  
  /**
   * Sets the value for the distribution spread
   *
   * @param spread the new distribution spread
   */
  public void setDistributionSpread(double spread) {

    m_DistributionSpread = spread;
  }

  /**
   * Gets the value for the distribution spread
   *
   * @return the distribution spread
   */    
  public double getDistributionSpread() {

    return m_DistributionSpread;
  }
    
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String maxCountTipText() {
    return "The maximum count for any class value (0 = unlimited).";
  }
  
  /**
   * Sets the value for the max count
   *
   * @param maxcount the new max count
   */
  public void setMaxCount(double maxcount) {

    m_MaxCount = (int)maxcount;
  }

  /**
   * Gets the value for the max count
   *
   * @return the max count
   */    
  public double getMaxCount() {

    return m_MaxCount;
  }
    
  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String randomSeedTipText() {
    return "Sets the random number seed for subsampling.";
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
   * @param instanceInfo an Instances object containing the input 
   * instance structure (any instances contained in the object are 
   * ignored - only the structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws UnassignedClassException if no class attribute has been set.
   * @throws UnsupportedClassTypeException if the class attribute
   * is not nominal. 
   */
  public boolean setInputFormat(Instances instanceInfo) 
       throws Exception {

    super.setInputFormat(instanceInfo);
    setOutputFormat(instanceInfo);
    return true;
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
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
   * @throws IllegalStateException if no input structure has been defined
   */
  public boolean batchFinished() {

    if (getInputFormat() == null) {
      throw new IllegalStateException("No input instance format defined");
    }

    if (!isFirstBatchDone()) {
      // Do the subsample, and clear the input instances.
      createSubsample();
    }

    flushInput();
    m_NewBatch = true;
    m_FirstBatchDone = true;
    return (numPendingOutput() != 0);
  }


  /**
   * Creates a subsample of the current set of input instances. The output
   * instances are pushed onto the output queue for collection.
   */
  private void createSubsample() {

    int classI = getInputFormat().classIndex();
    // Sort according to class attribute.
    getInputFormat().sort(classI);
    // Determine where each class starts in the sorted dataset
    int [] classIndices = getClassIndices();

    // Get the existing class distribution
    int [] counts = new int [getInputFormat().numClasses()];
    double [] weights = new double [getInputFormat().numClasses()];
    int min = -1;
    for (int i = 0; i < getInputFormat().numInstances(); i++) {
      Instance current = getInputFormat().instance(i);
      if (current.classIsMissing() == false) {
        counts[(int)current.classValue()]++;
        weights[(int)current.classValue()]+= current.weight();
      }
    }

    // Convert from total weight to average weight
    for (int i = 0; i < counts.length; i++) {
      if (counts[i] > 0) {
        weights[i] = weights[i] / counts[i];
      }
      /*
      System.err.println("Class:" + i + " " + getInputFormat().classAttribute().value(i)
                         + " Count:" + counts[i]
                         + " Total:" + weights[i] * counts[i]
                         + " Avg:" + weights[i]);
      */
    }
    
    // find the class with the minimum number of instances
    int minIndex = -1;
    for (int i = 0; i < counts.length; i++) {
      if ( (min < 0) && (counts[i] > 0) ) {
        min = counts[i];
        minIndex = i;
      } else if ((counts[i] < min) && (counts[i] > 0)) {
        min = counts[i];
        minIndex = i;
      }
    }

    if (min < 0) { 
	System.err.println("SpreadSubsample: *warning* none of the classes have any values in them.");
	return;
    }

    // determine the new distribution 
    int [] new_counts = new int [getInputFormat().numClasses()];
    for (int i = 0; i < counts.length; i++) {
      new_counts[i] = (int)Math.abs(Math.min(counts[i],
                                             min * m_DistributionSpread));
      if (i == minIndex) {
        if (m_DistributionSpread > 0 && m_DistributionSpread < 1.0) {
          // don't undersample the minority class!
          new_counts[i] = counts[i];
        }
      }
      if (m_DistributionSpread == 0) {
        new_counts[i] = counts[i];
      }

      if (m_MaxCount > 0) {
        new_counts[i] = Math.min(new_counts[i], m_MaxCount);
      }
    }

    // Sample without replacement
    Random random = new Random(m_RandomSeed);
    Hashtable t = new Hashtable();
    for (int j = 0; j < new_counts.length; j++) {
      double newWeight = 1.0;
      if (m_AdjustWeights && (new_counts[j] > 0)) {
        newWeight = weights[j] * counts[j] / new_counts[j];
        /*
        System.err.println("Class:" + j + " " + getInputFormat().classAttribute().value(j) 
                           + " Count:" + counts[j]
                           + " Total:" + weights[j] * counts[j]
                           + " Avg:" + weights[j]
                           + " NewCount:" + new_counts[j]
                           + " NewAvg:" + newWeight);
        */
      }
      for (int k = 0; k < new_counts[j]; k++) {
        boolean ok = false;
        do {
	  int index = classIndices[j] + (Math.abs(random.nextInt()) 
                                         % (classIndices[j + 1] - classIndices[j])) ;
	  // Have we used this instance before?
          if (t.get("" + index) == null) {
            // if not, add it to the hashtable and use it
            t.put("" + index, "");
            ok = true;
	    if(index >= 0) {
              Instance newInst = (Instance)getInputFormat().instance(index).copy();
              if (m_AdjustWeights) {
                newInst.setWeight(newWeight);
              }
              push(newInst);
            }
          }
        } while (!ok);
      }
    }
  }

  /**
   * Creates an index containing the position where each class starts in 
   * the getInputFormat(). m_InputFormat must be sorted on the class attribute.
   * 
   * @return the positions
   */
  private int[] getClassIndices() {

    // Create an index of where each class value starts
    int [] classIndices = new int [getInputFormat().numClasses() + 1];
    int currentClass = 0;
    classIndices[currentClass] = 0;
    for (int i = 0; i < getInputFormat().numInstances(); i++) {
      Instance current = getInputFormat().instance(i);
      if (current.classIsMissing()) {
        for (int j = currentClass + 1; j < classIndices.length; j++) {
          classIndices[j] = i;
        }
        break;
      } else if (current.classValue() != currentClass) {
        for (int j = currentClass + 1; j <= current.classValue(); j++) {
          classIndices[j] = i;
        }          
        currentClass = (int) current.classValue();
      }
    }
    if (currentClass <= getInputFormat().numClasses()) {
      for (int j = currentClass + 1; j < classIndices.length; j++) {
        classIndices[j] = getInputFormat().numInstances();
      }
    }
    return classIndices;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5492 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: 
   * use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new SpreadSubsample(), argv);
  }
}

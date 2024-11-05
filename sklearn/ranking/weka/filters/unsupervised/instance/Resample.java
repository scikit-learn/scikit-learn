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
 *    Resample.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.instance;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Produces a random subsample of a dataset using either sampling with replacement or without replacement. The original dataset must fit entirely in memory. The number of instances in the generated dataset may be specified. When used in batch mode, subsequent batches are NOT resampled.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -S &lt;num&gt;
 *  Specify the random number seed (default 1)</pre>
 * 
 * <pre> -Z &lt;num&gt;
 *  The size of the output dataset, as a percentage of
 *  the input dataset (default 100)</pre>
 * 
 * <pre> -no-replacement
 *  Disables replacement of instances
 *  (default: with replacement)</pre>
 * 
 * <pre> -V
 *  Inverts the selection - only available with '-no-replacement'.</pre>
 * 
 <!-- options-end -->
 *
 * @author Len Trigg (len@reeltwo.com)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5499 $ 
 */
public class Resample 
  extends Filter 
  implements UnsupervisedFilter, OptionHandler {
  
  /** for serialization */
  static final long serialVersionUID = 3119607037607101160L;

  /** The subsample size, percent of original set, default 100% */
  protected double m_SampleSizePercent = 100;
  
  /** The random number generator seed */
  protected int m_RandomSeed = 1;

  /** Whether to perform sampling with replacement or without */
  protected boolean m_NoReplacement = false;

  /** Whether to invert the selection (only if instances are drawn WITHOUT 
   * replacement)
   * @see #m_NoReplacement */
  protected boolean m_InvertSelection = false;
  
  /**
   * Returns a string describing this classifier
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Produces a random subsample of a dataset using either sampling with "
      + "replacement or without replacement. The original dataset must fit "
      + "entirely in memory. The number of instances in the generated dataset "
      + "may be specified. When used in batch mode, subsequent batches are "
      + "NOT resampled.";
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector result = new Vector();

    result.addElement(new Option(
	"\tSpecify the random number seed (default 1)",
	"S", 1, "-S <num>"));

    result.addElement(new Option(
	"\tThe size of the output dataset, as a percentage of\n"
	+"\tthe input dataset (default 100)",
	"Z", 1, "-Z <num>"));

    result.addElement(new Option(
	"\tDisables replacement of instances\n"
	+"\t(default: with replacement)",
	"no-replacement", 0, "-no-replacement"));

    result.addElement(new Option(
	"\tInverts the selection - only available with '-no-replacement'.",
	"V", 0, "-V"));

    return result.elements();
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
   * <pre> -Z &lt;num&gt;
   *  The size of the output dataset, as a percentage of
   *  the input dataset (default 100)</pre>
   * 
   * <pre> -no-replacement
   *  Disables replacement of instances
   *  (default: with replacement)</pre>
   * 
   * <pre> -V
   *  Inverts the selection - only available with '-no-replacement'.</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('S', options);
    if (tmpStr.length() != 0)
      setRandomSeed(Integer.parseInt(tmpStr));
    else
      setRandomSeed(1);

    tmpStr = Utils.getOption('Z', options);
    if (tmpStr.length() != 0)
      setSampleSizePercent(Double.parseDouble(tmpStr));
    else
      setSampleSizePercent(100);

    setNoReplacement(Utils.getFlag("no-replacement", options));

    if (getNoReplacement())
      setInvertSelection(Utils.getFlag('V', options));

    if (getInputFormat() != null)
      setInputFormat(getInputFormat());
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {
    Vector<String>	result;

    result = new Vector<String>();

    result.add("-S");
    result.add("" + getRandomSeed());

    result.add("-Z");
    result.add("" + getSampleSizePercent());

    if (getNoReplacement()) {
      result.add("-no-replacement");
      if (getInvertSelection())
	result.add("-V");
    }
    
    return result.toArray(new String[result.size()]);
  }

  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String randomSeedTipText() {
    return "The seed used for random sampling.";
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
  public String sampleSizePercentTipText() {
    return "Size of the subsample as a percentage of the original dataset.";
  }

  /**
   * Gets the subsample size as a percentage of the original set.
   *
   * @return the subsample size
   */
  public double getSampleSizePercent() {
    return m_SampleSizePercent;
  }
  
  /**
   * Sets the size of the subsample, as a percentage of the original set.
   *
   * @param newSampleSizePercent the subsample set size, between 0 and 100.
   */
  public void setSampleSizePercent(double newSampleSizePercent) {
    m_SampleSizePercent = newSampleSizePercent;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String noReplacementTipText() {
    return "Disables the replacement of instances.";
  }

  /**
   * Gets whether instances are drawn with or without replacement.
   * 
   * @return true if the replacement is disabled
   */
  public boolean getNoReplacement() {
    return m_NoReplacement;
  }
  
  /**
   * Sets whether instances are drawn with or with out replacement.
   * 
   * @param value if true then the replacement of instances is disabled
   */
  public void setNoReplacement(boolean value) {
    m_NoReplacement = value;
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String invertSelectionTipText() {
    return "Inverts the selection (only if instances are drawn WITHOUT replacement).";
  }

  /**
   * Gets whether selection is inverted (only if instances are drawn WIHTOUT 
   * replacement).
   * 
   * @return true if the replacement is disabled
   * @see #m_NoReplacement
   */
  public boolean getInvertSelection() {
    return m_InvertSelection;
  }
  
  /**
   * Sets whether the selection is inverted (only if instances are drawn WIHTOUT 
   * replacement).
   * 
   * @param value if true then selection is inverted
   */
  public void setInvertSelection(boolean value) {
    m_InvertSelection = value;
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
   * creates the subsample with replacement
   * 
   * @param random	the random number generator to use
   * @param origSize	the original size of the dataset
   * @param sampleSize	the size to generate
   */
  public void createSubsampleWithReplacement(Random random, int origSize, 
      int sampleSize) {
    
    for (int i = 0; i < sampleSize; i++) {
	int index = random.nextInt(origSize);
	push((Instance) getInputFormat().instance(index).copy());
    }
  }

  /**
   * creates the subsample without replacement
   * 
   * @param random	the random number generator to use
   * @param origSize	the original size of the dataset
   * @param sampleSize	the size to generate
   */
  public void createSubsampleWithoutReplacement(Random random, int origSize, 
      int sampleSize) {
    
    if (sampleSize > origSize) {
      sampleSize = origSize;
      System.err.println(
	  "Resampling with replacement can only use percentage <=100% - "
	  + "Using full dataset!");
    }

    Vector<Integer> indices = new Vector<Integer>(origSize);
    Vector<Integer> indicesNew = new Vector<Integer>(sampleSize);

    // generate list of all indices to draw from
    for (int i = 0; i < origSize; i++)
      indices.add(i);

    // draw X random indices (selected ones get removed before next draw)
    for (int i = 0; i < sampleSize; i++) {
      int index = random.nextInt(indices.size());
      indicesNew.add(indices.get(index));
      indices.remove(index);
    }

    if (getInvertSelection())
      indicesNew = indices;
    else
      Collections.sort(indicesNew);

    for (int i = 0; i < indicesNew.size(); i++)
      push((Instance) getInputFormat().instance(indicesNew.get(i)).copy());

    // clean up
    indices.clear();
    indicesNew.clear();
    indices = null;
    indicesNew = null;
  }
  
  /**
   * Creates a subsample of the current set of input instances. The output
   * instances are pushed onto the output queue for collection.
   */
  protected void createSubsample() {
    int origSize = getInputFormat().numInstances();
    int sampleSize = (int) (origSize * m_SampleSizePercent / 100);
    Random random = new Random(m_RandomSeed);
    
    if (getNoReplacement())
      createSubsampleWithoutReplacement(random, origSize, sampleSize);
    else
      createSubsampleWithReplacement(random, origSize, sampleSize);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5499 $");
  }
  
  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: 
   * use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new Resample(), argv);
  }
}

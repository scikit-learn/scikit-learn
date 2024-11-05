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
 * SimpleStreamFilter.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters;

import weka.core.Instance;
import weka.core.Instances;

/** 
 * This filter is a superclass for simple stream filters. <p/>
 * <p/>
 * 
 * <b>General notes:</b><br/>
 * <ul>
 *   <li>After the first call of batchFinished() the field m_FirstBatchDone is
 *   set to <code>true</code>. </li>
 * </ul>
 * <p/>
 * 
 * <b>Example:</b><br/>
 * The following code snippet uses the filter <code>SomeFilter</code> on a
 * dataset that is loaded from <code>filename</code>.
 * <pre>
 * import weka.core.*;
 * import weka.filters.*;
 * import java.io.*;
 * ...
 * SomeFilter filter = new SomeFilter();
 * // set necessary options for the filter
 * Instances data = new Instances(
 *                    new BufferedReader(
 *                      new FileReader(filename)));
 * Instances filteredData = Filter.useFilter(data, filter);
 * </pre>
 *
 * <b>Implementation:</b><br/>
 * Only the following abstract methods need to be implemented:
 * <ul>
 *   <li>globalInfo()</li>
 *   <li>determineOutputFormat(Instances)</li>
 *   <li>process(Instance)</li>
 * </ul>
 * <br/>
 * And the <b>getCapabilities()</b> method must return what kind of
 * attributes and classes the filter can handle.
 * <p/>
 * 
 * If more options are necessary, then the following methods need to be
 * overriden:
 * <ul>
 *   <li>listOptions()</li>
 *   <li>setOptions(String[])</li>
 *   <li>getOptions()</li>
 * </ul>
 * <p/>
 *
 * To make the filter available from commandline one must add the following
 * main method for correct execution (&lt;Filtername&gt; must be replaced
 * with the actual filter classname):
 * <pre>
 *  public static void main(String[] args) {
 *    runFilter(new &lt;Filtername&gt;(), args);
 *  }
 * </pre>
 * <p/>
 *
 * <b>Example implementation:</b><br/>
 * <pre>
 * import weka.core.*;
 * import weka.core.Capabilities.*;
 * import weka.filters.*;
 *
 * import java.util.Random;
 * 
 * public class SimpleStream
 *   extends SimpleStreamFilter {
 *   
 *   public String globalInfo() {
 *     return "A simple stream filter that adds an attribute 'bla' at the end containing a random number.";
 *   }
 *     
 *   public Capabilities getCapabilities() {
 *     Capabilities result = super.getCapabilities();
 *     result.enableAllAttributes();
 *     result.enableAllClasses();
 *     result.enable(Capability.NO_CLASS);  // filter doesn't need class to be set
 *     return result;
 *   }
 * 
 *   protected Instances determineOutputFormat(Instances inputFormat) {
 *     Instances result = new Instances(inputFormat, 0);
 *     result.insertAttributeAt(new Attribute("bla"), result.numAttributes());
 *     return result;
 *   }
 * 
 *   protected Instance process(Instance inst) {
 *     double[] values = new double[inst.numAttributes() + 1];
 *     for (int n = 0; n &lt; inst.numAttributes(); n++)
 *       values[n] = inst.value(n);
 *     values[values.length - 1] = new Random().nextInt();
 *     Instance result = new Instance(1, values);
 *     return result;
 *   }
 * 
 *   public static void main(String[] args) {
 *     runFilter(new SimpleStream(), args);
 *   }
 * }
 * </pre>
 * <p/>
 *
 * <b>Options:</b><br/>
 * Valid filter-specific options are:<p/>
 *
 * -D <br/>
 * Turns on output of debugging information.<p/>
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5047 $
 * @see     SimpleBatchFilter 
 * @see     #input(Instance)
 * @see     #batchFinished()
 * @see     #m_FirstBatchDone
 */
public abstract class SimpleStreamFilter
  extends SimpleFilter 
  implements StreamableFilter {

  /** for serialization */
  private static final long serialVersionUID = 2754882676192747091L;

  /**
   * Returns true if the output format is immediately available after the
   * input format has been set and not only after all the data has been
   * seen (see batchFinished()). This method should normally return true
   * for a stream filter, since the data will be processed in a batch
   * manner instead (or at least for the second batch of files, see
   * m_FirstBatchDone).
   *
   * @return      true if the output format is immediately available
   * @see         #batchFinished()
   * @see         #setInputFormat(Instances)
   * @see         #m_FirstBatchDone
   */
  protected boolean hasImmediateOutputFormat() {
    return true;
  }
  
  /**
   * Determines the output format based on the input format and returns 
   * this. In case the output format cannot be returned immediately, i.e.,
   * hasImmediateOutputFormat() returns false, then this method will called
   * from batchFinished() after the call of preprocess(Instances), in which,
   * e.g., statistics for the actual processing step can be gathered.
   *
   * @param inputFormat     the input format to base the output format on
   * @return                the output format
   * @throws Exception      in case the determination goes wrong
   * @see                   #hasImmediateOutputFormat()
   * @see                   #batchFinished()
   * @see                   #preprocess(Instances)
   */
  protected abstract Instances determineOutputFormat(Instances inputFormat) throws Exception;

  /**
   * processes the given instance (may change the provided instance) and
   * returns the modified version.
   *
   * @param instance    the instance to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   */
  protected abstract Instance process(Instance instance) throws Exception;

  /**
   * Processes the given data (may change the provided dataset) and returns
   * the modified version. This method is called in batchFinished().
   * This implementation only calls process(Instance) for each instance
   * in the given dataset.
   *
   * @param instances   the data to process
   * @return            the modified data
   * @throws Exception  in case the processing goes wrong
   * @see               #batchFinished()
   * @see               #process(Instance)
   */
  protected Instances process(Instances instances) throws Exception {
    Instances     result;
    int           i;

    result = new Instances(getOutputFormat(), 0);
    
    for (i = 0; i < instances.numInstances(); i++)
      result.add(process(instances.instance(i)));
    
    return result;
  }

  /**
   * In case the output format cannot be returned immediately, this method
   * is called before the actual processing of the instances. Derived classes
   * can implement specific behavior here.
   *
   * @param instances   the instances to work on
   * @see               #hasImmediateOutputFormat()
   * @see               #determineOutputFormat(Instances)
   */
  protected void preprocess(Instances instances) {
  }

  /**
   * Input an instance for filtering. Filter requires all
   * training instances be read before producing output.
   *
   * @param instance    the input instance
   * @return            true if the filtered instance may now be
   *                    collected with output().
   * @throws IllegalStateException 	if no input structure has been defined
   * @throws Exception	if something goes wrong
   */
  public boolean input(Instance instance) throws Exception {
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    if (m_NewBatch) {
      resetQueue();
      m_NewBatch = false;
    }
    
    try {
      if (hasImmediateOutputFormat() || isFirstBatchDone()) {
        push(process((Instance) instance.copy()));
        return true;
      }
      else {
        bufferInput(instance);
        return false;
      }
    }
    catch (Exception e) {
      return false;
    }
  }

  /**
   * Signify that this batch of input to the filter is finished. If
   * the filter requires all instances prior to filtering, output()
   * may now be called to retrieve the filtered instances. Any
   * subsequent instances filtered should be filtered based on setting
   * obtained from the first batch (unless the setInputFormat has been
   * re-assigned or new options have been set).
   *
   * @return 		true if there are instances pending output
   * @throws IllegalStateException 	if no input format has been set. 
   */
  public boolean batchFinished() throws Exception {
    int         i;
    Instances   inst;
    
    if (getInputFormat() == null)
      throw new IllegalStateException("No input instance format defined");

    inst = new Instances(getInputFormat());
    flushInput();

    if (!hasImmediateOutputFormat())
      preprocess(inst);

    // process data
    inst = process(inst);

    // if output format hasn't been set yet, do it now
    if (!hasImmediateOutputFormat() && !isFirstBatchDone())
      setOutputFormat(inst);

    // move data to the output
    for (i = 0; i < inst.numInstances(); i++)
      push(inst.instance(i));
    
    m_NewBatch       = true;
    m_FirstBatchDone = true;
    
    return (numPendingOutput() != 0);
  }
}


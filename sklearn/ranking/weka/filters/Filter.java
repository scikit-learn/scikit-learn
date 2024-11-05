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
 *    Filter.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters;

import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Queue;
import weka.core.RelationalLocator;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.SerializedObject;
import weka.core.StringLocator;
import weka.core.UnsupportedAttributeTypeException;
import weka.core.Utils;
import weka.core.Version;
import weka.core.Capabilities.Capability;
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.labelranking.PreferenceDenseInstance;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Date;
import java.util.Enumeration;
import java.util.Iterator;

/** 
 * An abstract class for instance filters: objects that take instances
 * as input, carry out some transformation on the instance and then
 * output the instance. The method implementations in this class
 * assume that most of the work will be done in the methods overridden
 * by subclasses.<p>
 *
 * A simple example of filter use. This example doesn't remove
 * instances from the output queue until all instances have been
 * input, so has higher memory consumption than an approach that
 * uses output instances as they are made available:<p>
 *
 * <code> <pre>
 *  Filter filter = ..some type of filter..
 *  Instances instances = ..some instances..
 *  for (int i = 0; i < data.numInstances(); i++) {
 *    filter.input(data.instance(i));
 *  }
 *  filter.batchFinished();
 *  Instances newData = filter.outputFormat();
 *  Instance processed;
 *  while ((processed = filter.output()) != null) {
 *    newData.add(processed);
 *  }
 *  ..do something with newData..
 * </pre> </code>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6797 $
 */
public abstract class Filter
  implements Serializable, CapabilitiesHandler, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -8835063755891851218L;

  /** The output format for instances */
  private Instances m_OutputFormat = null;

  /** The output instance queue */
  private Queue m_OutputQueue = null;

  /** Indices of string attributes in the output format */
  protected StringLocator m_OutputStringAtts = null;

  /** Indices of string attributes in the input format */
  protected StringLocator m_InputStringAtts = null;

  /** Indices of relational attributes in the output format */
  protected RelationalLocator m_OutputRelAtts = null;

  /** Indices of relational attributes in the input format */
  protected RelationalLocator m_InputRelAtts = null;

  /** The input format for instances */
  private Instances m_InputFormat = null;

  /** Record whether the filter is at the start of a batch */
  protected boolean m_NewBatch = true;

  /** True if the first batch has been done */
  protected boolean m_FirstBatchDone = false;

  /**
   * Returns true if the a new batch was started, either a new instance of the 
   * filter was created or the batchFinished() method got called.
   * 
   * @return true if a new batch has been initiated
   * @see #m_NewBatch
   * @see #batchFinished()
   */
  public boolean isNewBatch() {
    return m_NewBatch;
  }
  
  /**
   * Returns true if the first batch of instances got processed. Necessary for
   * supervised filters, which "learn" from the first batch and then shouldn't
   * get updated with subsequent calls of batchFinished().
   * 
   * @return true if the first batch has been processed
   * @see #m_FirstBatchDone
   * @see #batchFinished()
   */
  public boolean isFirstBatchDone() {
    return m_FirstBatchDone;
  }
  
  /**
   * Default implementation returns false. Some filters may not
   * necessarily be able to produce an instance for output for
   * every instance input after the first batch has been 
   * completed - such filters should override this method
   * and return true.
   * 
   * @return false by default
   */
  public boolean mayRemoveInstanceAfterFirstBatchDone() {
    return false;
  }

  /** 
   * Returns the Capabilities of this filter. Derived filters have to
   * override this method to enable capabilities.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities 	result;

    result = new Capabilities(this);
    result.enableAll();
    
    result.setMinimumNumberInstances(0);
    
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return            the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6797 $");
  }

  /** 
   * Returns the Capabilities of this filter, customized based on the data.
   * I.e., if removes all class capabilities, in case there's not class
   * attribute present or removes the NO_CLASS capability, in case that
   * there's a class present.
   *
   * @param data	the data to use for customization
   * @return            the capabilities of this object, based on the data
   * @see               #getCapabilities()
   */
  public Capabilities getCapabilities(Instances data) {
    Capabilities 	result;
    Capabilities 	classes;
    Iterator		iter;
    Capability		cap;

    result = getCapabilities();

    // no class? -> remove all class capabilites apart from NO_CLASS
    if (data.classIndex() == -1) {
      classes = result.getClassCapabilities();
      iter    = classes.capabilities();
      while (iter.hasNext()) {
	cap = (Capability) iter.next();
	if (cap != Capability.NO_CLASS) {
	  result.disable(cap);
	  result.disableDependency(cap);
	}
      }
    }
    // class? -> remove NO_CLASS
    else {
      result.disable(Capability.NO_CLASS);
      result.disableDependency(Capability.NO_CLASS);
    }
    
    return result;
  }

  /**
   * Sets the format of output instances. The derived class should use this
   * method once it has determined the outputformat. The 
   * output queue is cleared.
   *
   * @param outputFormat the new output format
   */
  protected void setOutputFormat(Instances outputFormat) {

    if (outputFormat != null) {
      m_OutputFormat = outputFormat.stringFreeStructure();
      initOutputLocators(m_OutputFormat, null);

      // Rename the relation
      String relationName = outputFormat.relationName() 
        + "-" + this.getClass().getName();
      if (this instanceof OptionHandler) {
        String [] options = ((OptionHandler)this).getOptions();
        for (int i = 0; i < options.length; i++) {
          relationName += options[i].trim();
        }
      }
      m_OutputFormat.setRelationName(relationName);
    } else {
      m_OutputFormat = null;
    }
    m_OutputQueue = new Queue();
  }

  /**
   * Gets the currently set inputformat instances. This dataset may contain
   * buffered instances.
   *
   * @return the input Instances.
   */
  protected Instances getInputFormat() {

    return m_InputFormat;
  }

  /**
   * Returns a reference to the current input format without
   * copying it.
   *
   * @return a reference to the current input format
   */
  protected Instances inputFormatPeek() {

    return m_InputFormat;
  }

  /**
   * Returns a reference to the current output format without
   * copying it.
   *
   * @return a reference to the current output format
   */
  protected Instances outputFormatPeek() {

    return m_OutputFormat;
  }

  /**
   * Adds an output instance to the queue. The derived class should use this
   * method for each output instance it makes available. 
   *
   * @param instance the instance to be added to the queue.
   */
  protected void push(Instance instance) {

    if (instance != null) {
      if (instance.dataset() != null)
	copyValues(instance, false);
      instance.setDataset(m_OutputFormat);
      m_OutputQueue.push(instance);
    }
  }

  /**
   * Clears the output queue.
   */
  protected void resetQueue() {

    m_OutputQueue = new Queue();
  }

  /**
   * Adds the supplied input instance to the inputformat dataset for
   * later processing.  Use this method rather than
   * getInputFormat().add(instance). Or else. Note that the provided
   * instance gets copied when buffered. 
   *
   * @param instance the <code>Instance</code> to buffer.  
   */
  protected void bufferInput(Instance instance) {
    if (instance != null) {
      copyValues(instance, true);
      m_InputFormat.add(instance);
    }
  }

  /**
   * Initializes the input attribute locators. If indices is null then all 
   * attributes of the data will be considered, otherwise only the ones
   * that were provided.
   * 
   * @param data		the data to initialize the locators with
   * @param indices		if not null, the indices to which to restrict
   * 				the locating
   */
  protected void initInputLocators(Instances data, int[] indices) {
    if (indices == null) {
      m_InputStringAtts = new StringLocator(data);
      m_InputRelAtts    = new RelationalLocator(data);
    }
    else {
      m_InputStringAtts = new StringLocator(data, indices);
      m_InputRelAtts    = new RelationalLocator(data, indices);
    }
  }

  /**
   * Initializes the output attribute locators. If indices is null then all 
   * attributes of the data will be considered, otherwise only the ones
   * that were provided.
   * 
   * @param data		the data to initialize the locators with
   * @param indices		if not null, the indices to which to restrict
   * 				the locating
   */
  protected void initOutputLocators(Instances data, int[] indices) {
    if (indices == null) {
      m_OutputStringAtts = new StringLocator(data);
      m_OutputRelAtts    = new RelationalLocator(data);
    }
    else {
      m_OutputStringAtts = new StringLocator(data, indices);
      m_OutputRelAtts    = new RelationalLocator(data, indices);
    }
  }
  
  /**
   * Copies string/relational values contained in the instance copied to a new
   * dataset. The Instance must already be assigned to a dataset. This
   * dataset and the destination dataset must have the same structure.
   *
   * @param instance		the Instance containing the string/relational 
   * 				values to copy.
   * @param isInput		if true the input format and input attribute 
   * 				locators are used otherwise the output format 
   * 				and output locators
   */
  protected void copyValues(Instance instance, boolean isInput) {

    RelationalLocator.copyRelationalValues(
	instance, 
	(isInput) ? m_InputFormat : m_OutputFormat, 
	(isInput) ? m_InputRelAtts : m_OutputRelAtts);

    StringLocator.copyStringValues(
	instance, 
	(isInput) ? m_InputFormat : m_OutputFormat, 
	(isInput) ? m_InputStringAtts : m_OutputStringAtts);
  }

  /**
   * Takes string/relational values referenced by an Instance and copies them 
   * from a source dataset to a destination dataset. The instance references are
   * updated to be valid for the destination dataset. The instance may have the 
   * structure (i.e. number and attribute position) of either dataset (this
   * affects where references are obtained from). Only works if the number
   * of string/relational attributes is the same in both indices (implicitly 
   * these string/relational attributes should be semantically same but just 
   * with shifted positions).
   *
   * @param instance 		the instance containing references to strings/
   * 				relational values in the source dataset that 
   * 				will have references updated to be valid for 
   * 				the destination dataset.
   * @param instSrcCompat 	true if the instance structure is the same as 
   * 				the source, or false if it is the same as the 
   * 				destination (i.e. which of the string/relational 
   * 				attribute indices contains the correct locations 
   * 				for this instance).
   * @param srcDataset 		the dataset for which the current instance 
   * 				string/relational value references are valid 
   * 				(after any position mapping if needed)
   * @param destDataset 	the dataset for which the current instance 
   * 				string/relational value references need to be 
   * 				inserted (after any position mapping if needed)
   */
  protected void copyValues(Instance instance, boolean instSrcCompat,
                         Instances srcDataset, Instances destDataset) {

    RelationalLocator.copyRelationalValues(
	instance, instSrcCompat, 
	srcDataset, m_InputRelAtts,
	destDataset, m_OutputRelAtts);

    StringLocator.copyStringValues(
	instance, instSrcCompat, 
	srcDataset, m_InputStringAtts,
	getOutputFormat(), m_OutputStringAtts);
  }

  /**
   * This will remove all buffered instances from the inputformat dataset.
   * Use this method rather than getInputFormat().delete();
   */
  protected void flushInput() {

    if (    (m_InputStringAtts.getAttributeIndices().length > 0) 
	 || (m_InputRelAtts.getAttributeIndices().length > 0) ) {
      m_InputFormat = m_InputFormat.stringFreeStructure();
    } else {
      // This more efficient than new Instances(m_InputFormat, 0);
      m_InputFormat.delete();
    }
  }
  
  /**
   * tests the data whether the filter can actually handle it
   * 
   * @param instanceInfo	the data to test
   * @throws Exception		if the test fails
   */
  protected void testInputFormat(Instances instanceInfo) throws Exception {
    getCapabilities(instanceInfo).testWithFail(instanceInfo);
  }

  /**
   * Sets the format of the input instances. If the filter is able to
   * determine the output format before seeing any input instances, it
   * does so here. This default implementation clears the output format
   * and output queue, and the new batch flag is set. Overriders should
   * call <code>super.setInputFormat(Instances)</code>
   *
   * @param instanceInfo an Instances object containing the input instance
   * structure (any instances contained in the object are ignored - only the
   * structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the inputFormat can't be set successfully 
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {
    testInputFormat(instanceInfo);
    
    m_InputFormat = instanceInfo.stringFreeStructure();
    m_OutputFormat = null;
    m_OutputQueue = new Queue();
    m_NewBatch = true;
    m_FirstBatchDone = false;
    initInputLocators(m_InputFormat, null);
    return false;
  }

  /**
   * Gets the format of the output instances. This should only be called
   * after input() or batchFinished() has returned true. The relation
   * name of the output instances should be changed to reflect the
   * action of the filter (eg: add the filter name and options).
   *
   * @return an Instances object containing the output instance
   * structure only.
   * @throws NullPointerException if no input structure has been
   * defined (or the output format hasn't been determined yet) 
   */
  public Instances getOutputFormat() {

    if (m_OutputFormat == null) {
      throw new NullPointerException("No output format defined.");
    }
    return new Instances(m_OutputFormat, 0);
  }

  /**
   * Input an instance for filtering. Ordinarily the instance is
   * processed and made available for output immediately. Some filters
   * require all instances be read before producing output, in which
   * case output instances should be collected after calling
   * batchFinished(). If the input marks the start of a new batch, the
   * output queue is cleared. This default implementation assumes all
   * instance conversion will occur when batchFinished() is called.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws NullPointerException if the input format has not been
   * defined.
   * @throws Exception if the input instance was not of the correct 
   * format or if there was a problem with the filtering.  
   */
  public boolean input(Instance instance) throws Exception {

    if (m_InputFormat == null) {
      throw new NullPointerException("No input instance format defined");
    }
    if (m_NewBatch) {
      m_OutputQueue = new Queue();
      m_NewBatch = false;
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
   * re-assigned or new options have been set). This default
   * implementation assumes all instance processing occurs during
   * inputFormat() and input().
   *
   * @return true if there are instances pending output
   * @throws NullPointerException if no input structure has been defined,
   * @throws Exception if there was a problem finishing the batch.
   */
  public boolean batchFinished() throws Exception {

    if (m_InputFormat == null) {
      throw new NullPointerException("No input instance format defined");
    }
    flushInput();
    m_NewBatch = true;
    m_FirstBatchDone = true;
    return (numPendingOutput() != 0);
  }


  /**
   * Output an instance after filtering and remove from the output queue.
   *
   * @return the instance that has most recently been filtered (or null if
   * the queue is empty).
   * @throws NullPointerException if no output structure has been defined
   */
  public Instance output() {

    if (m_OutputFormat == null) {
      throw new NullPointerException("No output instance format defined");
    }
    if (m_OutputQueue.empty()) {
      return null;
    }
    Instance result = (Instance)m_OutputQueue.pop();
  
    //result = new PreferenceDenseInstance((DenseInstance)result);
    
    // Clear out references to old strings/relationals occasionally
    /*if (m_OutputQueue.empty() && m_NewBatch) {
      if (    (m_OutputStringAtts.getAttributeIndices().length > 0)
	   || (m_OutputRelAtts.getAttributeIndices().length > 0) ) {
        m_OutputFormat = m_OutputFormat.stringFreeStructure();
      }
    }*/
    return result;
  }
  
  /**
   * Output an instance after filtering but do not remove from the
   * output queue.
   *
   * @return the instance that has most recently been filtered (or null if
   * the queue is empty).
   * @throws NullPointerException if no input structure has been defined 
   */
  public Instance outputPeek() {

    if (m_OutputFormat == null) {
      throw new NullPointerException("No output instance format defined");
    }
    if (m_OutputQueue.empty()) {
      return null;
    }
    Instance result = (Instance)m_OutputQueue.peek();
    return result;
  }

  /**
   * Returns the number of instances pending output
   *
   * @return the number of instances  pending output
   * @throws NullPointerException if no input structure has been defined
   */
  public int numPendingOutput() {

    if (m_OutputFormat == null) {
      throw new NullPointerException("No output instance format defined");
    }
    return m_OutputQueue.size();
  }

  /**
   * Returns whether the output format is ready to be collected
   *
   * @return true if the output format is set
   */
  public boolean isOutputFormatDefined() {

    return (m_OutputFormat != null);
  }

  /**
   * Creates a deep copy of the given filter using serialization.
   *
   * @param model 	the filter to copy
   * @return 		a deep copy of the filter
   * @throws Exception 	if an error occurs
   */
  public static Filter makeCopy(Filter model) throws Exception {
    return (Filter)new SerializedObject(model).getObject();
  }

  /**
   * Creates a given number of deep copies of the given filter using 
   * serialization.
   * 
   * @param model 	the filter to copy
   * @param num 	the number of filter copies to create.
   * @return 		an array of filters.
   * @throws Exception 	if an error occurs
   */
  public static Filter[] makeCopies(Filter model, int num) throws Exception {

    if (model == null) {
      throw new Exception("No model filter set");
    }
    Filter[] filters = new Filter[num];
    SerializedObject so = new SerializedObject(model);
    for (int i = 0; i < filters.length; i++) {
      filters[i] = (Filter) so.getObject();
    }
    return filters;
  }
  
  /**
   * Filters an entire set of instances through a filter and returns
   * the new set. 
   *
   * @param data the data to be filtered
   * @param filter the filter to be used
   * @return the filtered set of data
   * @throws Exception if the filter can't be used successfully
   */
  public static Instances useFilter(Instances data,
				    Filter filter) throws Exception {
    /*
    System.err.println(filter.getClass().getName() 
                       + " in:" + data.numInstances());
    */
    for (int i = 0; i < data.numInstances(); i++) {
      filter.input(data.instance(i));
    }
    filter.batchFinished();
    Instances newData = filter.getOutputFormat();
    Instance processed;
    while ((processed = filter.output()) != null) {
      newData.add(processed);
    }

    /*
    System.err.println(filter.getClass().getName() 
                       + " out:" + newData.numInstances());
    */
    return newData;
  }

  /**
   * Returns a description of the filter, by default only the classname.
   * 
   * @return a string describing the filter
   */
  public String toString() {
    return this.getClass().getName();
  }
  
  /**
   * generates source code from the filter
   * 
   * @param filter the filter to output as source
   * @param className the name of the generated class
   * @param input the input data the header is generated for
   * @param output the output data the header is generated for
   * @return the generated source code
   * @throws Exception if source code cannot be generated
   */
  public static String wekaStaticWrapper(
      Sourcable filter, String className, Instances input, Instances output) 
    throws Exception {
    
    StringBuffer	result;
    int			i;
    int			n;
    
    result = new StringBuffer();
    
    result.append("// Generated with Weka " + Version.VERSION + "\n");
    result.append("//\n");
    result.append("// This code is public domain and comes with no warranty.\n");
    result.append("//\n");
    result.append("// Timestamp: " + new Date() + "\n");
    result.append("// Relation: " + input.relationName() + "\n");
    result.append("\n");
    
    result.append("package weka.filters;\n");
    result.append("\n");
    result.append("import weka.core.Attribute;\n");
    result.append("import weka.core.Capabilities;\n");
    result.append("import weka.core.Capabilities.Capability;\n");
    result.append("import weka.core.FastVector;\n");
    result.append("import weka.core.Instance;\n");
    result.append("import weka.core.Instances;\n");
    result.append("import weka.filters.Filter;\n");
    result.append("\n");
    result.append("public class WekaWrapper\n");
    result.append("  extends Filter {\n");

    // globalInfo
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns only the toString() method.\n");
    result.append("   *\n");
    result.append("   * @return a string describing the filter\n");
    result.append("   */\n");
    result.append("  public String globalInfo() {\n");
    result.append("    return toString();\n");
    result.append("  }\n");
    
    // getCapabilities
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns the capabilities of this filter.\n");
    result.append("   *\n");
    result.append("   * @return the capabilities\n");
    result.append("   */\n");
    result.append("  public Capabilities getCapabilities() {\n");
    result.append(((Filter) filter).getCapabilities().toSource("result", 4));
    result.append("    return result;\n");
    result.append("  }\n");

    // objectsToInstance
    result.append("\n");
    result.append("  /**\n");
    result.append("   * turns array of Objects into an Instance object\n");
    result.append("   *\n");
    result.append("   * @param obj	the Object array to turn into an Instance\n");
    result.append("   * @param format	the data format to use\n");
    result.append("   * @return		the generated Instance object\n");
    result.append("   */\n");
    result.append("  protected Instance objectsToInstance(Object[] obj, Instances format) {\n");
    result.append("    Instance		result;\n");
    result.append("    double[]		values;\n");
    result.append("    int		i;\n");
    result.append("\n");  
    result.append("    values = new double[obj.length];\n");
    result.append("\n");
    result.append("    for (i = 0 ; i < obj.length; i++) {\n");
    result.append("      if (obj[i] == null)\n");
    result.append("        values[i] = Instance.missingValue();\n");
    result.append("      else if (format.attribute(i).isNumeric())\n");
    result.append("        values[i] = (Double) obj[i];\n");
    result.append("      else if (format.attribute(i).isNominal() || format.attribute(i).isRanking())\n");
    result.append("        values[i] = format.attribute(i).indexOfValue((String) obj[i]);\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    // create new instance\n");
    result.append("    result = new Instance(1.0, values);\n");
    result.append("    result.setDataset(format);\n");
    result.append("\n");
    result.append("    return result;\n");
    result.append("  }\n");

    // instanceToObjects
    result.append("\n");
    result.append("  /**\n");
    result.append("   * turns the Instance object into an array of Objects\n");
    result.append("   *\n");
    result.append("   * @param inst	the instance to turn into an array\n");
    result.append("   * @return		the Object array representing the instance\n");
    result.append("   */\n");
    result.append("  protected Object[] instanceToObjects(Instance inst) {\n");
    result.append("    Object[]	result;\n");
    result.append("    int		i;\n");
    result.append("\n");  
    result.append("    result = new Object[inst.numAttributes()];\n");
    result.append("\n");
    result.append("    for (i = 0 ; i < inst.numAttributes(); i++) {\n");
    result.append("      if (inst.isMissing(i))\n");
    result.append("  	result[i] = null;\n");
    result.append("      else if (inst.attribute(i).isNumeric())\n");
    result.append("  	result[i] = inst.value(i);\n");
    result.append("      else\n");
    result.append("  	result[i] = inst.stringValue(i);\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    return result;\n");
    result.append("  }\n");

    // instancesToObjects
    result.append("\n");
    result.append("  /**\n");
    result.append("   * turns the Instances object into an array of Objects\n");
    result.append("   *\n");
    result.append("   * @param data	the instances to turn into an array\n");
    result.append("   * @return		the Object array representing the instances\n");
    result.append("   */\n");
    result.append("  protected Object[][] instancesToObjects(Instances data) {\n");
    result.append("    Object[][]	result;\n");
    result.append("    int		i;\n");
    result.append("\n");  
    result.append("    result = new Object[data.numInstances()][];\n");
    result.append("\n");  
    result.append("    for (i = 0; i < data.numInstances(); i++)\n");
    result.append("      result[i] = instanceToObjects(data.instance(i));\n");
    result.append("\n");  
    result.append("    return result;\n");
    result.append("  }\n");
    
    // setInputFormat
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Only tests the input data.\n");
    result.append("   *\n");
    result.append("   * @param instanceInfo the format of the data to convert\n");
    result.append("   * @return always true, to indicate that the output format can \n");
    result.append("   *         be collected immediately.\n");
    result.append("   */\n");
    result.append("  public boolean setInputFormat(Instances instanceInfo) throws Exception {\n");
    result.append("    super.setInputFormat(instanceInfo);\n");
    result.append("    \n");
    result.append("    // generate output format\n");
    result.append("    FastVector atts = new FastVector();\n");
    result.append("    FastVector attValues;\n");
    for (i = 0; i < output.numAttributes(); i++) {
      result.append("    // " + output.attribute(i).name() + "\n");
      if (output.attribute(i).isNumeric()) {
	result.append("    atts.addElement(new Attribute(\"" 
	    + output.attribute(i).name() + "\"));\n");
      }
      else if (output.attribute(i).isNominal() || output.attribute(i).isRanking()) {
	result.append("    attValues = new FastVector();\n");
	for (n = 0; n < output.attribute(i).numValues(); n++) {
	  result.append("    attValues.addElement(\"" + output.attribute(i).value(n) + "\");\n");
	}
	result.append("    atts.addElement(new Attribute(\"" 
	    + output.attribute(i).name() + "\", attValues));\n");
      }
      else {
	throw new UnsupportedAttributeTypeException(
	    "Attribute type '" + output.attribute(i).type() + "' (position " 
	    + (i+1) + ") is not supported!");
      }
    }
    result.append("    \n");
    result.append("    Instances format = new Instances(\"" + output.relationName() + "\", atts, 0);\n");
    result.append("    format.setClassIndex(" + output.classIndex() + ");\n");
    result.append("    setOutputFormat(format);\n");
    result.append("    \n");
    result.append("    return true;\n");
    result.append("  }\n");
    
    // input
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Directly filters the instance.\n");
    result.append("   *\n");
    result.append("   * @param instance the instance to convert\n");
    result.append("   * @return always true, to indicate that the output can \n");
    result.append("   *         be collected immediately.\n");
    result.append("   */\n");
    result.append("  public boolean input(Instance instance) throws Exception {\n");
    result.append("    Object[] filtered = " + className + ".filter(instanceToObjects(instance));\n");
    result.append("    push(objectsToInstance(filtered, getOutputFormat()));\n");
    result.append("    return true;\n");
    result.append("  }\n");
    
    // batchFinished
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Performs a batch filtering of the buffered data, if any available.\n");
    result.append("   *\n");
    result.append("   * @return true if instances were filtered otherwise false\n");
    result.append("   */\n");
    result.append("  public boolean batchFinished() throws Exception {\n");
    result.append("    if (getInputFormat() == null)\n");
    result.append("      throw new NullPointerException(\"No input instance format defined\");;\n");
    result.append("\n");
    result.append("    Instances inst = getInputFormat();\n");
    result.append("    if (inst.numInstances() > 0) {\n");
    result.append("      Object[][] filtered = " + className + ".filter(instancesToObjects(inst));\n");
    result.append("      for (int i = 0; i < filtered.length; i++) {\n");
    result.append("        push(objectsToInstance(filtered[i], getOutputFormat()));\n");
    result.append("      }\n");
    result.append("    }\n");
    result.append("\n");
    result.append("    flushInput();\n");
    result.append("    m_NewBatch = true;\n");
    result.append("    m_FirstBatchDone = true;\n");
    result.append("\n");
    result.append("    return (inst.numInstances() > 0);\n");
    result.append("  }\n");

    // toString
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Returns only the classnames and what filter it is based on.\n");
    result.append("   *\n");
    result.append("   * @return a short description\n");
    result.append("   */\n");
    result.append("  public String toString() {\n");
    result.append("    return \"Auto-generated filter wrapper, based on " 
	+ filter.getClass().getName() + " (generated with Weka " + Version.VERSION + ").\\n" 
	+ "\" + this.getClass().getName() + \"/" + className + "\";\n");
    result.append("  }\n");
    
    // main
    result.append("\n");
    result.append("  /**\n");
    result.append("   * Runs the filter from commandline.\n");
    result.append("   *\n");
    result.append("   * @param args the commandline arguments\n");
    result.append("   */\n");
    result.append("  public static void main(String args[]) {\n");
    result.append("    runFilter(new WekaWrapper(), args);\n");
    result.append("  }\n");
    result.append("}\n");

    // actual filter code
    result.append("\n");
    result.append(filter.toSource(className, input));
    
    return result.toString();
  }
  
  /**
   * Method for testing filters.
   *
   * @param filter the filter to use
   * @param options should contain the following arguments: <br/>
   * -i input_file <br/>
   * -o output_file <br/>
   * -c class_index <br/>
   * -z classname (for filters implementing weka.filters.Sourcable) <br/>
   * or -h for help on options
   * @throws Exception if something goes wrong or the user requests help on
   * command options
   */
  public static void filterFile(Filter filter, String [] options) 
    throws Exception {

    boolean debug = false;
    Instances data = null;
    DataSource input = null;
    PrintWriter output = null;
    boolean helpRequest;
    String sourceCode = "";

    try {
       helpRequest = Utils.getFlag('h', options);

      if (Utils.getFlag('d', options)) {
	debug = true;
      }
      String infileName = Utils.getOption('i', options);
      String outfileName = Utils.getOption('o', options); 
      String classIndex = Utils.getOption('c', options);
      if (filter instanceof Sourcable)
	sourceCode = Utils.getOption('z', options);
      
      if (filter instanceof OptionHandler) {
	((OptionHandler)filter).setOptions(options);
      }

      Utils.checkForRemainingOptions(options);
      if (helpRequest) {
	throw new Exception("Help requested.\n");
      }
      if (infileName.length() != 0) {
	input = new DataSource(infileName);
      } else {
	input = new DataSource(System.in);
      }
      if (outfileName.length() != 0) {
	output = new PrintWriter(new FileOutputStream(outfileName));
      } else { 
	output = new PrintWriter(System.out);
      }

      data = input.getStructure();
      if (classIndex.length() != 0) {
	if (classIndex.equals("first")) {
	  data.setClassIndex(0);
	} else if (classIndex.equals("last")) {
	  data.setClassIndex(data.numAttributes() - 1);
	} else {
	  data.setClassIndex(Integer.parseInt(classIndex) - 1);
	}
      }
    } catch (Exception ex) {
      String filterOptions = "";
      // Output the error and also the valid options
      if (filter instanceof OptionHandler) {
	filterOptions += "\nFilter options:\n\n";
	Enumeration enu = ((OptionHandler)filter).listOptions();
	while (enu.hasMoreElements()) {
	  Option option = (Option) enu.nextElement();
	  filterOptions += option.synopsis() + '\n'
	    + option.description() + "\n";
	}
      }

      String genericOptions = "\nGeneral options:\n\n"
	+ "-h\n"
	+ "\tGet help on available options.\n"
	+ "\t(use -b -h for help on batch mode.)\n"
	+ "-i <file>\n"
	+ "\tThe name of the file containing input instances.\n"
	+ "\tIf not supplied then instances will be read from stdin.\n"
	+ "-o <file>\n"
	+ "\tThe name of the file output instances will be written to.\n"
	+ "\tIf not supplied then instances will be written to stdout.\n"
	+ "-c <class index>\n"
	+ "\tThe number of the attribute to use as the class.\n"
	+ "\t\"first\" and \"last\" are also valid entries.\n"
	+ "\tIf not supplied then no class is assigned.\n";

      if (filter instanceof Sourcable) {
	genericOptions +=
	  "-z <class name>\n"
	  + "\tOutputs the source code representing the trained filter.\n";
      }
      
      throw new Exception('\n' + ex.getMessage()
			  + filterOptions+genericOptions);
    }
    
    if (debug) {
      System.err.println("Setting input format");
    }
    boolean printedHeader = false;
    if (filter.setInputFormat(data)) {
      if (debug) {
	System.err.println("Getting output format");
      }
      output.println(filter.getOutputFormat().toString());
      printedHeader = true;
    }
    
    // Pass all the instances to the filter
    Instance inst;
    while (input.hasMoreElements(data)) {
      inst = input.nextElement(data);
      if (debug) {
	System.err.println("Input instance to filter");
      }
      if (filter.input(inst)) {
	if (debug) {
	  System.err.println("Filter said collect immediately");
	}
	if (!printedHeader) {
	  throw new Error("Filter didn't return true from setInputFormat() "
			  + "earlier!");
	}
	if (debug) {
	  System.err.println("Getting output instance");
	}
	output.println(filter.output().toString());
      }
    }

    // Say that input has finished, and print any pending output instances
    if (debug) {
      System.err.println("Setting end of batch");
    }
    if (filter.batchFinished()) {
      if (debug) {
	System.err.println("Filter said collect output");
      }
      if (!printedHeader) {
	if (debug) {
	  System.err.println("Getting output format");
	}
	output.println(filter.getOutputFormat().toString());
      }
      if (debug) {
	System.err.println("Getting output instance");
      }
      while (filter.numPendingOutput() > 0) {
	output.println(filter.output().toString());
	if (debug){
	  System.err.println("Getting output instance");
	}
      }
    }
    if (debug) {
      System.err.println("Done");
    }
    
    if (output != null) {
      output.close();
    }
    
    if (sourceCode.length() != 0)
      System.out.println(
	  wekaStaticWrapper(
	      (Sourcable) filter, sourceCode, data, filter.getOutputFormat()));
  }

  /**
   * Method for testing filters ability to process multiple batches.
   *
   * @param filter the filter to use
   * @param options should contain the following arguments: <br/>
   * -i (first) input file <br/>
   * -o (first) output file <br/>
   * -r (second) input file <br/>
   * -s (second) output file <br/>
   * -c class_index <br/>
   * -z classname (for filters implementing weka.filters.Sourcable) <br/>
   * or -h for help on options
   * @throws Exception if something goes wrong or the user requests help on
   * command options
   */
  public static void batchFilterFile(Filter filter, String [] options) 
    throws Exception {

    Instances firstData = null;
    Instances secondData = null;
    DataSource firstInput = null;
    DataSource secondInput = null;
    PrintWriter firstOutput = null;
    PrintWriter secondOutput = null;
    boolean helpRequest;
    String sourceCode = "";

    try {
      helpRequest = Utils.getFlag('h', options);

      String fileName = Utils.getOption('i', options); 
      if (fileName.length() != 0) {
	firstInput = new DataSource(fileName);
      } else {
	throw new Exception("No first input file given.\n");
      }

      fileName = Utils.getOption('r', options); 
      if (fileName.length() != 0) {
	secondInput = new DataSource(fileName);
      } else {
	throw new Exception("No second input file given.\n");
      }

      fileName = Utils.getOption('o', options); 
      if (fileName.length() != 0) {
	firstOutput = new PrintWriter(new FileOutputStream(fileName));
      } else {
	firstOutput = new PrintWriter(System.out);
      }
      
      fileName = Utils.getOption('s', options); 
      if (fileName.length() != 0) {
	secondOutput = new PrintWriter(new FileOutputStream(fileName));
      } else {
	secondOutput = new PrintWriter(System.out);
      }
      String classIndex = Utils.getOption('c', options);
      if (filter instanceof Sourcable)
	sourceCode = Utils.getOption('z', options);

      if (filter instanceof OptionHandler) {
	((OptionHandler)filter).setOptions(options);
      }
      Utils.checkForRemainingOptions(options);
      
      if (helpRequest) {
	throw new Exception("Help requested.\n");
      }
      firstData = firstInput.getStructure();
      secondData = secondInput.getStructure();
      if (!secondData.equalHeaders(firstData)) {
	throw new Exception("Input file formats differ.\n" + secondData.equalHeadersMsg(firstData) + "\n");
      }
      if (classIndex.length() != 0) {
	if (classIndex.equals("first")) {
	  firstData.setClassIndex(0);
	  secondData.setClassIndex(0);
	} else if (classIndex.equals("last")) {
	  firstData.setClassIndex(firstData.numAttributes() - 1);
	  secondData.setClassIndex(secondData.numAttributes() - 1);
	} else {
	  firstData.setClassIndex(Integer.parseInt(classIndex) - 1);
	  secondData.setClassIndex(Integer.parseInt(classIndex) - 1);
	}
      }
    } catch (Exception ex) {
      String filterOptions = "";
      // Output the error and also the valid options
      if (filter instanceof OptionHandler) {
	filterOptions += "\nFilter options:\n\n";
	Enumeration enu = ((OptionHandler)filter).listOptions();
	while (enu.hasMoreElements()) {
	  Option option = (Option) enu.nextElement();
	  filterOptions += option.synopsis() + '\n'
	    + option.description() + "\n";
	}
      }

      String genericOptions = "\nGeneral options:\n\n"
	+ "-h\n"
	+ "\tGet help on available options.\n"
	+ "-i <filename>\n"
	+ "\tThe file containing first input instances.\n"
	+ "-o <filename>\n"
	+ "\tThe file first output instances will be written to.\n"
	+ "-r <filename>\n"
	+ "\tThe file containing second input instances.\n"
	+ "-s <filename>\n"
	+ "\tThe file second output instances will be written to.\n"
	+ "-c <class index>\n"
	+ "\tThe number of the attribute to use as the class.\n"
	+ "\t\"first\" and \"last\" are also valid entries.\n"
	+ "\tIf not supplied then no class is assigned.\n";

      if (filter instanceof Sourcable) {
	genericOptions +=
	  "-z <class name>\n"
	  + "\tOutputs the source code representing the trained filter.\n";
      }
      
      throw new Exception('\n' + ex.getMessage()
			  + filterOptions+genericOptions);
    }
    boolean printedHeader = false;
    if (filter.setInputFormat(firstData)) {
      firstOutput.println(filter.getOutputFormat().toString());
      printedHeader = true;
    }
    
    // Pass all the instances to the filter
    Instance inst;
    while (firstInput.hasMoreElements(firstData)) {
      inst = firstInput.nextElement(firstData);
      if (filter.input(inst)) {
	if (!printedHeader) {
	  throw new Error("Filter didn't return true from setInputFormat() "
			  + "earlier!");
	}
	firstOutput.println(filter.output().toString());
      }
    }
    
    // Say that input has finished, and print any pending output instances
    if (filter.batchFinished()) {
      if (!printedHeader) {
	firstOutput.println(filter.getOutputFormat().toString());
      }
      while (filter.numPendingOutput() > 0) {
	firstOutput.println(filter.output().toString());
      }
    }
    
    if (firstOutput != null) {
      firstOutput.close();
    }    
    printedHeader = false;
    if (filter.isOutputFormatDefined()) {
      secondOutput.println(filter.getOutputFormat().toString());
      printedHeader = true;
    }
    // Pass all the second instances to the filter
    while (secondInput.hasMoreElements(secondData)) {
      inst = secondInput.nextElement(secondData);
      if (filter.input(inst)) {
	if (!printedHeader) {
	  throw new Error("Filter didn't return true from"
			  + " isOutputFormatDefined() earlier!");
	}
	secondOutput.println(filter.output().toString());
      }
    }
    
    // Say that input has finished, and print any pending output instances
    if (filter.batchFinished()) {
      if (!printedHeader) {
	secondOutput.println(filter.getOutputFormat().toString());
      }
      while (filter.numPendingOutput() > 0) {
	secondOutput.println(filter.output().toString());
      }
    }
    if (secondOutput != null) {
      secondOutput.close();
    }

    if (sourceCode.length() != 0)
      System.out.println(
	  wekaStaticWrapper(
	      (Sourcable) filter, sourceCode, firstData, filter.getOutputFormat()));
  }

  /**
   * runs the filter instance with the given options.
   * 
   * @param filter	the filter to run
   * @param options	the commandline options
   */
  public static void runFilter(Filter filter, String[] options) {
    try {
      if (Utils.getFlag('b', options)) {
	Filter.batchFilterFile(filter, options);
      } else {
	Filter.filterFile(filter, options);
      }
    } catch (Exception e) {
      if (    (e.toString().indexOf("Help requested") == -1) 
	   && (e.toString().indexOf("Filter options") == -1) )
	e.printStackTrace();
      else
	System.err.println(e.getMessage());
    }
  }
  
  /**
   * Main method for testing this class.
   *
   * @param args should contain arguments to the filter: use -h for help
   */
  public static void main(String [] args) {
    
    try {
      if (args.length == 0) {
        throw new Exception("First argument must be the class name of a Filter");
      }
      String fname = args[0];
      Filter f = (Filter)Class.forName(fname).newInstance();
      args[0] = "";
      runFilter(f, args);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

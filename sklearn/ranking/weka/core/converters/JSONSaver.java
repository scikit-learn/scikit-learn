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
 * JSONSaver.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.SingleIndex;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.json.JSONInstances;
import weka.core.json.JSONNode;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.Vector;
import java.util.zip.GZIPOutputStream;

/**
 <!-- globalinfo-start -->
 * Writes to a destination that is in JSON format.<br/>
 * The data can be compressed with gzip, in order to save space.<br/>
 * <br/>
 * For more information, see JSON homepage:<br/>
 * http://www.json.org/
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -i &lt;the input file&gt;
 *  The input file</pre>
 * 
 * <pre> -o &lt;the output file&gt;
 *  The output file</pre>
 * 
 * <pre> -C &lt;class index&gt;
 *  The class index (first and last are valid as well).
 *  (default: last)</pre>
 * 
 * <pre> -compress
 *  Compresses the data (uses '.json.gz' as extension instead of '.json')
 *  (default: off)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see Saver
 */
public class JSONSaver 
  extends AbstractFileSaver 
  implements BatchConverter {
  
  /** for serialization. */
  private static final long serialVersionUID = -1047134047244534557L;

  /** the class index. */
  protected SingleIndex m_ClassIndex = new SingleIndex(); 
  
  /** whether to compress the output. */
  protected boolean m_CompressOutput = false;
  
  /**
   * Constructor.
   */
  public JSONSaver(){
    resetOptions();
  }
  
  /**
   * Returns a string describing this Saver.
   * 
   * @return 		a description of the Saver suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Writes to a destination that is in JSON format.\n"
      + "The data can be compressed with gzip, in order to save space.\n\n"
      + "For more information, see JSON homepage:\n"
      + "http://www.json.org/";
  }
  
  /**
   * Returns an enumeration describing the available options.
   * 
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option>      result;
    
    result = new Vector<Option>();
    
    Enumeration en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement((Option)en.nextElement());
    
    result.addElement(
        new Option(
            "\tThe class index (first and last are valid as well).\n"
            + "\t(default: last)",
            "C", 1, "-C <class index>"));
    
    result.addElement(
        new Option(
            "\tCompresses the data (uses '" 
            + JSONLoader.FILE_EXTENSION_COMPRESSED 
            + "' as extension instead of '" 
            + JSONLoader.FILE_EXTENSION + "')\n"
            + "\t(default: off)",
            "compress", 0, "-compress"));
    
    return result.elements();
  }
  
  /**
   * returns the options of the current setup.
   *
   * @return		the current options
   */
  public String[] getOptions(){
    int       	i;
    Vector<String>    	result;
    String[]  	options;

    result = new Vector<String>();

    if (getClassIndex().length() != 0) {
      result.add("-C");
      result.add(getClassIndex());
    }

    if (getCompressOutput())
      result.add("-compress");
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    return (String[]) result.toArray(new String[result.size()]);	  
  }

  /**
   * Parses the options for this object. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -i &lt;the input file&gt;
   *  The input file</pre>
   * 
   * <pre> -o &lt;the output file&gt;
   *  The output file</pre>
   * 
   * <pre> -C &lt;class index&gt;
   *  The class index (first and last are valid as well).
   *  (default: last)</pre>
   * 
   * <pre> -compress
   *  Compresses the data (uses '.json.gz' as extension instead of '.json')
   *  (default: off)</pre>
   * 
   <!-- options-end -->
   *
   * @param options	the options to use
   * @throws Exception	if setting of options fails
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;

    tmpStr = Utils.getOption('C', options);
    if (tmpStr.length() != 0)
      setClassIndex(tmpStr);
    else
      setClassIndex("last");

    setCompressOutput(Utils.getFlag("compress", options));
    
    super.setOptions(options);
  }
  
  /**
   * Returns a description of the file type.
   *
   * @return a short file description
   */
  public String getFileDescription() {
    return "JSON data files";
  }

  /**
   * Gets all the file extensions used for this type of file.
   *
   * @return the file extensions
   */
  public String[] getFileExtensions() {
    return new String[]{JSONLoader.FILE_EXTENSION, JSONLoader.FILE_EXTENSION_COMPRESSED};
  }
  
  /** 
   * Sets the destination file.
   * 
   * @param outputFile the destination file.
   * @throws IOException throws an IOException if file cannot be set
   */
  public void setFile(File outputFile) throws IOException  {
    if (outputFile.getAbsolutePath().endsWith(JSONLoader.FILE_EXTENSION_COMPRESSED))
      setCompressOutput(true);
    
    super.setFile(outputFile);
  }
  
  /**
   * Resets the Saver.
   */
  public void resetOptions() {
    super.resetOptions();
    
    if (getCompressOutput())
      setFileExtension(JSONLoader.FILE_EXTENSION_COMPRESSED);
    else
      setFileExtension(JSONLoader.FILE_EXTENSION);
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String classIndexTipText() {
    return "Sets the class index (\"first\" and \"last\" are valid values)";
  }

  /**
   * Get the index of the class attribute.
   *
   * @return 		the index of the class attribute
   */
  public String getClassIndex() {
    return m_ClassIndex.getSingleIndex();
  }

  /**
   * Sets index of the class attribute.
   *
   * @param value 	the index of the class attribute
   */
  public void setClassIndex(String value) {
    m_ClassIndex.setSingleIndex(value);
  }

  /**
   * Returns the tip text for this property.
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String compressOutputTipText() {
    return "Optional compression of the output data";
  }

  /**
   * Gets whether the output data is compressed.
   *
   * @return 		true if the output data is compressed
   */
  public boolean getCompressOutput() {
    return m_CompressOutput;
  }

  /**
   * Sets whether to compress the output.
   *
   * @param value 	if truee the output will be compressed
   */
  public void setCompressOutput(boolean value) {
    m_CompressOutput = value;
  }

  /** 
   * Returns the Capabilities of this saver.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    
    // attributes
    result.enable(Capability.NOMINAL_ATTRIBUTES);
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    result.enable(Capability.DATE_ATTRIBUTES);
    result.enable(Capability.STRING_ATTRIBUTES);
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enable(Capability.NOMINAL_CLASS);
    result.enable(Capability.NUMERIC_CLASS);
    result.enable(Capability.DATE_CLASS);
    result.enable(Capability.STRING_CLASS);
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);
    
    return result;
  }
  
  /**
   * Sets instances that should be stored.
   *
   * @param instances 	the instances
   */
  public void setInstances(Instances instances) {
    if (m_ClassIndex.getSingleIndex().length() != 0) {
      m_ClassIndex.setUpper(instances.numAttributes() - 1);
      instances.setClassIndex(m_ClassIndex.getIndex());
    }
    
    super.setInstances(instances);
  }
  
  /** 
   * Sets the destination output stream.
   * 
   * @param output 		the output stream.
   * @throws IOException 	throws an IOException if destination cannot be set
   */
  public void setDestination(OutputStream output) throws IOException {
    if (getCompressOutput())
      super.setDestination(new GZIPOutputStream(output));
    else
      super.setDestination(output);
  }
  
  /**
   * Writes a Batch of instances.
   * 
   * @throws IOException 	throws IOException if saving in batch mode 
   * 				is not possible
   */
  public void writeBatch() throws IOException {
    if (getInstances() == null)
      throw new IOException("No instances to save");
    
    if (getRetrieval() == INCREMENTAL)
      throw new IOException("Batch and incremental saving cannot be mixed.");
    
    setRetrieval(BATCH);
    setWriteMode(WRITE);
    
    PrintWriter outW;
    if ((retrieveFile() == null) && (getWriter() == null))
      outW = new PrintWriter(System.out);
    else
      outW = new PrintWriter(getWriter());
    
    JSONNode json = JSONInstances.toJSON(getInstances());
    StringBuffer buffer = new StringBuffer();
    json.toString(buffer);
    outW.println(buffer.toString());
    outW.flush();
    
    if (getWriter() != null)
      outW.close();
    
    setWriteMode(WAIT);
    outW = null;
    resetWriter();
    setWriteMode(CANCEL);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5953 $");
  }
  
  /**
   * Main method.
   *
   * @param args 	should contain the options of a Saver.
   */
  public static void main(String[] args) {
    runFileSaver(new JSONSaver(), args);
  }
}

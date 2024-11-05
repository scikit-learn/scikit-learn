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
 * XRFFSaver.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
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
import weka.core.xml.XMLInstances;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.Vector;
import java.util.zip.GZIPOutputStream;

/**
 <!-- globalinfo-start -->
 * Writes to a destination that is in the XML version of the ARFF format. The data can be compressed with gzip, in order to save space.
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
 *  Compresses the data (uses '.xrff.gz' as extension instead of '.xrff')
 *  (default: off)</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see Saver
 */
public class XRFFSaver 
  extends AbstractFileSaver 
  implements BatchConverter {
  
  /** for serialization */
  private static final long serialVersionUID = -7226404765213522043L;

  /** the class index */
  protected SingleIndex m_ClassIndex = new SingleIndex(); 

  /** the generated XML document */
  protected XMLInstances m_XMLInstances;
  
  /** whether to compress the output */
  protected boolean m_CompressOutput = false;
  
  /**
   * Constructor
   */
  public XRFFSaver(){
    resetOptions();
  }
  
  /**
   * Returns a string describing this Saver
   * 
   * @return 		a description of the Saver suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Writes to a destination that is in the XML version of the ARFF format. "
      + "The data can be compressed with gzip, in order to save space.";
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
            + XRFFLoader.FILE_EXTENSION_COMPRESSED 
            + "' as extension instead of '" 
            + XRFFLoader.FILE_EXTENSION + "')\n"
            + "\t(default: off)",
            "compress", 0, "-compress"));
    
    return result.elements();
  }
  
  /**
   * returns the options of the current setup
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
   *  Compresses the data (uses '.xrff.gz' as extension instead of '.xrff')
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
    return "XRFF data files";
  }

  /**
   * Gets all the file extensions used for this type of file
   *
   * @return the file extensions
   */
  public String[] getFileExtensions() {
    return new String[]{XRFFLoader.FILE_EXTENSION, XRFFLoader.FILE_EXTENSION_COMPRESSED};
  }
  
  /** 
   * Sets the destination file.
   * 
   * @param outputFile the destination file.
   * @throws IOException throws an IOException if file cannot be set
   */
  public void setFile(File outputFile) throws IOException  {
    if (outputFile.getAbsolutePath().endsWith(XRFFLoader.FILE_EXTENSION_COMPRESSED))
      setCompressOutput(true);
    
    super.setFile(outputFile);
  }
  
  /**
   * Resets the Saver 
   */
  public void resetOptions() {
    super.resetOptions();
    
    if (getCompressOutput())
      setFileExtension(XRFFLoader.FILE_EXTENSION_COMPRESSED);
    else
      setFileExtension(XRFFLoader.FILE_EXTENSION);
    
    try {
      m_XMLInstances = new XMLInstances();
    }
    catch (Exception e) {
      m_XMLInstances = null;
    }
  }

  /**
   * Returns the tip text for this property
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
   * Returns the tip text for this property
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
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enableAllClasses();
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
   * Writes a Batch of instances
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

    // generate XML
    m_XMLInstances.setInstances(getInstances());
    
    if ((retrieveFile() == null) && (getWriter() == null)) {
      System.out.println(m_XMLInstances.toString());
      setWriteMode(WAIT);
    }
    else {
      PrintWriter outW = new PrintWriter(getWriter());
      outW.println(m_XMLInstances.toString());
      outW.flush();
      outW.close();
      setWriteMode(WAIT);
      outW = null;
      resetWriter();
      setWriteMode(CANCEL);
    }
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
    runFileSaver(new XRFFSaver(), args);
  }
}

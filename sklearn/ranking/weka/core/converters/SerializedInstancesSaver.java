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
 *    SerializedInstancesSaver.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Capabilities;
import weka.core.RevisionUtils;
import weka.core.Capabilities.Capability;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStream;

/**
 <!-- globalinfo-start -->
 * Serializes the instances to a file with extension bsi.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -i &lt;the input file&gt;
 * The input file</pre>
 * 
 * <pre> -o &lt;the output file&gt;
 * The output file</pre>
 * 
 <!-- options-end -->
 *
 * @author Stefan Mutter (mutter@cs.waikato.ac.nz)
 * @version $Revision: 4912 $
 * @see Saver
 */
public class SerializedInstancesSaver 
  extends AbstractFileSaver 
  implements BatchConverter {

  /** for serialization. */
  static final long serialVersionUID = -7717010648500658872L;
  
  /** the output stream. */
  protected ObjectOutputStream m_objectstream;
  
  /** Constructor. */  
  public SerializedInstancesSaver(){
      resetOptions();
  }
    
  /**
   * Returns a string describing this Saver.
   * 
   * @return a description of the Saver suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Serializes the instances to a file with extension bsi.";
  }
 
  /**
   * Returns a description of the file type.
   *
   * @return a short file description
   */
  public String getFileDescription() {
    return "Binary serialized instances";
  }

  /**
   * Resets the Saver.
   */
  public void resetOptions() {

    super.resetOptions();
    setFileExtension(".bsi");
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
   * Resets the writer, setting writer and objectstream to null.
   */  
  public void resetWriter() {
    super.resetWriter();
    
    m_objectstream = null;
  }
  
  /**
   * Sets the destination output stream.
   * 
   * @param output the output stream.
   * @throws IOException throws an IOException if destination cannot be set
   */
  public void setDestination(OutputStream output) throws IOException {
    super.setDestination(output);
    
    m_objectstream = new ObjectOutputStream(output);
  }
  
  /** 
   * Writes a Batch of instances.
   * 
   * @throws IOException throws IOException if saving in batch mode is not possible
   */
  public void writeBatch() throws IOException {
    if(getRetrieval() == INCREMENTAL)
      throw new IOException("Batch and incremental saving cannot be mixed.");
    
    if(getInstances() == null)
      throw new IOException("No instances to save");
    
    setRetrieval(BATCH);
    
    if (m_objectstream == null)
      throw new IOException("No output for serialization.");

    setWriteMode(WRITE);
    m_objectstream.writeObject(getInstances());
    m_objectstream.flush();
    m_objectstream.close();
    setWriteMode(WAIT);
    resetWriter();
    setWriteMode(CANCEL);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 4912 $");
  }

  /**
   * Main method.
   *
   * @param args should contain the options of a Saver.
   */
  public static void main(String[] args) {
    runFileSaver(new SerializedInstancesSaver(), args);
  }
}

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
 *    ArffSaver.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.Capabilities.Capability;

import java.io.IOException;
import java.io.PrintWriter;

/**
 * Writes to a destination in arff text format. <p/>
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
 * @version $Revision: 5953 $
 * @see Saver
 */
public class ArffSaver 
  extends AbstractFileSaver 
  implements BatchConverter, IncrementalConverter {

  /** for serialization */
  static final long serialVersionUID = 2223634248900042228L;    
  
  /** Constructor */  
  public ArffSaver(){
  
      resetOptions();
  }
   
   
  /**
   * Returns a string describing this Saver
   * @return a description of the Saver suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Writes to a destination that is in arff (attribute relation file format) "
      +"format. ";
  }

  
  /**
   * Returns a description of the file type.
   *
   * @return a short file description
   */
  public String getFileDescription() {
    return "Arff data files";
  }

  /**
   * Resets the Saver 
   */
  public void resetOptions() {

    super.resetOptions();
    setFileExtension(".arff");
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
  
  /** Saves an instances incrementally. Structure has to be set by using the
   * setStructure() method or setInstances() method.
   * @param inst the instance to save
   * @throws IOException throws IOEXception if an instance cannot be saved incrementally.
   */  
  public void writeIncremental(Instance inst) throws IOException{
  
      int writeMode = getWriteMode();
      Instances structure = getInstances();
      PrintWriter outW = null;
      
      if(getRetrieval() == BATCH || getRetrieval() == NONE)
          throw new IOException("Batch and incremental saving cannot be mixed.");
      if(getWriter() != null)
          outW = new PrintWriter(getWriter());
          
      if(writeMode == WAIT){
        if(structure == null){
            setWriteMode(CANCEL);
            if(inst != null)
                System.err.println("Structure(Header Information) has to be set in advance");
        }
        else
            setWriteMode(STRUCTURE_READY);
        writeMode = getWriteMode();
      }
      if(writeMode == CANCEL){
          if(outW != null)
              outW.close();
          cancel();
      }
      if(writeMode == STRUCTURE_READY){
          setWriteMode(WRITE);
          //write header
          Instances header = new Instances(structure,0);
          if(retrieveFile() == null || outW == null)
              System.out.println(header.toString());
          else{
              outW.print(header.toString());
              outW.print("\n");
              outW.flush();
          }
          writeMode = getWriteMode();
      }
      if(writeMode == WRITE){
          if(structure == null)
              throw new IOException("No instances information available.");
          if(inst != null){
          //write instance 
              if(retrieveFile() == null || outW == null)
                System.out.println(inst);
              else{
                outW.println(inst);
                m_incrementalCounter++;
                //flush every 100 instances
                if(m_incrementalCounter > 100){
                    m_incrementalCounter = 0;
                    outW.flush();
                }
              }
          }
          else{
          //close
              if(outW != null){
                outW.flush();
                outW.close();
              }
              m_incrementalCounter = 0;
              resetStructure();
              outW = null;
              resetWriter();
          }
      }
  }
  
  /** Writes a Batch of instances
   * @throws IOException throws IOException if saving in batch mode is not possible
   */
  public void writeBatch() throws IOException {
  
      if(getInstances() == null)
          throw new IOException("No instances to save");
      if(getRetrieval() == INCREMENTAL)
          throw new IOException("Batch and incremental saving cannot be mixed.");
      setRetrieval(BATCH);
      setWriteMode(WRITE);
      if(retrieveFile() == null && getWriter() == null){
          System.out.println((getInstances()).toString());
          setWriteMode(WAIT);
          return;
      }
      PrintWriter outW = new PrintWriter(getWriter());
      outW.print((getInstances()).toString());
      outW.flush();
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
   * @param args should contain the options of a Saver.
   */
  public static void main(String[] args) {
    runFileSaver(new ArffSaver(), args);
  }
}

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
 *    Saver.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;

/** 
 * Interface to something that can save Instances to an output destination in some
 * format.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author Stefan Mutter (mutter@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public interface Saver
  extends Serializable, RevisionHandler {
    
    /** The retrieval modes */
  static final int NONE = 0;
  static final int BATCH = 1;
  static final int INCREMENTAL = 2;
  


  /*@ public model instance boolean model_structureDetermined
    @   initially: model_structureDetermined == false;
    @*/

  /*@ public model instance boolean model_sourceSupplied
    @   initially: model_sourceSupplied == false;
    @*/

  /**
   * Resets the Saver object and sets the destination to be 
   * the supplied File object.
   *
   * @param file the File
   * @exception IOException if an error occurs
   * support loading from a File.
   *
   * <pre><jml>
   *    public_normal_behavior
   *      requires: file != null
   *                && (* file exists *);
   *      modifiable: model_sourceSupplied, model_structureDetermined;
   *      ensures: model_sourceSupplied == true 
   *               && model_structureDetermined == false;
   *  also
   *    public_exceptional_behavior
   *      requires: file == null
   *                || (* file does not exist *);
   *    signals: (IOException);
   * </jml></pre>
   */
  void setDestination(File file) throws IOException;

  /** Resets the Saver object and sets the destination to be
   * the supplied InputStream.
   * @param output the output stream
   * @exception IOException if this Loader doesn't
   * support loading from a File.
   */
  void setDestination(OutputStream output) throws IOException;
  
  /** Sets the retrieval mode
   * @param mode an integer representing a retrieval mode
   */  
  void setRetrieval(int mode);
  
  /** Gets the file extension
   * @return a string conatining the file extension (including the '.')
   * @throws Exception exception if a Saver not implementing FileSourcedConverter is used.
   */  
  String getFileExtension() throws Exception;
  
  /** Sets the output file
   * @param file the output file
   * @throws IOException exception if new output file cannot be set
   */  
  void setFile(File file)throws IOException;
  
  /** Sets the file prefix.
   * This method is used in the KnowledgeFlow GUI.
   * @param prefix the prefix of the file name
   * @throws Exception exception if a Saver not implementing FileSourcedConverter is used.
   */  
  void setFilePrefix(String prefix) throws Exception;
  
  /** Gets the file prefix
   * This method is used in the KnowledgeFlow GUI.
   * @return the prefix of the file name
   * @throws Exception exception if a Saver not implementing FileSourcedConverter is used.
   */  
  String filePrefix() throws Exception;
  
  /** Sets the directory of the output file.
   * This method is used in the KnowledgeFlow GUI.
   * @param dir a string containing the path and name of the directory
   * @throws IOException exception if a Saver not implementing FileSourcedConverter is used.
   */  
  void setDir(String dir) throws IOException;
  
  /** Sets the file prefix and the directory.
   * This method is used in the KnowledgeFlow GUI.
   * @param relationName the name of the realtion to be saved
   * @param add additional String for the file name
   * @throws IOException exception if a Saver not implementing FileSourcedConverter is used.
   */  
  public void setDirAndPrefix(String relationName, String add) throws IOException; 
  
  /** Gets the driectory of the output file
   * This method is used in the KnowledgeFlow GUI.
   * @return the directory as a string
   * @throws IOException exception if a Saver not implementing FileSourcedConverter is used.
   */  
  String retrieveDir() throws IOException;
  
  /** Sets the instances to be saved
   * @param instances the instances
   */  
  void setInstances(Instances instances);

  /** Writes to a destination in batch mode
   * @throws IOException throws exection if writting in batch mode is not possible
   */  
  void writeBatch() throws IOException;
  
  /** Writes to a destination in incremental mode.
   * If the instance is null, the outputfile will be closed.
   * @param inst the instance to write, if null the output file is closed
   * @throws IOException throws exception if incremental writting is not possible
   */  
  void writeIncremental(Instance inst) throws IOException;
  
  /** Gets the write mode
   * @return an integer representing the write mode
   */  
  public int getWriteMode();
  
}






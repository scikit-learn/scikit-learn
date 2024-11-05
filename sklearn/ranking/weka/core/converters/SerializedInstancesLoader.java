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
 *    SerializedInstancesLoader.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;

/**
 <!-- globalinfo-start -->
 * Reads a source that contains serialized Instances.
 * <p/>
 <!-- globalinfo-end -->
 * 
 * @author <a href="mailto:len@reeltwo.com">Len Trigg</a>
 * @version $Revision: 5953 $
 * @see Loader
 */
public class SerializedInstancesLoader 
  extends AbstractFileLoader 
  implements BatchConverter, IncrementalConverter {

  /** for serialization */
  static final long serialVersionUID = 2391085836269030715L;
  
  /** the file extension */
  public static String FILE_EXTENSION = 
    Instances.SERIALIZED_OBJ_FILE_EXTENSION;
  
  /** Holds the structure (header) of the data set. */
  protected Instances m_Dataset = null;

  /** The current index position for incremental reading */
  protected int m_IncrementalIndex = 0;
  
  /**
   * Returns a string describing this object
   * 
   * @return a description of the classifier suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Reads a source that contains serialized Instances.";
  }

  /** Resets the Loader ready to read a new data set */
  public void reset() {

    m_Dataset = null;
    m_IncrementalIndex = 0;
  }
  
  /**
   * Get the file extension used for arff files
   *
   * @return the file extension
   */
  public String getFileExtension() {
    return FILE_EXTENSION;
  }

  /**
   * Gets all the file extensions used for this type of file
   *
   * @return the file extensions
   */
  public String[] getFileExtensions() {
    return new String[]{getFileExtension()};
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
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied InputStream.
   *
   * @param in the source InputStream.
   * @throws IOException if there is a problem with IO
   */
  public void setSource(InputStream in) throws IOException {

    ObjectInputStream oi = new ObjectInputStream(new BufferedInputStream(in));
    try {
      m_Dataset = (Instances)oi.readObject();
    } catch (ClassNotFoundException ex) {
      throw new IOException("Could not deserialize instances from this source.");
    }

    // close the stream
    oi.close();
  }

  /**
   * Determines and returns (if possible) the structure (internally the 
   * header) of the data set as an empty set of instances.
   *
   * @return the structure of the data set as an empty set of Instances
   * @throws IOException if an error occurs
   */
  public Instances getStructure() throws IOException {

    if (m_Dataset == null) {
      throw new IOException("No source has been specified");
    }

    // We could cache a structure-only if getStructure is likely to be called
    // many times.
    return new Instances(m_Dataset, 0);
  }

  /**
   * Return the full data set. If the structure hasn't yet been determined
   * by a call to getStructure then method should do so before processing
   * the rest of the data set.
   *
   * @return the structure of the data set as an empty set of Instances
   * @throws IOException if there is no source or parsing fails
   */
  public Instances getDataSet() throws IOException {

    if (m_Dataset == null) {
      throw new IOException("No source has been specified");
    }

    return m_Dataset;
  }

  /**
   * Read the data set incrementally---get the next instance in the data 
   * set or returns null if there are no
   * more instances to get. If the structure hasn't yet been 
   * determined by a call to getStructure then method should do so before
   * returning the next instance in the data set.
   *
   * @param structure ignored
   * @return the next instance in the data set as an Instance object or null
   * if there are no more instances to be read
   * @throws IOException if there is an error during parsing
   */
  public Instance getNextInstance(Instances structure) throws IOException {

    if (m_Dataset == null) {
      throw new IOException("No source has been specified");
    }

    // We have to fake this method, since we can only deserialize an entire
    // dataset at a time.
    if (m_IncrementalIndex == m_Dataset.numInstances()) {
      return null;
    }
 
    return m_Dataset.instance(m_IncrementalIndex++);
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
   * @param args should contain the name of an input file.
   */
  public static void main(String[] args) {
    runFileLoader(new SerializedInstancesLoader(), args);
  }
}

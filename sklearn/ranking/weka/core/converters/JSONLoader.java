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
 * JSONLoader.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.json.JSONInstances;
import weka.core.json.JSONNode;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.util.zip.GZIPInputStream;

/**
 <!-- globalinfo-start -->
 * Reads a source that is in the JSON format.<br/>
 * It automatically decompresses the data if the extension is '.json.gz'.<br/>
 * <br/>
 * For more information, see JSON homepage:<br/>
 * http://www.json.org/
 * <p/>
 <!-- globalinfo-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5784 $
 * @see Loader
 */
public class JSONLoader 
  extends AbstractFileLoader 
  implements BatchConverter, URLSourcedLoader {

  /** for serialization. */
  private static final long serialVersionUID = 3764533621135196582L;

  /** the file extension. */
  public static String FILE_EXTENSION = ".json";

  /** the extension for compressed files. */
  public static String FILE_EXTENSION_COMPRESSED = FILE_EXTENSION + ".gz";

  /** the url. */
  protected String m_URL = "http://";

  /** The reader for the source file. */
  protected transient Reader m_sourceReader = null;

  /** the loaded JSON object. */
  protected JSONNode m_JSON;
  
  /**
   * Returns a string describing this Loader.
   * 
   * @return 		a description of the Loader suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Reads a source that is in the JSON format.\n"
      + "It automatically decompresses the data if the extension is '" 
      + FILE_EXTENSION_COMPRESSED + "'.\n\n"
      + "For more information, see JSON homepage:\n"
      + "http://www.json.org/";
  }

  /**
   * Get the file extension used for JSON files.
   *
   * @return 		the file extension
   */
  public String getFileExtension() {
    return FILE_EXTENSION;
  }

  /**
   * Gets all the file extensions used for this type of file.
   *
   * @return the file extensions
   */
  public String[] getFileExtensions() {
    return new String[]{FILE_EXTENSION, FILE_EXTENSION_COMPRESSED};
  }

  /**
   * Returns a description of the file type.
   *
   * @return 		a short file description
   */
  public String getFileDescription() {
    return "JSON Instances files";
  }

  /**
   * Resets the Loader ready to read a new data set.
   * 
   * @throws IOException 	if something goes wrong
   */
  public void reset() throws IOException {
    m_structure = null;
    m_JSON      = null;

    setRetrieval(NONE);
    
    if (m_File != null) {
      setFile(new File(m_File));
    }
    else if ((m_URL != null) && !m_URL.equals("http://")) {
      setURL(m_URL);
    }
  }

  /**
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied File object.
   *
   * @param file 		the source file.
   * @throws IOException 	if an error occurs
   */
  public void setSource(File file) throws IOException {
    m_structure = null;
    m_JSON      = null;
    
    setRetrieval(NONE);

    if (file == null)
      throw new IOException("Source file object is null!");

    try {
      if (file.getName().endsWith(FILE_EXTENSION_COMPRESSED))
	setSource(new GZIPInputStream(new FileInputStream(file)));
      else
	setSource(new FileInputStream(file));
    }
    catch (FileNotFoundException ex) {
      throw new IOException("File not found");
    }
    
    m_sourceFile = file;
    m_File       = file.getAbsolutePath();
  }

  /**
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied url.
   *
   * @param url 	the source url.
   * @throws IOException 	if an error occurs
   */
  public void setSource(URL url) throws IOException {
    m_structure = null;
    m_JSON      = null;
    
    setRetrieval(NONE);
    
    setSource(url.openStream());

    m_URL = url.toString();
  }
  
  /**
   * Set the url to load from.
   *
   * @param url 		the url to load from
   * @throws IOException 	if the url can't be set.
   */
  public void setURL(String url) throws IOException {
    m_URL = url;
    setSource(new URL(url));
  }

  /**
   * Return the current url.
   *
   * @return the current url
   */
  public String retrieveURL() {
    return m_URL;
  }

  /**
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied InputStream.
   *
   * @param in 			the source InputStream.
   * @throws IOException 	if initialization of reader fails.
   */
  public void setSource(InputStream in) throws IOException {
    m_File = (new File(System.getProperty("user.dir"))).getAbsolutePath();
    m_URL  = "http://";

    m_sourceReader = new BufferedReader(new InputStreamReader(in));
  }
  
  /**
   * Determines and returns (if possible) the structure (internally the 
   * header) of the data set as an empty set of instances.
   *
   * @return 			the structure of the data set as an empty set 
   * 				of Instances
   * @throws IOException 	if an error occurs
   */
  public Instances getStructure() throws IOException {
    if (m_sourceReader == null)
      throw new IOException("No source has been specified");

    if (m_structure == null) {
      try {
	m_JSON      = JSONNode.read(m_sourceReader);
	m_structure = new Instances(JSONInstances.toHeader(m_JSON), 0);
      }
      catch (IOException ioe) {
	// just re-throw it
	throw ioe;
      }
      catch (Exception e) {
	throw new RuntimeException(e);
      }
    }

    return new Instances(m_structure, 0);
  }
  
  /**
   * Return the full data set. If the structure hasn't yet been determined
   * by a call to getStructure then method should do so before processing
   * the rest of the data set.
   *
   * @return 			the structure of the data set as an empty 
   * 				set of Instances
   * @throws IOException 	if there is no source or parsing fails
   */
  public Instances getDataSet() throws IOException {
    if (m_sourceReader == null)
      throw new IOException("No source has been specified");
    
    if (getRetrieval() == INCREMENTAL)
      throw new IOException("Cannot mix getting Instances in both incremental and batch modes");

    setRetrieval(BATCH);
    if (m_structure == null)
      getStructure();

    try {
      // close the stream
      m_sourceReader.close();
    } catch (Exception ex) {
    }

    return JSONInstances.toInstances(m_JSON);
  }

  /**
   * JSONLoader is unable to process a data set incrementally.
   *
   * @param structure		ignored
   * @return 			never returns without throwing an exception
   * @throws IOException 	always. JSONLoader is unable to process a 
   * 				data set incrementally.
   */
  public Instance getNextInstance(Instances structure) throws IOException {
    throw new IOException("JSONLoader can't read data sets incrementally.");
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5784 $");
  }

  /**
   * Main method.
   *
   * @param args 	should contain the name of an input file.
   */
  public static void main(String[] args) {
    runFileLoader(new JSONLoader(), args);
  }
}

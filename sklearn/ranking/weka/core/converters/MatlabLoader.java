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
 * MatlabLoader.java
 * Copyright (C) 2009 University of Waikato, Hamilton, NZ
 *
 */

package weka.core.converters;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.RevisionUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.util.Vector;
import java.util.ArrayList;

/**
 <!-- globalinfo-start -->
 * Reads a Matlab file containing a single matrix in ASCII format.
 * <p/>
 <!-- globalinfo-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5987 $
 * @see Loader
 */
public class MatlabLoader 
  extends AbstractFileLoader 
  implements BatchConverter, URLSourcedLoader {

  /** for serialization. */
  private static final long serialVersionUID = -8861142318612875251L;

  /** the file extension. */
  public static String FILE_EXTENSION = ".m";

  /** the url. */
  protected String m_URL = "http://";

  /** The reader for the source file. */
  protected transient Reader m_sourceReader = null;

  /** the buffer of the rows read so far. */
  protected Vector<Vector<Double>> m_Buffer = null;
  
  /**
   * Returns a string describing this Loader.
   * 
   * @return 		a description of the Loader suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return "Reads a Matlab file containing a single matrix in ASCII format.";
  }

  /**
   * Get the file extension used for libsvm files.
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
    return new String[]{getFileExtension()};
  }

  /**
   * Returns a description of the file type.
   *
   * @return 		a short file description
   */
  public String getFileDescription() {
    return "Matlab ASCII files";
  }

  /**
   * Resets the Loader ready to read a new data set.
   * 
   * @throws IOException 	if something goes wrong
   */
  public void reset() throws IOException {
    m_structure = null;
    m_Buffer    = null;
    
    setRetrieval(NONE);
    
    if ((m_File != null) && (new File(m_File)).isFile()) {
      setFile(new File(m_File));
    }
    else if ((m_URL != null) && !m_URL.equals("http://")) {
      setURL(m_URL);
    }
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
    m_Buffer    = null;
    
    setRetrieval(NONE);
    
    setSource(url.openStream());

    m_URL = url.toString();
  }

  /**
   * Set the url to load from.
   *
   * @param url 		the url to load from
   * @throws IOException 		if the url can't be set.
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
    int			numAtt;
    ArrayList<Attribute>		atts;
    int			i;
    String		relName;
    Vector<Double>	row;
    int			c;
    char		chr;
    StringBuffer	str;
    boolean		isComment;
    
    if (m_sourceReader == null)
      throw new IOException("No source has been specified");

    if (m_structure == null) {
      numAtt    = 0;
      m_Buffer  = new Vector<Vector<Double>>();
      row       = new Vector<Double>();
      str       = new StringBuffer();
      isComment = false;
      m_Buffer.add(row);
      try {
	// determine number of attributes
	while ((c = m_sourceReader.read()) != -1) {
	  chr = (char) c;
	  
	  // comment found?
	  if (chr == '%')
	    isComment = true;
	  
	  // end of line reached
	  if ((chr == '\n') || (chr == '\r')) {
	    isComment = false;
	    if (str.length() > 0)
	      row.add(new Double(str.toString()));
	    if (numAtt == 0)
	      numAtt = row.size();
	    if (row.size() > 0) {
	      row = new Vector<Double>();
	      m_Buffer.add(row);
	    }
	    str = new StringBuffer();
	    continue;
	  }

	  // skip till end of comment line
	  if (isComment)
	    continue;
	  
	  // separator found?
	  if ((chr == '\t') || (chr == ' ')) {
	    if (str.length() > 0) {
	      row.add(new Double(str.toString()));
	      str = new StringBuffer();
	    }
	  }
	  else {
	    str.append(chr);
	  }
	}
	
	// last number?
	if (str.length() > 0)
	  row.add(new Double(str.toString()));
	
	// generate header
	atts = new ArrayList<Attribute>(numAtt);
	for (i = 0; i < numAtt; i++)
	  atts.add(new Attribute("att_" + (i+1)));
	
	if (!m_URL.equals("http://"))
	  relName = m_URL;
	else
	  relName = m_File;
	
	m_structure = new Instances(relName, atts, 0);
	m_structure.setClassIndex(m_structure.numAttributes() - 1);
      }
      catch (Exception ex) {
	ex.printStackTrace();
	throw new IOException("Unable to determine structure as Matlab ASCII file: " + ex);
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
    Instances 		result;
    Vector<Double>	row;
    double[]		data;
    int			i;
    int			n;

    if (m_sourceReader == null)
      throw new IOException("No source has been specified");
    
    if (getRetrieval() == INCREMENTAL)
      throw new IOException("Cannot mix getting Instances in both incremental and batch modes");

    setRetrieval(BATCH);
    if (m_structure == null)
      getStructure();

    result = new Instances(m_structure, 0);

    // create instances from buffered data
    for (i = 0; i < m_Buffer.size(); i++) {
      row = m_Buffer.get(i);
      if (row.size() == 0)
	continue;
      data = new double[row.size()];
      for (n = 0; n < row.size(); n++)
	data[n] = row.get(n);
      
      result.add(new DenseInstance(1.0, data));
    }

    // close the stream
    try {
      m_sourceReader.close();
    }
    catch (Exception ex) {
      // ignored
    }
    
    return result;
  }

  /**
   * MatlabLoader is unable to process a data set incrementally.
   *
   * @param structure 		ignored
   * @return 			never returns without throwing an exception
   * @throws IOException 	always. MatlabLoader is unable to process a 
   * 				data set incrementally.
   */
  public Instance getNextInstance(Instances structure) throws IOException {
    throw new IOException("MatlabLoader can't read data sets incrementally.");
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5987 $");
  }

  /**
   * Main method.
   *
   * @param args 	should contain the name of an input file.
   */
  public static void main(String[] args) {
    runFileLoader(new MatlabLoader(), args);
  }
}

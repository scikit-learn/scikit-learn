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
 * LibSVMLoader.java
 * Copyright (C) 2006 University of Waikato, Hamilton, NZ
 *
 */

package weka.core.converters;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionUtils;
import weka.core.SparseInstance;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.ArrayList;

/**
 <!-- globalinfo-start -->
 * Reads a source that is in libsvm format.<br/>
 * <br/>
 * For more information about libsvm see:<br/>
 * <br/>
 * http://www.csie.ntu.edu.tw/~cjlin/libsvm/
 * <p/>
 <!-- globalinfo-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see Loader
 */
public class LibSVMLoader 
  extends AbstractFileLoader 
  implements BatchConverter, URLSourcedLoader {

  /** for serialization. */
  private static final long serialVersionUID = 4988360125354664417L;

  /** the file extension. */
  public static String FILE_EXTENSION = ".libsvm";

  /** the url. */
  protected String m_URL = "http://";

  /** The reader for the source file. */
  protected transient Reader m_sourceReader = null;

  /** the buffer of the rows read so far. */
  protected Vector<double[]> m_Buffer = null;
  
  /**
   * Returns a string describing this Loader.
   * 
   * @return 		a description of the Loader suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Reads a source that is in libsvm format.\n\n"
      + "For more information about libsvm see:\n\n"
      + "http://www.csie.ntu.edu.tw/~cjlin/libsvm/";
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
    return "libsvm data files";
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
   * turns a libsvm row into a double array with the class as the last
   * entry.
   * 
   * @param row		the row to turn into a double array
   * @return		the corresponding double array
   */
  protected double[] libsvmToArray(String row) {
    double[]		result;
    StringTokenizer	tok;
    int			index;
    int			max;
    String		col;
    double		value;

    // determine max index
    max = 0;
    tok = new StringTokenizer(row, " \t");
    tok.nextToken();  // skip class
    while (tok.hasMoreTokens()) {
      col   = tok.nextToken();
      index = Integer.parseInt(col.substring(0, col.indexOf(":")));
      if (index > max)
	max = index;
    }

    // read values into array
    tok    = new StringTokenizer(row, " \t");
    result = new double[max + 1];
    
    // 1. class
    result[result.length - 1] = Double.parseDouble(tok.nextToken());
    
    // 2. attributes
    while (tok.hasMoreTokens()) {
      col   = tok.nextToken();
      index = Integer.parseInt(col.substring(0, col.indexOf(":")));
      value = Double.parseDouble(col.substring(col.indexOf(":") + 1));
      result[index - 1] = value;
    }
    
    return result;
  }
  
  /**
   * determines the number of attributes, if the number of attributes in the
   * given row is greater than the current amount then this number will be
   * returned, otherwise the current number.
   * 
   * @param row		row to determine the number of attributes from
   * @param num		the current number of attributes
   * @return 		the new number of attributes
   */
  protected int determineNumAttributes(String row, int num) {
    int		result;
    int		count;
    
    result = num;
    
    count = libsvmToArray(row).length;
    if (count > result)
      result = count;
    
    return result;
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
    StringBuffer	line;
    int			cInt;
    char		c;
    int			numAtt;
    ArrayList<Attribute>		atts;
    int			i;
    String		relName;
    
    if (m_sourceReader == null)
      throw new IOException("No source has been specified");

    if (m_structure == null) {
      m_Buffer = new Vector<double[]>();
      try {
	// determine number of attributes
	numAtt = 0;
	line   = new StringBuffer();
	while ((cInt = m_sourceReader.read()) != -1) {
	  c = (char) cInt;
	  if ((c == '\n') || (c == '\r')) {
	    if (line.length() > 0) {
	      m_Buffer.add(libsvmToArray(line.toString()));
	      numAtt = determineNumAttributes(line.toString(), numAtt);
	    }
	    line = new StringBuffer();
	  }
	  else {
	    line.append(c);
	  }
	}
	
	// last line?
	if (line.length() != 0) {
	  m_Buffer.add(libsvmToArray(line.toString()));
	  numAtt = determineNumAttributes(line.toString(), numAtt);
	}
	
	// generate header
	atts = new ArrayList<Attribute>(numAtt);
	for (i = 0; i < numAtt - 1; i++)
	  atts.add(new Attribute("att_" + (i+1)));
	atts.add(new Attribute("class"));
	
	if (!m_URL.equals("http://"))
	  relName = m_URL;
	else
	  relName = m_File;
	
	m_structure = new Instances(relName, atts, 0);
	m_structure.setClassIndex(m_structure.numAttributes() - 1);
      }
      catch (Exception ex) {
	ex.printStackTrace();
	throw new IOException("Unable to determine structure as libsvm: " + ex);
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
    Instances 	result;
    double[]	sparse;
    double[]	data;
    int		i;

    if (m_sourceReader == null)
      throw new IOException("No source has been specified");
    
    if (getRetrieval() == INCREMENTAL)
      throw new IOException("Cannot mix getting Instances in both incremental and batch modes");

    setRetrieval(BATCH);
    if (m_structure == null)
      getStructure();

    result = new Instances(m_structure, 0);

    // create instances from buffered arrays
    for (i = 0; i < m_Buffer.size(); i++) {
      sparse = (double[]) m_Buffer.get(i);
      
      if (sparse.length != m_structure.numAttributes()) {
	data = new double[m_structure.numAttributes()];
	// attributes
	System.arraycopy(sparse, 0, data, 0, sparse.length - 1);
	// class
	data[data.length - 1] = sparse[sparse.length - 1];
      }
      else {
	data = sparse;
      }
      
      result.add(new SparseInstance(1, data));
    }

    try {
      // close the stream
      m_sourceReader.close();
    } catch (Exception ex) {

    }
    
    return result;
  }

  /**
   * LibSVmLoader is unable to process a data set incrementally.
   *
   * @param structure 		ignored
   * @return 			never returns without throwing an exception
   * @throws IOException 	always. LibSVMLoader is unable to process a 
   * 				data set incrementally.
   */
  public Instance getNextInstance(Instances structure) throws IOException {
    throw new IOException("LibSVMLoader can't read data sets incrementally.");
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
   * @param args 	should contain the name of an input file.
   */
  public static void main(String[] args) {
    runFileLoader(new LibSVMLoader(), args);
  }
}

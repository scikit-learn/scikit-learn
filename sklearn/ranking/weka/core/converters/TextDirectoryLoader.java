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
 * TextDirectoryLoader.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.converters;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Vector;
import java.util.ArrayList;

/**
 <!-- globalinfo-start -->
 * Loads all text files in a directory and uses the subdirectory names as class labels. The content of the text files will be stored in a String attribute, the filename can be stored as well.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -D
 *  Enables debug output.
 *  (default: off)</pre>
 * 
 * <pre> -F
 *  Stores the filename in an additional attribute.
 *  (default: off)</pre>
 * 
 * <pre> -dir &lt;directory&gt;
 *  The directory to work on.
 *  (default: current directory)</pre>
 * 
 <!-- options-end -->
 *
 * Based on code from the TextDirectoryToArff tool:
 * <ul>
 *    <li><a href="https://list.scms.waikato.ac.nz/mailman/htdig/wekalist/2002-October/000685.html" target="_blank">Original tool</a></li>
 *    <li><a href="https://list.scms.waikato.ac.nz/mailman/htdig/wekalist/2004-January/002160.html" target="_blank">Current version</a></li>
 *    <li><a href="http://weka.wikispaces.com/ARFF+files+from+Text+Collections" target="_blank">Wiki article</a></li>
 * </ul>
 *
 * @author Ashraf M. Kibriya (amk14 at cs.waikato.ac.nz)
 * @author Richard Kirkby (rkirkby at cs.waikato.ac.nz)
 * @author fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6766 $
 * @see Loader
 */
public class TextDirectoryLoader
  extends AbstractLoader
  implements BatchConverter, OptionHandler {
  
  /** for serialization */
  private static final long serialVersionUID = 2592118773712247647L;
  
  /** Holds the determined structure (header) of the data set. */
  protected Instances m_structure = null;
  
  /** Holds the source of the data set. */
  protected File m_sourceFile = new File(System.getProperty("user.dir"));
  
  /** whether to print some debug information */
  protected boolean m_Debug = false;
  
  /** whether to include the filename as an extra attribute */
  protected boolean m_OutputFilename = false;
  
  /**
   * default constructor
   */
  public TextDirectoryLoader() {
    // No instances retrieved yet
    setRetrieval(NONE);
  }
  
  /**
   * Returns a string describing this loader
   * 
   * @return 		a description of the evaluator suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return 
        "Loads all text files in a directory and uses the subdirectory names "
      + "as class labels. The content of the text files will be stored in a "
      + "String attribute, the filename can be stored as well.";
  }
  
  /** 
   * Lists the available options
   * 
   * @return 		an enumeration of the available options
   */  
  public Enumeration listOptions() {
    
    Vector<Option> result = new Vector<Option>();
    
    result.add(new Option(
	"\tEnables debug output.\n"
	+ "\t(default: off)",
	"D", 0, "-D"));
    
    result.add(new Option(
	"\tStores the filename in an additional attribute.\n"
	+ "\t(default: off)",
	"F", 0, "-F"));
    
    result.add(new Option(
	"\tThe directory to work on.\n"
	+ "\t(default: current directory)",
	"dir", 0, "-dir <directory>"));
    
    return  result.elements();
  }
  
  /** 
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -D
   *  Enables debug output.
   *  (default: off)</pre>
   * 
   * <pre> -F
   *  Stores the filename in an additional attribute.
   *  (default: off)</pre>
   * 
   * <pre> -dir &lt;directory&gt;
   *  The directory to work on.
   *  (default: current directory)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the options
   * @throws Exception if options cannot be set
   */  
  public void setOptions(String[] options) throws Exception {
    setDebug(Utils.getFlag("D", options));
    
    setOutputFilename(Utils.getFlag("F", options));
    
    setDirectory(new File(Utils.getOption("dir", options)));
  }
  
  /** 
   * Gets the setting
   * 
   * @return the current setting
   */  
  public String[] getOptions() {
    Vector<String> options = new Vector<String>();
    
    if (getDebug())
      options.add("-D");
    
    if (getOutputFilename())
      options.add("-F");

    options.add("-dir");
    options.add(getDirectory().getAbsolutePath());
    
    return (String[]) options.toArray(new String[options.size()]);
  }
  
  /**
   * Sets whether to print some debug information.
   * 
   * @param value	if true additional debug information will be printed.
   */
  public void setDebug(boolean value) {
    m_Debug = value;
  }
  
  /**
   * Gets whether additional debug information is printed.
   * 
   * @return		true if additional debug information is printed
   */
  public boolean getDebug() {
    return m_Debug;
  }
  
  /**
   * the tip text for this property
   * 
   * @return 		the tip text
   */
  public String debugTipText(){
    return "Whether to print additional debug information to the console.";
  }
  
  /**
   * Sets whether the filename will be stored as an extra attribute.
   * 
   * @param value	if true the filename will be stored in an extra
   * 			attribute
   */
  public void setOutputFilename(boolean value) {
    m_OutputFilename = value;
    reset();
  }
  
  /**
   * Gets whether the filename will be stored as an extra attribute.
   * 
   * @return		true if the filename is stored in an extra attribute
   */
  public boolean getOutputFilename() {
    return m_OutputFilename;
  }
  
  /**
   * the tip text for this property
   * 
   * @return 		the tip text
   */
  public String outputFilenameTipText(){
    return "Whether to store the filename in an additional attribute.";
  }
  
  /**
   * Returns a description of the file type, actually it's directories.
   *
   * @return 		a short file description
   */
  public String getFileDescription() {
    return "Directories";
  }
  
  /**
   * get the Dir specified as the source
   *
   * @return 		the source directory
   */
  public File getDirectory() {
    return new File(m_sourceFile.getAbsolutePath());
  }
  
  /**
   * sets the source directory
   *
   * @param 			dir the source directory
   * @throws IOException 	if an error occurs
   */
  public void setDirectory(File dir) throws IOException {
    setSource(dir);
  }
  
  /**
   * Resets the loader ready to read a new data set
   */
  public void reset() {
    m_structure = null;
    setRetrieval(NONE);
  }
  
  /**
   * Resets the Loader object and sets the source of the data set to be 
   * the supplied File object.
   *
   * @param dir 		the source directory.
   * @throws IOException 	if an error occurs
   */
  public void setSource(File dir) throws IOException {
    reset();
    
    if (dir == null) {
      throw new IOException("Source directory object is null!");
    }
    
    m_sourceFile = dir;
    if (!dir.exists() || !dir.isDirectory())
      throw new IOException("Directory '" + dir + "' not found");
  }
  
  /**
   * Determines and returns (if possible) the structure (internally the 
   * header) of the data set as an empty set of instances.
   *
   * @return 			the structure of the data set as an empty 
   * 				set of Instances
   * @throws IOException 	if an error occurs
   */
  public Instances getStructure() throws IOException {
    if (getDirectory() == null) {
      throw new IOException("No directory/source has been specified");
    }
    
    // determine class labels, i.e., sub-dirs
    if (m_structure == null) {
      String directoryPath = getDirectory().getAbsolutePath();
      ArrayList<Attribute> atts = new ArrayList<Attribute>();
      ArrayList<String> classes = new ArrayList<String>();
      
      File dir = new File(directoryPath);
      String[] subdirs = dir.list();
      
      for (int i = 0; i < subdirs.length; i++) {
	File subdir = new File(directoryPath + File.separator + subdirs[i]);
	if (subdir.isDirectory())
	  classes.add(subdirs[i]);
      }
      
      atts.add(new Attribute("text", (ArrayList<String>) null));
      if (m_OutputFilename)
	atts.add(new Attribute("filename", (ArrayList<String>) null));
      // make sure that the name of the class attribute is unlikely to 
      // clash with any attribute created via the StringToWordVector filter
      atts.add(new Attribute("@@class@@", classes));
      
      String relName = directoryPath.replaceAll("/", "_");
      relName = relName.replaceAll("\\\\", "_").replaceAll(":", "_");
      m_structure = new Instances(relName, atts, 0);    
      m_structure.setClassIndex(m_structure.numAttributes() - 1);
    }
    
    return m_structure;
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
    if (getDirectory() == null)
      throw new IOException("No directory/source has been specified");
    
    String directoryPath = getDirectory().getAbsolutePath();
    ArrayList<String> classes = new ArrayList<String>();
    Enumeration enm = getStructure().classAttribute().enumerateValues();
    while (enm.hasMoreElements())
      classes.add((String)enm.nextElement());
    
    Instances data = getStructure();
    int fileCount = 0;
    for (int k = 0; k < classes.size(); k++) {
      String subdirPath = (String) classes.get(k);
      File subdir = new File(directoryPath + File.separator + subdirPath);
      String[] files = subdir.list();
      for (int j = 0; j < files.length; j++) {
	try {
	  fileCount++;
	  if (getDebug())
	    System.err.println(
		"processing " + fileCount + " : " + subdirPath + " : " + files[j]); 
	  
	  double[] newInst = null;
	  if (m_OutputFilename)
	    newInst = new double[3];
	  else
	    newInst = new double[2];		    
	  File txt = new File(directoryPath + File.separator + subdirPath + File.separator + files[j]);
	  BufferedInputStream is;
	  is = new BufferedInputStream(new FileInputStream(txt));
	  StringBuffer txtStr = new StringBuffer();
	  int c;
	  while ((c = is.read()) != -1) {
	    txtStr.append((char) c);
	  }
	  
	  newInst[0] = (double) data.attribute(0).addStringValue(txtStr.toString());
	  if (m_OutputFilename)
	    newInst[1] = (double) data.attribute(1).addStringValue(subdirPath + File.separator + files[j]);
	  newInst[data.classIndex()] = (double) k;
	  data.add(new DenseInstance(1.0, newInst));
          is.close();
	}
	catch (Exception e) {
	  System.err.println("failed to convert file: " + directoryPath + File.separator + subdirPath + File.separator + files[j]);
	}
      }
    }
    
    return data;
  }
  
  /**
   * TextDirectoryLoader is unable to process a data set incrementally.
   *
   * @param structure ignored
   * @return never returns without throwing an exception
   * @throws IOException always. TextDirectoryLoader is unable to process a data
   * set incrementally.
   */
  public Instance getNextInstance(Instances structure) throws IOException {
    throw new IOException("TextDirectoryLoader can't read data sets incrementally.");
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6766 $");
  }
  
  /**
   * Main method.
   *
   * @param args should contain the name of an input file.
   */
  public static void main(String[] args) {
    if (args.length > 0) {
      try {
	TextDirectoryLoader loader = new TextDirectoryLoader();
	loader.setOptions(args);
	System.out.println(loader.getDataSet());
      } 
      catch (Exception e) {
	e.printStackTrace();
      }
    } 
    else {
      System.err.println(
	  "\nUsage:\n" 
	  + "\tTextDirectoryLoader [options]\n"
	  + "\n"
	  + "Options:\n");

      Enumeration enm = ((OptionHandler) new TextDirectoryLoader()).listOptions();
      while (enm.hasMoreElements()) {
	Option option = (Option) enm.nextElement();
	System.err.println(option.synopsis());
	System.err.println(option.description());
      }
      
      System.err.println();
    }
  }
}

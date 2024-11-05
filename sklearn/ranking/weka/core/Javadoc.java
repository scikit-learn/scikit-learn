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
 * Javadoc.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Enumeration;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * Abstract superclass for classes that generate Javadoc comments and replace
 * the content between certain comment tags.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public abstract class Javadoc 
  implements OptionHandler, RevisionHandler {

  /** the start tag */
  protected String[] m_StartTag = null;

  /** the end tag */
  protected String[] m_EndTag = null;
  
  /** the classname */
  protected String m_Classname = Javadoc.class.getName();
  
  /** whether to include the stars in the Javadoc */
  protected boolean m_UseStars = true;

  /** the directory above the class to update */
  protected String m_Dir = "";
  
  /** whether to suppress error messages (no printout in the console) */
  protected boolean m_Silent = false;
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option> result = new Vector<Option>();

    result.addElement(new Option(
        "\tThe class to load.",
        "W", 1, "-W <classname>"));
    
    result.addElement(new Option(
        "\tSuppresses the '*' in the Javadoc.",
        "nostars", 0, "-nostars"));
    
    result.addElement(new Option(
        "\tThe directory above the package hierarchy of the class.",
        "dir", 1, "-dir <dir>"));
    
    result.addElement(new Option(
        "\tSuppresses printing in the console.",
        "silent", 0, "-silent"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. 
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    String      		tmpStr;
    
    tmpStr = Utils.getOption('W', options);
    if (tmpStr.length() > 0)
      setClassname(tmpStr);
    else
      setClassname(this.getClass().getName());

    setUseStars(!Utils.getFlag("nostars", options));

    setDir(Utils.getOption("dir", options));

    setSilent(Utils.getFlag("silent", options));
  }
  
  /**
   * Gets the current settings of this object.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String> 	result;

    result = new Vector<String>();
    
    result.add("-W");
    result.add(getClassname());
    
    if (!getUseStars())
      result.add("-nostars");
    
    if (getDir().length() != 0) {
      result.add("-dir");
      result.add(getDir());
    }
    
    if (getSilent())
      result.add("-silent");
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * sets the classname of the class to generate the Javadoc for
   * 
   * @param value	the new classname
   */
  public void setClassname(String value) {
    m_Classname = value;
  }
  
  /**
   * returns the current classname
   * 
   * @return	the current classname
   */
  public String getClassname() {
    return m_Classname;
  }
  
  /**
   * sets whether to prefix the Javadoc with "*"
   * 
   * @param value	true if stars are used
   */
  public void setUseStars(boolean value) {
    m_UseStars = value;
  }
  
  /**
   * whether the Javadoc is prefixed with "*"
   * 
   * @return 		whether stars are used
   */
  public boolean getUseStars() {
    return m_UseStars;
  }
  
  /**
   * sets the dir containing the file that is to be updated. It is the dir
   * above the package hierarchy of the class.
   * 
   * @param value	the directory containing the classes
   */
  public void setDir(String value) {
    m_Dir = value;
  }
  
  /**
   * returns the current dir containing the class to update. It is the dir
   * above the package name of the class.
   * 
   * @return		the  current directory
   */
  public String getDir() {
    return m_Dir;
  }
  
  /**
   * sets whether to suppress output in the console
   * 
   * @param value	true if output is to be suppressed
   */
  public void setSilent(boolean value) {
    m_Silent = value;
  }
  
  /**
   * whether output in the console is suppressed
   * 
   * @return 		true if output is suppressed
   */
  public boolean getSilent() {
    return m_Silent;
  }
  
  /**
   * prints the given object to System.err
   * 
   * @param o		the object to print
   */
  protected void println(Object o) {
    if (!getSilent())
      System.err.println(o.toString());
  }

  /**
   * returns true if the class can be instantiated, i.e., has a default
   * constructor.
   * 
   * @return true if the class can be instantiated
   */
  protected boolean canInstantiateClass() {
    boolean	result;
    Class	cls;

    result = true;
    cls    = null;

    try {
      cls = Class.forName(getClassname());
    }
    catch (Exception e) {
      result = false;
      println("Cannot instantiate '" + getClassname() + "'! Class in CLASSPATH?");
    }

    if (result) {
      try {
	cls.newInstance();
      }
      catch (Exception e) {
	result = false;
	println("Cannot instantiate '" + getClassname() + "'! Missing default constructor?");
      }
    }
    
    return result;
  }
  
  /**
   * Returns a new instance of the class
   * 
   * @return a new instance of the class
   */
  protected Object getInstance() {
    Object	result;
    Class	cls;

    result = null;
    
    try {
      cls    = Class.forName(getClassname());
      result = cls.newInstance();
    }
    catch (Exception e) {
      result = null;
    }
    
    return result;
  }
  
  /**
   * converts the given String into HTML, i.e., replacing some char entities
   * with HTML entities.
   * 
   * @param s		the string to convert
   * @return		the HTML conform string
   */
  protected String toHTML(String s) {
    String	result;
    
    result = s;
    
    result = result.replaceAll("&", "&amp;");
    result = result.replaceAll("<", "&lt;");
    result = result.replaceAll(">", "&gt;");
    result = result.replaceAll("@", "&#64;");
    result = result.replaceAll("\n", "<br/>\n");
    
    return result;
  }

  /**
   * indents the given string by a given number of indention strings
   * 
   * @param content	the string to indent
   * @param count	the number of times to indent one line
   * @param indentStr	the indention string
   * @return		the indented content
   */
  protected String indent(String content, int count, String indentStr) {
    String		result;
    StringTokenizer	tok;
    int			i;
    
    tok = new StringTokenizer(content, "\n", true);
    result = "";
    while (tok.hasMoreTokens()) {
      if (result.endsWith("\n") || (result.length() == 0)) {
	for (i = 0; i < count; i++)
	  result += indentStr;
      }
      result += tok.nextToken();
    }
    
    return result;
  }
  
  /**
   * generates and returns the Javadoc for the specified start/end tag pair.
   * 
   * @param index	the index in the start/end tag array
   * @return		the generated Javadoc
   * @throws Exception 	in case the generation fails
   */
  protected abstract String generateJavadoc(int index) throws Exception;
  
  /**
   * generates and returns the Javadoc
   * 
   * @return		the generated Javadoc
   * @throws Exception 	in case the generation fails
   */
  protected String generateJavadoc() throws Exception {
    String	result;
    int		i;
    
    result = "";
    
    for (i = 0; i < m_StartTag.length; i++) {
      if (i > 0)
	result += "\n\n";
      result += generateJavadoc(i).trim();
    }
    
    return result;
  }

  /**
   * determines the base string of the given indention string, whether it's
   * either only spaces (one space will be retured) or mixed mode (tabs and 
   * spaces, in that case the same string will be returned)
   * 
   * @param str		the string to analyze
   * @return 		the indention string
   */
  protected String getIndentionString(String str) {
    String	result;
    
    // only spaces?
    if (str.replaceAll(" ", "").length() == 0)
      result = " ";
    // only tabs?
    else if (str.replaceAll("\t", "").length() == 0)
      result = "\t";
    else
      result = str;
      
    return result;
  }
  
  /**
   * determines the number of indention strings that have to be inserted to
   * generated the given indention string.
   * 
   * @param str 	the string to analyze
   * @return		the number of base indention strings to insert
   */
  protected int getIndentionLength(String str) {
    int		result;
    
    // only spaces?
    if (str.replaceAll(" ", "").length() == 0)
      result = str.length();
    // only tabs?
    else if (str.replaceAll("\t", "").length() == 0)
      result = str.length();
    else
      result = 1;
    
    return result;
  }
  
  /**
   * generates and returns the Javadoc for the specified start/end tag pair
   * 
   * @param content	the current source code
   * @param index	the index in the start/end tag array
   * @return		the generated Javadoc
   * @throws Exception 	in case the generation fails
   */
  protected String updateJavadoc(String content, int index) throws Exception {
    StringBuffer	resultBuf;
    int			indentionLen;
    String		indentionStr;
    String		part;
    String		tmpStr;

    // start and end tag?
    if (    (content.indexOf(m_StartTag[index]) == -1)
	   || (content.indexOf(m_EndTag[index]) == -1) ) {
	println(
	    "No start and/or end tags found: " 
	    + m_StartTag[index] + "/" + m_EndTag[index]);
	return content;
    }

    // replace option-tags
    resultBuf = new StringBuffer();
    while (content.length() > 0) {
      if (content.indexOf(m_StartTag[index]) > -1) {
	part = content.substring(0, content.indexOf(m_StartTag[index]));
	// is it a Java constant? -> skip
	if (part.endsWith("\"")) {
	  resultBuf.append(part);
	  resultBuf.append(m_StartTag[index]);
	  content = content.substring(part.length() + m_StartTag[index].length());
	}
	else {
	  tmpStr       = part.substring(part.lastIndexOf("\n") + 1);
	  indentionLen = getIndentionLength(tmpStr);
	  indentionStr = getIndentionString(tmpStr);
	  part         = part.substring(0, part.lastIndexOf("\n") + 1);
	  resultBuf.append(part);
	  resultBuf.append(indent(m_StartTag[index], indentionLen, indentionStr) + "\n");
	  resultBuf.append(indent(generateJavadoc(index), indentionLen, indentionStr));
	  resultBuf.append(indent(m_EndTag[index], indentionLen, indentionStr));
	  content = content.substring(content.indexOf(m_EndTag[index]));
	  content = content.substring(m_EndTag[index].length());
	}
      }
      else {
	resultBuf.append(content);
	content = "";
      }
    }
    
    return resultBuf.toString().trim();
  }
  
  /**
   * updates the Javadoc in the given source code.
   * 
   * @param content	the source code
   * @return		the updated source code
   * @throws Exception 	in case the generation fails
   */
  protected String updateJavadoc(String content) throws Exception {
    String	result;
    int		i;
    
    result = content;
    
    for (i = 0; i < m_StartTag.length; i++) {
      result = updateJavadoc(result, i);
    }
    
    return result;
  }
  
  /**
   * generates the Javadoc and returns it applied to the source file if one
   * was provided, otherwise an empty string.
   * 
   * @return		the generated Javadoc
   * @throws Exception 	in case the generation fails
   */
  public String updateJavadoc() throws Exception {
    StringBuffer	contentBuf;
    BufferedReader	reader;
    String		line;
    String		result;
    File		file;
    
    result = "";
    
    // non-existing?
    file = new File(getDir() + "/" + getClassname().replaceAll("\\.", "/") + ".java");
    if (!file.exists()) {
      println("File '" + file.getAbsolutePath() + "' doesn't exist!");
      return result;
    }
    
    try {
      // load file
      reader     = new BufferedReader(new FileReader(file));
      contentBuf = new StringBuffer();
      while ((line = reader.readLine()) != null) {
	contentBuf.append(line + "\n");
      }
      reader.close();
      result = updateJavadoc(contentBuf.toString());
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    
    return result.trim();
  }
  
  /**
   * generates either the plain Javadoc (if no filename specified) or the
   * updated file (if a filename is specified). The start and end tag for
   * the global info have to be specified in the file in the latter case.
   * 
   * @return 		either the plain Javadoc or the modified file
   * @throws Exception 	in case the generation fails
   */
  public String generate() throws Exception {
    if (getDir().length() == 0)
      return generateJavadoc();
    else
      return updateJavadoc();
  }
  
  /**
   * generates a string to print as help on the console
   * 
   * @return 	the generated help
   */
  public String generateHelp() {
    String 		result;
    Enumeration 	enm;
    Option 		option;
    
    result = getClass().getName().replaceAll(".*\\.", "") + " Options:\n\n";
    enm = listOptions();
    while (enm.hasMoreElements()) {
      option = (Option) enm.nextElement();
      result += option.synopsis() + "\n" + option.description() + "\n";
    }
    
    return result;
  }
  
  /**
   * runs the javadoc producer with the given commandline options
   * 
   * @param javadoc	the javadoc producer to execute
   * @param options	the commandline options
   */
  protected static void runJavadoc(Javadoc javadoc, String[] options) {
    try {
      try {
	if (Utils.getFlag('h', options))
	  throw new Exception("Help requested");

	javadoc.setOptions(options);
        Utils.checkForRemainingOptions(options);

        // directory is necessary!
        if (javadoc.getDir().length() == 0)
          throw new Exception("No directory provided!");
      } 
      catch (Exception ex) {
        String result = "\n" + ex.getMessage() + "\n\n" + javadoc.generateHelp();
        throw new Exception(result);
      }

      System.out.println(javadoc.generate() + "\n");
    } 
    catch (Exception ex) {
      System.err.println(ex.getMessage());
    }
  }
}

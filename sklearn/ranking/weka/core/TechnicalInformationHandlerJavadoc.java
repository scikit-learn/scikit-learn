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
 * TechnicalInformationHandlerJavadoc.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.util.Enumeration;
import java.util.Vector;

/**
 * Generates Javadoc comments from the TechnicalInformationHandler's data.
 * Update the BibTex references and the plaintext techincal information. 
 * <p/>
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -W &lt;classname&gt;
 *  The class to load.</pre>
 * 
 * <pre> -nostars
 *  Suppresses the '*' in the Javadoc.</pre>
 * 
 * <pre> -dir &lt;dir&gt;
 *  The directory above the package hierarchy of the class.</pre>
 * 
 * <pre> -silent
 *  Suppresses printing in the console.</pre>
 * 
 * <pre> -noprolog
 *  Suppresses the 'BibTex:' prolog in the Javadoc.</pre>
 * 
 <!-- options-end -->
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see #PLAINTEXT_STARTTAG
 * @see #PLAINTEXT_ENDTAG
 * @see #BIBTEX_STARTTAG
 * @see #BIBTEX_ENDTAG
 */
public class TechnicalInformationHandlerJavadoc 
  extends Javadoc {
  
  /** the start comment tag for inserting the generated BibTex */
  public final static String PLAINTEXT_STARTTAG = "<!-- technical-plaintext-start -->";
  
  /** the end comment tag for inserting the generated BibTex */
  public final static String PLAINTEXT_ENDTAG = "<!-- technical-plaintext-end -->";
  
  /** the start comment tag for inserting the generated BibTex */
  public final static String BIBTEX_STARTTAG = "<!-- technical-bibtex-start -->";
  
  /** the end comment tag for inserting the generated BibTex */
  public final static String BIBTEX_ENDTAG = "<!-- technical-bibtex-end -->";
  
  /** whether to include the "Valid options..." prolog in the Javadoc */
  protected boolean m_Prolog = true;
  
  /**
   * default constructor 
   */
  public TechnicalInformationHandlerJavadoc() {
    super();
    
    m_StartTag    = new String[2];
    m_EndTag      = new String[2];
    m_StartTag[0] = PLAINTEXT_STARTTAG;
    m_EndTag[0]   = PLAINTEXT_ENDTAG;
    m_StartTag[1] = BIBTEX_STARTTAG;
    m_EndTag[1]   = BIBTEX_ENDTAG;
  }
  
  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {
    Vector<Option>        result;
    Enumeration   en;
    
    result = new Vector<Option>();
    
    en = super.listOptions();
    while (en.hasMoreElements())
      result.addElement((Option)en.nextElement());

    result.addElement(new Option(
        "\tSuppresses the 'BibTex:' prolog in the Javadoc.",
        "noprolog", 0, "-noprolog"));
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. 
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
    super.setOptions(options);

    setProlog(!Utils.getFlag("noprolog", options));
  }
  
  /**
   * Gets the current settings of this object.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    Vector<String>        result;
    String[]      options;
    int           i;
    
    result  = new Vector<String>();
    
    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);

    if (!getProlog())
      result.add("-noprolog");
    
    return (String[]) result.toArray(new String[result.size()]);
  }
  
  /**
   * sets whether to add the "Valid options are..." prolog
   * 
   * @param value	true if the prolog is to be used
   */
  public void setProlog(boolean value) {
    m_Prolog = value;
  }
  
  /**
   * whether "Valid options are..." prolog is included in the Javadoc
   * 
   * @return		whether the prolog is currently used
   */
  public boolean getProlog() {
    return m_Prolog;
  }
  
  /**
   * generates and returns the Javadoc for the specified start/end tag pair.
   * 
   * @param index	the index in the start/end tag array
   * @return		the generated Javadoc
   * @throws Exception 	in case the generation fails
   */
  protected String generateJavadoc(int index) throws Exception {
    String			result;
    TechnicalInformationHandler	handler;
    
    result = "";
    
    if (!canInstantiateClass())
      return result;
    
    if (!ClassDiscovery.hasInterface(TechnicalInformationHandler.class, getInstance().getClass()))
	throw new Exception("Class '" + getClassname() + "' is not a TechnicalInformationHandler!");

    handler = (TechnicalInformationHandler) getInstance();
    
    switch (index) {
      case 0:  // plain text
	result = toHTML(handler.getTechnicalInformation().toString()) + "\n";
	break;
	
      case 1:  // bibtex
	// prolog?
	if (getProlog())
	  result = "BibTeX:\n";
	result += "<pre>\n";
	result += toHTML(handler.getTechnicalInformation().toBibTex()).replaceAll("<br/>", "") + "\n";
	result += "</pre>\n<p/>\n";
	break;
    }
    
    // stars?
    if (getUseStars()) 
      result = indent(result, 1, "* ");
    
    return result;
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
   * Parses the given commandline parameters and generates the Javadoc.
   * 
   * @param args	the commandline parameters for the object
   */
  public static void main(String[] args) {
    runJavadoc(new TechnicalInformationHandlerJavadoc(), args);
  }
}

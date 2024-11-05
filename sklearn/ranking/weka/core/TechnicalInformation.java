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
 * TechnicalInformation.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

/**
 * Used for paper references in the Javadoc and for BibTex generation.
 * Based on documentation found here: <p/>
 * <a href="http://www.ecst.csuchico.edu/~jacobsd/bib/formats/bibtex.html" target="_blank">http://www.ecst.csuchico.edu/~jacobsd/bib/formats/bibtex.html</a>
 * <p/>
 * BibTex examples can be found here: <p/>
 * <a href="http://bib2web.djvuzone.org/bibtex.html" target="_blank">http://bib2web.djvuzone.org/bibtex.html</a>
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 * @see TechnicalInformationHandler
 */
public class TechnicalInformation
  implements RevisionHandler {

  /** the different types of information */
  public enum Type {
    /** An article from a journal or magazine. */
    ARTICLE("article", "An article from a journal or magazine."),  
    /** A book with an explicit publisher. */
    BOOK("book", "A book with an explicit publisher."),
    /** A work that is printed and bound, but without a named publisher or sponsoring institution. */
    BOOKLET("booklet", "A work that is printed and bound, but without a named publisher or sponsoring institution."),
    /** The same as inproceedings. */
    CONFERENCE("conference", "The same as inproceedings."),
    /** A part of a book, which may be a chapter (or section or whatever) and/or a range of pages. */
    INBOOK("inbook", "A part of a book, which may be a chapter (or section or whatever) and/or a range of pages."),
    /** A part of a book having its own title. */
    INCOLLECTION("incollection", "A part of a book having its own title."),
    /** An article in a conference proceedings. */
    INPROCEEDINGS("inproceedings", "An article in a conference proceedings."),
    /** Technical documentation. */
    MANUAL("manual", "Technical documentation."),
    /** A Master's thesis. */
    MASTERSTHESIS("mastersthesis", "A Master's thesis."),
    /** Use this type when nothing else fits. */
    MISC("misc", "Use this type when nothing else fits."),
    /** A PhD thesis. */
    PHDTHESIS("phdthesis", "A PhD thesis."),
    /** The proceedings of a conference. */
    PROCEEDINGS("proceedings", "The proceedings of a conference."),
    /** A report published by a school or other institution, usually numbered within a series. */
    TECHREPORT("techreport", "A report published by a school or other institution, usually numbered within a series."),
    /** A document having an author and title, but not formally published. */
    UNPUBLISHED("unpublished", "A document having an author and title, but not formally published.");
    
    /** the string used in toString()  */
    protected String m_Display;
    
    /** the comment about this type */
    protected String m_Comment;
    
    /**
     * the constructor
     * 
     * @param display	the string to return in toString()
     * @param comment	a comment about the type
     */
    private Type(String display, String comment) {
      m_Display = display;
      m_Comment = comment;
    }
    
    /**
     * returns the display string
     * 
     * @return 		the display string
     */
    public String getDisplay() {
      return m_Display;
    }
    
    /**
     * returns the comment string
     * 
     * @return		the comment string
     */
    public String getComment() {
      return m_Comment;
    }
    
    /**
     * returns the display string of the Type
     * 
     * @return		the display string
     */
    public String toString() {
      return m_Display;
    }
  }
  
  /** the possible fields */
  public enum Field {
    /** Usually the address of the publisher or other type of institution. For major publishing houses, van Leunen recommends omitting the information entirely. For small publishers, on the other hand, you can help the reader by giving the complete address. */
    ADDRESS("address", "Usually the address of the publisher or other type of institution. For major publishing houses, van Leunen recommends omitting the information entirely. For small publishers, on the other hand, you can help the reader by giving the complete address."),
    /** An annotation. It is not used by the standard bibliography styles, but may be used by others that produce an annotated bibliography. */
    ANNOTE("annote", "An annotation. It is not used by the standard bibliography styles, but may be used by others that produce an annotated bibliography."),
    /** The name(s) of the author(s), in the format described in the LaTeX book. */
    AUTHOR("author", "The name(s) of the author(s), in the format described in the LaTeX book."),
    /** Title of a book, part of which is being cited. See the LaTeX book for how to type titles. For book entries, use the title field instead. */
    BOOKTITLE("booktitle", "Title of a book, part of which is being cited. See the LaTeX book for how to type titles. For book entries, use the title field instead."),
    /** A chapter (or section or whatever) number. */
    CHAPTER("chapter", "A chapter (or section or whatever) number."),
    /** The database key of the entry being cross referenced. Any fields that are missing from the current record are inherited from the field being cross referenced. */
    CROSSREF("crossref", "The database key of the entry being cross referenced. Any fields that are missing from the current record are inherited from the field being cross referenced."),
    /** The edition of a book---for example, ``Second''. This should be an ordinal, and should have the first letter capitalized, as shown here; the standard styles convert to lower case when necessary. */
    EDITION("edition", "The edition of a book---for example, ``Second''. This should be an ordinal, and should have the first letter capitalized, as shown here; the standard styles convert to lower case when necessary."),
    /** Name(s) of editor(s), typed as indicated in the LaTeX book. If there is also an author field, then the editor field gives the editor of the book or collection in which the reference appears. */
    EDITOR("editor", "Name(s) of editor(s), typed as indicated in the LaTeX book. If there is also an author field, then the editor field gives the editor of the book or collection in which the reference appears."),
    /** How something strange has been published. The first word should be capitalized. */
    HOWPUBLISHED("howpublished", "How something strange has been published. The first word should be capitalized."),
    /** The sponsoring institution of a technical report. */
    INSTITUTION("institution", "The sponsoring institution of a technical report."),
    /** A journal name. Abbreviations are provided for many journals. */
    JOURNAL("journal", "A journal name. Abbreviations are provided for many journals."),
    /** Used for alphabetizing, cross referencing, and creating a label when the ``author'' information is missing. This field should not be confused with the key that appears in the cite command and at the beginning of the database entry. */
    KEY("key", "Used for alphabetizing, cross referencing, and creating a label when the ``author'' information is missing. This field should not be confused with the key that appears in the cite command and at the beginning of the database entry."),
    /** The month in which the work was published or, for an unpublished work, in which it was written. You should use the standard three-letter abbreviation, as described in Appendix B.1.3 of the LaTeX book. */
    MONTH("month", "The month in which the work was published or, for an unpublished work, in which it was written. You should use the standard three-letter abbreviation, as described in Appendix B.1.3 of the LaTeX book."),
    /** Any additional information that can help the reader. The first word should be capitalized. */
    NOTE("note", "Any additional information that can help the reader. The first word should be capitalized."),
    /** The number of a journal, magazine, technical report, or of a work in a series. An issue of a journal or magazine is usually identified by its volume and number; the organization that issues a technical report usually gives it a number; and sometimes books are given numbers in a named series. */
    NUMBER("number", "The number of a journal, magazine, technical report, or of a work in a series. An issue of a journal or magazine is usually identified by its volume and number; the organization that issues a technical report usually gives it a number; and sometimes books are given numbers in a named series."),
    /** The organization that sponsors a conference or that publishes a manual. */
    ORGANIZATION("organization", "The organization that sponsors a conference or that publishes a manual."),
    /** One or more page numbers or range of numbers, such as 42--111 or 7,41,73--97 or 43+ (the `+' in this last example indicates pages following that don't form a simple range). To make it easier to maintain Scribe-compatible databases, the standard styles convert a single dash (as in 7-33) to the double dash used in TeX to denote number ranges (as in 7--33). */
    PAGES("pages", "One or more page numbers or range of numbers, such as 42--111 or 7,41,73--97 or 43+ (the `+' in this last example indicates pages following that don't form a simple range). To make it easier to maintain Scribe-compatible databases, the standard styles convert a single dash (as in 7-33) to the double dash used in TeX to denote number ranges (as in 7--33)."),
    /** The publisher's name. */
    PUBLISHER("publisher", "The publisher's name."),
    /** The name of the school where a thesis was written. */
    SCHOOL("school", "The name of the school where a thesis was written."),
    /** The name of a series or set of books. When citing an entire book, the the title field gives its title and an optional series field gives the name of a series or multi-volume set in which the book is published. */
    SERIES("series", "The name of a series or set of books. When citing an entire book, the the title field gives its title and an optional series field gives the name of a series or multi-volume set in which the book is published."),
    /** The work's title, typed as explained in the LaTeX book. */
    TITLE("title", "The work's title, typed as explained in the LaTeX book."),
    /** The type of a technical report---for example, ``Research Note''. */
    TYPE("type", "The type of a technical report---for example, ``Research Note''."),
    /** The volume of a journal or multi-volume book. */
    VOLUME("volume", "The volume of a journal or multi-volume book."),
    /** The year of publication or, for an unpublished work, the year it was written. Generally it should consist of four numerals, such as 1984, although the standard styles can handle any year whose last four nonpunctuation characters are numerals, such as `\\hbox{(about 1984)}'. */
    YEAR("year", "The year of publication or, for an unpublished work, the year it was written. Generally it should consist of four numerals, such as 1984, although the standard styles can handle any year whose last four nonpunctuation characters are numerals, such as `\\hbox{(about 1984)}'."),
    // other fields
    /** The authors affiliation. */
    AFFILIATION("affiliation", "The authors affiliation."),
    /** An abstract of the work. */
    ABSTRACT("abstract", "An abstract of the work."),
    /** A Table of Contents.  */
    CONTENTS("contents", "A Table of Contents "),
    /** Copyright information. */
    COPYRIGHT("copyright", "Copyright information."),
    /** The International Standard Book Number (10 digits). */
    ISBN("ISBN", "The International Standard Book Number (10 digits)."),
    /** The International Standard Book Number (13 digits). */
    ISBN13("ISBN-13", "The International Standard Book Number (13 digits)."),
    /** The International Standard Serial Number. Used to identify a journal. */
    ISSN("ISSN", "The International Standard Serial Number. Used to identify a journal."),
    /** Key words used for searching or possibly for annotation. */
    KEYWORDS("keywords", "Key words used for searching or possibly for annotation."),
    /** The language the document is in. */
    LANGUAGE("language", "The language the document is in."),
    /** A location associated with the entry, such as the city in which a conference took place. */
    LOCATION("location", "A location associated with the entry, such as the city in which a conference took place."),
    /** The Library of Congress Call Number. I've also seen this as lib-congress. */
    LCCN("LCCN", "The Library of Congress Call Number. I've also seen this as lib-congress."),
    /** The Mathematical Reviews number. */
    MRNUMBER("mrnumber", "The Mathematical Reviews number."),
    /** The price of the document. */
    PRICE("price", "The price of the document."),
    /** The physical dimensions of a work. */
    SIZE("size", "The physical dimensions of a work."),
    /** The WWW Universal Resource Locator that points to the item being referenced. This often is used for technical reports to point to the ftp site where the postscript source of the report is located. */
    URL("URL", "The WWW Universal Resource Locator that points to the item being referenced. This often is used for technical reports to point to the ftp site where the postscript source of the report is located."),
    // additional fields
    /** A link to a postscript file. */
    PS("PS", "A link to a postscript file."),
    /** A link to a postscript file. */
    PDF("PDF", "A link to a PDF file."),
    /** A link to a postscript file. */
    HTTP("HTTP", "A hyperlink to a resource.");
    
    /** the string used in toString()  */
    protected String m_Display;
    
    /** the comment about this type */
    protected String m_Comment;
    
    /**
     * the constructor
     * 
     * @param display	the string to return in toString()
     * @param comment	a comment about the type
     */
    private Field(String display, String comment) {
      m_Display = display;
      m_Comment = comment;
    }
    
    /**
     * returns the display string
     * 
     * @return		the display string
     */
    public String getDisplay() {
      return m_Display;
    }
    
    /**
     * returns the comment string
     * 
     * @return		the comment string
     */
    public String getComment() {
      return m_Comment;
    }
    
    /**
     * returns the display string of the Type
     * 
     * @return		the display string
     */
    public String toString() {
      return m_Display;
    }
  }

  /** will be returned if no ID can be generated */
  protected final static String MISSING_ID = "missing_id";
  
  /** the type of this technical information */
  protected Type m_Type = null;
  
  /** the unique identifier of this information, will be generated 
   * automatically if left empty */
  protected String m_ID = "";
  
  /** stores all the values associated with the fields (FIELD - String) */
  protected Hashtable<Field,String> m_Values = new Hashtable<Field,String>();

  /** additional technical informations */
  protected Vector<TechnicalInformation> m_Additional = new Vector<TechnicalInformation>();
  
  /**
   * Initializes the information with the given type
   * 
   * @param type	the type of this information
   * @see Type
   */
  public TechnicalInformation(Type type) {
    this(type, "");
  }

  /**
   * Initializes the information with the given type
   * 
   * @param type	the type of this information
   * @param id		the unique ID (for BibTex), can be empty
   * @see Type
   */
  public TechnicalInformation(Type type, String id) {
    m_Type = type;
    m_ID   = id;
  }
  
  /**
   * returns the type of this technical information
   * 
   * @return		the type of this information
   */
  public Type getType() {
    return m_Type;
  }
  
  /**
   * splits the authors on the " and " and returns a vector with the names
   * 
   * @return		the authors in an array
   */
  protected String[] getAuthors() {
    return getValue(Field.AUTHOR).split(" and ");
  }
  
  /**
   * Generates an ID based on the current settings and returns it. If nothing
   * can be generated the MISSING_ID string will be returned
   * 
   * @return		the ID
   * @see #MISSING_ID
   */
  protected String generateID() {
    String 	result;
    String[]	authors;
    String[]	parts;
    
    result = m_ID;
    
    // try surname of first author and year
    if (result.length() == 0) {
      if (exists(Field.AUTHOR) && exists(Field.YEAR)) {
	authors = getAuthors();
	if (authors[0].indexOf(",") > -1) {
	  parts = authors[0].split(",");
	  result = parts[0];
	}
	else {
	  parts = authors[0].split(" ");
	  if (parts.length == 1)
	    result = parts[0];
	  else
	    result = parts[parts.length - 1];
	}
	result += getValue(Field.YEAR);
	result = result.replaceAll(" ", "");
      }
    }
    
    // still nothing?
    if (result.length() == 0)
      result = MISSING_ID;
    
    return result;
  }
  
  /**
   * returns the unique ID (either the one used in creating this instance
   * or the automatically generated one)
   * 
   * @return		the ID for this information
   */
  public String getID() {
    return generateID();
  }
  
  /**
   * sets the value for the given field, overwrites any previously existing one.
   * 
   * @param field	the field to set the value for
   * @param value	the value of the field
   */
  public void setValue(Field field, String value) {
    m_Values.put(field, value);
  }
  
  /**
   * returns the value associated with the given field, or empty if field is
   * not currently stored.
   * 
   * @param field	the field to retrieve the value for
   * @return		the value associated with this field, empty if not existing
   */
  public String getValue(Field field) {
    if (m_Values.containsKey(field))
      return (String) m_Values.get(field);
    else
      return "";
  }
  
  /**
   * returns TRUE if the field is stored and has a value different from the
   * empty string.
   * 
   * @param field	the field to check
   * @return 		true if a value is stored for the field and non-empty
   */
  public boolean exists(Field field) {
    return (    m_Values.containsKey(field) 
	     && (((String) m_Values.get(field)).length() != 0) );
  }
  
  /**
   * returns an enumeration over all the stored fields
   * 
   * @return		all currently stored fields
   */
  public Enumeration<Field> fields() {
    return m_Values.keys();
  }
  
  /**
   * returns true if there are further technical informations stored in this
   * 
   * @return true if there are further technical informations available
   */
  public boolean hasAdditional() {
    return (m_Additional.size() > 0);
  }
  
  /**
   * returns an enumeration of all the additional technical informations (if
   * there are any)
   * 
   * @return an enumeration over all additional technical informations
   */
  public Enumeration<TechnicalInformation> additional() {
    return m_Additional.elements();
  }
  
  /**
   * adds the given information to the list of additional technical 
   * informations
   * 
   * @param value the information to add
   */
  public void add(TechnicalInformation value) {
    if (value == this)
      throw new IllegalArgumentException("Can't add object to itself!");
    m_Additional.add(value);
  }
  /**
   * Adds an empty technical information with the given type to the list
   * of additional informations and returns the instance.
   * 
   * @param type	the type of the new information to add
   * @return 		the generated information
   */
  public TechnicalInformation add(Type type) {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(type);
    add(result);
    
    return result;
  }
  
  /**
   * Returns a plain-text string representing this technical information.
   * Note: it only returns a string based on some fields. At least AUTHOR,
   * YEAR and TITLE are necessary.
   * 
   * @return		a string representation of this information
   */
  public String toString() {
    String	result;
    String[]	authors;
    int		i;
    Enumeration<TechnicalInformation>	enm;
    
    result  = "";
    authors = getAuthors();
    
    // BOOK
    if (getType() == Type.BOOK) {
      for (i = 0; i < authors.length; i++) {
	if (i > 0)
	  result += ", ";
	result += authors[i];
      }
      if (exists(Field.YEAR))
	result += " (" + getValue(Field.YEAR) + ").";
      else
	result += ".";
      result += " " + getValue(Field.TITLE) + ".";
      result += " " + getValue(Field.PUBLISHER);
      if (exists(Field.ADDRESS))
	result += ", " + getValue(Field.ADDRESS);
      result += ".";
    }
    // ARTICLE
    else if (getType() == Type.ARTICLE) {
      for (i = 0; i < authors.length; i++) {
	if (i > 0)
	  result += ", ";
	result += authors[i];
      }
      if (exists(Field.YEAR))
	result += " (" + getValue(Field.YEAR) + ").";
      else
	result += ".";
      result += " " + getValue(Field.TITLE) + ".";
      
      // journal
      if (exists(Field.JOURNAL)) {
	result += " " + getValue(Field.JOURNAL) + ".";

	if (exists(Field.VOLUME))
	  result += " " + getValue(Field.VOLUME);
	if (exists(Field.NUMBER))
	  result += "(" + getValue(Field.NUMBER) + ")";
	if (exists(Field.PAGES))
	  result += ":" + getValue(Field.PAGES);
	
	result += ".";
      }
      
      // other than JOURNAL???

      // URL
      if (exists(Field.URL))
	result += " URL " + getValue(Field.URL) + ".";
    }
    // CONFERENCE/INPROCEEDINGS
    else if ( (getType() == Type.CONFERENCE) || (getType() == Type.INPROCEEDINGS) ) {
      for (i = 0; i < authors.length; i++) {
	if (i > 0)
	  result += ", ";
	result += authors[i];
      }
      result += ": " + getValue(Field.TITLE) + ".";
      result += " In: " + getValue(Field.BOOKTITLE);
      
      if (exists(Field.ADDRESS))
	result += ", " + getValue(Field.ADDRESS);
      if (exists(Field.PAGES))
	result += ", " + getValue(Field.PAGES);
	
      if (exists(Field.YEAR))
	result += ", " + getValue(Field.YEAR) + ".";
      else
	result += ".";
    }
    // INCOLLECTION
    else if (getType() == Type.INCOLLECTION) {
      for (i = 0; i < authors.length; i++) {
	if (i > 0)
	  result += ", ";
	result += authors[i];
      }
      result += ": " + getValue(Field.TITLE) + ".";
      result += " In ";
      if (exists(Field.EDITOR))
	result += getValue(Field.EDITOR) + ", editors, ";
      result += getValue(Field.BOOKTITLE);
      
      if (exists(Field.ADDRESS))
	result += ", " + getValue(Field.ADDRESS);
      if (exists(Field.PAGES))
	result += ", " + getValue(Field.PAGES);
	
      if (exists(Field.YEAR))
	result += ", " + getValue(Field.YEAR) + ".";
      else
	result += ".";
    }
    // default
    else {
      for (i = 0; i < authors.length; i++) {
	if (i > 0)
	  result += ", ";
	result += authors[i];
      }
      if (exists(Field.YEAR))
	result += " (" + getValue(Field.YEAR) + ").";
      else
	result += ".";
      result += " " + getValue(Field.TITLE) + ".";
      if (exists(Field.ADDRESS))
	result += " " + getValue(Field.ADDRESS) + ".";
      if (exists(Field.URL))
	result += " URL " + getValue(Field.URL) + ".";
    }
    
    // additional informations?
    enm = additional();
    while (enm.hasMoreElements()) {
      result += "\n\n" + enm.nextElement().toString();
    }
    
    return result;
  }
  
  /**
   * Returns a BibTex string representing this technical information. 
   * Note: this is just a very raw implementation, special characters need to
   * be escaped manually for LaTeX.
   * 
   * @return		the BibTeX representation of this information
   */
  public String toBibTex() {
    String		result;
    Field		field;
    Vector<Field>		list;
    int			i;
    String		value;
    
    result = "@" + getType() + "{" + getID() + "";
    
    // sort the fields
    list = new Vector<Field>();
    Enumeration<Field> enm  = fields();
    while (enm.hasMoreElements())
      list.add(enm.nextElement());
    Collections.sort(list);
    
    // list field=value pairs
    for (i = 0; i < list.size(); i++) {
      field = (Field) list.get(i);
      if (!exists(field))
	continue;
      value = getValue(field);
      value = value.replaceAll("\\~", "\\\\~");
      result += ",\n   " + field + " = {" + value +  "}";
    }
    
    result += "\n}";
    
    // additional informations?
    Enumeration<TechnicalInformation> enm2 = additional();
    while (enm2.hasMoreElements()) {
      result += "\n\n" + ((TechnicalInformation) enm2.nextElement()).toBibTex();
    }
    
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
   * Prints some examples of technical informations if there are no 
   * commandline options given. Otherwise the information of a given 
   * TechnicalInformationHandler can be printed. <p/>
   * 
   * Valid options are: <p/>
   * 
   * -W classname <br/>
   *  The classname of the TechnicalInformationHandler to print the 
   *  information for <p/>
   *  
   * -bibtex <br/>
   *  Print the information in BibTeX format <p/>
   *  
   * -plaintext <br/>
   *  Print the information in plain text format <p/>
   * 
   * @param args 	the commandline options
   * @throws Exception	if the option parsing fails
   */
  public static void main(String[] args) throws Exception {
    TechnicalInformation	info;
    TechnicalInformation	additional;
    String			tmpStr;
    Class			cls;
    TechnicalInformationHandler	handler;
    
    // example from command line
    if (args.length != 0) {
      info = null;
      
      tmpStr = Utils.getOption('W', args);
      if (tmpStr.length() != 0) {
	cls     = Class.forName(tmpStr);
	handler = (TechnicalInformationHandler) cls.newInstance();
	info    = handler.getTechnicalInformation();
      }
      else {
	throw new IllegalArgumentException("A classname has to be provided with the -W option!");
      }
      
      if (Utils.getFlag("bibtex", args))
        System.out.println("\n" + handler.getClass().getName() + ":\n" + info.toBibTex());

      if (Utils.getFlag("plaintext", args))
        System.out.println("\n" + handler.getClass().getName() + ":\n" + info.toString());
    }
    else {
      // first example
      info = new TechnicalInformation(Type.BOOK);
      info.setValue(Field.AUTHOR, "Ross Quinlan");
      info.setValue(Field.YEAR, "1993");
      info.setValue(Field.TITLE, "C4.5: Programs for Machine Learning");
      info.setValue(Field.PUBLISHER, "Morgan Kaufmann Publishers");
      info.setValue(Field.ADDRESS, "San Mateo, CA");
      additional = info;
      
      System.out.println("\ntoString():\n" + info.toString());
      System.out.println("\ntoBibTex():\n" + info.toBibTex());
      
      // second example
      info = new TechnicalInformation(Type.INPROCEEDINGS);
      info.setValue(Field.AUTHOR, "Freund, Y. and Mason, L.");
      info.setValue(Field.YEAR, "1999");
      info.setValue(Field.TITLE, "The alternating decision tree learning algorithm");
      info.setValue(Field.BOOKTITLE, "Proceeding of the Sixteenth International Conference on Machine Learning");
      info.setValue(Field.ADDRESS, "Bled, Slovenia");
      info.setValue(Field.PAGES, "124-133");
      
      System.out.println("\ntoString():\n" + info.toString());
      System.out.println("\ntoBibTex():\n" + info.toBibTex());
      
      // third example
      info = new TechnicalInformation(Type.ARTICLE);
      info.setValue(Field.AUTHOR, "R. Quinlan");
      info.setValue(Field.YEAR, "1986");
      info.setValue(Field.TITLE, "Induction of decision trees");
      info.setValue(Field.JOURNAL, "Machine Learning");
      info.setValue(Field.VOLUME, "1");
      info.setValue(Field.NUMBER, "1");
      info.setValue(Field.PAGES, "81-106");
      
      additional = new TechnicalInformation(Type.BOOK);
      additional.setValue(Field.AUTHOR, "Ross Quinlan");
      additional.setValue(Field.YEAR, "1993");
      additional.setValue(Field.TITLE, "C4.5: Programs for Machine Learning");
      additional.setValue(Field.PUBLISHER, "Morgan Kaufmann Publishers");
      additional.setValue(Field.ADDRESS, "San Mateo, CA");
      info.add(additional);
      
      System.out.println("\ntoString():\n" + info.toString());
      System.out.println("\ntoBibTex():\n" + info.toBibTex());
    }
  }
}

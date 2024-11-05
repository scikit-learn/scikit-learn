// Filename: SyntaxDocument.java
// http://forums.sun.com/thread.jspa?threadID=743407&messageID=9497681
// http://www.dound.com/src/MultiSyntaxDocument.java

/**
 * (C) camickr (primary author; java sun forums user)
 * (C) David Underhill
 * (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import weka.gui.visualize.VisualizeUtils;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Toolkit;
import java.util.HashMap;
import java.util.Properties;

import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultEditorKit;
import javax.swing.text.DefaultStyledDocument;
import javax.swing.text.Element;
import javax.swing.text.MutableAttributeSet;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.TabSet;
import javax.swing.text.TabStop;

/**
 * Highlights syntax in a DefaultStyledDocument. Allows any number of keywords
 * to be formatted in any number of user-defined styles.
 * 
 * @author camickr (primary author; java sun forums user)
 * @author David Underhill
 * @author FracPete (fracpete at waikato dot ac dot nz) - use of a properties file to setup syntax highlighting instead of hard-coded
 */
public class SyntaxDocument
  extends DefaultStyledDocument {

  /** for serialization. */
  protected static final long serialVersionUID = -3642426465631271381L;

  /** the maximum number of tabs. */
  public static final int MAX_TABS = 35;
  
  /** the font family. */
  public static final String DEFAULT_FONT_FAMILY = "monospaced";

  /** the font size. */
  public static final int DEFAULT_FONT_SIZE = 12;

  /** the attribute set for normal code. */
  public static final SimpleAttributeSet DEFAULT_NORMAL;

  /** the attribute set for comments. */
  public static final SimpleAttributeSet DEFAULT_COMMENT;

  /** the attribute set for strings. */
  public static final SimpleAttributeSet DEFAULT_STRING;

  /** the attribute set for keywords. */
  public static final SimpleAttributeSet DEFAULT_KEYWORD;

  static {
    DEFAULT_NORMAL = new SimpleAttributeSet();
    StyleConstants.setForeground(DEFAULT_NORMAL, Color.BLACK);
    StyleConstants.setFontFamily(DEFAULT_NORMAL, DEFAULT_FONT_FAMILY);
    StyleConstants.setFontSize(DEFAULT_NORMAL, DEFAULT_FONT_SIZE);

    DEFAULT_COMMENT = new SimpleAttributeSet();
    StyleConstants.setForeground(DEFAULT_COMMENT, Color.GRAY);
    StyleConstants.setFontFamily(DEFAULT_COMMENT, DEFAULT_FONT_FAMILY);
    StyleConstants.setFontSize(DEFAULT_COMMENT, DEFAULT_FONT_SIZE);

    DEFAULT_STRING = new SimpleAttributeSet();
    StyleConstants.setForeground(DEFAULT_STRING, Color.RED);
    StyleConstants.setFontFamily(DEFAULT_STRING, DEFAULT_FONT_FAMILY);
    StyleConstants.setFontSize(DEFAULT_STRING, DEFAULT_FONT_SIZE);

    // default style for new keyword types
    DEFAULT_KEYWORD = new SimpleAttributeSet();
    StyleConstants.setForeground(DEFAULT_KEYWORD, Color.BLUE);
    StyleConstants.setBold(DEFAULT_KEYWORD, false);
    StyleConstants.setFontFamily(DEFAULT_KEYWORD, DEFAULT_FONT_FAMILY);
    StyleConstants.setFontSize(DEFAULT_KEYWORD, DEFAULT_FONT_SIZE);
  }

  /**
   * The attribute type.
   */
  public enum ATTR_TYPE {
    /** normal string. */
    Normal,
    /** a comment. */
    Comment,
    /** a quoted string. */
    Quote;
  }

  /** the document. */
  protected DefaultStyledDocument m_Self;

  /** the root element. */
  protected Element m_RootElement;

  /** whether we're currently in a multi-line comment. */
  protected boolean m_InsideMultiLineComment;

  /** the keywords. */
  protected HashMap<String, MutableAttributeSet> m_Keywords;
  
  /** the delimiters. */
  protected String m_Delimiters;
  
  /** the quote delimiter. */
  protected String m_QuoteDelimiters;
  
  /** the quote escape. */
  protected String m_QuoteEscape;
  
  /** the multi-line comment start. */
  protected String m_MultiLineCommentStart;
  
  /** the multi-line comment end. */
  protected String m_MultiLineCommentEnd;
  
  /** the single-line comment start. */
  protected String m_SingleLineCommentStart;

  /** the start of a block. */
  protected String m_BlockStart;

  /** the end of a block. */
  protected String m_BlockEnd;
  
  /** the font size. */
  protected int m_FontSize;

  /** the font name. */
  protected String m_FontName;
  
  /** the background color. */
  protected Color m_BackgroundColor;
  
  /** the number of spaces used for indentation. */
  protected String m_Indentation;

  /** whether to add matching brackets. */
  protected boolean m_AddMatchingEndBlocks;
  
  /** whether to use blanks instead of tabs. */
  protected boolean m_UseBlanks;
  
  /** whether multi-line comments are enabled. */
  protected boolean m_MultiLineComment;
  
  /** whether keywords are case-sensitive. */
  protected boolean m_CaseSensitive;
  
  /**
   * Initializes the document.
   * 
   * @param props	the properties to obtain the setup from
   */
  public SyntaxDocument(Properties props) {
    m_Self        = this;
    m_RootElement = m_Self.getDefaultRootElement();
    m_Keywords    = new HashMap<String, MutableAttributeSet>();
    m_FontSize    = DEFAULT_FONT_SIZE;
    m_FontName    = DEFAULT_FONT_FAMILY;
    putProperty(DefaultEditorKit.EndOfLineStringProperty, "\n");
    
    setup(props);
  }
  
  /**
   * Sets up the document according to the properties.
   * 
   * @param props	the properties to use
   */
  protected void setup(Properties props) {
    setDelimiters(props.getProperty("Delimiters", ";:{}()[]+-/%<=>!&|^~*"));
    setQuoteDelimiters(props.getProperty("QuoteDelimiters", "\"\'"));
    setQuoteEscape(props.getProperty("QuoteEscape", "\\"));
    setSingleLineCommentStart(props.getProperty("SingleLineCommentStart", "//"));
    setMultiLineComment(props.getProperty("MultiLineComment", "false").equals("true"));
    setMultiLineCommentStart(props.getProperty("MultiLineCommentStart", "/*"));
    setMultiLineCommentEnd(props.getProperty("MultiLineCommentEnd", "*/"));
    setBlockStart(props.getProperty("BlockStart", "{"));
    setBlockEnd(props.getProperty("BlockEnd", "}"));
    setAddMatchingEndBlocks(props.getProperty("AddMatchingBlockEnd", "false").equals("true"));
    setUseBlanks(props.getProperty("UseBlanks", "false").equals("true"));
    setCaseSensitive(props.getProperty("CaseSensitive", "true").equals("true"));
    addKeywords(props.getProperty("Keywords", "").trim().replaceAll(" ", "").split(","), DEFAULT_KEYWORD);
    setTabs(Integer.parseInt(props.getProperty("Tabs", "2")));
    setAttributeColor(DEFAULT_NORMAL, VisualizeUtils.processColour(props.getProperty("ForegroundColor", "black"), Color.BLACK));
    setAttributeColor(DEFAULT_COMMENT, VisualizeUtils.processColour(props.getProperty("CommentColor", "gray"), Color.GRAY));
    setAttributeColor(DEFAULT_STRING, VisualizeUtils.processColour(props.getProperty("StringColor", "red"), Color.RED));
    setAttributeColor(DEFAULT_KEYWORD, VisualizeUtils.processColour(props.getProperty("KeywordColor", "blue"), Color.BLUE));
    setBackgroundColor(VisualizeUtils.processColour(props.getProperty("BackgroundColor", "white"), Color.WHITE));
    setFontName(props.getProperty("FontName", "monospaced"));
    setFontSize(Integer.parseInt(props.getProperty("FontSize", "12")));
    setIndentationSize(Integer.parseInt(props.getProperty("Indentation", "2")));
  }

  /**
   * Sets the font of the specified attribute.
   * 
   * @param attr
   *          the attribute to apply this font to (normal, comment, string)
   * @param style
   *          font style (Font.BOLD, Font.ITALIC, Font.PLAIN)
   */
  public void setAttributeFont(ATTR_TYPE attr, int style) {
    Font f = new Font(m_FontName, style, m_FontSize);

    if (attr == ATTR_TYPE.Comment)
      setAttributeFont(DEFAULT_COMMENT, f);
    else if (attr == ATTR_TYPE.Quote)
      setAttributeFont(DEFAULT_STRING, f);
    else
      setAttributeFont(DEFAULT_NORMAL, f);
  }

  /**
   * Sets the font of the specified attribute.
   * 
   * @param attr
   *          attribute to apply this font to
   * @param f
   *          the font to use
   */
  public static void setAttributeFont(MutableAttributeSet attr, Font f) {
    StyleConstants.setBold(attr, f.isBold());
    StyleConstants.setItalic(attr, f.isItalic());
    StyleConstants.setFontFamily(attr, f.getFamily());
    StyleConstants.setFontSize(attr, f.getSize());
  }

  /**
   * Sets the foreground (font) color of the specified attribute.
   * 
   * @param attr
   *          the attribute to apply this font to (normal, comment, string)
   * @param c
   *          the color to use
   */
  public void setAttributeColor(ATTR_TYPE attr, Color c) {
    if (attr == ATTR_TYPE.Comment)
      setAttributeColor(DEFAULT_COMMENT, c);
    else if (attr == ATTR_TYPE.Quote)
      setAttributeColor(DEFAULT_STRING, c);
    else
      setAttributeColor(DEFAULT_NORMAL, c);
  }

  /**
   * Sets the foreground (font) color of the specified attribute.
   * 
   * @param attr
   *          attribute to apply this color to
   * @param c
   *          the color to use
   */
  public static void setAttributeColor(MutableAttributeSet attr, Color c) {
    StyleConstants.setForeground(attr, c);
  }

  /**
   * Associates the keywords with a particular formatting style.
   * 
   * @param keywords
   *          the tokens or words to format
   * @param attr
   *          how to format the keywords
   */
  public void addKeywords(String[] keywords, MutableAttributeSet attr) {
    int		i;
    
    for (i = 0; i < keywords.length; i++)
      addKeyword(keywords[i], attr);
  }
  
  /**
   * Associates a keyword with a particular formatting style.
   * 
   * @param keyword
   *          the token or word to format
   * @param attr
   *          how to format keyword
   */
  public void addKeyword(String keyword, MutableAttributeSet attr) {
    if (m_CaseSensitive)
      m_Keywords.put(keyword, attr);
    else
      m_Keywords.put(keyword.toLowerCase(), attr);
  }

  /**
   * Gets the formatting for a keyword.
   * 
   * @param keyword
   *          the token or word to stop formatting
   * 
   * @return how keyword is formatted, or null if no formatting is applied to it
   */
  public MutableAttributeSet getKeywordFormatting(String keyword) {
    if (m_CaseSensitive)
      return m_Keywords.get(keyword);
    else
      return m_Keywords.get(keyword.toLowerCase());
  }

  /**
   * Removes an association between a keyword with a particular formatting style.
   * 
   * @param keyword
   *          the token or word to stop formatting
   */
  public void removeKeyword(String keyword) {
    if (m_CaseSensitive)
      m_Keywords.remove(keyword);
    else
      m_Keywords.remove(keyword.toLowerCase());
  }

  /** 
   * sets the number of characters per tab.
   * 
   * @param charactersPerTab
   * 		the characters per tab
   */
  public void setTabs(int charactersPerTab) {
    Font f = new Font(m_FontName, Font.PLAIN, m_FontSize);

    FontMetrics fm = Toolkit.getDefaultToolkit().getFontMetrics(f);
    int charWidth = fm.charWidth('w');
    int tabWidth = charWidth * charactersPerTab;

    TabStop[] tabs = new TabStop[MAX_TABS];

    for (int j = 0; j < tabs.length; j++)
      tabs[j] = new TabStop((j+1) * tabWidth);

    TabSet tabSet = new TabSet(tabs);
    SimpleAttributeSet attributes = new SimpleAttributeSet();
    StyleConstants.setTabSet(attributes, tabSet);
    int length = getLength();
    setParagraphAttributes(0, length, attributes, false);
  }

  /**
   * Override to apply syntax highlighting after the document has been updated.
   * 
   * @param offset
   * 		the offset
   * @param str
   * 		the string to insert
   * @param a
   * 		the attribute set, can be null
   * @throws BadLocationException
   * 		if offset is invalid
   */
  public void insertString(int offset, String str, AttributeSet a)
      throws BadLocationException {
    if (m_AddMatchingEndBlocks && (m_BlockStart.length() > 0) && str.equals(m_BlockStart))
      str = addMatchingBlockEnd(offset);
    else if (m_UseBlanks && str.equals("\t"))
      str = m_Indentation;
    
    super.insertString(offset, str, a);
    processChangedLines(offset, str.length());
  }

  /**
   * Applies syntax highlighting after the document has been updated.
   * 
   * @param offset
   * 		the offset of the deletion
   * @param length
   * 		the length of the deletion
   * @throws BadLocationException
   * 		if offsets are invalid
   */
  public void remove(int offset, int length) throws BadLocationException {
    super.remove(offset, length);
    processChangedLines(offset, 0);
  }

  /**
   * Determine how many lines have been changed, then apply highlighting to each
   * line.
   * 
   * @param offset
   * 		the offset of the changed lines
   * @param length
   * 		the length of the change
   * @throws BadLocationException
   * 		if offset is invalid
   */
  public void processChangedLines(int offset, int length)
      throws BadLocationException {
    String content = m_Self.getText(0, m_Self.getLength());

    // The lines affected by the latest document update

    int startLine = m_RootElement.getElementIndex(offset);
    int endLine = m_RootElement.getElementIndex(offset + length);

    // Make sure all comment lines prior to the start line are commented
    // and determine if the start line is still in a multi line comment

    if (getMultiLineComment())
      setInsideMultiLineComment(commentLinesBefore(content, startLine));

    // Do the actual highlighting

    for (int i = startLine; i <= endLine; i++) {
      applyHighlighting(content, i);
    }

    // Resolve highlighting to the next end multi line delimiter

    if (isMultiLineComment())
      commentLinesAfter(content, endLine);
    else
      highlightLinesAfter(content, endLine);
  }

  /**
   * Highlight lines when a multi line comment is still 'open' (ie. matching end
   * delimiter has not yet been encountered).
   * 
   * @param content
   * 		the content to check
   * @param line
   * 		the line number
   * @return
   * 		true if there are comment lines before
   */
  protected boolean commentLinesBefore(String content, int line) {
    int offset = m_RootElement.getElement(line).getStartOffset();

    // Start of comment not found, nothing to do

    int startDelimiter = -1;
    if (getMultiLineComment())
      startDelimiter = lastIndexOf(content, getMultiLineCommentStart(), offset - 2);

    if (startDelimiter < 0)
      return false;

    // Matching start/end of comment found, nothing to do

    int endDelimiter = indexOf(content, getMultiLineCommentEnd(), startDelimiter);

    if (endDelimiter < offset & endDelimiter != -1)
      return false;

    // End of comment not found, highlight the lines

    m_Self.setCharacterAttributes(startDelimiter, offset - startDelimiter + 1,
	DEFAULT_COMMENT, false);
    return true;
  }

  /**
   * Highlight comment lines to matching end delimiter.
   * 
   * @param content
   * 		the content to parse
   * @param line
   * 		the line number
   */
  protected void commentLinesAfter(String content, int line) {
    int offset = m_RootElement.getElement(line).getEndOffset();

    // End of comment not found, nothing to do

    int endDelimiter = -1;
    if (getMultiLineComment())
      endDelimiter = indexOf(content, getMultiLineCommentEnd(), offset);

    if (endDelimiter < 0)
      return;

    // Matching start/end of comment found, comment the lines

    int startDelimiter = lastIndexOf(content, getMultiLineCommentStart(), endDelimiter);

    if (startDelimiter < 0 || startDelimiter <= offset) {
      m_Self.setCharacterAttributes(offset, endDelimiter - offset + 1, DEFAULT_COMMENT,
	  false);
    }
  }

  /**
   * Highlight lines to start or end delimiter.
   * 
   * @param content
   * 		the content to parse
   * @param line
   * 		the line number
   * @throws BadLocationException
   * 		if offsets are wrong
   */
  protected void highlightLinesAfter(String content, int line)
      throws BadLocationException {
    int offset = m_RootElement.getElement(line).getEndOffset();

    // Start/End delimiter not found, nothing to do

    int startDelimiter = -1;
    int endDelimiter = -1;
    if (getMultiLineComment()) {
      startDelimiter = indexOf(content, getMultiLineCommentStart(), offset);
      endDelimiter = indexOf(content, getMultiLineCommentEnd(), offset);
    }

    if (startDelimiter < 0)
      startDelimiter = content.length();

    if (endDelimiter < 0)
      endDelimiter = content.length();

    int delimiter = Math.min(startDelimiter, endDelimiter);

    if (delimiter < offset)
      return;

    // Start/End delimiter found, reapply highlighting

    int endLine = m_RootElement.getElementIndex(delimiter);

    for (int i = line + 1; i < endLine; i++) {
      Element branch = m_RootElement.getElement(i);
      Element leaf = m_Self.getCharacterElement(branch.getStartOffset());
      AttributeSet as = leaf.getAttributes();

      if (as.isEqual(DEFAULT_COMMENT))
	applyHighlighting(content, i);
    }
  }

  /**
   * Parse the line to determine the appropriate highlighting.
   * 
   * @param content
   * 		the content to parse
   * @param line
   * 		the line number
   * @throws BadLocationException
   * 		if offsets are invalid
   */
  protected void applyHighlighting(String content, int line)
      throws BadLocationException {
    int startOffset = m_RootElement.getElement(line).getStartOffset();
    int endOffset = m_RootElement.getElement(line).getEndOffset() - 1;

    int lineLength = endOffset - startOffset;
    int contentLength = content.length();

    if (endOffset >= contentLength)
      endOffset = contentLength - 1;

    // check for multi line comments
    // (always set the comment attribute for the entire line)

    if (getMultiLineComment()) {
      if (endingMultiLineComment(content, startOffset, endOffset)
	  || isMultiLineComment()
	  || startingMultiLineComment(content, startOffset, endOffset)) {
	m_Self.setCharacterAttributes(startOffset, endOffset - startOffset + 1,
	    DEFAULT_COMMENT, false);
	return;
      }
    }

    // set normal attributes for the line

    m_Self.setCharacterAttributes(startOffset, lineLength, DEFAULT_NORMAL, true);

    // check for single line comment

    int index = content.indexOf(getSingleLineCommentStart(), startOffset);

    if ((index > -1) && (index < endOffset)) {
      m_Self.setCharacterAttributes(index, endOffset - index + 1, DEFAULT_COMMENT, false);
      endOffset = index - 1;
    }

    // check for tokens

    checkForTokens(content, startOffset, endOffset);
  }

  /**
   * Does this line contain the start of a multi-line comment.
   * 
   * @param content
   * 		the content to search
   * @param startOffset
   * 		the start of the search
   * @param endOffset
   * 		the end of the search
   * @return
   * 		true if it contains the start delimiter
   * @throws BadLocationException
   * 		if offsets are invalid
   */
  protected boolean startingMultiLineComment(String content, int startOffset,
      int endOffset) throws BadLocationException {
    if (!getMultiLineComment())
      return false;
    
    int index = indexOf(content, getMultiLineCommentStart(), startOffset);

    if ((index < 0) || (index > endOffset))
      return false;
    else {
      setInsideMultiLineComment(true);
      return true;
    }
  }

  /**
   * Does this line contain the end delimiter of a multi-line comment.
   * 
   * @param content
   * 		the content to search
   * @param startOffset
   * 		the start of the search
   * @param endOffset
   * 		the end of the search
   * @return
   * 		true if the line contains the end delimiter
   * @throws BadLocationException
   * 		if offsets are invalid
   */
  protected boolean endingMultiLineComment(String content, int startOffset,
      int endOffset) throws BadLocationException {
    if (!getMultiLineComment())
      return false;
    
    int index = indexOf(content, getMultiLineCommentEnd(), startOffset);

    if ((index < 0) || (index > endOffset))
      return false;
    else {
      setInsideMultiLineComment(false);
      return true;
    }
  }

  /**
   * We have found a start delimiter and are still searching for the end
   * delimiter.
   * 
   * @return
   * 		true if currently within a multi-line comment
   */
  protected boolean isMultiLineComment() {
    return m_InsideMultiLineComment;
  }

  /**
   * Sets whether we're currently within a multi-line comment or not.
   * 
   * @param value
   * 		true if currently within a multi-line comment
   */
  protected void setInsideMultiLineComment(boolean value) {
    m_InsideMultiLineComment = value;
  }

  /**
   * Parse the line for tokens to highlight.
   * 
   * @param content
   * 		the content to parse
   * @param startOffset
   * 		the start position
   * @param endOffset
   * 		the end position
   */
  protected void checkForTokens(String content, int startOffset, int endOffset) {
    while (startOffset <= endOffset) {
      // skip the delimiters to find the start of a new token

      while (isDelimiter(content.substring(startOffset, startOffset + 1))) {
	if (startOffset < endOffset)
	  startOffset++;
	else
	  return;
      }

      // Extract and process the entire token

      if (isQuoteDelimiter(content.substring(startOffset, startOffset + 1)))
	startOffset = getQuoteToken(content, startOffset, endOffset);
      else
	startOffset = getOtherToken(content, startOffset, endOffset);
    }
  }

  /**
   * Searches for a quote token.
   *
   * @param content
   * 		the content to search
   * @param startOffset
   * 		the start of the search
   * @param endOffset
   * 		the end of the search
   * @return
   * 		the new position
   */
  protected int getQuoteToken(String content, int startOffset, int endOffset) {
    String quoteDelimiter = content.substring(startOffset, startOffset + 1);
    String escapeString = escapeQuote(quoteDelimiter);

    int index;
    int endOfQuote = startOffset;

    // skip over the escape quotes in this quote

    index = content.indexOf(escapeString, endOfQuote + 1);

    while ((index > -1) && (index < endOffset)) {
      endOfQuote = index + 1;
      index = content.indexOf(escapeString, endOfQuote);
    }

    // now find the matching delimiter

    index = content.indexOf(quoteDelimiter, endOfQuote + 1);

    if ((index < 0) || (index > endOffset))
      endOfQuote = endOffset;
    else
      endOfQuote = index;

    m_Self.setCharacterAttributes(startOffset, endOfQuote - startOffset + 1,
	DEFAULT_STRING, false);

    return endOfQuote + 1;
  }

  /**
   * Searches for a keyword token.
   *
   * @param content
   * 		the content to search in
   * @param startOffset
   * 		the position to start the search fromm
   * @param endOffset
   * 		the position to end the search
   * @return
   * 		the new position
   */
  protected int getOtherToken(String content, int startOffset, int endOffset) {
    int endOfToken = startOffset + 1;

    while (endOfToken <= endOffset) {
      if (isDelimiter(content.substring(endOfToken, endOfToken + 1)))
	break;
      endOfToken++;
    }

    String token = content.substring(startOffset, endOfToken);

    // see if this token has a highlighting format associated with it
    MutableAttributeSet attr = getKeywordFormatting(token);
    if (attr != null)
      m_Self.setCharacterAttributes(startOffset, endOfToken - startOffset, attr, false);

    return endOfToken + 1;
  }

  /**
   * Assume the needle will the found at the start/end of the line.
   * 
   * @param content
   * 		the content to search
   * @param needle
   * 		the string to look for
   * @param offset
   * 		the offset to start at
   * @return
   * 		the index
   */
  protected int indexOf(String content, String needle, int offset) {
    int index;

    while ((index = content.indexOf(needle, offset)) != -1) {
      String text = getLine(content, index).trim();

      if (text.startsWith(needle) || text.endsWith(needle))
	break;
      else
	offset = index + 1;
    }

    return index;
  }

  /**
   * Assume the needle will the found at the start/end of the line.
   * 
   * @param content
   * 		the content search
   * @param needle
   * 		what to look for
   * @param offset
   * 		the offset to start
   * @return
   * 		the index
   */
  protected int lastIndexOf(String content, String needle, int offset) {
    int index;

    while ((index = content.lastIndexOf(needle, offset)) != -1) {
      String text = getLine(content, index).trim();

      if (text.startsWith(needle) || text.endsWith(needle))
	break;
      else
	offset = index - 1;
    }

    return index;
  }

  /**
   * Returns the line.
   * 
   * @param content
   * 		the content
   * @param offset
   * 		the offset to start at
   * @return
   * 		the line
   */
  protected String getLine(String content, int offset) {
    int line = m_RootElement.getElementIndex(offset);
    Element lineElement = m_RootElement.getElement(line);
    int start = lineElement.getStartOffset();
    int end = lineElement.getEndOffset();
    return content.substring(start, end - 1);
  }

  /**
   * Checks whether the character is a delimiter.
   * 
   * @param character
   * 		the character to check
   * @return
   * 		true if a delimiter
   */
  public boolean isDelimiter(String character) {
    return Character.isWhitespace(character.charAt(0)) || (m_Delimiters.indexOf(character.charAt(0)) > -1);
  }

  /**
   * Checks whether the character is quote delimiter.
   * 
   * @param character
   * 		the character to check
   * @return
   * 		true if a quote delimiter
   */
  public boolean isQuoteDelimiter(String character) {
    return (m_QuoteDelimiters.indexOf(character.charAt(0)) > -1);
  }

  /**
   * Escapes the quote delimiter.
   * 
   * @param quoteDelimiter
   * 		the string to escape
   * @return
   * 		the escaped string
   */
  public String escapeQuote(String quoteDelimiter) {
    return m_QuoteEscape + quoteDelimiter;
  }

  /**
   * Adds the matching block end.
   * 
   * @param offset
   * 		the offset
   * @return
   * 		the string after adding the matching block end
   * @throws BadLocationException
   * 		if the offset is invalid
   */
  protected String addMatchingBlockEnd(int offset) throws BadLocationException {
    StringBuffer result;
    StringBuffer whiteSpace = new StringBuffer();
    int line = m_RootElement.getElementIndex(offset);
    int i = m_RootElement.getElement(line).getStartOffset();

    while (true) {
      String temp = m_Self.getText(i, 1);

      if (temp.equals(" ") || temp.equals("\t")) {
	whiteSpace.append(temp);
	i++;
      }
      else {
	break;
      }
    }

    // assemble string
    result = new StringBuffer();
    result.append(m_BlockStart);
    result.append("\n");
    result.append(whiteSpace.toString());
    if (m_UseBlanks)
      result.append(m_Indentation);
    else
      result.append("\t");
    result.append("\n");
    result.append(whiteSpace.toString());
    result.append(m_BlockEnd);
    
    return result.toString();
  }

  /** 
   * gets the current font size.
   * 
   * @return
   * 		the font size
   */
  public int getFontSize() {
    return m_FontSize;
  }

  /** 
   * sets the current font size (affects all built-in styles).
   * 
   * @param fontSize
   * 		the size
   */
  public void setFontSize(int fontSize) {
    m_FontSize = fontSize;
    StyleConstants.setFontSize(DEFAULT_NORMAL, fontSize);
    StyleConstants.setFontSize(DEFAULT_STRING, fontSize);
    StyleConstants.setFontSize(DEFAULT_COMMENT, fontSize);
  }

  /**
   * gets the current font family.
   * 
   * @return
   * 		the font name
   */
  public String getFontName() {
    return m_FontName;
  }

  /**
   * sets the current font family (affects all built-in styles).
   * 
   * @param fontName
   * 		the font name
   */
  public void setFontName(String fontName) {
    m_FontName = fontName;
    StyleConstants.setFontFamily(DEFAULT_NORMAL, fontName);
    StyleConstants.setFontFamily(DEFAULT_STRING, fontName);
    StyleConstants.setFontFamily(DEFAULT_COMMENT, fontName);
  }
  
  /**
   * Sets the number of blanks to use for indentation.
   * 
   * @param value	
   * 		the number of blanks
   */
  public void setIndentationSize(int value) {
    int		i;
    
    m_Indentation = "";
    for (i = 0; i < value; i++)
      m_Indentation += " ";
  }
  
  /**
   * Returns the number of blanks used for indentation.
   * 
   * @return		
   * 		the number of blanks
   */
  public int getIndentationSize() {
    return m_Indentation.length();
  }
  
  /**
   * Sets the delimiter characters to use.
   * 
   * @param value	
   * 		the characters
   */
  public void setDelimiters(String value) {
    m_Delimiters = value;
  }
  
  /**
   * Returns the delimiter characters to use.
   * 
   * @return		
   * 		the characters
   */
  public String getDelimiters() {
    return m_Delimiters;
  }
  
  /**
   * Sets the quote delimiter characters to use.
   * 
   * @param value	
   * 		the characters
   */
  public void setQuoteDelimiters(String value) {
    m_QuoteDelimiters = value;
  }
  
  /**
   * Returns the quote delimiter characters to use.
   * 
   * @return		
   * 		the characters
   */
  public String getQuoteDelimiters() {
    return m_QuoteDelimiters;
  }
  
  /**
   * Sets the character to use for escaping a quote character.
   * 
   * @param value	
   * 		the character
   */
  public void setQuoteEscape(String value) {
    m_QuoteEscape = value;
  }
  
  /**
   * Returns the character for escaping a quote delimiter.
   * 
   * @return		
   * 		the character
   */
  public String getQuoteEscape() {
    return m_QuoteEscape;
  }
  
  /**
   * Sets the string that is the start of a single-line comment.
   * 
   * @param value	
   * 		the string
   */
  public void setSingleLineCommentStart(String value) {
    m_SingleLineCommentStart = value;
  }

  /**
   * Retrusn the single line comment start string.
   * 
   * @return
   * 		the start string
   */
  public String getSingleLineCommentStart() {
    return m_SingleLineCommentStart;
  }
  
  /**
   * Sets the string that is the start of a multi-line comment.
   * 
   * @param value	
   * 		the string
   */
  public void setMultiLineCommentStart(String value) {
    m_MultiLineCommentStart = value;
  }
  
  /**
   * Returns the string that is the start of a multi-line comment.
   * 
   * @return		
   * 		the string
   */
  public String getMultiLineCommentStart() {
    return m_MultiLineCommentStart;
  }
  
  /**
   * Sets the string that is the end of a multi-line comment.
   * 
   * @param value	the string
   */
  public void setMultiLineCommentEnd(String value) {
    m_MultiLineCommentEnd = value;
  }

  /**
   * Returns the end of a multi-line comment.
   * 
   * @return
   * 		the end string
   */
  public String getMultiLineCommentEnd() {
    return m_MultiLineCommentEnd;
  }
  
  /**
   * Sets the string that is the start of a block.
   * 
   * @param value
   * 		the string
   */
  public void setBlockStart(String value) {
    m_BlockStart = value;
  }

  /**
   * Returns the start of a block.
   * 
   * @return
   * 		the end string
   */
  public String getBlockStart() {
    return m_BlockStart;
  }
  
  /**
   * Sets the string that is the end of a block.
   * 
   * @param value	
   * 		the string
   */
  public void setBlockEnd(String value) {
    m_BlockEnd = value;
  }

  /**
   * Returns the end of a block.
   * 
   * @return
   * 		the end string
   */
  public String getBlockEnd() {
    return m_BlockEnd;
  }
  
  /**
   * Sets whether matching block ends are inserted or not.
   * 
   * @param value	
   * 		if true then matching block ends are inserted
   */
  public void setAddMatchingEndBlocks(boolean value) {
    m_AddMatchingEndBlocks = value;
  }

  /**
   * Returns whether matching block ends are inserted or not.
   * 
   * @return
   * 		true if matching block ends are inserted
   */
  public boolean getAddMatchingEndBlocks() {
    return m_AddMatchingEndBlocks;
  }
  
  /**
   * Sets whether to use blanks instead of tabs.
   * 
   * @param value	
   * 		if true then blanks are used instead of tabs
   */
  public void setUseBlanks(boolean value) {
    m_UseBlanks = value;
  }

  /**
   * Returns whether blanks are used instead of tabs.
   * 
   * @return
   * 		true if blanks are used instead of tabs
   */
  public boolean getUseBlanks() {
    return m_UseBlanks;
  }
  
  /**
   * Sets the background color.
   * 
   * @param value	
   * 		the background color
   */
  public void setBackgroundColor(Color value) {
    m_BackgroundColor = value;
  }
  
  /**
   * Returns the background color.
   * 
   * @return
   * 		the background color
   */
  public Color getBackgroundColor() {
    return m_BackgroundColor;
  }
  
  /**
   * Sets whether to enable multi-line comments.
   * 
   * @param value	
   * 		if true then multi-line comments are enabled
   */
  public void setMultiLineComment(boolean value) {
    m_MultiLineComment = value;
  }

  /**
   * Returns whether multi-line comments are enabled.
   * 
   * @return
   * 		true if multi-line comments are enabled
   */
  public boolean getMultiLineComment() {
    return m_MultiLineComment;
  }
  
  /**
   * Sets whether the keywords are case-sensitive or not.
   * 
   * @param value	
   * 		if true then keywords are treated case-sensitive
   */
  public void setCaseSensitive(boolean value) {
    m_CaseSensitive = value;
  }

  /**
   * Returns whether blanks are used instead of tabs.
   * 
   * @return
   * 		true if keywords are case-sensitive
   */
  public boolean getCaseSensitive() {
    return m_CaseSensitive;
  }
}

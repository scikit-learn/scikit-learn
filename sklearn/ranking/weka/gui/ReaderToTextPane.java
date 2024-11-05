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
 * ReaderToTextArea.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui;

import java.awt.Color;
import java.io.LineNumberReader;
import java.io.Reader;

import javax.swing.JTextPane;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

/**
 * A class that sends all lines from a reader to a JTextPane component.
 * 
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public class ReaderToTextPane
  extends Thread {

  /** The reader being monitored. */
  protected LineNumberReader m_Input;

  /** The output text component. */
  protected JTextPane m_Output;
  
  /** the color to use. */
  protected Color m_Color;

  /**
   * Sets up the thread. Using black as color for displaying the text.
   *
   * @param input 	the Reader to monitor
   * @param output 	the TextArea to send output to
   */
  public ReaderToTextPane(Reader input, JTextPane output) {
    this(input, output, Color.BLACK);
  }

  /**
   * Sets up the thread.
   *
   * @param input 	the Reader to monitor
   * @param output 	the TextArea to send output to
   * @param color	the color to use
   */
  public ReaderToTextPane(Reader input, JTextPane output, Color color) {
    StyledDocument      doc;
    Style               style;

    setDaemon(true);
    
    m_Color  = color;
    m_Input  = new LineNumberReader(input);
    m_Output = output;
    
    doc   = m_Output.getStyledDocument();
    style = StyleContext.getDefaultStyleContext()
                        .getStyle(StyleContext.DEFAULT_STYLE);
    style = doc.addStyle(getStyleName(), style);
    StyleConstants.setFontFamily(style, "monospaced");
    StyleConstants.setForeground(style, m_Color);
  }
  
  /**
   * Returns the color in use.
   * 
   * @return		the color
   */
  public Color getColor() {
    return m_Color;
  }
  
  /**
   * Returns the style name.
   * 
   * @return		the style name
   */
  protected String getStyleName() {
    return "" + m_Color.hashCode();
  }

  /**
   * Sit here listening for lines of input and appending them straight
   * to the text component.
   */
  public void run() {
    while (true) {
      try {
	StyledDocument doc = m_Output.getStyledDocument();
	doc.insertString(
	    doc.getLength(), 
	    m_Input.readLine() + '\n', 
	    doc.getStyle(getStyleName()));
	m_Output.setCaretPosition(doc.getLength());
      }
      catch (Exception ex) {
	try {
	  sleep(100);
	}
	catch (Exception e) {
	  // ignored
	}
      }
    }
  }
}

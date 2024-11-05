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
 *    TextEvent.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

import java.util.EventObject;

/**
 * Event that encapsulates some textual information
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.4 $
 */
public class TextEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 4196810607402973744L;
  
  /**
   * The text
   */
  protected String m_text;

  /**
   * The title for the text. Could be used in a list component
   */
  protected String m_textTitle;

  /**
   * Creates a new <code>TextEvent</code> instance.
   *
   * @param source an <code>Object</code> value
   * @param text a <code>String</code> value
   */
  public TextEvent(Object source, String text, String textTitle) {
    super(source);
    
    m_text = text;
    m_textTitle = textTitle;
  }

  /**
   * Describe <code>getText</code> method here.
   *
   * @return a <code>String</code> value
   */
  public String getText() {
    return m_text;
  }

  /**
   * Describe <code>getTextTitle</code> method here.
   *
   * @return a <code>String</code> value
   */
  public String getTextTitle() {
    return m_textTitle;
  }
}

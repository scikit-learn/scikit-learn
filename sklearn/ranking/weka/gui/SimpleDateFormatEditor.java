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
 * SimpleDateFormatEditor.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.beans.PropertyEditor;
import java.text.SimpleDateFormat;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

/**
 * Class for editing SimpleDateFormat strings. 
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.2 $
 * @see SimpleDateFormat
 */
public class SimpleDateFormatEditor 
  implements PropertyEditor {

  /** the default format */
  public final static String DEFAULT_FORMAT = "yyyy-MM-dd'T'HH:mm:ss";
  
  /** The date format being edited */
  private SimpleDateFormat m_Format;

  /** A helper class for notifying listeners */
  private PropertyChangeSupport m_propSupport;

  /** An instance of the custom editor */
  private CustomEditor m_customEditor;

  /**
   * This class presents a GUI for editing the cost matrix, and saving and 
   * loading from files.
   */
  private class CustomEditor 
    extends JPanel 
    implements ActionListener, DocumentListener {

    /** for serialization */
    private static final long serialVersionUID = -4018834274636309987L;
    
    /** The text field for the format */
    private JTextField m_FormatText;
    
    /** The button for setting a default date format */
    private JButton m_DefaultButton;

    /** The button for applying the format */
    private JButton m_ApplyButton;

    /**
     * Constructs a new CustomEditor.
     */
    public CustomEditor() {
      m_FormatText    = new JTextField(20);
      m_DefaultButton = new JButton("Default");
      m_ApplyButton   = new JButton("Apply");

      m_DefaultButton.setMnemonic('D');
      m_ApplyButton.setMnemonic('A');

      m_FormatText.getDocument().addDocumentListener(this);
      m_DefaultButton.addActionListener(this);
      m_ApplyButton.addActionListener(this);

      setLayout(new FlowLayout());
      add(new JLabel("ISO 8601 Date format"));
      add(m_FormatText);
      add(m_DefaultButton);
      add(m_ApplyButton);
    }

    /**
     * Responds to the user perfoming an action.
     *
     * @param e the action event that occured
     */
    public void actionPerformed(ActionEvent e) {
      if (e.getSource() == m_DefaultButton)
	defaultFormat();
      else if (e.getSource() == m_ApplyButton)
	applyFormat();
    }

    /**
     * sets the format to default 
     */
    public void defaultFormat() {
      m_FormatText.setText(DEFAULT_FORMAT);
      formatChanged();
    }

    /**
     * returns true if the current format is a valid format
     */
    protected boolean isValidFormat() {
      boolean 	result;
      
      result = false;
      
      try {
	new SimpleDateFormat(m_FormatText.getText());
	result = true;
      }
      catch (Exception e) {
	// we can ignore this exception
      }
      
      return result;
    }
    
    /**
     * sets the format, but only if it's a valid one
     */
    public void applyFormat() {
      if (isValidFormat()) {
	m_Format = new SimpleDateFormat(m_FormatText.getText());
	m_propSupport.firePropertyChange(null, null, null);
      }
      else {
	throw new IllegalArgumentException(
	    "Date format '" 
	    + m_FormatText.getText() 
	    + "' is invalid! Cannot execute applyFormat!");
      }
    }
    
    /**
     * Responds to a change of the text field.
     */
    public void formatChanged() {
      m_FormatText.setText(m_Format.toPattern());
      m_propSupport.firePropertyChange(null, null, null);
    }
    
    /**
     * Gives notification that an attribute or set of attributes changed.
     */
    public void changedUpdate(DocumentEvent e) {
      m_ApplyButton.setEnabled(isValidFormat());
    }
    
    /**
     * Gives notification that there was an insert into the document.
     */
    public void insertUpdate(DocumentEvent e) {
      m_ApplyButton.setEnabled(isValidFormat());
    }
    
    /**
     * Gives notification that a portion of the document has been removed.
     */
    public void removeUpdate(DocumentEvent e) {
      m_ApplyButton.setEnabled(isValidFormat());
    }
  }

  /**
   * Constructs a new SimpleDateFormatEditor.
   *
   */
  public SimpleDateFormatEditor() {
    m_propSupport = new PropertyChangeSupport(this);
    m_customEditor = new CustomEditor();
  }

  /**
   * Sets the value of the date format to be edited.
   *
   * @param value a SimpleDateFormat object to be edited
   */
  public void setValue(Object value) {
    m_Format = (SimpleDateFormat) value;
    m_customEditor.formatChanged();
  }

  /**
   * Gets the date format that is being edited.
   *
   * @return the edited SimpleDateFormat object
   */  
  public Object getValue() {
    return m_Format;
  }

  /**
   * Indicates whether the object can be represented graphically. In this case
   * it can.
   *
   * @return true
   */  
  public boolean isPaintable() {
    return true;
  }

  /**
   * Paints a graphical representation of the object. It just prints the 
   * format.
   *
   * @param gfx the graphics context to draw the representation to
   * @param box the bounds within which the representation should fit.
   */    
  public void paintValue(Graphics gfx,
			 Rectangle box) {
    gfx.drawString(m_Format.toPattern(), box.x, box.y + box.height);
  }

  /**
   * Returns the Java code that generates an object the same as the one being edited.
   *
   * @return the initialization string
   */   
  public String getJavaInitializationString() {
    return ("new SimpleDateFormat(" + m_Format.toPattern() + ")");
  }

  /**
   * Returns the date format string.
   *
   * @return the date format string
   */   
  public String getAsText() {
    return m_Format.toPattern();
  }

  /**
   * Sets the date format string.
   *
   * @param text the date format string
   */   
  public void setAsText(String text) {
    m_Format = new SimpleDateFormat(text);
  }

  /**
   * Some objects can return tags, but a date format cannot.
   *
   * @return null
   */  
  public String[] getTags() {
    return null;
  }

  /**
   * Gets a GUI component with which the user can edit the date format.
   *
   * @return an editor GUI component
   */    
  public Component getCustomEditor() {
    return m_customEditor;
  }

  /**
   * Indicates whether the date format can be edited in a GUI, which it can.
   *
   * @return true
   */     
  public boolean supportsCustomEditor() {
    return true;
  }

  /**
   * Adds an object to the list of those that wish to be informed when the
   * date format changes.
   *
   * @param listener a new listener to add to the list
   */   
  public void addPropertyChangeListener(PropertyChangeListener listener) {
    m_propSupport.addPropertyChangeListener(listener);
  }

  /**
   * Removes an object from the list of those that wish to be informed when the
   * date format changes.
   *
   * @param listener the listener to remove from the list
   */  
  public void removePropertyChangeListener(PropertyChangeListener listener) {
    m_propSupport.removePropertyChangeListener(listener);
  }
}

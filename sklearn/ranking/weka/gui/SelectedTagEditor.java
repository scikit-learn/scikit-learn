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
 *    SelectedTagEditor.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;


import weka.core.Tag;
import weka.core.SelectedTag;

import java.awt.BorderLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyEditor;
import java.beans.PropertyEditorSupport;
import javax.swing.JFrame;
import javax.swing.JLabel;

/** 
 * A PropertyEditor that uses tags, where the tags are obtained from a
 * weka.core.SelectedTag object.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public class SelectedTagEditor extends PropertyEditorSupport {

  /**
   * Returns a description of the property value as java source.
   *
   * @return a value of type 'String'
   */
  public String getJavaInitializationString() {

    SelectedTag s = (SelectedTag)getValue();
    Tag [] tags = s.getTags();
    String result = "new SelectedTag("
      + s.getSelectedTag().getID()
      + ", {\n";
    for (int i = 0; i < tags.length; i++) {
      result += "new Tag(" + tags[i].getID()
	+ ",\"" + tags[i].getReadable()
	+ "\")";
      if (i < tags.length - 1) {
	result += ',';
      }
      result += '\n';
    }
    return result + "})";
  }

  /**
   * Gets the current value as text.
   *
   * @return a value of type 'String'
   */
  public String getAsText() {

    SelectedTag s = (SelectedTag)getValue();
    return s.getSelectedTag().getReadable();
  }

  /**
   * Sets the current property value as text.
   *
   * @param text the text of the selected tag.
   * @exception java.lang.IllegalArgumentException if an error occurs
   */
  public void setAsText(String text)
    {

    SelectedTag s = (SelectedTag)getValue();
    Tag [] tags = s.getTags();
    try {
      for (int i = 0; i < tags.length; i++) {
	if (text.equals(tags[i].getReadable())) {
	  setValue(new SelectedTag(tags[i].getID(), tags));
	  return;
	}
      }
    } catch (Exception ex) {
      throw new java.lang.IllegalArgumentException(text);
    }
  }

  /**
   * Gets the list of tags that can be selected from.
   *
   * @return an array of string tags.
   */
  public String[] getTags() {

    SelectedTag s = (SelectedTag)getValue();
    Tag [] tags = s.getTags();
    String [] result = new String [tags.length];
    for (int i = 0; i < tags.length; i++) {
      result[i] = tags[i].getReadable();
    }
    return result;
  }
  
  /**
   * Tests out the selectedtag editor from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      GenericObjectEditor.registerEditors();
      Tag [] tags =  {
	new Tag(1, "First option"),
	new Tag(2, "Second option"),
	new Tag(3, "Third option"),
	new Tag(4, "Fourth option"),
	new Tag(5, "Fifth option"),
      };
      SelectedTag initial = new SelectedTag(1, tags);
      SelectedTagEditor ce = new SelectedTagEditor();
      ce.setValue(initial);
      PropertyValueSelector ps = new PropertyValueSelector(ce);
      JFrame f = new JFrame(); 
      f.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  System.exit(0);
	}
      });
      f.getContentPane().setLayout(new BorderLayout());
      f.getContentPane().add(ps, BorderLayout.CENTER);
      f.pack();
      f.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}


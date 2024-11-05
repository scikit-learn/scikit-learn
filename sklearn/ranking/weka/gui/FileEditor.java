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
 *    FileEditor.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import java.awt.Container;
import java.awt.Dialog;
import java.awt.FontMetrics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyEditorSupport;
import java.io.File;

import javax.swing.JFileChooser;


/** 
 * A PropertyEditor for File objects that lets the user select a file.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5204 $
 */
public class FileEditor extends PropertyEditorSupport {

  /** The file chooser used for selecting files. */
  protected JFileChooser m_FileChooser;
  
  /**
   * Returns a representation of the current property value as java source.
   *
   * @return a value of type 'String'
   */
  public String getJavaInitializationString() {

    File f = (File) getValue();
    if (f == null) {
      return "null";
    }
    return "new File(\"" + f.getName() + "\")";
  }

  /**
   * Returns true because we do support a custom editor.
   *
   * @return true
   */
  public boolean supportsCustomEditor() {
    return true;
  }
  
  /**
   * Gets the custom editor component.
   *
   * @return a value of type 'java.awt.Component'
   */
  public java.awt.Component getCustomEditor() {

    if (m_FileChooser == null) {
      File currentFile = (File) getValue();
      if (currentFile != null) {
	m_FileChooser 
	  = new JFileChooser();
	m_FileChooser.setSelectedFile(currentFile);
      } else {
	m_FileChooser 
	  = new JFileChooser(new File(System.getProperty("user.dir")));
      }
      m_FileChooser.setApproveButtonText("Select");
      m_FileChooser.setApproveButtonMnemonic('S');
      m_FileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
      m_FileChooser.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  String cmdString = e.getActionCommand();
	  if (cmdString.equals(JFileChooser.APPROVE_SELECTION)) {
	    File newVal = m_FileChooser.getSelectedFile();
	    setValue(newVal);
	  }
	  closeDialog();
	}
      });
    }
    return m_FileChooser;
  }

  /**
   * Returns true since this editor is paintable.
   *
   * @return true.
   */
  public boolean isPaintable() {
    return true;
  }

  /**
   * Paints a representation of the current Object.
   *
   * @param gfx the graphics context to use
   * @param box the area we are allowed to paint into
   */
  public void paintValue(java.awt.Graphics gfx, java.awt.Rectangle box) {

    FontMetrics fm = gfx.getFontMetrics();
    int vpad = (box.height - fm.getHeight()) / 2 ;
    File f = (File) getValue();
    String val = "No file";
    if (f != null) {
      val = f.getName();
    }
    gfx.drawString(val, 2, fm.getHeight() + vpad);
  }  
  
  /**
   * Closes the dialog.
   */
  protected void closeDialog() {
    if (m_FileChooser instanceof Container) {
      Dialog dlg = PropertyDialog.getParentDialog((Container) m_FileChooser);
      if (dlg != null)
	dlg.setVisible(false);
    }
  }
}


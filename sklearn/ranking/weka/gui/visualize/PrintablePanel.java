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
  *    PrintablePanel.java
  *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
  *
  */

package weka.gui.visualize;

import java.util.Hashtable;

import javax.swing.JPanel;

/** 
 * This Panel enables the user to print the panel to various file formats.
 * The Print dialog is accessible via Ctrl-Shft-Left Mouse Click. <p>
 * The individual JComponentWriter-descendants can be accessed by the
 * <code>getWriter(String)</code> method, if the parameters need to be changed.
 *
 * @see #getWriters()
 * @see #getWriter(String)
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.3 $
 */
public class PrintablePanel
  extends JPanel
  implements PrintableHandler {

  /** for serialization */
  private static final long serialVersionUID = 6281532227633417538L;
  
  /** the class responsible for printing */
  protected PrintableComponent m_Printer = null;
  
  /**
   * initializes the panel
   */
  public PrintablePanel() {
    super();
    m_Printer = new PrintableComponent(this);
  }
  
  /**
   * returns a Hashtable with the current available JComponentWriters in the 
   * save dialog. the key of the Hashtable is the description of the writer.
   * 
   * @return all currently available JComponentWriters 
   * @see JComponentWriter#getDescription()
   */
  public Hashtable getWriters() {
    return m_Printer.getWriters();
  }
  
  /**
   * returns the JComponentWriter associated with the given name, is 
   * <code>null</code> if not found
   * 
   * @return the writer associated with the given name
   * @see JComponentWriter#getDescription()
   */
  public JComponentWriter getWriter(String name) {
    return m_Printer.getWriter(name);
  }

  /**
   * sets the title for the save dialog
   */
  public void setSaveDialogTitle(String title) {
    m_Printer.setSaveDialogTitle(title);
  }
  
  /**
   * returns the title for the save dialog
   */
  public String getSaveDialogTitle() {
    return m_Printer.getSaveDialogTitle();
  }
  
  /**
   * sets the scale factor
   * @param x the scale factor for the x-axis 
   * @param y the scale factor for the y-axis 
   */
  public void setScale(double x, double y) {
    m_Printer.setScale(x, y);
  }
  
  /**
   * returns the scale factor for the x-axis
   */
  public double getXScale() {
    return m_Printer.getXScale();
  }
  
  /**
   * returns the scale factor for the y-axis
   */
  public double getYScale() {
    return m_Printer.getYScale();
  }
  
  /**
   * displays a save dialog for saving the panel to a file.  
   */
  public void saveComponent() {
    m_Printer.saveComponent();
  }
}

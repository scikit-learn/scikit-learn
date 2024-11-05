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
 *    ListSelectorDialog.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.regex.Pattern;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;

/** 
 * A dialog to present the user with a list of items, that the user can
 * make a selection from, or cancel the selection.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.9 $
 */
public class ListSelectorDialog
  extends JDialog {

  /** for serialization */
  private static final long serialVersionUID = 906147926840288895L;
  
  /** Click to choose the currently selected property */
  protected JButton m_SelectBut = new JButton("Select");

  /** Click to cancel the property selection */
  protected JButton m_CancelBut = new JButton("Cancel");

  /** Click to enter a regex pattern for selection */
  protected JButton m_PatternBut = new JButton("Pattern");

  /** The list component */
  protected JList m_List;
  
  /** Whether the selection was made or cancelled */
  protected int m_Result;

  /** Signifies an OK property selection */
  public static final int APPROVE_OPTION = 0;

  /** Signifies a cancelled property selection */
  public static final int CANCEL_OPTION = 1;

  /** The current regular expression. */
  protected String m_PatternRegEx = ".*";
  
  /**
   * Create the list selection dialog.
   *
   * @param parentFrame the parent frame of the dialog
   * @param userList the JList component the user will select from
   */
  public ListSelectorDialog(Frame parentFrame, JList userList) {
    
    super(parentFrame, "Select items", true);
    m_List = userList;
    m_CancelBut.setMnemonic('C');
    m_CancelBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	m_Result = CANCEL_OPTION;
	setVisible(false);
      }
    });
    m_SelectBut.setMnemonic('S');
    m_SelectBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	m_Result = APPROVE_OPTION;
	setVisible(false);
      }
    });
    m_PatternBut.setMnemonic('P');
    m_PatternBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        selectPattern();
      }
    });
    
    Container c = getContentPane();
    c.setLayout(new BorderLayout());
    //    setBorder(BorderFactory.createTitledBorder("Select a property"));
    Box b1 = new Box(BoxLayout.X_AXIS);
    b1.add(m_SelectBut);
    b1.add(Box.createHorizontalStrut(10));
    b1.add(m_PatternBut);
    b1.add(Box.createHorizontalStrut(10));
    b1.add(m_CancelBut);
    c.add(b1, BorderLayout.SOUTH);
    c.add(new JScrollPane(m_List), BorderLayout.CENTER);

    getRootPane().setDefaultButton(m_SelectBut);
    
    pack();

    // make sure, it's not bigger than the screen
    Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int width  = getWidth() > screen.getWidth() 
                    ? (int) screen.getWidth() : getWidth();
    int height = getHeight() > screen.getHeight() 
                    ? (int) screen.getHeight() : getHeight();
    setSize(width, height);
  }

  /**
   * Pops up the modal dialog and waits for cancel or a selection.
   *
   * @return either APPROVE_OPTION, or CANCEL_OPTION
   */
  public int showDialog() {

    m_Result = CANCEL_OPTION;
    int [] origSelected = m_List.getSelectedIndices();
    setVisible(true);
    if (m_Result == CANCEL_OPTION) {
      m_List.setSelectedIndices(origSelected);
    }
    return m_Result;
  }

  /**
   * opens a separate dialog for entering a regex pattern for selecting
   * elements from the provided list
   */
  protected void selectPattern() {
    String pattern = JOptionPane.showInputDialog(
                        m_PatternBut.getParent(),
                        "Enter a Perl regular expression ('.*' for all)",
                        m_PatternRegEx);
    if (pattern != null) {
      try {
        Pattern.compile(pattern);
        m_PatternRegEx = pattern;
        m_List.clearSelection();
        for (int i = 0; i < m_List.getModel().getSize(); i++) {
          if (Pattern.matches(
                pattern, m_List.getModel().getElementAt(i).toString()))
            m_List.addSelectionInterval(i, i);
        }
      }
      catch (Exception ex) {
        JOptionPane.showMessageDialog(
          m_PatternBut.getParent(),
          "'" + pattern + "' is not a valid Perl regular expression!\n" 
          + "Error: " + ex, 
          "Error in Pattern...", 
          JOptionPane.ERROR_MESSAGE);
      }
    }
  }
  
  /**
   * Tests out the list selector from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      DefaultListModel lm = new DefaultListModel();      
      lm.addElement("one");
      lm.addElement("two");
      lm.addElement("three");
      lm.addElement("four");
      lm.addElement("five");
      JList jl = new JList(lm);
      final ListSelectorDialog jd = new ListSelectorDialog(null, jl);
      int result = jd.showDialog();
      if (result == ListSelectorDialog.APPROVE_OPTION) {
	System.err.println("Fields Selected");
	int [] selected = jl.getSelectedIndices();
	for (int i = 0; i < selected.length; i++) {
	  System.err.println("" + selected[i]
			     + " " + lm.elementAt(selected[i]));
	}
      } else {
	System.err.println("Cancelled");
      }
      System.exit(0);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

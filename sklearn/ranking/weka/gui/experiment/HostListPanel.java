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
 *    HostListPanel.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.experiment;

import weka.experiment.RemoteExperiment;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

/** 
 * This panel controls setting a list of hosts for a RemoteExperiment to
 * use.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class HostListPanel
  extends JPanel
  implements ActionListener {

  /** for serialization */
  private static final long serialVersionUID = 7182791134585882197L;
  
  /** The remote experiment to set the host list of */
  protected RemoteExperiment m_Exp;

 /** The component displaying the host list */
  protected JList m_List;

  /** Click to remove the selected host from the list */
  protected JButton m_DeleteBut = new JButton("Delete selected");

  /** The field with which to enter host names */
  protected JTextField m_HostField = new JTextField(25);

  /**
   * Creates the host list panel with the given experiment.
   *
   * @param exp a value of type 'Experiment'
   */
  public HostListPanel(RemoteExperiment exp) {
    this();
    setExperiment(exp);
  }

  /**
   * Create the host list panel initially disabled.
   */
  public HostListPanel() {
    m_List = new JList();
    m_List.setModel(new DefaultListModel());
    m_DeleteBut.setEnabled(false);
    m_DeleteBut.addActionListener(this);
    m_HostField.addActionListener(this);
    setLayout(new BorderLayout());
    setBorder(BorderFactory.createTitledBorder("Hosts"));

    JPanel topLab = new JPanel();
    GridBagLayout gb = new GridBagLayout();
    GridBagConstraints constraints = new GridBagConstraints();
    topLab.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));
    //    topLab.setLayout(new GridLayout(1,2,5,5));
    topLab.setLayout(gb);
   
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;
    constraints.insets = new Insets(0,2,0,2);
    topLab.add(m_DeleteBut,constraints);
    constraints.gridx=1;constraints.gridy=0;constraints.weightx=5;
    constraints.gridwidth=1;constraints.gridheight=1;
    topLab.add(m_HostField,constraints);
    
    add(topLab, BorderLayout.NORTH);
    add(new JScrollPane(m_List), BorderLayout.CENTER);
  }

  /**
   * Tells the panel to act on a new experiment.
   *
   * @param exp a value of type 'Experiment'
   */
  public void setExperiment(RemoteExperiment exp) {
    m_Exp = exp;
    m_List.setModel(m_Exp.getRemoteHosts());
    if (((DefaultListModel)m_List.getModel()).size() > 0) {
      m_DeleteBut.setEnabled(true);
    }
  }

  /**
   * Handle actions when text is entered into the host field or the
   * delete button is pressed.
   *
   * @param e a value of type 'ActionEvent'
   */
  public void actionPerformed(ActionEvent e) {
    if (e.getSource() == m_HostField) {
      ((DefaultListModel)m_List.getModel())
	.addElement(m_HostField.getText());
      m_DeleteBut.setEnabled(true);
    } else if (e.getSource() == m_DeleteBut) {
      int [] selected = m_List.getSelectedIndices();
      if (selected != null) {
	for (int i = selected.length - 1; i >= 0; i--) {
	  int current = selected[i];
	  ((DefaultListModel)m_List.getModel()).removeElementAt(current);
	  if (((DefaultListModel)m_List.getModel()).size() > current) {
	    m_List.setSelectedIndex(current);
	  } else {
	    m_List.setSelectedIndex(current - 1);
	  }
	}
      }
      if (((DefaultListModel)m_List.getModel()).size() == 0) {
	m_DeleteBut.setEnabled(false);
      }
    }
  }

  /**
   * Tests out the host list panel from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      final JFrame jf = new JFrame("Host List Editor");
      jf.getContentPane().setLayout(new BorderLayout());
      HostListPanel dp = new HostListPanel();
      jf.getContentPane().add(dp,
			      BorderLayout.CENTER);
      jf.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      /* System.err.println("Short nap");
      Thread.currentThread().sleep(3000);
      System.err.println("Done"); */
      //      dp.setExperiment(new Experiment());
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

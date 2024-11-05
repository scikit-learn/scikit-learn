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
 *    DistributeExperimentPanel.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.experiment;

import weka.experiment.Experiment;
import weka.experiment.RemoteExperiment;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/** 
 * This panel enables an experiment to be distributed to multiple hosts;
 * it also allows remote host names to be specified.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.6 $
 */
public class DistributeExperimentPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = 5206721431754800278L;

  /**
   * The experiment to configure.
   */
  RemoteExperiment m_Exp = null;

  /** Distribute the current experiment to remote hosts */
  protected JCheckBox m_enableDistributedExperiment = 
    new JCheckBox();

  /** Popup the HostListPanel */
  protected JButton m_configureHostNames = new JButton("Hosts");

  /** The host list panel */
  protected HostListPanel m_hostList = new HostListPanel();

  /**
   * Split experiment up by data set.
   */
  protected JRadioButton m_splitByDataSet = new JRadioButton("By data set");

  /**
   * Split experiment up by run number.
   */
  protected JRadioButton m_splitByRun = new JRadioButton("By run");

  /** Handle radio buttons */
  ActionListener m_radioListener = new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	updateRadioLinks();
      }
    };
  
  /**
   * Constructor
   */
  public DistributeExperimentPanel() {
    m_enableDistributedExperiment.setSelected(false);
    m_enableDistributedExperiment.
      setToolTipText("Allow this experiment to be distributed to remote hosts");
    m_enableDistributedExperiment.setEnabled(false);
    m_configureHostNames.setEnabled(false);
    m_configureHostNames.setToolTipText("Edit the list of remote hosts");
    
    m_enableDistributedExperiment.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  m_configureHostNames.setEnabled(m_enableDistributedExperiment.
					  isSelected());
	  m_splitByDataSet.setEnabled(m_enableDistributedExperiment.
					  isSelected());
	  m_splitByRun.setEnabled(m_enableDistributedExperiment.
					  isSelected());
	  
	}
      });

    m_configureHostNames.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  popupHostPanel();
	}
      });

    m_splitByDataSet.setToolTipText("Distribute experiment by data set");
    m_splitByRun.setToolTipText("Distribute experiment by run number");
    m_splitByDataSet.setSelected(true);
    m_splitByDataSet.setEnabled(false);
    m_splitByRun.setEnabled(false);
    m_splitByDataSet.addActionListener(m_radioListener);
    m_splitByRun.addActionListener(m_radioListener);

    ButtonGroup bg = new ButtonGroup();
    bg.add(m_splitByDataSet);
    bg.add(m_splitByRun);

    JPanel rbuts = new JPanel();
    rbuts.setLayout(new GridLayout(1, 2));
    rbuts.add(m_splitByDataSet);
    rbuts.add(m_splitByRun);

    setLayout(new BorderLayout());
    setBorder(BorderFactory.createTitledBorder("Distribute experiment"));
    add(m_enableDistributedExperiment, BorderLayout.WEST);
    add(m_configureHostNames, BorderLayout.CENTER);
    add(rbuts, BorderLayout.SOUTH);
  }

  /**
   * Creates the panel with the supplied initial experiment.
   *
   * @param exp a value of type 'Experiment'
   */
  public DistributeExperimentPanel(Experiment exp) {
    this();
    setExperiment(exp);
  }

  /**
   * Sets the experiment to be configured.
   *
   * @param exp a value of type 'Experiment'
   */
  public void setExperiment(Experiment exp) {
    
    m_enableDistributedExperiment.setEnabled(true);
    m_Exp = null;
    if (exp instanceof RemoteExperiment) {
      m_Exp = (RemoteExperiment)exp;
      m_enableDistributedExperiment.setSelected(true);
      m_configureHostNames.setEnabled(true);
      m_hostList.setExperiment(m_Exp);
      m_splitByDataSet.setEnabled(true);
      m_splitByRun.setEnabled(true);
      m_splitByDataSet.setSelected(m_Exp.getSplitByDataSet());
      m_splitByRun.setSelected(!m_Exp.getSplitByDataSet());
    }
  }

  /**
   * Pop up the host list panel
   */
  private void popupHostPanel() {
    try {
      final JFrame jf = new JFrame("Edit host names");
      
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(m_hostList,
			      BorderLayout.CENTER);
      jf.addWindowListener(new WindowAdapter() {
	  public void windowClosing(WindowEvent e) {
	    jf.dispose();
	  }
	});
      jf.pack();
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }

  /**
   * Returns true if the distribute experiment checkbox is selected
   * @return true if the user has opted for distributing the experiment
   */
  public boolean distributedExperimentSelected() {
    return m_enableDistributedExperiment.isSelected();
  }

  /**
   * Enable objects to listen for changes to the check box
   * @param al an ActionListener
   */
  public void addCheckBoxActionListener(ActionListener al) {
    m_enableDistributedExperiment.addActionListener(al);
  }

  /**
   * Updates the remote experiment when a radio button is clicked
   */
  private void updateRadioLinks() {
    if (m_Exp != null) {
      m_Exp.setSplitByDataSet(m_splitByDataSet.isSelected());
    }
  }

  /**
   * Tests out the panel from the command line.
   *
   * @param args ignored.
   */
  public static void main(String [] args) {
    try {
      final JFrame jf = new JFrame("DistributeExperiment");
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(new DistributeExperimentPanel(new Experiment()),
			      BorderLayout.CENTER);
      jf.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

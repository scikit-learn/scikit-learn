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
 *    SetupModePanel.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.experiment;

import weka.experiment.Experiment;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeListener;

import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/** 
 * This panel switches between simple and advanced experiment setup panels.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class SetupModePanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -3758035565520727822L;

  /** The button for choosing simple setup mode */
  protected JRadioButton m_SimpleSetupRBut = 
    new JRadioButton("Simple");

  /** The button for choosing advanced setup mode */
  protected JRadioButton m_AdvancedSetupRBut = 
    new JRadioButton("Advanced");  

  /** The simple setup panel */
  protected SimpleSetupPanel m_simplePanel = new SimpleSetupPanel();

  /** The advanced setup panel */
  protected SetupPanel m_advancedPanel = new SetupPanel();

  /**
   * Creates the setup panel with no initial experiment.
   */
  public SetupModePanel() {

    m_simplePanel.setModePanel(this);

    m_SimpleSetupRBut.setMnemonic('S');
    m_SimpleSetupRBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  switchToSimple(null);
	}
      });

    m_AdvancedSetupRBut.setMnemonic('A');
    m_AdvancedSetupRBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  switchToAdvanced(null);
	}
      });

    ButtonGroup modeBG = new ButtonGroup();
    modeBG.add(m_SimpleSetupRBut);
    modeBG.add(m_AdvancedSetupRBut);
    m_SimpleSetupRBut.setSelected(true);

    JPanel modeButtons = new JPanel();
    modeButtons.setLayout(new GridLayout(1,0));
    modeButtons.add(m_SimpleSetupRBut);
    modeButtons.add(m_AdvancedSetupRBut);

    JPanel switchPanel = new JPanel();
    switchPanel.setLayout(new GridLayout(1,0));
    switchPanel.add(new JLabel("Experiment Configuration Mode:"));
    switchPanel.add(modeButtons);

    setLayout(new BorderLayout());
    add(switchPanel, BorderLayout.NORTH);
    add(m_simplePanel, BorderLayout.CENTER);
  }

  /**
   * Switches to the advanced setup mode.
   *
   * @param exp the experiment to configure
   */
  public void switchToAdvanced(Experiment exp) {
 
    if (exp == null) {
      exp = m_simplePanel.getExperiment();
    }
    if (exp != null) {
      m_AdvancedSetupRBut.setSelected(true);
      m_advancedPanel.setExperiment(exp);
    }
    remove(m_simplePanel);
    m_simplePanel.removeNotesFrame();
    add(m_advancedPanel, BorderLayout.CENTER);
    validate();
    repaint();
  }
  
  /**
   * Switches to the simple setup mode only if allowed to.
   *
   * @param exp the experiment to configure
   */
  public void switchToSimple(Experiment exp) {
    
    if (exp == null) {
      exp = m_advancedPanel.getExperiment();
    }
    if (exp != null && !m_simplePanel.setExperiment(exp)) {
      m_AdvancedSetupRBut.setSelected(true);
      switchToAdvanced(exp);
    } else {
      remove(m_advancedPanel);
      m_advancedPanel.removeNotesFrame();
      add(m_simplePanel, BorderLayout.CENTER);
      validate();
      repaint();
    }
  }

  /**
   * Adds a PropertyChangeListener who will be notified of value changes.
   *
   * @param l a value of type 'PropertyChangeListener'
   */
  public void addPropertyChangeListener(PropertyChangeListener l) {

    m_simplePanel.addPropertyChangeListener(l);
    m_advancedPanel.addPropertyChangeListener(l);
  }

  /**
   * Gets the currently configured experiment.
   *
   * @return the currently configured experiment.
   */
  public Experiment getExperiment() {

    if (m_SimpleSetupRBut.isSelected()) return m_simplePanel.getExperiment();
    else return m_advancedPanel.getExperiment();
  }
}

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
 *    RunNumberPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.experiment;

import weka.experiment.Experiment;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;



/** 
 * This panel controls configuration of lower and upper run numbers
 * in an experiment.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.6 $
 */
public class RunNumberPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -1644336658426067852L;

  /** Configures the lower run number */
  protected JTextField m_LowerText = new JTextField("1");

  /** Configures the upper run number */
  protected JTextField m_UpperText = new JTextField("10");

  /** The experiment being configured */
  protected Experiment m_Exp;
  
  /**
   * Creates the panel with no initial experiment.
   */
  public RunNumberPanel() {
    
    // Updates occur to the values in exp whenever enter is pressed
    // or the component loses focus
    m_LowerText.addKeyListener(new KeyAdapter() {
      public void keyReleased(KeyEvent e) {
	m_Exp.setRunLower(getLower());
      }
    });
    m_LowerText.addFocusListener(new FocusAdapter() {
      public void focusLost(FocusEvent e) {
	m_Exp.setRunLower(getLower());
      }
    });
    m_UpperText.addKeyListener(new KeyAdapter() {
      public void keyReleased(KeyEvent e) {
	m_Exp.setRunUpper(getUpper());
      }
    });
    m_UpperText.addFocusListener(new FocusAdapter() {
      public void focusLost(FocusEvent e) {
	m_Exp.setRunUpper(getUpper());
      }
    });
    m_LowerText.setEnabled(false);
    m_UpperText.setEnabled(false);

    // Set the GUI layout
    setLayout(new GridLayout(1,2));
    setBorder(BorderFactory.createTitledBorder("Runs"));
    Box b1 = new Box(BoxLayout.X_AXIS);
    b1.add(Box.createHorizontalStrut(10));
    b1.add(new JLabel("From:", SwingConstants.RIGHT));
    b1.add(Box.createHorizontalStrut(5));
    b1.add(m_LowerText);
    add(b1);
    Box b2 = new Box(BoxLayout.X_AXIS);
    b2.add(Box.createHorizontalStrut(10));
    b2.add(new JLabel("To:", SwingConstants.RIGHT));
    b2.add(Box.createHorizontalStrut(5));
    b2.add(m_UpperText);
    add(b2);
  }

  /**
   * Creates the panel with the supplied initial experiment.
   *
   * @param exp a value of type 'Experiment'
   */
  public RunNumberPanel(Experiment exp) {

    this();
    setExperiment(exp);
  }

  /**
   * Sets the experiment to be configured.
   *
   * @param exp a value of type 'Experiment'
   */
  public void setExperiment(Experiment exp) {
    
    m_Exp = exp;
    m_LowerText.setText("" + m_Exp.getRunLower());
    m_UpperText.setText("" + m_Exp.getRunUpper());
    m_LowerText.setEnabled(true);
    m_UpperText.setEnabled(true);
  }
  
  /**
   * Gets the current lower run number.
   *
   * @return the lower run number.
   */
  public int getLower() {

    int result = 1;
    try {
      result = Integer.parseInt(m_LowerText.getText());
    } catch (Exception ex) {
    }
    return Math.max(1, result);
  }
  
  /**
   * Gets the current upper run number.
   *
   * @return the upper run number.
   */
  public int getUpper() {

    int result = 1;
    try {
      result = Integer.parseInt(m_UpperText.getText());
    } catch (Exception ex) {
    }
    return Math.max(1, result);
  }
  
  /**
   * Tests out the panel from the command line.
   *
   * @param args ignored.
   */
  public static void main(String [] args) {

    try {
      final JFrame jf = new JFrame("Dataset List Editor");
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(new RunNumberPanel(new Experiment()),
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

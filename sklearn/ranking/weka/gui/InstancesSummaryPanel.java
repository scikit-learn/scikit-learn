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
 *    InstancesSummaryPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import weka.core.Instances;
import weka.core.Utils;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

/** 
 * This panel just displays relation name, number of instances, and number of
 * attributes.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5297 $
 */
public class InstancesSummaryPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -5243579535296681063L;

  /** Message shown when no instances have been loaded */
  protected static final String NO_SOURCE = "None";

  /** Displays the name of the relation */
  protected JLabel m_RelationNameLab = new JLabel(NO_SOURCE);
  
  /** Displays the number of instances */
  protected JLabel m_NumInstancesLab = new JLabel(NO_SOURCE);
  
  /** Displays the number of attributes */
  protected JLabel m_NumAttributesLab = new JLabel(NO_SOURCE);
  
  /** Displays the sum of instance weights */
  protected JLabel m_sumOfWeightsLab = new JLabel(NO_SOURCE);
    
  /** The instances we're playing with */
  protected Instances m_Instances;
  
  /** 
   * Whether to display 0 or ? for the number of instances
   * in cases where a dataset has only structure. Depending
   * on where this panel is used from, the user may have
   * loaded a dataset with no instances or a Loader that
   * can read incrementally may be being used (in which case
   * we don't know how many instances are in the dataset... 
   * yet).
   */
  protected boolean m_showZeroInstancesAsUnknown = false;

  /**
   * Creates the instances panel with no initial instances.
   */
  public InstancesSummaryPanel() {

    GridBagLayout gbLayout = new GridBagLayout();
    setLayout(gbLayout);
    JLabel lab = new JLabel("Relation:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
    GridBagConstraints gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0;     gbConstraints.gridx = 0;
    gbLayout.setConstraints(lab, gbConstraints);
    add(lab);
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.WEST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0;     gbConstraints.gridx = 1;
    gbConstraints.weightx = 100; gbConstraints.gridwidth = 3;
    gbLayout.setConstraints(m_RelationNameLab, gbConstraints);
    add(m_RelationNameLab);
    m_RelationNameLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 10));
    
    lab = new JLabel("Instances:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1;     gbConstraints.gridx = 0;
    gbLayout.setConstraints(lab, gbConstraints);
    add(lab);
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.WEST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1;     gbConstraints.gridx = 1;
    gbConstraints.weightx = 100;
    gbLayout.setConstraints(m_NumInstancesLab, gbConstraints);
    add(m_NumInstancesLab);
    m_NumInstancesLab.setBorder(BorderFactory.createEmptyBorder(0, 5, 0, 10));

    lab = new JLabel("Attributes:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0;     gbConstraints.gridx = 2;
    gbLayout.setConstraints(lab, gbConstraints);
    add(lab);
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.WEST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 0;     gbConstraints.gridx = 3;
    gbConstraints.weightx = 100;
    gbLayout.setConstraints(m_NumAttributesLab, gbConstraints);
    add(m_NumAttributesLab);
    m_NumAttributesLab.setBorder(BorderFactory.createEmptyBorder(0, 5,
								 0, 10));
    
    lab = new JLabel("Sum of weights:", SwingConstants.RIGHT);
    lab.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.EAST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1;     gbConstraints.gridx = 2;
    gbLayout.setConstraints(lab, gbConstraints);
    add(lab);
    gbConstraints = new GridBagConstraints();
    gbConstraints.anchor = GridBagConstraints.WEST;
    gbConstraints.fill = GridBagConstraints.HORIZONTAL;
    gbConstraints.gridy = 1;     gbConstraints.gridx = 3;
    gbConstraints.weightx = 100;
    gbLayout.setConstraints(m_sumOfWeightsLab, gbConstraints);
    add(m_sumOfWeightsLab);
    m_sumOfWeightsLab.setBorder(BorderFactory.createEmptyBorder(0, 5,
                                                                 0, 10));
    
  }
  
  public void setShowZeroInstancesAsUnknown(boolean zeroAsUnknown) {
    m_showZeroInstancesAsUnknown = zeroAsUnknown;
  }
  
  public boolean getShowZeroInstancesAsUnknown() {
    return m_showZeroInstancesAsUnknown;
  }

  /**
   * Tells the panel to use a new set of instances.
   *
   * @param inst a set of Instances
   */
  public void setInstances(Instances inst) {
    m_Instances = inst;
    m_RelationNameLab.setText(m_Instances.relationName());
    m_NumInstancesLab.setText("" + 
        ((m_showZeroInstancesAsUnknown && m_Instances.numInstances() == 0) 
            ? "?" 
            : "" + m_Instances.numInstances()));
    m_NumAttributesLab.setText("" + m_Instances.numAttributes());
    m_sumOfWeightsLab.setText("" + 
        ((m_showZeroInstancesAsUnknown && m_Instances.numInstances() == 0) 
            ? "?" 
            : "" + Utils.doubleToString(m_Instances.sumOfWeights(), 3)));
  }

  /**
   * Tests out the instance summary panel from the command line.
   *
   * @param args optional name of dataset to load
   */
  public static void main(String [] args) {

    try {
      final javax.swing.JFrame jf = new javax.swing.JFrame("Instances Panel");
      jf.getContentPane().setLayout(new BorderLayout());
      final InstancesSummaryPanel p = new InstancesSummaryPanel();
      p.setBorder(BorderFactory.createTitledBorder("Relation"));
      jf.getContentPane().add(p, BorderLayout.CENTER);
      jf.addWindowListener(new java.awt.event.WindowAdapter() {
	public void windowClosing(java.awt.event.WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      if (args.length == 1) {
	java.io.Reader r = new java.io.BufferedReader(
			   new java.io.FileReader(args[0]));
	Instances i = new Instances(r);
	p.setInstances(i);
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
  
}

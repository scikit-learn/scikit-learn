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
 *    OutputFormatDialog.java
 *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.experiment;

import weka.experiment.ResultMatrix;
import weka.experiment.ResultMatrixPlainText;
import weka.gui.GenericObjectEditor;
import weka.gui.PropertyPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/** 
 * A dialog for setting various output format parameters.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5346 $
 */
public class OutputFormatDialog
  extends JDialog {

  /** for serialization. */
  private static final long serialVersionUID = 2169792738187807378L;

  /** Signifies an OK property selection. */
  public static final int APPROVE_OPTION = 0;

  /** Signifies a cancelled property selection. */
  public static final int CANCEL_OPTION = 1;

  /** the result of the user's action, either OK or CANCEL. */
  protected int m_Result;
  
  /** the different classes for outputting the comparison tables. */
  protected static Vector<Class> m_OutputFormatClasses;
  
  /** the different names of matrices for outputting the comparison tables. */
  protected static Vector<String> m_OutputFormatNames;
  
  /** Lets the user configure the result matrix. */
  protected GenericObjectEditor m_ResultMatrixEditor;

  /** the panel for the GOE. */
  protected PropertyPanel m_ResultMatrixPanel;

  /** the label for the GOE. */
  protected JLabel m_ResultMatrixLabel;
  
  /** the current result matrix. */
  protected ResultMatrix m_ResultMatrix;

  /** lets the user choose the format for the output. */
  protected JComboBox m_OutputFormatComboBox;

  /** the label for the format. */
  protected JLabel m_OutputFormatLabel;

  /** the spinner to choose the precision for the mean from. */
  protected JSpinner m_MeanPrecSpinner;

  /** the label for the mean precision. */
  protected JLabel m_MeanPrecLabel;

  /** the spinner to choose the precision for the std. deviation from */
  protected JSpinner m_StdDevPrecSpinner;

  /** the label for the std dev precision. */
  protected JLabel m_StdDevPrecLabel;

  /** the checkbox for outputting the average. */
  protected JCheckBox m_ShowAverageCheckBox;

  /** the label for showing the average. */
  protected JLabel m_ShowAverageLabel;

  /** the checkbox for the removing of filter classnames. */
  protected JCheckBox m_RemoveFilterNameCheckBox;

  /** the label for the removing the filter classnames. */
  protected JLabel m_RemoveFilterNameLabel;
  
  /** Click to activate the current set parameters. */
  protected JButton m_OkButton;

  /** Click to cancel the dialog. */
  protected JButton m_CancelButton;

  /** whether to ignore updates in the GUI momentarily. */
  protected boolean m_IgnoreChanges;
  
  /**
   * initializes the dialog with the given parent frame.
   * 
   * @param parent the parent of this dialog
   */
  public OutputFormatDialog(Frame parent) {
    super(parent, "Output Format...", true);

    m_IgnoreChanges = true;
    
    initialize();
    initGUI();

    m_IgnoreChanges = false;
  }
  
  /**
   * initializes the member variables.
   */
  protected void initialize() {
    Vector 		classes;
    int			i;
    Class 		cls;
    ResultMatrix 	matrix;
    
    m_Result = CANCEL_OPTION;

    if (m_OutputFormatClasses == null) {
      classes = GenericObjectEditor.getClassnames(ResultMatrix.class.getName());

      // set names and classes
      m_OutputFormatClasses = new Vector<Class>();
      m_OutputFormatNames   = new Vector<String>();
      for (i = 0; i < classes.size(); i++) {
        try {
          cls    = Class.forName(classes.get(i).toString());
          matrix = (ResultMatrix) cls.newInstance();
          m_OutputFormatClasses.add(cls);
          m_OutputFormatNames.add(matrix.getDisplayName());
        }
        catch (Exception e) {
          e.printStackTrace();
        }
      }
    }
  }
  
  /**
   * performs the creation of the dialog and all its components.
   */
  protected void initGUI() {
    JPanel              panel;
    SpinnerNumberModel  model;
    JPanel		panel2;
    
    getContentPane().setLayout(new BorderLayout());
    
    panel = new JPanel(new GridLayout(6, 1));
    panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    getContentPane().add(panel, BorderLayout.CENTER);
    
    // mean precision
    m_MeanPrecSpinner = new JSpinner();
    m_MeanPrecSpinner.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
        getData();
      }
    });
    model = (SpinnerNumberModel) m_MeanPrecSpinner.getModel();
    model.setMaximum(new Integer(20));
    model.setMinimum(new Integer(0));
    m_MeanPrecLabel = new JLabel("Mean Precision");
    m_MeanPrecLabel.setDisplayedMnemonic('M');
    m_MeanPrecLabel.setLabelFor(m_MeanPrecSpinner);
    panel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panel2.add(m_MeanPrecLabel);
    panel2.add(m_MeanPrecSpinner);
    panel.add(panel2);
    
    // stddev precision
    m_StdDevPrecSpinner = new JSpinner();
    m_StdDevPrecSpinner.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
        getData();
      }
    });
    model = (SpinnerNumberModel) m_StdDevPrecSpinner.getModel();
    model.setMaximum(new Integer(20));
    model.setMinimum(new Integer(0));
    m_StdDevPrecLabel = new JLabel("StdDev. Precision");
    m_StdDevPrecLabel.setDisplayedMnemonic('S');
    m_StdDevPrecLabel.setLabelFor(m_StdDevPrecSpinner);
    panel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panel2.add(m_StdDevPrecLabel);
    panel2.add(m_StdDevPrecSpinner);
    panel.add(panel2);
    
    // Format
    m_OutputFormatComboBox = new JComboBox(m_OutputFormatNames);
    m_OutputFormatComboBox.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	getData();
      }
    });
    m_OutputFormatLabel = new JLabel("Output Format");
    m_OutputFormatLabel.setDisplayedMnemonic('F');
    m_OutputFormatLabel.setLabelFor(m_OutputFormatComboBox);
    panel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panel2.add(m_OutputFormatLabel);
    panel2.add(m_OutputFormatComboBox);
    panel.add(panel2);

    // Average
    m_ShowAverageCheckBox = new JCheckBox("");
    m_ShowAverageCheckBox.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
	getData();
      }
    });
    m_ShowAverageLabel = new JLabel("Show Average");
    m_ShowAverageLabel.setDisplayedMnemonic('A');
    m_ShowAverageLabel.setLabelFor(m_ShowAverageCheckBox);
    panel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panel2.add(m_ShowAverageLabel);
    panel2.add(m_ShowAverageCheckBox);
    panel.add(panel2);

    // Remove filter classname
    m_RemoveFilterNameCheckBox = new JCheckBox("");
    m_RemoveFilterNameCheckBox.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
	getData();
      }
    });
    m_RemoveFilterNameLabel = new JLabel("Remove filter classnames");
    m_RemoveFilterNameLabel.setDisplayedMnemonic('R');
    m_RemoveFilterNameLabel.setLabelFor(m_RemoveFilterNameCheckBox);
    panel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panel2.add(m_RemoveFilterNameLabel);
    panel2.add(m_RemoveFilterNameCheckBox);
    panel.add(panel2);

    // Advanced setup
    m_ResultMatrix       = ExperimenterDefaults.getOutputFormat();
    m_ResultMatrixEditor = new GenericObjectEditor(true);
    m_ResultMatrixEditor.setClassType(ResultMatrix.class);
    m_ResultMatrixEditor.setValue(m_ResultMatrix);
    m_ResultMatrixEditor.addPropertyChangeListener(new PropertyChangeListener() {
	public void propertyChange(PropertyChangeEvent e) {
	  // user selected different class?
	  if (!m_ResultMatrix.getClass().equals(m_ResultMatrixEditor.getValue().getClass())) {
	    // if it's the preferred class, then automaticallly use the Experimenter defaults
	    if (m_ResultMatrixEditor.getValue().getClass().equals(ExperimenterDefaults.getOutputFormat().getClass())) {
	      m_ResultMatrix = ExperimenterDefaults.getOutputFormat();
	      m_ResultMatrixEditor.setValue(ExperimenterDefaults.getOutputFormat());
	    }
	    else {
	      m_ResultMatrix = (ResultMatrix) m_ResultMatrixEditor.getValue();
	    }
	    setData();
	  }
	  repaint();
	}
      });
    ((GenericObjectEditor.GOEPanel) m_ResultMatrixEditor.getCustomEditor()).addOkListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
          m_ResultMatrix = (ResultMatrix) m_ResultMatrixEditor.getValue();
          setData();
	}
      });
    m_ResultMatrixPanel = new PropertyPanel(m_ResultMatrixEditor, false);
    m_ResultMatrixLabel = new JLabel("Advanced setup");
    panel2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
    panel2.add(m_ResultMatrixLabel);
    panel2.add(m_ResultMatrixPanel);
    panel.add(panel2);
    
    // Buttons
    panel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
    getContentPane().add(panel, BorderLayout.SOUTH);
    m_CancelButton = new JButton("Cancel");
    m_CancelButton.setMnemonic('C');
    m_CancelButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        m_Result = CANCEL_OPTION;
        setVisible(false);
      }
    });
    m_OkButton = new JButton("OK");
    m_OkButton.setMnemonic('O');
    m_OkButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        getData();
        m_Result = APPROVE_OPTION;
        setVisible(false);
      }
    });
    panel.add(m_OkButton);
    panel.add(m_CancelButton);

    // default button
    getRootPane().setDefaultButton(m_OkButton);
    
    // initial layout (to get widths and heights)
    pack();
    
    // adjust dimensions
    m_MeanPrecLabel.setPreferredSize(new Dimension(m_RemoveFilterNameLabel.getWidth(), m_MeanPrecLabel.getHeight()));
    m_MeanPrecSpinner.setPreferredSize(new Dimension(m_MeanPrecSpinner.getWidth() * 3, m_MeanPrecSpinner.getHeight()));
    m_StdDevPrecLabel.setPreferredSize(new Dimension(m_RemoveFilterNameLabel.getWidth(), m_StdDevPrecLabel.getHeight()));
    m_StdDevPrecSpinner.setPreferredSize(new Dimension(m_StdDevPrecSpinner.getWidth() * 3, m_StdDevPrecSpinner.getHeight()));
    m_OutputFormatLabel.setPreferredSize(new Dimension(m_RemoveFilterNameLabel.getWidth(), m_OutputFormatLabel.getHeight()));
    m_ShowAverageLabel.setPreferredSize(new Dimension(m_RemoveFilterNameLabel.getWidth(), m_ShowAverageLabel.getHeight()));
    m_ResultMatrixLabel.setPreferredSize(new Dimension(m_RemoveFilterNameLabel.getWidth(), m_ResultMatrixLabel.getHeight()));
    m_ResultMatrixPanel.setPreferredSize(new Dimension((int) (m_ResultMatrixPanel.getWidth() * 1.5), m_ResultMatrixPanel.getHeight()));
    
    // final layout
    pack();
  }
  
  /**
   * initializes the GUI components with the data.
   */
  private void setData() {
    m_IgnoreChanges = true;
    
    // Precision
    m_MeanPrecSpinner.setValue(m_ResultMatrix.getMeanPrec());
    m_StdDevPrecSpinner.setValue(m_ResultMatrix.getStdDevPrec());
    
    // format
    for (int i = 0; i < m_OutputFormatClasses.size(); i++) {
      if (m_OutputFormatClasses.get(i).equals(m_ResultMatrix.getClass())) {
	m_OutputFormatComboBox.setSelectedItem(m_OutputFormatNames.get(i));
        break;
      }
    }

    // average
    m_ShowAverageCheckBox.setSelected(m_ResultMatrix.getShowAverage());

    // filter names
    m_RemoveFilterNameCheckBox.setSelected(m_ResultMatrix.getRemoveFilterName());

    // GOE
    m_ResultMatrixEditor.setValue(m_ResultMatrix);
    
    m_IgnoreChanges = false;
  }    
  
  /**
   *  gets the data from GUI components.
   */
  private void getData() {
    if (m_IgnoreChanges)
      return;
    
    // format
    try {
      if (!m_ResultMatrix.getClass().equals(m_OutputFormatClasses.get(m_OutputFormatComboBox.getSelectedIndex()))) {
	if (m_OutputFormatClasses.get(m_OutputFormatComboBox.getSelectedIndex()).equals(ExperimenterDefaults.getOutputFormat().getClass()))
	  m_ResultMatrix = ExperimenterDefaults.getOutputFormat();
	else
	  m_ResultMatrix = (ResultMatrix) ((Class) m_OutputFormatClasses.get(m_OutputFormatComboBox.getSelectedIndex())).newInstance();
      }
    }
    catch (Exception e) {
      e.printStackTrace();
      m_ResultMatrix = new ResultMatrixPlainText();
    }
    
    // Precision
    m_ResultMatrix.setMeanPrec(Integer.parseInt(m_MeanPrecSpinner.getValue().toString()));
    m_ResultMatrix.setStdDevPrec(Integer.parseInt(m_StdDevPrecSpinner.getValue().toString()));

    // average
    m_ResultMatrix.setShowAverage(m_ShowAverageCheckBox.isSelected());

    // filter names
    m_ResultMatrix.setRemoveFilterName(m_RemoveFilterNameCheckBox.isSelected());
    
    // update GOE
    m_ResultMatrixEditor.setValue(m_ResultMatrix);
  }

  /**
   * Sets the matrix to use as initial selected output format.
   * 
   * @param matrix the matrix to use as initial selected output format
   */
  public void setResultMatrix(ResultMatrix matrix) {
    m_ResultMatrix = matrix;
    setData();
  }

  /**
   * Gets the currently selected output format result matrix.
   * 
   * @return the currently selected matrix to use as output
   */
  public ResultMatrix getResultMatrix() {
    return m_ResultMatrix;
  }

  /**
   * sets the class of the chosen result matrix.
   */
  protected void setFormat() {
    for (int i = 0; i < m_OutputFormatClasses.size(); i++) {
      if (m_OutputFormatNames.get(i).equals(
            m_OutputFormatComboBox.getItemAt(i).toString())) {
        m_OutputFormatComboBox.setSelectedIndex(i);
        break;
      }
    }
  }
  
  /**
   * the result from the last display of the dialog, the same is returned
   * from <code>showDialog</code>.
   * 
   * @return the result from the last display of the dialog; 
   *         either APPROVE_OPTION, or CANCEL_OPTION
   * @see #showDialog()
   */
  public int getResult() {
    return m_Result;
  }

  /**
   * Pops up the modal dialog and waits for cancel or a selection.
   *
   * @return either APPROVE_OPTION, or CANCEL_OPTION
   */
  public int showDialog() {
    m_Result = CANCEL_OPTION;
    setData();
    setVisible(true);
    return m_Result;
  }

  /**
   * for testing only.
   * 
   * @param args ignored
   */
  public static void main(String[] args) {
    OutputFormatDialog      dialog;
    
    dialog = new OutputFormatDialog(null);
    if (dialog.showDialog() == APPROVE_OPTION)
      System.out.println("Accepted");
    else
      System.out.println("Aborted");
  }
}

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
 *    SetupPanel.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.experiment;

import weka.core.Utils;
import weka.core.xml.KOML;
import weka.experiment.Experiment;
import weka.experiment.PropertyNode;
import weka.experiment.RemoteExperiment;
import weka.experiment.ResultListener;
import weka.experiment.ResultProducer;
import weka.gui.ExtensionFileFilter;
import weka.gui.GenericObjectEditor;
import weka.gui.PropertyPanel;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.filechooser.FileFilter;

/** 
 * This panel controls the configuration of an experiment.
 * <p>
 * If <a href="http://koala.ilog.fr/XML/serialization/" target="_blank">KOML</a>
 * is in the classpath the experiments can also be saved to XML instead of a
 * binary format.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @author FracPete (fracpete at waikato dot ac dot nz) 
 * @version $Revision: 5841 $
 */
public class SetupPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = 6552671886903170033L;

  /** The experiment being configured */
  protected Experiment m_Exp;

  /** Click to load an experiment */
  protected JButton m_OpenBut = new JButton("Open...");

  /** Click to save an experiment */
  protected JButton m_SaveBut = new JButton("Save...");

  /** Click to create a new experiment with default settings */
  protected JButton m_NewBut = new JButton("New");

  /** A filter to ensure only experiment files get shown in the chooser */
  protected FileFilter m_ExpFilter = 
    new ExtensionFileFilter(Experiment.FILE_EXTENSION, 
                            "Experiment configuration files (*" + Experiment.FILE_EXTENSION + ")");

  /** A filter to ensure only experiment (in KOML format) files get shown in the chooser */
  protected FileFilter m_KOMLFilter = 
    new ExtensionFileFilter(KOML.FILE_EXTENSION, 
                            "Experiment configuration files (*" + KOML.FILE_EXTENSION + ")");

  /** A filter to ensure only experiment (in XML format) files get shown in the chooser */
  protected FileFilter m_XMLFilter = 
    new ExtensionFileFilter(".xml", 
                            "Experiment configuration files (*.xml)");

  /** The file chooser for selecting experiments */
  protected JFileChooser m_FileChooser = new JFileChooser(new File(System.getProperty("user.dir")));

  /** The ResultProducer editor */
  protected GenericObjectEditor m_RPEditor = new GenericObjectEditor();

  /** The panel to contain the ResultProducer editor */
  protected PropertyPanel m_RPEditorPanel = new PropertyPanel(m_RPEditor);

  /** The ResultListener editor */
  protected GenericObjectEditor m_RLEditor = new GenericObjectEditor();

  /** The panel to contain the ResultListener editor */
  protected PropertyPanel m_RLEditorPanel = new PropertyPanel(m_RLEditor);

  /** The panel that configures iteration on custom resultproducer property */
  protected GeneratorPropertyIteratorPanel m_GeneratorPropertyPanel
    = new GeneratorPropertyIteratorPanel();

  /** The panel for configuring run numbers */
  protected RunNumberPanel m_RunNumberPanel = new RunNumberPanel();

  /** The panel for enabling a distributed experiment */
  protected DistributeExperimentPanel m_DistributeExperimentPanel = 
    new DistributeExperimentPanel();

  /** The panel for configuring selected datasets */
  protected DatasetListPanel m_DatasetListPanel = new DatasetListPanel();

  /** A button for bringing up the notes */
  protected JButton m_NotesButton =  new JButton("Notes");

  /** Frame for the notes */
  protected JFrame m_NotesFrame = new JFrame("Notes");

  /** Area for user notes Default of 10 rows */
  protected JTextArea m_NotesText = new JTextArea(null, 10, 0);

  /**
   * Manages sending notifications to people when we change the experiment,
   * at this stage, only the resultlistener so the resultpanel can update.
   */
  protected PropertyChangeSupport m_Support = new PropertyChangeSupport(this);

  /** Click to advacne data set before custom generator */
  protected JRadioButton m_advanceDataSetFirst = 
    new JRadioButton("Data sets first");

  /** Click to advance custom generator before data set */
  protected JRadioButton m_advanceIteratorFirst = 
    new JRadioButton("Custom generator first");

  /** Handle radio buttons */
  ActionListener m_RadioListener = new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	updateRadioLinks();
      }
    };
  
  // Registers the appropriate property editors
  static {
    GenericObjectEditor.registerEditors();
  }
  
  /**
   * Creates the setup panel with the supplied initial experiment.
   *
   * @param exp a value of type 'Experiment'
   */
  public SetupPanel(Experiment exp) {

    this();
    setExperiment(exp);
  }
  
  /**
   * Creates the setup panel with no initial experiment.
   */
  public SetupPanel() {

    m_DistributeExperimentPanel.
      addCheckBoxActionListener(new ActionListener() {
	  public void actionPerformed(ActionEvent e) {
	    if (m_DistributeExperimentPanel.distributedExperimentSelected()) {
	      if (!(m_Exp instanceof RemoteExperiment)) {
		try {
		  RemoteExperiment re = new RemoteExperiment(m_Exp);
		  setExperiment(re);
		} catch (Exception ex) {
		  ex.printStackTrace();
		}
	      }
	    } else {
	      if (m_Exp instanceof RemoteExperiment) {
		setExperiment(((RemoteExperiment)m_Exp).getBaseExperiment());
	      }
	    }
	  }
	});	      

    m_NewBut.setMnemonic('N');
    m_NewBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	setExperiment(new Experiment());
      }
    });
    m_SaveBut.setMnemonic('S');
    m_SaveBut.setEnabled(false);
    m_SaveBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	saveExperiment();
      }
    });
    m_OpenBut.setMnemonic('O');
    m_OpenBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	openExperiment();
      }
    });
    m_FileChooser.addChoosableFileFilter(m_ExpFilter);
    if (KOML.isPresent())
       m_FileChooser.addChoosableFileFilter(m_KOMLFilter);
    m_FileChooser.addChoosableFileFilter(m_XMLFilter);
    m_FileChooser.setFileFilter(m_ExpFilter);
    m_FileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
    
    m_GeneratorPropertyPanel.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  updateRadioLinks();
	}
      });

    m_RPEditor.setClassType(ResultProducer.class);
    m_RPEditor.setEnabled(false);
    m_RPEditor.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
	m_Exp.setResultProducer((ResultProducer) m_RPEditor.getValue());
	m_Exp.setUsePropertyIterator(false);
	m_Exp.setPropertyArray(null);
	m_Exp.setPropertyPath(null);
	m_GeneratorPropertyPanel.setExperiment(m_Exp);
	repaint();
      }
    });

    m_RLEditor.setClassType(ResultListener.class);
    m_RLEditor.setEnabled(false);
    m_RLEditor.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
	m_Exp.setResultListener((ResultListener) m_RLEditor.getValue());
	m_Support.firePropertyChange("", null, null);
	repaint();
      }
    });

    m_NotesFrame.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  m_NotesButton.setEnabled(true);
	}
      });
    m_NotesFrame.getContentPane().add(new JScrollPane(m_NotesText));
    m_NotesFrame.setSize(600, 400);

    m_NotesButton.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  m_NotesButton.setEnabled(false);
	  m_NotesFrame.setVisible(true);
	}
      });
    m_NotesButton.setEnabled(false);

    m_NotesText.setEditable(true);
    //m_NotesText.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    m_NotesText.addKeyListener(new KeyAdapter() {
	public void keyReleased(KeyEvent e) {
	  m_Exp.setNotes(m_NotesText.getText());
	}
      });
    m_NotesText.addFocusListener(new FocusAdapter() {
	public void focusLost(FocusEvent e) {
	  m_Exp.setNotes(m_NotesText.getText());
	}
      });

    // Set up the GUI layout
    JPanel buttons = new JPanel();
    GridBagLayout gb = new GridBagLayout();
    GridBagConstraints constraints = new GridBagConstraints();
    buttons.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));
    //    buttons.setLayout(new GridLayout(1, 3, 5, 5));
    buttons.setLayout(gb);
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;
    constraints.insets = new Insets(0,2,0,2);
    buttons.add(m_OpenBut,constraints);
    constraints.gridx=1;constraints.gridy=0;constraints.weightx=5;
    constraints.gridwidth=1;constraints.gridheight=1;
    buttons.add(m_SaveBut,constraints);
    constraints.gridx=2;constraints.gridy=0;constraints.weightx=5;
    constraints.gridwidth=1;constraints.gridheight=1;
    buttons.add(m_NewBut,constraints);
    
    JPanel src = new JPanel();
    src.setLayout(new BorderLayout());
    src.setBorder(BorderFactory.createCompoundBorder(
		  BorderFactory.createTitledBorder("Result generator"),
		  BorderFactory.createEmptyBorder(0, 5, 5, 5)
		  ));
    src.add(m_RPEditorPanel, BorderLayout.NORTH);
    m_RPEditorPanel.setEnabled(false);

    JPanel dest = new JPanel();
    dest.setLayout(new BorderLayout());
    dest.setBorder(BorderFactory.createCompoundBorder(
		   BorderFactory.createTitledBorder("Destination"),
		   BorderFactory.createEmptyBorder(0, 5, 5, 5)
		   ));
    dest.add(m_RLEditorPanel, BorderLayout.NORTH);
    m_RLEditorPanel.setEnabled(false);

    m_advanceDataSetFirst.setEnabled(false);
    m_advanceIteratorFirst.setEnabled(false);
    m_advanceDataSetFirst.
      setToolTipText("Advance data set before custom generator");
    m_advanceIteratorFirst.
      setToolTipText("Advance custom generator before data set");
    m_advanceDataSetFirst.setSelected(true);
    ButtonGroup bg = new ButtonGroup();
    bg.add(m_advanceDataSetFirst);
    bg.add(m_advanceIteratorFirst);
    m_advanceDataSetFirst.addActionListener(m_RadioListener);
    m_advanceIteratorFirst.addActionListener(m_RadioListener);

    JPanel radioButs = new JPanel();
    radioButs.setBorder(BorderFactory.
			createTitledBorder("Iteration control"));
    radioButs.setLayout(new GridLayout(1, 2));
    radioButs.add(m_advanceDataSetFirst);
    radioButs.add(m_advanceIteratorFirst);

    JPanel simpleIterators = new JPanel();
    simpleIterators.setLayout(new BorderLayout());
    
    JPanel tmp = new JPanel();
    tmp.setLayout(new GridBagLayout());
    //    tmp.setLayout(new GridLayout(1, 2));
    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;
    constraints.insets = new Insets(0,2,0,2);
    tmp.add(m_RunNumberPanel,constraints);
    
    constraints.gridx=1;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=2;
    tmp.add(m_DistributeExperimentPanel, constraints);

    JPanel tmp2 = new JPanel();
    //    tmp2.setLayout(new GridLayout(2, 1));
    tmp2.setLayout(new GridBagLayout());

    constraints.gridx=0;constraints.gridy=0;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;
    constraints.insets = new Insets(0,2,0,2);
    tmp2.add(tmp,constraints);

    constraints.gridx=0;constraints.gridy=1;constraints.weightx=5;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridwidth=1;constraints.gridheight=1;
    tmp2.add(radioButs,constraints);
		   

    simpleIterators.add(tmp2, BorderLayout.NORTH);
    simpleIterators.add(m_DatasetListPanel, BorderLayout.CENTER);
    JPanel iterators = new JPanel();
    iterators.setLayout(new GridLayout(1, 2));
    iterators.add(simpleIterators);
    iterators.add(m_GeneratorPropertyPanel);

    JPanel top = new JPanel();
    top.setLayout(new GridLayout(2, 1));
    top.add(dest);
    top.add(src);

    JPanel notes = new JPanel();
    notes.setLayout(new BorderLayout());
    notes.add(m_NotesButton, BorderLayout.CENTER);
    
    JPanel p2 = new JPanel();
    //    p2.setLayout(new GridLayout(2, 1));
    p2.setLayout(new BorderLayout());
    p2.add(iterators, BorderLayout.CENTER);
    p2.add(notes, BorderLayout.SOUTH);

    JPanel p3 = new JPanel();
    p3.setLayout(new BorderLayout());
    p3.add(buttons, BorderLayout.NORTH);
    p3.add(top, BorderLayout.SOUTH);
    setLayout(new BorderLayout());
    add(p3, BorderLayout.NORTH);
    add(p2, BorderLayout.CENTER);
  }
    
  /**
   * Deletes the notes frame.
   */
  protected void removeNotesFrame() {
    m_NotesFrame.setVisible(false);
  }

  /**
   * Sets the experiment to configure.
   *
   * @param exp a value of type 'Experiment'
   */
  public void setExperiment(Experiment exp) {

    boolean iteratorOn = exp.getUsePropertyIterator();
    Object propArray = exp.getPropertyArray();
    PropertyNode [] propPath = exp.getPropertyPath();

    m_Exp = exp;
    m_SaveBut.setEnabled(true);
    m_RPEditor.setValue(m_Exp.getResultProducer());
    m_RPEditor.setEnabled(true);
    m_RPEditorPanel.setEnabled(true);
    m_RPEditorPanel.repaint();
    m_RLEditor.setValue(m_Exp.getResultListener());
    m_RLEditor.setEnabled(true);
    m_RLEditorPanel.setEnabled(true);
    m_RLEditorPanel.repaint();

    m_NotesText.setText(exp.getNotes());
    m_NotesButton.setEnabled(true);

    m_advanceDataSetFirst.setSelected(m_Exp.getAdvanceDataSetFirst());
    m_advanceIteratorFirst.setSelected(!m_Exp.getAdvanceDataSetFirst());
    m_advanceDataSetFirst.setEnabled(true);
    m_advanceIteratorFirst.setEnabled(true);

    exp.setPropertyPath(propPath);
    exp.setPropertyArray(propArray);
    exp.setUsePropertyIterator(iteratorOn);

    m_GeneratorPropertyPanel.setExperiment(m_Exp);   
    m_RunNumberPanel.setExperiment(m_Exp);
    m_DatasetListPanel.setExperiment(m_Exp);
    m_DistributeExperimentPanel.setExperiment(m_Exp);
    m_Support.firePropertyChange("", null, null);
  }

  /**
   * Gets the currently configured experiment.
   *
   * @return the currently configured experiment.
   */
  public Experiment getExperiment() {

    return m_Exp;
  }
  
  /**
   * Prompts the user to select an experiment file and loads it.
   */
  private void openExperiment() {
    
    int returnVal = m_FileChooser.showOpenDialog(this);
    if (returnVal != JFileChooser.APPROVE_OPTION) {
      return;
    }
    File expFile = m_FileChooser.getSelectedFile();
    
    // add extension if necessary
    if (m_FileChooser.getFileFilter() == m_ExpFilter) {
      if (!expFile.getName().toLowerCase().endsWith(Experiment.FILE_EXTENSION))
        expFile = new File(expFile.getParent(), expFile.getName() + Experiment.FILE_EXTENSION);
    }
    else if (m_FileChooser.getFileFilter() == m_KOMLFilter) {
      if (!expFile.getName().toLowerCase().endsWith(KOML.FILE_EXTENSION))
        expFile = new File(expFile.getParent(), expFile.getName() + KOML.FILE_EXTENSION);
    }
    else if (m_FileChooser.getFileFilter() == m_XMLFilter) {
      if (!expFile.getName().toLowerCase().endsWith(".xml"))
        expFile = new File(expFile.getParent(), expFile.getName() + ".xml");
    }
    
    try {
      Experiment exp = Experiment.read(expFile.getAbsolutePath());
      setExperiment(exp);
      System.err.println("Opened experiment:\n" + m_Exp);
    } catch (Exception ex) {
      ex.printStackTrace();
      JOptionPane.showMessageDialog(this, "Couldn't open experiment file:\n"
				    + expFile
				    + "\nReason:\n" + ex.getMessage(),
				    "Open Experiment",
				    JOptionPane.ERROR_MESSAGE);
      // Pop up error dialog
    }
  }

  /**
   * Prompts the user for a filename to save the experiment to, then saves
   * the experiment.
   */
  private void saveExperiment() {

    int returnVal = m_FileChooser.showSaveDialog(this);
    if (returnVal != JFileChooser.APPROVE_OPTION) {
      return;
    }
    File expFile = m_FileChooser.getSelectedFile();
    
    // add extension if necessary
    if (m_FileChooser.getFileFilter() == m_ExpFilter) {
      if (!expFile.getName().toLowerCase().endsWith(Experiment.FILE_EXTENSION))
        expFile = new File(expFile.getParent(), expFile.getName() + Experiment.FILE_EXTENSION);
    }
    else if (m_FileChooser.getFileFilter() == m_KOMLFilter) {
      if (!expFile.getName().toLowerCase().endsWith(KOML.FILE_EXTENSION))
        expFile = new File(expFile.getParent(), expFile.getName() + KOML.FILE_EXTENSION);
    }
    else if (m_FileChooser.getFileFilter() == m_XMLFilter) {
      if (!expFile.getName().toLowerCase().endsWith(".xml"))
        expFile = new File(expFile.getParent(), expFile.getName() + ".xml");
    }
    
    try {
      Experiment.write(expFile.getAbsolutePath(), m_Exp);
      System.err.println("Saved experiment:\n" + m_Exp);
    } catch (Exception ex) {
      ex.printStackTrace();
      JOptionPane.showMessageDialog(this, "Couldn't save experiment file:\n"
				    + expFile
				    + "\nReason:\n" + ex.getMessage(),
				    "Save Experiment",
				    JOptionPane.ERROR_MESSAGE);
    }
  }

  /**
   * Adds a PropertyChangeListener who will be notified of value changes.
   *
   * @param l a value of type 'PropertyChangeListener'
   */
  public void addPropertyChangeListener(PropertyChangeListener l) {
    m_Support.addPropertyChangeListener(l);
  }

  /**
   * Removes a PropertyChangeListener.
   *
   * @param l a value of type 'PropertyChangeListener'
   */
  public void removePropertyChangeListener(PropertyChangeListener l) {
    m_Support.removePropertyChangeListener(l);
  }

  /**
   * Updates the primary loop iteration control of the experiment
   */
  private void updateRadioLinks() {
    
    m_advanceDataSetFirst.
      setEnabled(m_GeneratorPropertyPanel.getEditorActive());
    m_advanceIteratorFirst.
      setEnabled(m_GeneratorPropertyPanel.getEditorActive());

    if (m_Exp != null) {
      if (!m_GeneratorPropertyPanel.getEditorActive()) {
	m_Exp.setAdvanceDataSetFirst(true);
      } else {
	m_Exp.setAdvanceDataSetFirst(m_advanceDataSetFirst.isSelected());
      }
    }
  }

  /**
   * Tests out the experiment setup from the command line.
   *
   * @param args arguments to the program.
   */
  public static void main(String [] args) {

    try {
      boolean readExp = Utils.getFlag('l', args);
      final boolean writeExp = Utils.getFlag('s', args);
      final String expFile = Utils.getOption('f', args);
      if ((readExp || writeExp) && (expFile.length() == 0)) {
	throw new Exception("A filename must be given with the -f option");
      }
      Experiment exp = null;
      if (readExp) {
	FileInputStream fi = new FileInputStream(expFile);
	ObjectInputStream oi = new ObjectInputStream(
			       new BufferedInputStream(fi));
	exp = (Experiment)oi.readObject();
	oi.close();
      } else {
	exp = new Experiment();
      }
      System.err.println("Initial Experiment:\n" + exp.toString());
      final JFrame jf = new JFrame("Weka Experiment Setup");
      jf.getContentPane().setLayout(new BorderLayout());
      final SetupPanel sp = new SetupPanel();
      //sp.setBorder(BorderFactory.createTitledBorder("Setup"));
      jf.getContentPane().add(sp, BorderLayout.CENTER);
      jf.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  System.err.println("\nFinal Experiment:\n"
			     + sp.m_Exp.toString());
	  // Save the experiment to a file
	  if (writeExp) {
	    try {
	      FileOutputStream fo = new FileOutputStream(expFile);
	      ObjectOutputStream oo = new ObjectOutputStream(
				      new BufferedOutputStream(fo));
	      oo.writeObject(sp.m_Exp);
	      oo.close();
	    } catch (Exception ex) {
              ex.printStackTrace();
	      System.err.println("Couldn't write experiment to: " + expFile
				 + '\n' + ex.getMessage());
	    }
	  }
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setVisible(true);
      System.err.println("Short nap");
      Thread.currentThread().sleep(3000);
      System.err.println("Done");
      sp.setExperiment(exp);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

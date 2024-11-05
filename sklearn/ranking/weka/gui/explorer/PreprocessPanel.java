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
 *    PreprocessPanel.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.explorer;

import weka.core.Capabilities;
import weka.core.CapabilitiesHandler;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.converters.AbstractFileLoader;
import weka.core.converters.AbstractFileSaver;
import weka.core.converters.ConverterUtils;
import weka.core.converters.Loader;
import weka.core.converters.SerializedInstancesLoader;
import weka.core.converters.URLSourcedLoader;
import weka.datagenerators.DataGenerator;
import weka.experiment.InstanceQuery;
import weka.filters.Filter;
import weka.filters.SupervisedFilter;
import weka.filters.unsupervised.attribute.Remove;
import weka.gui.AttributeSelectionPanel;
import weka.gui.AttributeSummaryPanel;
import weka.gui.AttributeVisualizationPanel;
import weka.gui.ConverterFileChooser;
import weka.gui.GenericObjectEditor;
import weka.gui.InstancesSummaryPanel;
import weka.gui.Logger;
import weka.gui.PropertyDialog;
import weka.gui.PropertyPanel;
import weka.gui.SysErrLog;
import weka.gui.TaskLogger;
import weka.gui.ViewerDialog;
import weka.gui.explorer.Explorer.CapabilitiesFilterChangeEvent;
import weka.gui.explorer.Explorer.CapabilitiesFilterChangeListener;
import weka.gui.explorer.Explorer.ExplorerPanel;
import weka.gui.explorer.Explorer.LogHandler;
import weka.gui.sql.SqlViewerDialog;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectOutputStream;
import java.net.URL;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileFilter;

/** 
 * This panel controls simple preprocessing of instances. Summary
 * information on instances and attributes is shown. Filters may be
 * configured to alter the set of instances. Altered instances may
 * also be saved.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 6624 $
 */
public class PreprocessPanel
  extends JPanel 
  implements CapabilitiesFilterChangeListener, ExplorerPanel, LogHandler {

  /** for serialization */
  private static final long serialVersionUID = 6764850273874813049L;
  
  /** Displays simple stats on the working instances */
  protected InstancesSummaryPanel m_InstSummaryPanel =
    new InstancesSummaryPanel();

  /** Click to load base instances from a file */
  protected JButton m_OpenFileBut = new JButton("Open file...");

  /** Click to load base instances from a URL */
  protected JButton m_OpenURLBut = new JButton("Open URL...");

  /** Click to load base instances from a Database */
  protected JButton m_OpenDBBut = new JButton("Open DB...");

  /** Click to generate artificial data */
  protected JButton m_GenerateBut = new JButton("Generate...");

  /** Click to revert back to the last saved point */
  protected JButton m_UndoBut = new JButton("Undo");

  /** Click to open the current instances in a viewer */
  protected JButton m_EditBut = new JButton("Edit...");

  /** Click to apply filters and save the results */
  protected JButton m_SaveBut = new JButton("Save...");
  
  /** Panel to let the user toggle attributes */
  protected AttributeSelectionPanel m_AttPanel = new AttributeSelectionPanel();

  /** Button for removing attributes */
  protected JButton m_RemoveButton = new JButton("Remove");

  /** Displays summary stats on the selected attribute */
  protected AttributeSummaryPanel m_AttSummaryPanel =
    new AttributeSummaryPanel();

  /** Lets the user configure the filter */
  protected GenericObjectEditor m_FilterEditor =
    new GenericObjectEditor();

  /** Filter configuration */
  protected PropertyPanel m_FilterPanel = new PropertyPanel(m_FilterEditor);

  /** Click to apply filters and save the results */
  protected JButton m_ApplyFilterBut = new JButton("Apply");

  /** The file chooser for selecting data files */
  protected ConverterFileChooser m_FileChooser 
    = new ConverterFileChooser(new File(ExplorerDefaults.getInitialDirectory()));

  /** Stores the last URL that instances were loaded from */
  protected String m_LastURL = "http://";
  
  /** Stores the last sql query executed */
  protected String m_SQLQ = new String("SELECT * FROM ?");
 
  /** The working instances */
  protected Instances m_Instances;

  /** The last generator that was selected */
  protected DataGenerator m_DataGenerator = null;

  /** The visualization of the attribute values */
  protected AttributeVisualizationPanel m_AttVisualizePanel = 
    new AttributeVisualizationPanel();

  /** Keeps track of undo points */
  protected File[] m_tempUndoFiles = new File[20]; // set number of undo ops here

  /** The next available slot for an undo point */
  protected int m_tempUndoIndex = 0;
  
  /**
   * Manages sending notifications to people when we change the set of
   * working instances.
   */
  protected PropertyChangeSupport m_Support = new PropertyChangeSupport(this);

  /** A thread for loading/saving instances from a file or URL */
  protected Thread m_IOThread;

  /** The message logger */
  protected Logger m_Log = new SysErrLog();

  /** the parent frame */
  protected Explorer m_Explorer = null;
  
  static {
     weka.core.WekaPackageManager.loadPackages(false);
     GenericObjectEditor.registerEditors();
  }
  
  /**
   * Creates the instances panel with no initial instances.
   */
  public PreprocessPanel() {

    // Create/Configure/Connect components
    m_FilterEditor.setClassType(weka.filters.Filter.class);
    if (ExplorerDefaults.getFilter() != null)
      m_FilterEditor.setValue(ExplorerDefaults.getFilter());
    
    m_FilterEditor.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
        m_ApplyFilterBut.setEnabled(true);
        Capabilities currentCapabilitiesFilter = m_FilterEditor.getCapabilitiesFilter();
        Filter filter = (Filter) m_FilterEditor.getValue();
        Capabilities currentFilterCapabilities = null;
        if (filter != null && currentCapabilitiesFilter != null &&
            (filter instanceof CapabilitiesHandler)) {
          currentFilterCapabilities = ((CapabilitiesHandler)filter).getCapabilities();
          
          if (!currentFilterCapabilities.supportsMaybe(currentCapabilitiesFilter) &&
              !currentFilterCapabilities.supports(currentCapabilitiesFilter)) {
            m_ApplyFilterBut.setEnabled(false);
          }
        }
      }
    });
    m_OpenFileBut.setToolTipText("Open a set of instances from a file");
    m_OpenURLBut.setToolTipText("Open a set of instances from a URL");
    m_OpenDBBut.setToolTipText("Open a set of instances from a database");
    m_GenerateBut.setToolTipText("Generates artificial data");
    m_UndoBut.setToolTipText("Undo the last change to the dataset");
    m_EditBut.setToolTipText("Open the current dataset in a Viewer for editing");
    m_SaveBut.setToolTipText("Save the working relation to a file");
    m_ApplyFilterBut.setToolTipText("Apply the current filter to the data");

    m_FileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
    m_OpenURLBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	setInstancesFromURLQ();
      }
    });
    m_OpenDBBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        SqlViewerDialog dialog = new SqlViewerDialog(null);
        dialog.setVisible(true);
        if (dialog.getReturnValue() == JOptionPane.OK_OPTION)
          setInstancesFromDBQ(dialog.getURL(), dialog.getUser(),
                              dialog.getPassword(), dialog.getQuery(),
                              dialog.getGenerateSparseData());
      }
    });
    m_OpenFileBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	setInstancesFromFileQ();
      }
    });
    m_GenerateBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	generateInstances();
      }
    });
    m_UndoBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	undo();
      }
    });
    m_EditBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        edit();
      }
    });
    m_SaveBut.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	saveWorkingInstancesToFileQ();
      }
    });
    m_ApplyFilterBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  applyFilter((Filter) m_FilterEditor.getValue());
	}
      });
    m_AttPanel.getSelectionModel()
      .addListSelectionListener(new ListSelectionListener() {
	public void valueChanged(ListSelectionEvent e) {
	  if (!e.getValueIsAdjusting()) {	  
	    ListSelectionModel lm = (ListSelectionModel) e.getSource();
	    for (int i = e.getFirstIndex(); i <= e.getLastIndex(); i++) {
	      if (lm.isSelectedIndex(i)) {
		m_AttSummaryPanel.setAttribute(i);
		m_AttVisualizePanel.setAttribute(i);
		break;
	      }
	    }
	  }
	}
    });


    m_InstSummaryPanel.setBorder(BorderFactory
				 .createTitledBorder("Current relation"));
    JPanel attStuffHolderPanel = new JPanel();
    attStuffHolderPanel.setBorder(BorderFactory
				  .createTitledBorder("Attributes"));
    attStuffHolderPanel.setLayout(new BorderLayout());
    attStuffHolderPanel.add(m_AttPanel, BorderLayout.CENTER);
    m_RemoveButton.setEnabled(false);
    m_RemoveButton.setToolTipText("Remove selected attributes.");
    m_RemoveButton.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  try {
	    Remove r = new Remove();
	    int [] selected = m_AttPanel.getSelectedAttributes();
	    if (selected.length == 0) {
	      return;
	    }
	    if (selected.length == m_Instances.numAttributes()) {
	      // Pop up an error optionpane
	      JOptionPane.showMessageDialog(PreprocessPanel.this,
					    "Can't remove all attributes from data!\n",
					    "Remove Attributes",
					    JOptionPane.ERROR_MESSAGE);
	      m_Log.logMessage("Can't remove all attributes from data!");
	      m_Log.statusMessage("Problem removing attributes");
	      return;
	    }
	    r.setAttributeIndicesArray(selected);
	    applyFilter(r);
	  } catch (Exception ex) {
	    if (m_Log instanceof TaskLogger) {
	      ((TaskLogger)m_Log).taskFinished();
	    }
	    // Pop up an error optionpane
	    JOptionPane.showMessageDialog(PreprocessPanel.this,
					  "Problem filtering instances:\n"
					  + ex.getMessage(),
					  "Remove Attributes",
					  JOptionPane.ERROR_MESSAGE);
	    m_Log.logMessage("Problem removing attributes: " + ex.getMessage());
	    m_Log.statusMessage("Problem removing attributes");
	  }
	}
      });

    JPanel p1 = new JPanel();
    p1.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));
    p1.setLayout(new BorderLayout());
    p1.add(m_RemoveButton, BorderLayout.CENTER);
    attStuffHolderPanel.add(p1, BorderLayout.SOUTH);
    m_AttSummaryPanel.setBorder(BorderFactory
		    .createTitledBorder("Selected attribute"));
    m_UndoBut.setEnabled(false);
    m_EditBut.setEnabled(false);
    m_SaveBut.setEnabled(false);
    m_ApplyFilterBut.setEnabled(false);
    
    // Set up the GUI layout
    JPanel buttons = new JPanel();
    buttons.setBorder(BorderFactory.createEmptyBorder(10, 5, 10, 5));
    buttons.setLayout(new GridLayout(1, 6, 5, 5));
    buttons.add(m_OpenFileBut);
    buttons.add(m_OpenURLBut);
    buttons.add(m_OpenDBBut);
    buttons.add(m_GenerateBut);
    buttons.add(m_UndoBut);
    buttons.add(m_EditBut);
    buttons.add(m_SaveBut);

    JPanel attInfo = new JPanel();

    attInfo.setLayout(new BorderLayout());
    attInfo.add(attStuffHolderPanel, BorderLayout.CENTER);

    JPanel filter = new JPanel();
    filter.setBorder(BorderFactory
		    .createTitledBorder("Filter"));
    filter.setLayout(new BorderLayout());
    filter.add(m_FilterPanel, BorderLayout.CENTER);
    filter.add(m_ApplyFilterBut, BorderLayout.EAST); 

    JPanel attVis = new JPanel();
    attVis.setLayout( new GridLayout(2,1) );
    attVis.add(m_AttSummaryPanel);

    JComboBox colorBox = m_AttVisualizePanel.getColorBox();
    colorBox.setToolTipText("The chosen attribute will also be used as the " +
			    "class attribute when a filter is applied.");
    colorBox.addItemListener(new ItemListener() {
      public void itemStateChanged(ItemEvent ie) {
	if (ie.getStateChange() == ItemEvent.SELECTED) {
	  updateCapabilitiesFilter(m_FilterEditor.getCapabilitiesFilter());
	}
      }
    });
    final JButton visAllBut = new JButton("Visualize All");
    visAllBut.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent ae) {
	  if (m_Instances != null) {
	    try {
	      final weka.gui.beans.AttributeSummarizer as = 
		new weka.gui.beans.AttributeSummarizer();
	      as.setColoringIndex(m_AttVisualizePanel.getColoringIndex());
	      as.setInstances(m_Instances);
	      
	      final javax.swing.JFrame jf = new javax.swing.JFrame();
	      jf.getContentPane().setLayout(new java.awt.BorderLayout());
	      
	      jf.getContentPane().add(as, java.awt.BorderLayout.CENTER);
	      jf.addWindowListener(new java.awt.event.WindowAdapter() {
		  public void windowClosing(java.awt.event.WindowEvent e) {
		    visAllBut.setEnabled(true);
		    jf.dispose();
		  }
		});
	      jf.setSize(830,600);
	      jf.setVisible(true);
	    } catch (Exception ex) {
	      ex.printStackTrace();
	    }
	  }
	}
      });
    JPanel histoHolder = new JPanel();
    histoHolder.setLayout(new BorderLayout());
    histoHolder.add(m_AttVisualizePanel, BorderLayout.CENTER);
    JPanel histoControls = new JPanel();
    histoControls.setLayout(new BorderLayout());
    histoControls.add(colorBox, BorderLayout.CENTER);
    histoControls.add(visAllBut, BorderLayout.EAST);
    histoHolder.add(histoControls, BorderLayout.NORTH);
    attVis.add(histoHolder);

    JPanel lhs = new JPanel();
    lhs.setLayout(new BorderLayout());
    lhs.add(m_InstSummaryPanel, BorderLayout.NORTH);
    lhs.add(attInfo, BorderLayout.CENTER);

    JPanel rhs = new JPanel();
    rhs.setLayout(new BorderLayout());
    rhs.add(attVis, BorderLayout.CENTER);

    JPanel relation = new JPanel();
    relation.setLayout(new GridLayout(1, 2));
    relation.add(lhs);
    relation.add(rhs);

    JPanel middle = new JPanel();
    middle.setLayout(new BorderLayout());
    middle.add(filter, BorderLayout.NORTH);
    middle.add(relation, BorderLayout.CENTER);

    setLayout(new BorderLayout());
    add(buttons, BorderLayout.NORTH);
    add(middle, BorderLayout.CENTER);
  }

  /**
   * Sets the Logger to receive informational messages
   *
   * @param newLog the Logger that will now get info messages
   */
  public void setLog(Logger newLog) {

    m_Log = newLog;
  }
  
  /**
   * Tells the panel to use a new base set of instances.
   *
   * @param inst a set of Instances
   */
  public void setInstances(Instances inst) {

    m_Instances = inst;
    try {
      Runnable r = new Runnable() {
	public void run() {
	  m_InstSummaryPanel.setInstances(m_Instances);
	  m_AttPanel.setInstances(m_Instances);
	  m_RemoveButton.setEnabled(true);
	  m_AttSummaryPanel.setInstances(m_Instances);
	  m_AttVisualizePanel.setInstances(m_Instances);

	  // select the first attribute in the list
	  m_AttPanel.getSelectionModel().setSelectionInterval(0, 0);
	  m_AttSummaryPanel.setAttribute(0);
	  m_AttVisualizePanel.setAttribute(0);

	  m_ApplyFilterBut.setEnabled(true);

	  m_Log.logMessage("Base relation is now "
			   + m_Instances.relationName()
			   + " (" + m_Instances.numInstances()
			   + " instances)");
	  m_SaveBut.setEnabled(true);
	  m_EditBut.setEnabled(true);
	  m_Log.statusMessage("OK");
	  // Fire a propertychange event
	  m_Support.firePropertyChange("", null, null);
	  
	  // notify GOEs about change
	  try {
	    // get rid of old filter settings
	    getExplorer().notifyCapabilitiesFilterListener(null);

	    int oldIndex = m_Instances.classIndex();
	    m_Instances.setClassIndex(m_AttVisualizePanel.getColorBox().getSelectedIndex() - 1);
	    
	    // send new ones
	    if (ExplorerDefaults.getInitGenericObjectEditorFilter())
	      getExplorer().notifyCapabilitiesFilterListener(
		  Capabilities.forInstances(m_Instances));
	    else
	      getExplorer().notifyCapabilitiesFilterListener(
		  Capabilities.forInstances(new Instances(m_Instances, 0)));

	    m_Instances.setClassIndex(oldIndex);
	  }
	  catch (Exception e) {
	    e.printStackTrace();
	    m_Log.logMessage(e.toString());
	  }
	}
      };
      if (SwingUtilities.isEventDispatchThread()) {
	r.run();
      } else {
	SwingUtilities.invokeAndWait(r);
      }
    } catch (Exception ex) {
      ex.printStackTrace();
      JOptionPane.showMessageDialog(this,
				    "Problem setting base instances:\n"
				    + ex,
				    "Instances",
				    JOptionPane.ERROR_MESSAGE);
    }
  }

  /**
   * Gets the working set of instances.
   *
   * @return the working instances
   */
  public Instances getInstances() {

    return m_Instances;
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
   * Passes the dataset through the filter that has been configured for use.
   * 
   * @param filter	the filter to apply
   */
  protected void applyFilter(final Filter filter) {

    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	public void run() {
	  try {

	    if (filter != null) {
	    
	      if (m_Log instanceof TaskLogger) {
		((TaskLogger)m_Log).taskStarted();
	      }
	      m_Log.statusMessage("Passing dataset through filter "
		  + filter.getClass().getName());
	      String cmd = filter.getClass().getName();
	      if (filter instanceof OptionHandler)
		cmd += " " + Utils.joinOptions(((OptionHandler) filter).getOptions());
	      m_Log.logMessage("Command: " + cmd);
	      int classIndex = m_AttVisualizePanel.getColoringIndex();
	      if ((classIndex < 0) && (filter instanceof SupervisedFilter)) {
		throw new IllegalArgumentException("Class (colour) needs to " +
						   "be set for supervised " +
						   "filter.");
	      }
	      Instances copy = new Instances(m_Instances);
	      copy.setClassIndex(classIndex);
	      filter.setInputFormat(copy);
	      Instances newInstances = Filter.useFilter(copy, filter);
	      if (newInstances == null || newInstances.numAttributes() < 1) {
		throw new Exception("Dataset is empty.");
	      }
	      m_Log.statusMessage("Saving undo information");
	      addUndoPoint();
	      m_AttVisualizePanel.setColoringIndex(copy.classIndex());
	      // if class was not set before, reset it again after use of filter
	      if (m_Instances.classIndex() < 0)
		newInstances.setClassIndex(-1);
	      m_Instances = newInstances;
	      setInstances(m_Instances);
	      if (m_Log instanceof TaskLogger) {
		((TaskLogger)m_Log).taskFinished();
	      }
	    }
	    
	  } catch (Exception ex) {
	
	    if (m_Log instanceof TaskLogger) {
	      ((TaskLogger)m_Log).taskFinished();
	    }
	    // Pop up an error optionpane
	    JOptionPane.showMessageDialog(PreprocessPanel.this,
					  "Problem filtering instances:\n"
					  + ex.getMessage(),
					  "Apply Filter",
					  JOptionPane.ERROR_MESSAGE);
	    m_Log.logMessage("Problem filtering instances: " + ex.getMessage());
	    m_Log.statusMessage("Problem filtering instances");
	  }
	  m_IOThread = null;
	}
      };
      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    } else {
      JOptionPane.showMessageDialog(this,
				    "Can't apply filter at this time,\n"
				    + "currently busy with other IO",
				    "Apply Filter",
				    JOptionPane.WARNING_MESSAGE);
    }
  }

  /**
   * Queries the user for a file to save instances as, then saves the
   * instances in a background process. This is done in the IO
   * thread, and an error message is popped up if the IO thread is busy.
   */
  public void saveWorkingInstancesToFileQ() {
    
    if (m_IOThread == null) {
      m_FileChooser.setCapabilitiesFilter(m_FilterEditor.getCapabilitiesFilter());
      m_FileChooser.setAcceptAllFileFilterUsed(false);
      int returnVal = m_FileChooser.showSaveDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
	Instances inst = new Instances(m_Instances);
	inst.setClassIndex(m_AttVisualizePanel.getColoringIndex());
	saveInstancesToFile(m_FileChooser.getSaver(), inst);
      }
      FileFilter temp = m_FileChooser.getFileFilter();
      m_FileChooser.setAcceptAllFileFilterUsed(true);
      m_FileChooser.setFileFilter(temp);
    }
    else {
      JOptionPane.showMessageDialog(this,
				    "Can't save at this time,\n"
				    + "currently busy with other IO",
				    "Save Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }
  
  /**
   * saves the data with the specified saver
   * 
   * @param saver	the saver to use for storing the data
   * @param inst	the data to save
   */
  public void saveInstancesToFile(final AbstractFileSaver saver, final Instances inst) {
    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	  public void run() {
	    try {
	      m_Log.statusMessage("Saving to file...");

	      saver.setInstances(inst);
	      saver.writeBatch();
	      
	      m_Log.statusMessage("OK");
	    }
	    catch (Exception ex) {
	      ex.printStackTrace();
	      m_Log.logMessage(ex.getMessage());
	    }
	    m_IOThread = null;
	  }
	};
      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    }
    else {
      JOptionPane.showMessageDialog(this,
				    "Can't save at this time,\n"
				    + "currently busy with other IO",
				    "Saving instances",
				    JOptionPane.WARNING_MESSAGE);
    } 
  }
  
  /**
   * Queries the user for a file to load instances from, then loads the
   * instances in a background process. This is done in the IO
   * thread, and an error message is popped up if the IO thread is busy.
   */
  public void setInstancesFromFileQ() {
    
    if (m_IOThread == null) {
      int returnVal = m_FileChooser.showOpenDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
	try {
	  addUndoPoint();
	}
	catch (Exception ignored) {
	  // ignored
	}

	if (m_FileChooser.getLoader() == null) {
	  JOptionPane.showMessageDialog(this,
	      "Cannot determine file loader automatically, please choose one.",
	      "Load Instances",
	      JOptionPane.ERROR_MESSAGE);
	  converterQuery(m_FileChooser.getSelectedFile());
	}
	else {
	  setInstancesFromFile(m_FileChooser.getLoader());
	}
	    
      }
    } else {
      JOptionPane.showMessageDialog(this,
				    "Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }

  /**
   * Loads (non-sparse) instances from an SQL query the user provided with the
   * SqlViewerDialog, then loads the instances in a background process. This is
   * done in the IO thread, and an error message is popped up if the IO thread
   * is busy.
   * 
   * @param url           the database URL
   * @param user          the user to connect as
   * @param pw            the password of the user
   * @param query         the query for retrieving instances from
   */
  public void setInstancesFromDBQ(String url, String user, 
                                  String pw, String query) {
    setInstancesFromDBQ(url, user, pw, query, false);
  }

  /**
   * Loads instances from an SQL query the user provided with the
   * SqlViewerDialog, then loads the instances in a background process. This is
   * done in the IO thread, and an error message is popped up if the IO thread
   * is busy.
   * 
   * @param url		the database URL
   * @param user	the user to connect as
   * @param pw		the password of the user
   * @param query	the query for retrieving instances from
   * @param sparse	whether to create sparse or non-sparse instances
   */
  public void setInstancesFromDBQ(String url, String user, 
                                  String pw, String query,
                                  boolean sparse) {
    if (m_IOThread == null) {
      try {
	InstanceQuery InstQ = new InstanceQuery();
        InstQ.setDatabaseURL(url);
        InstQ.setUsername(user);
        InstQ.setPassword(pw);
        InstQ.setQuery(query);
        InstQ.setSparseData(sparse);
	
        // we have to disconnect, otherwise we can't change the DB!
        if (InstQ.isConnected())
          InstQ.disconnectFromDatabase();

	InstQ.connectToDatabase();      
	try {
	  addUndoPoint();
	} catch (Exception ignored) {}
	setInstancesFromDB(InstQ);
      } catch (Exception ex) {
	JOptionPane.showMessageDialog(this,
				      "Problem connecting to database:\n"
				      + ex.getMessage(),
				      "Load Instances",
				      JOptionPane.ERROR_MESSAGE);
      }
      
    } else {
      JOptionPane.showMessageDialog(this,
				     "Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }
    
  /**
   * Queries the user for a URL to load instances from, then loads the
   * instances in a background process. This is done in the IO
   * thread, and an error message is popped up if the IO thread is busy.
   */
  public void setInstancesFromURLQ() {
    
    if (m_IOThread == null) {
      try {
	String urlName = (String) JOptionPane.showInputDialog(this,
			"Enter the source URL",
			"Load Instances",
			JOptionPane.QUESTION_MESSAGE,
			null,
			null,
			m_LastURL);
	if (urlName != null) {
	  m_LastURL = urlName;
	  URL url = new URL(urlName);
	  try {
	    addUndoPoint();
	  } catch (Exception ignored) {}
	  setInstancesFromURL(url);
	}
      } catch (Exception ex) {
	ex.printStackTrace();
	JOptionPane.showMessageDialog(this,
				      "Problem with URL:\n"
				      + ex.getMessage(),
				      "Load Instances",
				      JOptionPane.ERROR_MESSAGE);
      }
    } else {
      JOptionPane.showMessageDialog(this,
				    "Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }
  
  /**
   * sets Instances generated via DataGenerators (pops up a Dialog)
   */
  public void generateInstances() {
    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	  public void run() {
	    try {
              // create dialog
              final DataGeneratorPanel generatorPanel = new DataGeneratorPanel();
              final JDialog dialog = new JDialog();
              final JButton generateButton = new JButton("Generate");
              final JCheckBox showOutputCheckBox = 
                                  new JCheckBox("Show generated data as text, incl. comments");

              showOutputCheckBox.setMnemonic('S');
              generatorPanel.setLog(m_Log);
              generatorPanel.setGenerator(m_DataGenerator);
              generatorPanel.setPreferredSize(
                  new Dimension(
                        300, 
                        (int) generatorPanel.getPreferredSize().getHeight()));
              generateButton.setMnemonic('G');
              generateButton.setToolTipText("Generates the dataset according the settings.");
              generateButton.addActionListener(new ActionListener(){
                  public void actionPerformed(ActionEvent evt){
                    // generate
                    generatorPanel.execute();
                    boolean generated = (generatorPanel.getInstances() != null);
                    if (generated)
                      setInstances(generatorPanel.getInstances());

                    // close dialog
                    dialog.dispose();

                    // get last generator
                    m_DataGenerator = generatorPanel.getGenerator();

                    // display output?
                    if ( (generated) && (showOutputCheckBox.isSelected()) )
                      showGeneratedInstances(generatorPanel.getOutput());
                }
              });
              dialog.setTitle("DataGenerator");
              dialog.getContentPane().add(generatorPanel, BorderLayout.CENTER);
              dialog.getContentPane().add(generateButton, BorderLayout.EAST);
              dialog.getContentPane().add(showOutputCheckBox, BorderLayout.SOUTH);
              dialog.pack();
              
              // display dialog
              dialog.setVisible(true);
	    } 
            catch (Exception ex) {
	      ex.printStackTrace();
	      m_Log.logMessage(ex.getMessage());
	    }
	    m_IOThread = null;
	  }
	};
      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    } 
    else {
      JOptionPane.showMessageDialog(this,
				    "Can't generate data at this time,\n"
				    + "currently busy with other IO",
				    "Generate Data",
				    JOptionPane.WARNING_MESSAGE);
    }
  }
  
  /**
   * displays a dialog with the generated instances from the DataGenerator
   * 
   * @param data	the data to display
   */
  protected void showGeneratedInstances(String data) {
    final JDialog dialog = new JDialog();
    final JButton saveButton = new JButton("Save");
    final JButton closeButton = new JButton("Close");
    final JTextArea textData = new JTextArea(data);
    final JPanel panel = new JPanel();
    panel.setLayout(new FlowLayout(FlowLayout.RIGHT));
    textData.setEditable(false);
    textData.setFont(
        new Font("Monospaced", Font.PLAIN, textData.getFont().getSize()));

    saveButton.setMnemonic('S');
    saveButton.setToolTipText("Saves the output to a file");
    saveButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent evt){
        JFileChooser filechooser = new JFileChooser();
        int result = filechooser.showSaveDialog(dialog);
        if (result == JFileChooser.APPROVE_OPTION) {
          try {
            BufferedWriter writer = new BufferedWriter(
                                      new FileWriter(
                                        filechooser.getSelectedFile()));
            writer.write(textData.getText());
            writer.flush();
            writer.close();
            JOptionPane.showMessageDialog(
              dialog, 
              "Output successfully saved to file '" 
              + filechooser.getSelectedFile() + "'!",
              "Information",
              JOptionPane.INFORMATION_MESSAGE);
          }
          catch (Exception e) {
            e.printStackTrace();
          }
          dialog.dispose();
        }
      }
    });
    closeButton.setMnemonic('C');
    closeButton.setToolTipText("Closes the dialog");
    closeButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent evt){
        dialog.dispose();
      }
    });
    panel.add(saveButton);
    panel.add(closeButton);
    dialog.setTitle("Generated Instances (incl. comments)");
    dialog.getContentPane().add(new JScrollPane(textData), BorderLayout.CENTER);
    dialog.getContentPane().add(panel, BorderLayout.SOUTH);
    dialog.pack();

    // make sure, it's not bigger than 80% of the screen
    Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int width  = dialog.getWidth() > screen.getWidth()*0.8
                    ? (int) (screen.getWidth()*0.8) : dialog.getWidth();
    int height = dialog.getHeight() > screen.getHeight()*0.8 
                    ? (int) (screen.getHeight()*0.8) : dialog.getHeight();
    dialog.setSize(width, height);
    
    // display dialog
    dialog.setVisible(true);
  }

  /**
   * Pops up generic object editor with list of conversion filters
   *
   * @param f the File
   */
  private void converterQuery(final File f) {
    final GenericObjectEditor convEd = new GenericObjectEditor(true);

    try {
      convEd.setClassType(weka.core.converters.Loader.class);
      convEd.setValue(new weka.core.converters.CSVLoader());
      ((GenericObjectEditor.GOEPanel)convEd.getCustomEditor())
	.addOkListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
	      tryConverter((Loader)convEd.getValue(), f);
	    }
	  });
    } catch (Exception ex) {
    }

    PropertyDialog pd;
    if (PropertyDialog.getParentDialog(this) != null)
      pd = new PropertyDialog(PropertyDialog.getParentDialog(this), convEd, 100, 100);
    else
      pd = new PropertyDialog(PropertyDialog.getParentFrame(this), convEd, 100, 100);
    pd.setVisible(true);
  }

  /**
   * Applies the selected converter
   *
   * @param cnv the converter to apply to the input file
   * @param f the input file
   */
  private void tryConverter(final Loader cnv, final File f) {

    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	  public void run() {
	    try {
	      cnv.setSource(f);
	      Instances inst = cnv.getDataSet();
	      setInstances(inst);
	    } catch (Exception ex) {
	      m_Log.statusMessage(cnv.getClass().getName()+" failed to load "
				 +f.getName());
	      JOptionPane.showMessageDialog(PreprocessPanel.this,
					    cnv.getClass().getName()+" failed to load '"
					    + f.getName() + "'.\n"
					    + "Reason:\n" + ex.getMessage(),
					    "Convert File",
					    JOptionPane.ERROR_MESSAGE);
	      m_IOThread = null;
	      converterQuery(f);
	    }
	    m_IOThread = null;
	  }
	};
      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    }
  }

  /**
   * Loads results from a set of instances retrieved with the supplied loader. 
   * This is started in the IO thread, and a dialog is popped up
   * if there's a problem.
   *
   * @param loader	the loader to use
   */
  public void setInstancesFromFile(final AbstractFileLoader loader) {
      
    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	public void run() {
	  try {
	    m_Log.statusMessage("Reading from file...");
	    Instances inst = loader.getDataSet();
	    setInstances(inst);
	  }
	  catch (Exception ex) {
		  ex.printStackTrace();
	    m_Log.statusMessage(
		"File '" + loader.retrieveFile() + "' not recognised as an '"
		+ loader.getFileDescription() + "' file.");
	    m_IOThread = null;
	    if (JOptionPane.showOptionDialog(PreprocessPanel.this,
					     "File '" + loader.retrieveFile()
					     + "' not recognised as an '"
					     + loader.getFileDescription() 
					     + "' file.\n"
					     + "Reason:\n" + ex.getMessage(),
					     "Load Instances",
					     0,
					     JOptionPane.ERROR_MESSAGE,
					     null,
					     new String[] {"OK", "Use Converter"},
					     null) == 1) {
	    
	      converterQuery(loader.retrieveFile());
	    }
	  }
	  m_IOThread = null;
	}
      };
      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    } else {
      JOptionPane.showMessageDialog(this,
				    "Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }
  
  /**
   * Loads instances from a database
   *
   * @param iq the InstanceQuery object to load from (this is assumed
   * to have been already connected to a valid database).
   */
  public void setInstancesFromDB(final InstanceQuery iq) {
    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	public void run() {
	  
	  try {
	    m_Log.statusMessage("Reading from database...");
	    final Instances i = iq.retrieveInstances();
	    SwingUtilities.invokeAndWait(new Runnable() {
	      public void run() {
		setInstances(new Instances(i));
	      }
	    });
	    iq.disconnectFromDatabase();
	  } catch (Exception ex) {
	    m_Log.statusMessage("Problem executing DB query "+m_SQLQ);
	    JOptionPane.showMessageDialog(PreprocessPanel.this,
					  "Couldn't read from database:\n"
					  + ex.getMessage(),
					  "Load Instances",
					  JOptionPane.ERROR_MESSAGE);
	  }

	   m_IOThread = null;
	}
      };

      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    } else {
       JOptionPane.showMessageDialog(this,
				    "Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }

  /**
   * Loads instances from a URL.
   *
   * @param u the URL to load from.
   */
  public void setInstancesFromURL(final URL u) {

    if (m_IOThread == null) {
      m_IOThread = new Thread() {
	public void run() {

	  try {
	    m_Log.statusMessage("Reading from URL...");
	    AbstractFileLoader loader = ConverterUtils.getURLLoaderForFile(u.toString());
	    if (loader == null)
	      throw new Exception("No suitable URLSourcedLoader found for URL!\n" + u);
	    ((URLSourcedLoader) loader).setURL(u.toString());
	    setInstances(loader.getDataSet());
	  } catch (Exception ex) {
	    ex.printStackTrace();
	    m_Log.statusMessage("Problem reading " + u);
	    JOptionPane.showMessageDialog(PreprocessPanel.this,
					  "Couldn't read from URL:\n"
					  + u + "\n"
					  + ex.getMessage(),
					  "Load Instances",
					  JOptionPane.ERROR_MESSAGE);
	  }

	  m_IOThread = null;
	}
      };
      m_IOThread.setPriority(Thread.MIN_PRIORITY); // UI has most priority
      m_IOThread.start();
    } else {
      JOptionPane.showMessageDialog(this,
				    "Can't load at this time,\n"
				    + "currently busy with other IO",
				    "Load Instances",
				    JOptionPane.WARNING_MESSAGE);
    }
  }

  /**
   * Backs up the current state of the dataset, so the changes can be undone.
   * 
   * @throws Exception 	if an error occurs
   */
  public void addUndoPoint() throws Exception {
    
    if (m_Instances != null) {
      // create temporary file
      File tempFile = File.createTempFile("weka", SerializedInstancesLoader.FILE_EXTENSION);
      tempFile.deleteOnExit();

      ObjectOutputStream oos = 
	new ObjectOutputStream(
	new BufferedOutputStream(
	new FileOutputStream(tempFile)));
    
      oos.writeObject(m_Instances);
      oos.flush();
      oos.close();

      // update undo file list
      if (m_tempUndoFiles[m_tempUndoIndex] != null) {
	// remove undo points that are too old
	m_tempUndoFiles[m_tempUndoIndex].delete();
      }
      m_tempUndoFiles[m_tempUndoIndex] = tempFile;
      if (++m_tempUndoIndex >= m_tempUndoFiles.length) {
	// wrap pointer around
	m_tempUndoIndex = 0;
      }

      m_UndoBut.setEnabled(true);
    }
  }

  /**
   * Reverts to the last backed up version of the dataset.
   */
  public void undo() {

    if (--m_tempUndoIndex < 0) {
      // wrap pointer around
      m_tempUndoIndex = m_tempUndoFiles.length-1;
    }
    
    if (m_tempUndoFiles[m_tempUndoIndex] != null) {
      // load instances from the temporary file
      AbstractFileLoader loader = ConverterUtils.getLoaderForFile(m_tempUndoFiles[m_tempUndoIndex]);
      try {
	loader.setFile(m_tempUndoFiles[m_tempUndoIndex]);
	setInstancesFromFile(loader);
      }
      catch (Exception e) {
	e.printStackTrace();
	m_Log.logMessage(e.toString());
	JOptionPane.showMessageDialog(PreprocessPanel.this,
	    "Cannot perform undo operation!\n" + e.toString(),
	    "Undo",
	    JOptionPane.ERROR_MESSAGE);
      }

      // update undo file list
      m_tempUndoFiles[m_tempUndoIndex] = null;
    }
    
    // update undo button
    int temp = m_tempUndoIndex-1;
    if (temp < 0) {
      temp = m_tempUndoFiles.length-1;
    }
    m_UndoBut.setEnabled(m_tempUndoFiles[temp] != null);
  }
  
  /**
   * edits the current instances object in the viewer 
   */
  public void edit() {
    ViewerDialog        dialog;
    int                 result;
    Instances           copy;
    Instances           newInstances;
    
    final int classIndex = m_AttVisualizePanel.getColoringIndex();
    copy   = new Instances(m_Instances);
    copy.setClassIndex(classIndex);
    dialog = new ViewerDialog(null);
    result = dialog.showDialog(copy);
    if (result == ViewerDialog.APPROVE_OPTION) {
      try {
        addUndoPoint();
      }
      catch (Exception e) {
        e.printStackTrace();
      }
      // if class was not set before, reset it again after use of filter
      newInstances = dialog.getInstances();
      if (m_Instances.classIndex() < 0)
        newInstances.setClassIndex(-1);
      setInstances(newInstances);
    }
  }

  /**
   * Sets the Explorer to use as parent frame (used for sending notifications
   * about changes in the data)
   * 
   * @param parent	the parent frame
   */
  public void setExplorer(Explorer parent) {
    m_Explorer = parent;
  }
  
  /**
   * returns the parent Explorer frame
   * 
   * @return		the parent
   */
  public Explorer getExplorer() {
    return m_Explorer;
  }
  
  /**
   * updates the capabilities filter of the GOE
   * 
   * @param filter	the new filter to use
   */
  protected void updateCapabilitiesFilter(Capabilities filter) {
    Instances 		tempInst;
    Capabilities 	filterClass;

    if (filter == null) {
      m_FilterEditor.setCapabilitiesFilter(new Capabilities(null));
      return;
    }
    
    if (!ExplorerDefaults.getInitGenericObjectEditorFilter())
      tempInst = new Instances(m_Instances, 0);
    else
      tempInst = new Instances(m_Instances);
    tempInst.setClassIndex(m_AttVisualizePanel.getColorBox().getSelectedIndex() - 1);

    try {
      filterClass = Capabilities.forInstances(tempInst);
    }
    catch (Exception e) {
      filterClass = new Capabilities(null);
    }
    
    // set new filter
    m_FilterEditor.setCapabilitiesFilter(filterClass);
    
    // check capabilities
    m_ApplyFilterBut.setEnabled(true);
    Capabilities currentCapabilitiesFilter = m_FilterEditor.getCapabilitiesFilter();
    Filter currentFilter = (Filter) m_FilterEditor.getValue();
    Capabilities currentFilterCapabilities = null;
    if (currentFilter != null && currentCapabilitiesFilter != null &&
        (currentFilter instanceof CapabilitiesHandler)) {
      currentFilterCapabilities = ((CapabilitiesHandler)currentFilter).getCapabilities();
      
      if (!currentFilterCapabilities.supportsMaybe(currentCapabilitiesFilter) &&
          !currentFilterCapabilities.supports(currentCapabilitiesFilter)) {
        m_ApplyFilterBut.setEnabled(false);
      }
    }
  }
  
  /**
   * method gets called in case of a change event
   * 
   * @param e		the associated change event
   */
  public void capabilitiesFilterChanged(CapabilitiesFilterChangeEvent e) {
    if (e.getFilter() == null)
      updateCapabilitiesFilter(null);
    else
      updateCapabilitiesFilter((Capabilities) e.getFilter().clone());
  }
  
  /**
   * Returns the title for the tab in the Explorer
   * 
   * @return 		the title of this tab
   */
  public String getTabTitle() {
    return "Preprocess";
  }
  
  /**
   * Returns the tooltip for the tab in the Explorer
   * 
   * @return 		the tooltip of this tab
   */
  public String getTabTitleToolTip() {
    return "Open/Edit/Save instances";
  }
  
  /**
   * Tests out the instance-preprocessing panel from the command line.
   *
   * @param args ignored
   */
  public static void main(String [] args) {

    try {
      final JFrame jf = new JFrame("Weka Explorer: Preprocess");
      jf.getContentPane().setLayout(new BorderLayout());
      final PreprocessPanel sp = new PreprocessPanel();
      jf.getContentPane().add(sp, BorderLayout.CENTER);
      weka.gui.LogPanel lp = new weka.gui.LogPanel();
      sp.setLog(lp);
      jf.getContentPane().add(lp, BorderLayout.SOUTH);
      jf.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  jf.dispose();
	  System.exit(0);
	}
      });
      jf.pack();
      jf.setSize(800, 600);
      jf.setVisible(true);
    } catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

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
 * ConverterFileChooser.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.gui;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.converters.AbstractFileLoader;
import weka.core.converters.AbstractFileSaver;
import weka.core.converters.AbstractLoader;
import weka.core.converters.AbstractSaver;
import weka.core.converters.ConverterUtils;
import weka.core.converters.FileSourcedConverter;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileFilter;

/**
 * A specialized JFileChooser that lists all available file Loaders and Savers.
 * To list only savers that can handle the data that is about to be saved, one
 * can set a Capabilities filter.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4915 $
 * @see	    #setCapabilitiesFilter(Capabilities)
 */
public class ConverterFileChooser
  extends JFileChooser {

  /** for serialization. */
  private static final long serialVersionUID = -5373058011025481738L;
  
  /** unhandled type of dialog. */
  public final static int UNHANDLED_DIALOG = 0;
  
  /** the loader dialog. */
  public final static int LOADER_DIALOG = 1;
  
  /** the saver dialog. */
  public final static int SAVER_DIALOG = 2;
  
  /** the file chooser itself. */
  protected ConverterFileChooser m_Self;

  /** the file filters for the loaders. */
  protected static Vector<ExtensionFileFilter> m_LoaderFileFilters;

  /** the file filters for the savers. */
  protected static Vector<ExtensionFileFilter> m_SaverFileFilters;
  
  /** the type of dialog to display. */
  protected int m_DialogType;

  /** the converter that was chosen by the user. */
  protected Object m_CurrentConverter;
  
  /** the configure button. */
  protected JButton m_ConfigureButton;

  /** the propertychangelistener. */
  protected PropertyChangeListener m_Listener;
  
  /** the last filter that was used for opening/saving. */
  protected FileFilter m_LastFilter;
  
  /** the Capabilities filter for the savers. */
  protected Capabilities m_CapabilitiesFilter;
  
  /** whether to popup a dialog in case the file already exists (only save
   * dialog). */
  protected boolean m_OverwriteWarning = true;

  /** whether the file to be opened must exist (only open dialog). */
  protected boolean m_FileMustExist = true;

  /** the checkbox for bringing up the GenericObjectEditor. */
  protected JCheckBox m_CheckBoxOptions;
  
  /** the note about the options dialog. */
  protected JLabel m_LabelOptions;

  /** the GOE for displaying the options of a loader/saver. */
  protected GenericObjectEditor m_Editor = null;
  
  /** whether the GOE was OKed or Canceled. */
  protected int m_EditorResult;
  
  /** whether to display only core converters (hardcoded in ConverterUtils).
   * Necessary for RMI/Remote Experiments for instance.
   * 
   * @see ConverterUtils#CORE_FILE_LOADERS
   * @see ConverterUtils#CORE_FILE_SAVERS */
  protected boolean m_CoreConvertersOnly = false;
  
  static {
    initFilters(true, ConverterUtils.getFileLoaders());
    initFilters(false, ConverterUtils.getFileSavers());
  }
  
  /**
   * onstructs a FileChooser pointing to the user's default directory.
   */
  public ConverterFileChooser() {
    super();
    initialize();
  }

  /**
   * Constructs a FileChooser using the given File as the path.
   * 
   * @param currentDirectory	the path to start in
   */
  public ConverterFileChooser(File currentDirectory) {
    super(currentDirectory);
    initialize();
  }
  
  /**
   * Constructs a FileChooser using the given path.
   * 
   * @param currentDirectory	the path to start in
   */
  public ConverterFileChooser(String currentDirectory) {
    super(currentDirectory);
    initialize();
  }
  
  /**
   * Further initializations.
   */
  protected void initialize() {
    JPanel	panel;
    JPanel	panel2;
    
    m_Self = this;
    
    m_CheckBoxOptions = new JCheckBox("Invoke options dialog");
    m_CheckBoxOptions.setMnemonic('I');
    m_LabelOptions = new JLabel("<html><br>Note:<br><br>Some file formats offer additional<br>options which can be customized<br>when invoking the options dialog.</html>");
    panel = new JPanel(new BorderLayout());
    panel.add(m_CheckBoxOptions, BorderLayout.NORTH);
    panel2 = new JPanel(new BorderLayout());
    panel2.add(m_LabelOptions, BorderLayout.NORTH);
    panel2.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    panel.add(panel2, BorderLayout.CENTER);
    setAccessory(panel);

    m_Editor = new GenericObjectEditor(false);
    ((GenericObjectEditor.GOEPanel) m_Editor.getCustomEditor()).addOkListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	m_EditorResult     = JFileChooser.APPROVE_OPTION;
	m_CurrentConverter = m_Editor.getValue();
	// thanks to serialization and transient readers/streams, we have
	// to set the file again to initialize the converter again
	try {
	  ((FileSourcedConverter) m_CurrentConverter).setFile(((FileSourcedConverter) m_CurrentConverter).retrieveFile());
	}
	catch (Exception ex) {
	  // ignored
	}
      }
    });
    ((GenericObjectEditor.GOEPanel) m_Editor.getCustomEditor()).addCancelListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
	m_EditorResult = JFileChooser.CANCEL_OPTION;
      }
    });
  }

  /**
   * filters out all non-core loaders if only those should be displayed.
   * 
   * @param list	the list of filters to check
   * @return		the filtered list of filters
   * @see		#m_CoreConvertersOnly
   */
  protected Vector<ExtensionFileFilter> filterNonCoreLoaderFileFilters(Vector<ExtensionFileFilter> list) {
    Vector<ExtensionFileFilter> result;
    int				i;
    ExtensionFileFilter		filter;
    AbstractLoader		loader;
    
    if (!getCoreConvertersOnly()) {
      result = list;
    }
    else {
      result = new Vector<ExtensionFileFilter>();
      for (i = 0; i < list.size(); i++) {
	filter = list.get(i);
	loader = ConverterUtils.getLoaderForExtension(filter.getExtensions()[0]);
	if (ConverterUtils.isCoreFileLoader(loader.getClass().getName()))
	  result.add(filter);
      }
    }
    
    return result;
  }

  /**
   * filters out all non-core savers if only those should be displayed.
   * 
   * @param list	the list of filters to check
   * @return		the filtered list of filters
   * @see		#m_CoreConvertersOnly
   */
  protected Vector<ExtensionFileFilter> filterNonCoreSaverFileFilters(Vector<ExtensionFileFilter> list) {
    Vector<ExtensionFileFilter> result;
    int				i;
    ExtensionFileFilter		filter;
    AbstractSaver		saver;
    
    if (!getCoreConvertersOnly()) {
      result = list;
    }
    else {
      result = new Vector<ExtensionFileFilter>();
      for (i = 0; i < list.size(); i++) {
	filter = list.get(i);
	saver = ConverterUtils.getSaverForExtension(filter.getExtensions()[0]);
	if (ConverterUtils.isCoreFileSaver(saver.getClass().getName()))
	  result.add(filter);
      }
    }
    
    return result;
  }
  
  /**
   * filters the list of file filters according to the currently set.
   * Capabilities
   * 
   * @param list	the filters to check
   * @return		the filtered list of filters
   */
  protected Vector<ExtensionFileFilter> filterSaverFileFilters(Vector<ExtensionFileFilter> list) {
    Vector<ExtensionFileFilter>	result;
    int				i;
    ExtensionFileFilter		filter;
    AbstractSaver		saver;
    
    if (m_CapabilitiesFilter == null) {
      result = list;
    }
    else {
      result = new Vector<ExtensionFileFilter>();
      
      for (i = 0; i < list.size(); i++) {
	filter = list.get(i);
	saver  = ConverterUtils.getSaverForExtension(filter.getExtensions()[0]);
	if (saver.getCapabilities().supports(m_CapabilitiesFilter))
	  result.add(filter);
      }
    }

    return result;
  }
  
  /**
   * initializes the ExtensionFileFilters.
   * 
   * @param loader	if true then the loader filter are initialized
   * @param classnames	the classnames of the converters
   */
  protected static void initFilters(boolean loader, Vector<String> classnames) {
    int				i;
    int				n;
    String 			classname;
    Class 			cls;
    String[] 			ext;
    String 			desc;
    FileSourcedConverter 	converter;
    ExtensionFileFilter 	filter;
    
    if (loader)
      m_LoaderFileFilters = new Vector<ExtensionFileFilter>();
    else
      m_SaverFileFilters  = new Vector<ExtensionFileFilter>();

    for (i = 0; i < classnames.size(); i++) {
      classname = (String) classnames.get(i);
      
      // get data from converter
      try {
	cls       = Class.forName(classname);
	converter = (FileSourcedConverter) cls.newInstance();
	ext       = converter.getFileExtensions();
	desc      = converter.getFileDescription();
      }
      catch (Exception e) {
	cls       = null;
	converter = null;
	ext       = new String[0];
	desc      = "";
      }
      
      if (converter == null)
	continue;

      // loader?
      if (loader) {
	for (n = 0; n < ext.length; n++) {
	  filter = new ExtensionFileFilter(ext[n], desc + " (*" + ext[n] + ")");
	  m_LoaderFileFilters.add(filter);
	}
      }
      else {
	for (n = 0; n < ext.length; n++) {
	  filter = new ExtensionFileFilter(ext[n], desc + " (*" + ext[n] + ")");
	  m_SaverFileFilters.add(filter);
	}
      }
    }
  }
  
  /**
   * initializes the GUI.
   * 
   * @param dialogType		the type of dialog to setup the GUI for
   */
  protected void initGUI(int dialogType) {
    Vector<ExtensionFileFilter>	list;
    int				i;
    boolean 			acceptAll;

    // backup current state
    acceptAll = isAcceptAllFileFilterUsed();

    // setup filters
    resetChoosableFileFilters();
    setAcceptAllFileFilterUsed(acceptAll);
    if (dialogType == LOADER_DIALOG)
      list = filterNonCoreLoaderFileFilters(m_LoaderFileFilters);
    else
      list = filterSaverFileFilters(filterNonCoreSaverFileFilters(m_SaverFileFilters));
    for (i = 0; i < list.size(); i++) {
      addChoosableFileFilter(list.get(i));
    }
    if (list.size() > 0) {
      if ( (m_LastFilter == null) || (!list.contains(m_LastFilter)) )
	setFileFilter(list.get(0));
      else
	setFileFilter(m_LastFilter);
    }

    // listener
    if (m_Listener != null)
      removePropertyChangeListener(m_Listener);
    m_Listener = new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent evt) {
	// filter changed
	if (evt.getPropertyName().equals(FILE_FILTER_CHANGED_PROPERTY)) {
	  updateCurrentConverter();
	}
      }
    };
    addPropertyChangeListener(m_Listener);
    
    // initial setup
    if (dialogType == LOADER_DIALOG) {
      m_Editor.setClassType(AbstractFileLoader.class);
      m_Editor.setValue(new weka.core.converters.ArffLoader());
    }
    else {
      m_Editor.setClassType(AbstractFileSaver.class);
      m_Editor.setValue(new weka.core.converters.ArffSaver());
    }
    
    updateCurrentConverter();
  }

  /**
   * sets the capabilities that the savers must have. use null if all should be
   * listed.
   * 
   * @param value	the minimum Capabilities the savers must have
   */
  public void setCapabilitiesFilter(Capabilities value) {
    m_CapabilitiesFilter = (Capabilities) value.clone();
  }
  
  /**
   * returns the capabilities filter for the savers, can be null if all are
   * listed.
   * 
   * @return		the minimum Capabilities the savers must have
   */
  public Capabilities getCapabilitiesFilter() {
    if (m_CapabilitiesFilter != null)
      return (Capabilities) m_CapabilitiesFilter.clone();
    else
      return null;
  }

  /**
   * Whether a warning is popped up if the file that is to be saved already
   * exists (only save dialog).
   * 
   * @param value	if true a warning will be popup
   */
  public void setOverwriteWarning(boolean value) {
    m_OverwriteWarning = value;
  }
  
  /**
   * Returns whether a popup appears with a warning that the file already 
   * exists (only save dialog).
   * 
   * @return		true if a warning pops up
   */
  public boolean getOverwriteWarning() {
    return m_OverwriteWarning;
  }

  /**
   * Whether the selected file must exst (only open dialog).
   * 
   * @param value	if true the file must exist
   */
  public void setFileMustExist(boolean value) {
    m_FileMustExist = value;
  }
  
  /**
   * Returns whether the selected file must exist (only open dialog).
   * 
   * @return		true if the file must exist
   */
  public boolean getFileMustExist() {
    return m_FileMustExist;
  }

  /**
   * Whether to display only the hardocded core converters. Necessary for
   * RMI/Remote Experiments (dynamic class discovery doesn't work there!).
   * 
   * @param value	if true only the core converters will be displayed
   * @see 		#m_CoreConvertersOnly
   */
  public void setCoreConvertersOnly(boolean value) {
    m_CoreConvertersOnly = value;
  }
  
  /**
   * Returns whether only the hardcoded core converters are displayed. 
   * Necessary for RMI/REmote Experiments (dynamic class discovery doesn't 
   * work there!).
   * 
   * @return		true if the file must exist
   * @see 		#m_CoreConvertersOnly
   */
  public boolean getCoreConvertersOnly() {
    return m_CoreConvertersOnly;
  }
  
  /**
   * Pops a custom file chooser dialog with a custom approve button. Throws
   * an exception, if the dialog type is UNHANDLED_DIALOG.
   * 
   * @param parent		the parent of this dialog
   * @param approveButtonText	the text for the OK button
   * @return			the user's action
   */
  public int showDialog(Component parent, String approveButtonText) {
    if (m_DialogType == UNHANDLED_DIALOG)
      throw new IllegalStateException("Either use showOpenDialog or showSaveDialog!");
    else
      return super.showDialog(parent, approveButtonText);
  }
  
  /**
   * Pops up an "Open File" file chooser dialog.
   * 
   * @param parent		the parent of this file chooser
   * @return			the result of the user's action
   */
  public int showOpenDialog(Component parent) {
    m_DialogType       = LOADER_DIALOG;
    m_CurrentConverter = null;
    
    initGUI(LOADER_DIALOG);
    
    int result = super.showOpenDialog(parent);
    
    m_DialogType = UNHANDLED_DIALOG;
    removePropertyChangeListener(m_Listener);

    // do we have to add the extension?
    if ( (result == APPROVE_OPTION) && (getSelectedFile().isFile()) ) {
      if (getFileFilter() instanceof ExtensionFileFilter) {
	String filename = getSelectedFile().getAbsolutePath();
	String[] extensions = ((ExtensionFileFilter) getFileFilter()).getExtensions();
	if (!filename.endsWith(extensions[0])) {
	  filename += extensions[0];
	  setSelectedFile(new File(filename));
	}
      }
    }
    
    // does file exist?
    if (    (result == APPROVE_OPTION) 
	 && (getFileMustExist()) 
	 && (getSelectedFile().isFile())
	 && (!getSelectedFile().exists()) ) {
      int retVal = JOptionPane.showConfirmDialog(
	  parent, 
	  "The file '" 
	  + getSelectedFile() 
	  + "' does not exist - please select again!");
      if (retVal == JOptionPane.OK_OPTION)
	result = showOpenDialog(parent);
      else
	result = CANCEL_OPTION;
    }

    if (result == APPROVE_OPTION) {
      m_LastFilter = getFileFilter();
      configureCurrentConverter(LOADER_DIALOG);

      // bring up options dialog?
      if (m_CheckBoxOptions.isSelected()) {
	m_EditorResult = JFileChooser.CANCEL_OPTION;
	m_Editor.setValue(m_CurrentConverter);
	PropertyDialog pd;
	if (PropertyDialog.getParentDialog(this) != null)
	  pd = new PropertyDialog(PropertyDialog.getParentDialog(this), m_Editor);
	else
	  pd = new PropertyDialog(PropertyDialog.getParentFrame(this), m_Editor);
	pd.setVisible(true);
	result = m_EditorResult;
      }
    }
    
    return result;
  }
  
  /**
   * Pops up an "Save File" file chooser dialog.
   * 
   * @param parent		the parent of this file chooser
   * @return			the result of the user's action
   */
  public int showSaveDialog(Component parent) {
    m_DialogType       = SAVER_DIALOG;
    m_CurrentConverter = null;
    
    initGUI(SAVER_DIALOG);
    
    boolean acceptAll = isAcceptAllFileFilterUsed();

    // using "setAcceptAllFileFilterUsed" messes up the currently selected 
    // file filter/file, hence backup/restore of currently selected 
    // file filter/file
    FileFilter currentFilter = getFileFilter();
    File currentFile = getSelectedFile();
    setAcceptAllFileFilterUsed(false);
    setFileFilter(currentFilter);
    setSelectedFile(currentFile);
    
    int result = super.showSaveDialog(parent);
    
    // do we have to add the extension?
    if (result == APPROVE_OPTION) {
      if (getFileFilter() instanceof ExtensionFileFilter) {
	String filename = getSelectedFile().getAbsolutePath();
	String[] extensions = ((ExtensionFileFilter) getFileFilter()).getExtensions();
	if (!filename.endsWith(extensions[0])) {
	  filename += extensions[0];
	  setSelectedFile(new File(filename));
	}
      }
    }
    
    // using "setAcceptAllFileFilterUsed" messes up the currently selected 
    // file filter/file, hence backup/restore of currently selected 
    // file filter/file
    currentFilter = getFileFilter();
    currentFile = getSelectedFile();
    setAcceptAllFileFilterUsed(acceptAll);
    setFileFilter(currentFilter);
    setSelectedFile(currentFile);
    
    m_DialogType = UNHANDLED_DIALOG;
    removePropertyChangeListener(m_Listener);

    // overwrite the file?
    if (    (result == APPROVE_OPTION) 
	 && (getOverwriteWarning()) 
	 && (getSelectedFile().exists()) ) {
      int retVal = JOptionPane.showConfirmDialog(
	  	  parent, 
	  	  "The file '" 
	  	  + getSelectedFile() 
	  	  + "' already exists - overwrite it?");
      if (retVal == JOptionPane.OK_OPTION)
	result = APPROVE_OPTION;
      else if (retVal == JOptionPane.NO_OPTION)
	result = showSaveDialog(parent);
      else
	result = CANCEL_OPTION;
    }
    
    if (result == APPROVE_OPTION) {
      m_LastFilter = getFileFilter();
      configureCurrentConverter(SAVER_DIALOG);

      // bring up options dialog?
      if (m_CheckBoxOptions.isSelected()) {
	m_EditorResult = JFileChooser.CANCEL_OPTION;
	m_Editor.setValue(m_CurrentConverter);
	PropertyDialog pd;
	if (PropertyDialog.getParentDialog(this) != null)
	  pd = new PropertyDialog(PropertyDialog.getParentDialog(this), m_Editor);
	else
	  pd = new PropertyDialog(PropertyDialog.getParentFrame(this), m_Editor);
	pd.setVisible(true);
	result = m_EditorResult;
      }
    }
    
    return result;
  }
  
  /**
   * returns the loader that was chosen by the user, can be null in case the
   * user aborted the dialog or the save dialog was shown.
   * 
   * @return		the chosen loader, if any
   */
  public AbstractFileLoader getLoader() {
    configureCurrentConverter(LOADER_DIALOG);
    
    if (m_CurrentConverter instanceof AbstractFileSaver)
      return null;
    else
      return (AbstractFileLoader) m_CurrentConverter;
  }
  
  /**
   * returns the saver that was chosen by the user, can be null in case the
   * user aborted the dialog or the open dialog was shown.
   * 
   * @return		the chosen saver, if any
   */
  public AbstractFileSaver getSaver() {
    configureCurrentConverter(SAVER_DIALOG);
    
    if (m_CurrentConverter instanceof AbstractFileLoader)
      return null;
    else
      return (AbstractFileSaver) m_CurrentConverter;
  }
  
  /**
   * sets the current converter according to the current filefilter.
   */
  protected void updateCurrentConverter() {
    String[]	extensions;
    Object	newConverter;
    
    if (getFileFilter() == null)
      return;
    
    if (!isAcceptAllFileFilterUsed()) {
      // determine current converter
      extensions = ((ExtensionFileFilter) getFileFilter()).getExtensions();
      if (m_DialogType == LOADER_DIALOG)
	newConverter = ConverterUtils.getLoaderForExtension(extensions[0]);
      else
	newConverter = ConverterUtils.getSaverForExtension(extensions[0]);

      try {
	if (m_CurrentConverter == null) {
	  m_CurrentConverter = newConverter;
	}
	else {
	  if (!m_CurrentConverter.getClass().equals(newConverter.getClass()))
	    m_CurrentConverter = newConverter;
	}
      }
      catch (Exception e) {
	m_CurrentConverter = null;
	e.printStackTrace();
      }
    }
    else {
      m_CurrentConverter = null;
    }
  }
  
  /**
   * configures the current converter.
   * 
   * @param dialogType		the type of dialog to configure for
   */
  protected void configureCurrentConverter(int dialogType) {
    String	filename;
    File	currFile;
    
    if ( (getSelectedFile() == null) || (getSelectedFile().isDirectory()) )
      return;
    
    filename = getSelectedFile().getAbsolutePath();
    
    if (m_CurrentConverter == null) {
      if (dialogType == LOADER_DIALOG)
	m_CurrentConverter = ConverterUtils.getLoaderForFile(filename);
      else if (dialogType == SAVER_DIALOG)
	m_CurrentConverter = ConverterUtils.getSaverForFile(filename);
      else
	throw new IllegalStateException("Cannot determine loader/saver!");
      
      // none found?
      if (m_CurrentConverter == null)
	return;
    }
    
    try {
      currFile = ((FileSourcedConverter) m_CurrentConverter).retrieveFile();
      if ((currFile == null) || (!currFile.getAbsolutePath().equals(filename)))
	((FileSourcedConverter) m_CurrentConverter).setFile(new File(filename));
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /**
   * For testing the file chooser.
   * 
   * @param args	the commandline options - ignored
   * @throws Exception	if something goes wrong with loading/saving
   */
  public static void main(String[] args) throws Exception {
    ConverterFileChooser	fc;
    int				retVal;
    AbstractFileLoader		loader;
    AbstractFileSaver		saver;
    Instances			data;
    
    fc     = new ConverterFileChooser();
    retVal = fc.showOpenDialog(null);
    
    // load file
    if (retVal == ConverterFileChooser.APPROVE_OPTION) {
      loader = fc.getLoader();
      data   = loader.getDataSet();
      retVal = fc.showSaveDialog(null);

      // save file
      if (retVal == ConverterFileChooser.APPROVE_OPTION) {
	saver = fc.getSaver();
	saver.setInstances(data);
	saver.writeBatch();
      }
      else {
	System.out.println("Saving aborted!");
      }
    }
    else {
      System.out.println("Loading aborted!");
    }
  }
}

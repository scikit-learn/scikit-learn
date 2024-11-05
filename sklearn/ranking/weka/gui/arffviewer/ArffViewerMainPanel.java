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
 * ArffViewerMainPanel.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.arffviewer;

import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.converters.AbstractSaver;
import weka.gui.ComponentHelper;
import weka.gui.ConverterFileChooser;
import weka.gui.JTableHelper;
import weka.gui.ListSelectorDialog;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * The main panel of the ArffViewer. It has a reference to the menu, that an
 * implementing JFrame only needs to add via the setJMenuBar(JMenuBar) method.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.7 $ 
 */

public class ArffViewerMainPanel 
  extends JPanel 
  implements ActionListener, ChangeListener {
  
  /** for serialization */
  static final long serialVersionUID = -8763161167586738753L;
  
  /** the default for width */
  public final static int    DEFAULT_WIDTH     = -1;
  /** the default for height */
  public final static int    DEFAULT_HEIGHT    = -1;
  /** the default for left */
  public final static int    DEFAULT_LEFT      = -1;
  /** the default for top */
  public final static int    DEFAULT_TOP       = -1;
  /** default width */
  public final static int    WIDTH             = 800;
  /** default height */
  public final static int    HEIGHT            = 600;

  protected Container             parent;
  protected JTabbedPane           tabbedPane;
  protected JMenuBar              menuBar;
  protected JMenu                 menuFile;
  protected JMenuItem             menuFileOpen;
  protected JMenuItem             menuFileSave;
  protected JMenuItem             menuFileSaveAs;
  protected JMenuItem             menuFileClose;
  protected JMenuItem             menuFileCloseAll;
  protected JMenuItem             menuFileProperties;
  protected JMenuItem             menuFileExit;
  protected JMenu                 menuEdit;
  protected JMenuItem             menuEditUndo;
  protected JMenuItem             menuEditCopy;
  protected JMenuItem             menuEditSearch;
  protected JMenuItem             menuEditClearSearch;
  protected JMenuItem             menuEditDeleteAttribute;
  protected JMenuItem             menuEditDeleteAttributes;
  protected JMenuItem             menuEditRenameAttribute;
  protected JMenuItem             menuEditAttributeAsClass;
  protected JMenuItem             menuEditDeleteInstance;
  protected JMenuItem             menuEditDeleteInstances;
  protected JMenuItem             menuEditSortInstances;
  protected JMenu                 menuView;
  protected JMenuItem             menuViewAttributes;
  protected JMenuItem             menuViewValues;
  protected JMenuItem             menuViewOptimalColWidths;
  
  protected ConverterFileChooser  fileChooser;
  protected String                frameTitle;
  protected boolean               confirmExit;
  protected int                   width;
  protected int                   height;
  protected int                   top;
  protected int                   left;
  protected boolean               exitOnClose;
  
  /**
   * initializes the object
   * 
   * @param parentFrame		the parent frame (JFrame or JInternalFrame)
   */
  public ArffViewerMainPanel(Container parentFrame) {
    parent     = parentFrame;
    frameTitle = "ARFF-Viewer"; 
    createPanel();
  }
  
  /**
   * creates all the components in the panel
   */
  protected void createPanel() {
    // basic setup
    setSize(WIDTH, HEIGHT);
    
    setConfirmExit(false);
    setLayout(new BorderLayout());
    
    // file dialog
    fileChooser = new ConverterFileChooser(new File(System.getProperty("user.dir")));
    fileChooser.setMultiSelectionEnabled(true);
    
    // menu
    menuBar        = new JMenuBar();
    menuFile       = new JMenu("File");
    menuFileOpen   = new JMenuItem("Open...", ComponentHelper.getImageIcon("open.gif"));
    menuFileOpen.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, KeyEvent.CTRL_MASK));
    menuFileOpen.addActionListener(this);
    menuFileSave   = new JMenuItem("Save", ComponentHelper.getImageIcon("save.gif"));
    menuFileSave.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_MASK));
    menuFileSave.addActionListener(this);
    menuFileSaveAs = new JMenuItem("Save as...", ComponentHelper.getImageIcon("empty.gif"));
    menuFileSaveAs.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK));
    menuFileSaveAs.addActionListener(this);
    menuFileClose  = new JMenuItem("Close", ComponentHelper.getImageIcon("empty.gif"));
    menuFileClose.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_W, KeyEvent.CTRL_MASK));
    menuFileClose.addActionListener(this);
    menuFileCloseAll = new JMenuItem("Close all", ComponentHelper.getImageIcon("empty.gif"));
    menuFileCloseAll.addActionListener(this);
    menuFileProperties  = new JMenuItem("Properties", ComponentHelper.getImageIcon("empty.gif"));
    menuFileProperties.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_ENTER, KeyEvent.CTRL_MASK));
    menuFileProperties.addActionListener(this);
    menuFileExit   = new JMenuItem("Exit", ComponentHelper.getImageIcon("forward.gif"));
    menuFileExit.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X, KeyEvent.ALT_MASK));
    menuFileExit.addActionListener(this);
    menuFile.add(menuFileOpen);
    menuFile.add(menuFileSave);
    menuFile.add(menuFileSaveAs);
    menuFile.add(menuFileClose);
    menuFile.add(menuFileCloseAll);
    menuFile.addSeparator();
    menuFile.add(menuFileProperties);
    menuFile.addSeparator();
    menuFile.add(menuFileExit);
    menuBar.add(menuFile);
    
    menuEdit       = new JMenu("Edit");
    menuEditUndo   = new JMenuItem("Undo", ComponentHelper.getImageIcon("undo.gif"));
    menuEditUndo.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Z, KeyEvent.CTRL_MASK));
    menuEditUndo.addActionListener(this);
    menuEditCopy   = new JMenuItem("Copy", ComponentHelper.getImageIcon("copy.gif"));
    menuEditCopy.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_INSERT, KeyEvent.CTRL_MASK));
    menuEditCopy.addActionListener(this);
    menuEditSearch   = new JMenuItem("Search...", ComponentHelper.getImageIcon("find.gif"));
    menuEditSearch.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, KeyEvent.CTRL_MASK));
    menuEditSearch.addActionListener(this);
    menuEditClearSearch   = new JMenuItem("Clear search", ComponentHelper.getImageIcon("empty.gif"));
    menuEditClearSearch.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK));
    menuEditClearSearch.addActionListener(this);
    menuEditRenameAttribute = new JMenuItem("Rename attribute", ComponentHelper.getImageIcon("empty.gif"));
    menuEditRenameAttribute.addActionListener(this);
    menuEditAttributeAsClass = new JMenuItem("Attribute as class", ComponentHelper.getImageIcon("empty.gif"));
    menuEditAttributeAsClass.addActionListener(this);
    menuEditDeleteAttribute = new JMenuItem("Delete attribute", ComponentHelper.getImageIcon("empty.gif"));
    menuEditDeleteAttribute.addActionListener(this);
    menuEditDeleteAttributes = new JMenuItem("Delete attributes", ComponentHelper.getImageIcon("empty.gif"));
    menuEditDeleteAttributes.addActionListener(this);
    menuEditDeleteInstance = new JMenuItem("Delete instance", ComponentHelper.getImageIcon("empty.gif"));
    menuEditDeleteInstance.addActionListener(this);
    menuEditDeleteInstances = new JMenuItem("Delete instances", ComponentHelper.getImageIcon("empty.gif"));
    menuEditDeleteInstances.addActionListener(this);
    menuEditSortInstances = new JMenuItem("Sort data (ascending)", ComponentHelper.getImageIcon("sort.gif"));
    menuEditSortInstances.addActionListener(this);
    menuEdit.add(menuEditUndo);
    menuEdit.addSeparator();
    menuEdit.add(menuEditCopy);
    menuEdit.addSeparator();
    menuEdit.add(menuEditSearch);
    menuEdit.add(menuEditClearSearch);
    menuEdit.addSeparator();
    menuEdit.add(menuEditRenameAttribute);
    menuEdit.add(menuEditAttributeAsClass);
    menuEdit.add(menuEditDeleteAttribute);
    menuEdit.add(menuEditDeleteAttributes);
    menuEdit.addSeparator();
    menuEdit.add(menuEditDeleteInstance);
    menuEdit.add(menuEditDeleteInstances);
    menuEdit.add(menuEditSortInstances);
    menuBar.add(menuEdit);
    
    menuView       = new JMenu("View");
    menuViewAttributes   = new JMenuItem("Attributes...", ComponentHelper.getImageIcon("objects.gif"));
    menuViewAttributes.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK));
    menuViewAttributes.addActionListener(this);
    menuViewValues   = new JMenuItem("Values...", ComponentHelper.getImageIcon("properties.gif"));
    menuViewValues.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK));
    menuViewValues.addActionListener(this);
    menuViewOptimalColWidths = new JMenuItem("Optimal column width (all)", ComponentHelper.getImageIcon("resize.gif"));
    menuViewOptimalColWidths.addActionListener(this);
    menuView.add(menuViewAttributes);
    menuView.add(menuViewValues);
    menuView.addSeparator();
    menuView.add(menuViewOptimalColWidths);
    menuBar.add(menuView);
    
    // tabbed pane
    tabbedPane  = new JTabbedPane();
    tabbedPane.addChangeListener(this);
    add(tabbedPane, BorderLayout.CENTER);
    
    updateMenu();
    updateFrameTitle();
  }
  
  /**
   * returns the parent frame, if it's a JFrame, otherwise null
   * 
   * @return		the parent frame
   */
  public JFrame getParentFrame() {
    if (parent instanceof JFrame)
      return (JFrame) parent;
    else
      return null;
  }
  
  /**
   * returns the parent frame, if it's a JInternalFrame, otherwise null
   * 
   * @return		the parent frame
   */
  public JInternalFrame getParentInternalFrame() {
    if (parent instanceof JInternalFrame)
      return (JInternalFrame) parent;
    else
      return null;
  }
  
  /**
   * sets the new parent frame
   * 
   * @param value	the parent frame
   */
  public void setParent(Container value) {
    parent = value;
  }
  
  /**
   * returns the menu bar to be added in a frame
   * 
   * @return		the menu bar
   */
  public JMenuBar getMenu() {
    return menuBar;
  }
  
  /**
   * returns the tabbedpane instance
   * 
   * @return		the tabbed pane
   */
  public JTabbedPane getTabbedPane() {
    return tabbedPane;
  }
  
  /**
   * whether to present a MessageBox on Exit or not
   * @param confirm           whether a MessageBox pops up or not to confirm
   *                          exit
   */
  public void setConfirmExit(boolean confirm) {
    confirmExit = confirm;
  }
  
  /**
   * returns the setting of whether to display a confirm messagebox or not
   * on exit
   * @return                  whether a messagebox is displayed or not
   */
  public boolean getConfirmExit() {
    return confirmExit;
  }

  /**
   * whether to do a System.exit(0) on close
   * 
   * @param value	enables/disables a System.exit(0) on close
   */
  public void setExitOnClose(boolean value) {
    exitOnClose = value;
  }

  /**
   * returns TRUE if a System.exit(0) is done on a close
   * 
   * @return		true if a System.exit(0) is done on close
   */
  public boolean getExitOnClose() {
    return exitOnClose;
  }
  
  /**
   * validates and repaints the frame
   */
  public void refresh() {
    validate();
    repaint();
  }
  
  /**
   * returns the title (incl. filename) for the frame
   * 
   * @return		the frame title
   */
  public String getFrameTitle() {
    if (getCurrentFilename().equals(""))
      return frameTitle;
    else
      return frameTitle + " - " + getCurrentFilename();
  }
  
  /**
   * sets the title of the parent frame, if one was provided
   */
  public void updateFrameTitle() {
    if (getParentFrame() != null)
      getParentFrame().setTitle(getFrameTitle());
    if (getParentInternalFrame() != null)
      getParentInternalFrame().setTitle(getFrameTitle());
  }
  
  /**
   * sets the enabled/disabled state of the menu 
   */
  protected void updateMenu() {
    boolean       fileOpen;
    boolean       isChanged;
    boolean       canUndo;
    
    fileOpen  = (getCurrentPanel() != null);
    isChanged = fileOpen && (getCurrentPanel().isChanged());
    canUndo   = fileOpen && (getCurrentPanel().canUndo());
    
    // File
    menuFileOpen.setEnabled(true);
    menuFileSave.setEnabled(isChanged);
    menuFileSaveAs.setEnabled(fileOpen);
    menuFileClose.setEnabled(fileOpen);
    menuFileCloseAll.setEnabled(fileOpen);
    menuFileProperties.setEnabled(fileOpen);
    menuFileExit.setEnabled(true);
    // Edit
    menuEditUndo.setEnabled(canUndo);
    menuEditCopy.setEnabled(fileOpen);
    menuEditSearch.setEnabled(fileOpen);
    menuEditClearSearch.setEnabled(fileOpen);
    menuEditAttributeAsClass.setEnabled(fileOpen);
    menuEditRenameAttribute.setEnabled(fileOpen);
    menuEditDeleteAttribute.setEnabled(fileOpen);
    menuEditDeleteAttributes.setEnabled(fileOpen);
    menuEditDeleteInstance.setEnabled(fileOpen);
    menuEditDeleteInstances.setEnabled(fileOpen);
    menuEditSortInstances.setEnabled(fileOpen);
    // View
    menuViewAttributes.setEnabled(fileOpen);
    menuViewValues.setEnabled(fileOpen);
    menuViewOptimalColWidths.setEnabled(fileOpen);
  }
  
  /**
   * sets the title of the tab that contains the given component
   * 
   * @param component		the component to set the title for
   */
  protected void setTabTitle(JComponent component) {
    int            index;
    
    if (!(component instanceof ArffPanel))
      return;
    
    index = tabbedPane.indexOfComponent(component);
    if (index == -1)
      return;
    
    tabbedPane.setTitleAt(index, ((ArffPanel) component).getTitle());
    updateFrameTitle();
  }
  
  /**
   * returns the number of panels currently open
   * 
   * @return		the number of open panels
   */
  public int getPanelCount() {
    return tabbedPane.getTabCount();
  }
  
  /**
   * returns the specified panel, <code>null</code> if index is out of bounds  
   * 
   * @param index	the index of the panel
   * @return		the panel
   */
  public ArffPanel getPanel(int index) {
    if ((index >= 0) && (index < getPanelCount()))
      return (ArffPanel) tabbedPane.getComponentAt(index);
    else
      return null;
  }
  
  /**
   * returns the currently selected tab index
   * 
   * @return		the index of the currently selected tab
   */
  public int getCurrentIndex() {
    return tabbedPane.getSelectedIndex();
  }
  
  /**
   * returns the currently selected panel
   * 
   * @return		the currently selected panel
   */
  public ArffPanel getCurrentPanel() {
    return getPanel(getCurrentIndex());
  }
  
  /**
   * checks whether a panel is currently selected
   * 
   * @return		true if a panel is currently selected
   */
  public boolean isPanelSelected() {
    return (getCurrentPanel() != null);
  }
  
  /**
   * returns the filename of the specified panel 
   * 
   * @param index	the index of the panel
   * @return		the filename for the panel
   */
  public String getFilename(int index) {
    String            result;
    ArffPanel         panel;
    
    result = "";
    panel  = getPanel(index);
    
    if (panel != null)
      result = panel.getFilename();
    
    return result;
  }
  
  /**
   * returns the filename of the current tab
   * 
   * @return		the current filename
   */
  public String getCurrentFilename() {
    return getFilename(getCurrentIndex());
  }
  
  /**
   * sets the filename of the specified panel
   * 
   * @param index	the index of the panel
   * @param filename	the new filename
   */
  public void setFilename(int index, String filename) {
    ArffPanel         panel;
    
    panel = getPanel(index);
    
    if (panel != null) {
      panel.setFilename(filename);
      setTabTitle(panel);
    }
  }
  
  /**
   * sets the filename of the current tab
   * 
   * @param filename	the new filename
   */
  public void setCurrentFilename(String filename) {
    setFilename(getCurrentIndex(), filename);
  }
  
  /**
   * if the file is changed it pops up a dialog whether to change the
   * settings. if the project wasn't changed or saved it returns TRUE
   * 
   * @return 		true if project wasn't changed or saved
   */
  protected boolean saveChanges() {
    return saveChanges(true);
  }
  
  /**
   * if the file is changed it pops up a dialog whether to change the
   * settings. if the project wasn't changed or saved it returns TRUE
   * 
   * @param showCancel	whether we have YES/NO/CANCEL or only YES/NO
   * @return 		true if project wasn't changed or saved
   */
  protected boolean saveChanges(boolean showCancel) {
    int            button;
    boolean        result;
    
    if (!isPanelSelected())
      return true;
    
    result = !getCurrentPanel().isChanged();
    
    if (getCurrentPanel().isChanged()) {
      try {
        if (showCancel)
          button = ComponentHelper.showMessageBox(
              this,
              "Changed",
              "The file is not saved - Do you want to save it?",
              JOptionPane.YES_NO_CANCEL_OPTION,
              JOptionPane.QUESTION_MESSAGE );
        else
          button = ComponentHelper.showMessageBox(
              this,
              "Changed",
              "The file is not saved - Do you want to save it?",
              JOptionPane.YES_NO_OPTION,
              JOptionPane.QUESTION_MESSAGE );
      }
      catch (Exception e) {
        button = JOptionPane.CANCEL_OPTION; 
      }
      
      switch (button) {
        case JOptionPane.YES_OPTION:
          saveFile();
          result = !getCurrentPanel().isChanged();
          break;
        case JOptionPane.NO_OPTION:
          result = true;
          break;
        case JOptionPane.CANCEL_OPTION: 
          result = false;
          break;
      }
    }
    
    return result;
  }

  /**
   * loads the specified file
   * 
   * @param filename	the file to load
   */
  public void loadFile(String filename) {
    ArffPanel         panel;

    panel    = new ArffPanel(filename);
    panel.addChangeListener(this);
    tabbedPane.addTab(panel.getTitle(), panel);
    tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
  }
  
  /**
   * loads the specified file into the table
   */
  public void loadFile() {
    int               retVal;
    int               i;
    String            filename;
    
    retVal = fileChooser.showOpenDialog(this);
    if (retVal != ConverterFileChooser.APPROVE_OPTION)
      return;
    
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    
    for (i = 0; i< fileChooser.getSelectedFiles().length; i++) {
      filename = fileChooser.getSelectedFiles()[i].getAbsolutePath();
      loadFile(filename);
    }
    
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
  }
  
  /**
   * saves the current data into a file
   */
  public void saveFile() {
    ArffPanel           panel;
    String              filename;
    AbstractSaver       saver;
    
    // no panel? -> exit
    panel = getCurrentPanel();
    if (panel == null)
      return;
    
    filename = panel.getFilename();
    if (filename.equals(ArffPanel.TAB_INSTANCES)) {
      saveFileAs();
    }
    else {
      saver = fileChooser.getSaver();
      try {
	saver.setInstances(panel.getInstances());
	saver.writeBatch();
	panel.setChanged(false);
	setCurrentFilename(filename);
      }
      catch (Exception e) {
	e.printStackTrace();
      }
    }
  }
  
  /**
   * saves the current data into a new file
   */
  public void saveFileAs() {
    int                  retVal;
    ArffPanel            panel;
    
    // no panel? -> exit
    panel = getCurrentPanel();
    if (panel == null) {
      System.out.println("nothing selected!");
      return;
    }
    
    if (!getCurrentFilename().equals("")) {
      try {
        fileChooser.setSelectedFile(new File(getCurrentFilename()));
      }
      catch (Exception e) {
        // ignore
      }
    }
    
    // set filter for savers
    try {
      fileChooser.setCapabilitiesFilter(Capabilities.forInstances(panel.getInstances()));
    }
    catch (Exception e) {
      fileChooser.setCapabilitiesFilter(null);
    }

    retVal = fileChooser.showSaveDialog(this);
    if (retVal != ConverterFileChooser.APPROVE_OPTION)
      return;
    
    panel.setChanged(false);
    setCurrentFilename(fileChooser.getSelectedFile().getAbsolutePath());
    saveFile();
  }
  
  /**
   * closes the current tab
   */
  public void closeFile() {
    closeFile(true);
  }
  
  /**
   * closes the current tab
   * @param showCancel           whether to show an additional CANCEL button
   *                             in the "Want to save changes"-dialog
   * @see                        #saveChanges(boolean)
   */
  public void closeFile(boolean showCancel) {
    if (getCurrentIndex() == -1)
      return;
    
    if (!saveChanges(showCancel))
      return;
    
    tabbedPane.removeTabAt(getCurrentIndex());
    updateFrameTitle();
    System.gc();
  }
  
  /**
   * closes all open files
   */
  public void closeAllFiles() {
    while (tabbedPane.getTabCount() > 0) {
      if (!saveChanges(true))
        return;
      
      tabbedPane.removeTabAt(getCurrentIndex());
      updateFrameTitle();
      System.gc();
    }
  }
  
  /**
   * displays some properties of the instances
   */
  public void showProperties() {
    ArffPanel             panel;
    ListSelectorDialog    dialog;
    Vector                props;
    Instances             inst;
    
    panel = getCurrentPanel();
    if (panel == null)
      return;
    
    inst  = panel.getInstances();
    if (inst == null)
      return;
    if (inst.classIndex() < 0)
      inst.setClassIndex(inst.numAttributes() - 1);
    
    // get some data
    props = new Vector();
    props.add("Filename: " + panel.getFilename());
    props.add("Relation name: " + inst.relationName());
    props.add("# of instances: " + inst.numInstances());
    props.add("# of attributes: " + inst.numAttributes());
    props.add("Class attribute: " + inst.classAttribute().name());
    props.add("# of class labels: " + inst.numClasses());
    
    dialog = new ListSelectorDialog(getParentFrame(), new JList(props));
    dialog.showDialog();
  }
  
  /**
   * closes the window, i.e., if the parent is not null and implements
   * the WindowListener interface it calls the windowClosing method
   */
  public void close() {
    if (getParentInternalFrame() != null)
      getParentInternalFrame().doDefaultCloseAction();
    else if (getParentFrame() != null)
      ((Window) getParentFrame()).dispatchEvent(
	  new WindowEvent(
	      (Window) getParentFrame(), WindowEvent.WINDOW_CLOSING));
  }
  
  /**
   * undoes the last action 
   */
  public void undo() {
    if (!isPanelSelected())
      return;
    
    getCurrentPanel().undo();
  }
  
  /**
   * copies the content of the selection to the clipboard
   */
  public void copyContent() {
    if (!isPanelSelected())
      return;
    
    getCurrentPanel().copyContent();
  }
  
  /**
   * searches for a string in the cells
   */
  public void search() {
    if (!isPanelSelected())
      return;

    getCurrentPanel().search();
  }
  
  /**
   * clears the search, i.e. resets the found cells
   */
  public void clearSearch() {
    if (!isPanelSelected())
      return;

    getCurrentPanel().clearSearch();
  }
  
  /**
   * renames the current selected Attribute
   */
  public void renameAttribute() {
    if (!isPanelSelected())
      return;
    
    getCurrentPanel().renameAttribute();
  }
  
  /**
   * sets the current selected Attribute as class attribute, i.e. it moves it
   * to the end of the attributes
   */
  public void attributeAsClass() {
    if (!isPanelSelected())
      return;
    
    getCurrentPanel().attributeAsClass();
  }
  
  /**
   * deletes the current selected Attribute or several chosen ones
   * 
   * @param multiple	whether to delete myultiple attributes
   */
  public void deleteAttribute(boolean multiple) {
    if (!isPanelSelected())
      return;
    
    if (multiple)
      getCurrentPanel().deleteAttributes();
    else
      getCurrentPanel().deleteAttribute();
  }
  
  /**
   * deletes the current selected Instance or several chosen ones
   * 
   * @param multiple		whether to delete multiple instances
   */
  public void deleteInstance(boolean multiple) {
    if (!isPanelSelected())
      return;
    
    if (multiple)
      getCurrentPanel().deleteInstances();
    else
      getCurrentPanel().deleteInstance();
  }
  
  /**
   * sorts the current selected attribute
   */
  public void sortInstances() {
    if (!isPanelSelected())
      return;
    
    getCurrentPanel().sortInstances();
  }
  
  /**
   * displays all the attributes, returns the selected item or NULL if canceled
   * 
   * @return		the name of the selected attribute
   */
  public String showAttributes() {
    ArffSortedTableModel     model;
    ListSelectorDialog  dialog;
    int                 i;
    JList               list;
    String              name;
    int                 result;
    
    if (!isPanelSelected())
      return null;
    
    list   = new JList(getCurrentPanel().getAttributes());
    dialog = new ListSelectorDialog(getParentFrame(), list);
    result = dialog.showDialog();
    
    if (result == ListSelectorDialog.APPROVE_OPTION) {
      model = (ArffSortedTableModel) getCurrentPanel().getTable().getModel();
      name  = list.getSelectedValue().toString();
      i     = model.getAttributeColumn(name);
      JTableHelper.scrollToVisible(getCurrentPanel().getTable(), 0, i);
      getCurrentPanel().getTable().setSelectedColumn(i);
      return name;
    }
    else {
      return null;
    }
  }
  
  /**
   * displays all the distinct values for an attribute
   */
  public void showValues() {
    String                attribute;
    ArffSortedTableModel       model;
    ArffTable             table;
    HashSet               values;
    Vector                items;
    Iterator              iter;
    ListSelectorDialog    dialog;
    int                   i;
    int                   col;
    
    // choose attribute to retrieve values for
    attribute = showAttributes();
    if (attribute == null)
      return;
    
    table  = (ArffTable) getCurrentPanel().getTable();
    model  = (ArffSortedTableModel) table.getModel();
    
    // get column index
    col    = -1;
    for (i = 0; i < table.getColumnCount(); i++) {
      if (table.getPlainColumnName(i).equals(attribute)) {
        col = i;
        break;
      }
    }
    // not found?
    if (col == -1)
      return;
    
    // get values
    values = new HashSet();
    items  = new Vector();
    for (i = 0; i < model.getRowCount(); i++)
      values.add(model.getValueAt(i, col).toString());
    if (values.isEmpty())
      return;
    iter = values.iterator();
    while (iter.hasNext())
      items.add(iter.next());
    Collections.sort(items);
    
    dialog = new ListSelectorDialog(getParentFrame(), new JList(items));
    dialog.showDialog();
  }
  
  /**
   * sets the optimal column width for all columns
   */
  public void setOptimalColWidths() {
    if (!isPanelSelected())
      return;
    
    getCurrentPanel().setOptimalColWidths();
  }
  
  /**
   * invoked when an action occurs
   * 
   * @param e		the action event
   */
  public void actionPerformed(ActionEvent e) {
    Object          o;
    
    o = e.getSource();
    
    if (o == menuFileOpen)
      loadFile();
    else if (o == menuFileSave)
      saveFile();
    else if (o == menuFileSaveAs)
      saveFileAs();
    else if (o == menuFileClose)
      closeFile();
    else if (o == menuFileCloseAll)
      closeAllFiles();
    else if (o == menuFileProperties)
      showProperties();
    else if (o == menuFileExit)
      close();
    else if (o == menuEditUndo)
      undo();
    else if (o == menuEditCopy)
      copyContent();
    else if (o == menuEditSearch)
      search();
    else if (o == menuEditClearSearch)
      clearSearch();
    else if (o == menuEditDeleteAttribute)
      deleteAttribute(false);
    else if (o == menuEditDeleteAttributes)
      deleteAttribute(true);
    else if (o == menuEditRenameAttribute)
      renameAttribute();
    else if (o == menuEditAttributeAsClass)
      attributeAsClass();
    else if (o == menuEditDeleteInstance)
      deleteInstance(false);
    else if (o == menuEditDeleteInstances)
      deleteInstance(true);
    else if (o == menuEditSortInstances)
      sortInstances();
    else if (o == menuViewAttributes)
      showAttributes();
    else if (o == menuViewValues)
      showValues();
    else if (o == menuViewOptimalColWidths)
      setOptimalColWidths();
    
    updateMenu();
  }
  
  /**
   * Invoked when the target of the listener has changed its state.
   * 
   * @param e		the change event
   */
  public void stateChanged(ChangeEvent e) {
    updateFrameTitle();
    updateMenu();
    
    // did the content of panel change? -> change title of tab
    if (e.getSource() instanceof JComponent)
      setTabTitle((JComponent) e.getSource());
  }
  
  /**
   * returns only the classname
   * 
   * @return		the classname
   */
  public String toString() {
    return this.getClass().getName();
  }
}

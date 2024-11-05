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
 * ArffPanel.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.arffviewer;

import weka.core.Instances;
import weka.core.Undoable;
import weka.core.Utils;
import weka.gui.ComponentHelper;
import weka.gui.JTableHelper;
import weka.gui.ListSelectorDialog;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TableModelEvent;

/**
 * A Panel representing an ARFF-Table and the associated filename.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6455 $ 
 */

public class ArffPanel 
  extends JPanel
  implements ActionListener, ChangeListener, MouseListener, Undoable {
  
  /** for serialization */
  static final long serialVersionUID = -4697041150989513939L;
  
  /** the name of the tab for instances that were set directly */
  public final static String TAB_INSTANCES = "Instances";

  /** the underlying table */
  private ArffTable m_TableArff;
  /** the popup menu for the header row */
  private JPopupMenu m_PopupHeader;
  /** the popup menu for the data rows */
  private JPopupMenu m_PopupRows;
  /** displays the relation name */
  private JLabel m_LabelName;
  /** whether to display the attribute index in the table header. */
  private boolean m_ShowAttributeIndex;

  // menu items
  private JMenuItem             menuItemMean;
  private JMenuItem             menuItemSetAllValues;
  private JMenuItem             menuItemSetMissingValues;
  private JMenuItem             menuItemReplaceValues;
  private JMenuItem             menuItemRenameAttribute;
  private JMenuItem             menuItemAttributeAsClass;
  private JMenuItem             menuItemDeleteAttribute;
  private JMenuItem             menuItemDeleteAttributes;
  private JMenuItem             menuItemSortInstances;
  private JMenuItem             menuItemDeleteSelectedInstance;
  private JMenuItem             menuItemDeleteAllSelectedInstances;
  private JMenuItem             menuItemSearch;
  private JMenuItem             menuItemClearSearch;
  private JMenuItem             menuItemUndo;
  private JMenuItem             menuItemCopy;
  private JMenuItem             menuItemOptimalColWidth;
  private JMenuItem             menuItemOptimalColWidths;
  
  /** the filename used in the title */
  private String m_Filename;
  /** the title prefix */
  private String m_Title;
  /** the currently selected column */
  private int m_CurrentCol;
  /** flag for whether data got changed */
  private boolean m_Changed;
  /** the listeners that listen for modifications */
  private HashSet m_ChangeListeners;
  /** the string used in the last search */
  private String m_LastSearch;
  /** the string used in the last replace */
  private String m_LastReplace;
  
  /**
   * initializes the panel with no data
   */
  public ArffPanel() {
    super();
    
    initialize();
    createPanel();
  }
  
  /**
   * initializes the panel and loads the specified file
   * 
   * @param filename	the file to load
   */
  public ArffPanel(String filename) {
    this();
    
    loadFile(filename);
  }
  
  /**
   * initializes the panel with the given data
   * 
   * @param data	the data to use
   */
  public ArffPanel(Instances data) {
    this();
    
    m_Filename = "";
    
    setInstances(data);
  }
  
  /**
   * any member variables are initialized here
   */
  protected void initialize() {
    m_Filename           = "";
    m_Title              = "";
    m_CurrentCol         = -1;
    m_LastSearch         = "";
    m_LastReplace        = "";
    m_ShowAttributeIndex = true;
    m_Changed            = false;
    m_ChangeListeners    = new HashSet();
  }
  
  /**
   * creates all the components in the frame
   */
  protected void createPanel() {
    JScrollPane                pane;
    
    setLayout(new BorderLayout());
    
    menuItemMean = new JMenuItem("Get mean...");
    menuItemMean.addActionListener(this);
    menuItemSetAllValues = new JMenuItem("Set all values to...");
    menuItemSetAllValues.addActionListener(this);
    menuItemSetMissingValues = new JMenuItem("Set missing values to...");
    menuItemSetMissingValues.addActionListener(this);
    menuItemReplaceValues = new JMenuItem("Replace values with...");
    menuItemReplaceValues.addActionListener(this);
    menuItemRenameAttribute = new JMenuItem("Rename attribute...");
    menuItemRenameAttribute.addActionListener(this);
    menuItemAttributeAsClass = new JMenuItem("Attribute as class");
    menuItemAttributeAsClass.addActionListener(this);
    menuItemDeleteAttribute = new JMenuItem("Delete attribute");
    menuItemDeleteAttribute.addActionListener(this);
    menuItemDeleteAttributes = new JMenuItem("Delete attributes...");
    menuItemDeleteAttributes.addActionListener(this);
    menuItemSortInstances = new JMenuItem("Sort data (ascending)");
    menuItemSortInstances.addActionListener(this);
    menuItemOptimalColWidth = new JMenuItem("Optimal column width (current)");
    menuItemOptimalColWidth.addActionListener(this);
    menuItemOptimalColWidths = new JMenuItem("Optimal column width (all)");
    menuItemOptimalColWidths.addActionListener(this);

    // row popup
    menuItemUndo = new JMenuItem("Undo");
    menuItemUndo.addActionListener(this);
    menuItemCopy = new JMenuItem("Copy");
    menuItemCopy.addActionListener(this);
    menuItemSearch = new JMenuItem("Search...");
    menuItemSearch.addActionListener(this);
    menuItemClearSearch = new JMenuItem("Clear search");
    menuItemClearSearch.addActionListener(this);
    menuItemDeleteSelectedInstance = new JMenuItem("Delete selected instance");
    menuItemDeleteSelectedInstance.addActionListener(this);
    menuItemDeleteAllSelectedInstances = new JMenuItem("Delete ALL selected instances");
    menuItemDeleteAllSelectedInstances.addActionListener(this);
    
    // table
    m_TableArff = new ArffTable();
    m_TableArff.setToolTipText("Right click (or left+alt) for context menu");
    m_TableArff.getTableHeader().addMouseListener(this);
    m_TableArff.getTableHeader().setToolTipText("<html><b>Sort view:</b> left click = ascending / Shift + left click = descending<br><b>Menu:</b> right click (or left+alt)</html>");
    m_TableArff.getTableHeader().setDefaultRenderer(new ArffTableCellRenderer());
    m_TableArff.addChangeListener(this);
    m_TableArff.addMouseListener(this);
    pane = new JScrollPane(m_TableArff);
    add(pane, BorderLayout.CENTER);
    
    // relation name
    m_LabelName   = new JLabel();
    add(m_LabelName, BorderLayout.NORTH);
  }

  /**
   * initializes the popup menus
   */
  private void initPopupMenus() {
    // header popup
    m_PopupHeader  = new JPopupMenu();
    m_PopupHeader.addMouseListener(this);
    m_PopupHeader.add(menuItemMean);
    if (!isReadOnly()) {
      m_PopupHeader.addSeparator();
      m_PopupHeader.add(menuItemSetAllValues);
      m_PopupHeader.add(menuItemSetMissingValues);
      m_PopupHeader.add(menuItemReplaceValues);
      m_PopupHeader.addSeparator();
      m_PopupHeader.add(menuItemRenameAttribute);
      m_PopupHeader.add(menuItemAttributeAsClass);
      m_PopupHeader.add(menuItemDeleteAttribute);
      m_PopupHeader.add(menuItemDeleteAttributes);
      m_PopupHeader.add(menuItemSortInstances);
    }
    m_PopupHeader.addSeparator();
    m_PopupHeader.add(menuItemOptimalColWidth);
    m_PopupHeader.add(menuItemOptimalColWidths);
    
    // row popup
    m_PopupRows = new JPopupMenu();
    m_PopupRows.addMouseListener(this);
    if (!isReadOnly()) {
      m_PopupRows.add(menuItemUndo);
      m_PopupRows.addSeparator();
    }
    m_PopupRows.add(menuItemCopy);
    m_PopupRows.addSeparator();
    m_PopupRows.add(menuItemSearch);
    m_PopupRows.add(menuItemClearSearch);
    if (!isReadOnly()) {
      m_PopupRows.addSeparator();
      m_PopupRows.add(menuItemDeleteSelectedInstance);
      m_PopupRows.add(menuItemDeleteAllSelectedInstances);
    }
  }
  
  /**
   * sets the enabled/disabled state of the menu items 
   */
  private void setMenu() {
    boolean			isNumeric;
    boolean			hasColumns;
    boolean			hasRows;
    boolean			attSelected;
    ArffSortedTableModel	model;
    boolean			isNull;
    
    model       = (ArffSortedTableModel) m_TableArff.getModel();
    isNull      = (model.getInstances() == null);
    hasColumns  = !isNull && (model.getInstances().numAttributes() > 0);
    hasRows     = !isNull && (model.getInstances().numInstances() > 0);
    attSelected = hasColumns && (m_CurrentCol > 0);
    isNumeric   = attSelected && (model.getAttributeAt(m_CurrentCol).isNumeric());
    
    menuItemUndo.setEnabled(canUndo());
    menuItemCopy.setEnabled(true);
    menuItemSearch.setEnabled(true);
    menuItemClearSearch.setEnabled(true);
    menuItemMean.setEnabled(isNumeric);
    menuItemSetAllValues.setEnabled(attSelected);
    menuItemSetMissingValues.setEnabled(attSelected);
    menuItemReplaceValues.setEnabled(attSelected);
    menuItemRenameAttribute.setEnabled(attSelected);
    menuItemDeleteAttribute.setEnabled(attSelected);
    menuItemDeleteAttributes.setEnabled(attSelected);
    menuItemSortInstances.setEnabled(hasRows && attSelected);
    menuItemDeleteSelectedInstance.setEnabled(hasRows && m_TableArff.getSelectedRow() > -1);
    menuItemDeleteAllSelectedInstances.setEnabled(hasRows && (m_TableArff.getSelectedRows().length > 0));
  }
  
  /**
   * returns the table component
   * 
   * @return 		the table
   */
  public ArffTable getTable() {
    return m_TableArff;
  }
  
  /**
   * returns the title for the Tab, i.e. the filename
   * 
   * @return 		the title for the tab
   */
  public String getTitle() {
    return m_Title;
  }
  
  /**
   * returns the filename
   * 
   * @return		the filename
   */
  public String getFilename() {
    return m_Filename;
  }
  
  /**
   * sets the filename
   * 
   * @param filename	the new filename
   */
  public void setFilename(String filename) {
    m_Filename = filename;
    createTitle();
  }
  
  /**
   * returns the instances of the panel, if none then NULL
   * 
   * @return		the instances of the panel
   */
  public Instances getInstances() {
    Instances            result;
    
    result = null;
    
    if (m_TableArff.getModel() != null)
      result = ((ArffSortedTableModel) m_TableArff.getModel()).getInstances();
    
    return result;
  }
  
  /**
   * displays the given instances, i.e. creates a tab with the title 
   * TAB_INSTANCES. if one already exists it closes it.<br>
   * if a different instances object is used here, don't forget to clear
   * the undo-history by calling <code>clearUndo()</code>
   * 
   * @param data	the instances to display
   * @see               #TAB_INSTANCES
   * @see               #clearUndo()
   */
  public void setInstances(Instances data) {
    ArffSortedTableModel         model;
    
    m_Filename = TAB_INSTANCES;
    
    createTitle();
    model = new ArffSortedTableModel(data);
    model.setShowAttributeIndex(m_ShowAttributeIndex);
    
    m_TableArff.setModel(model);
    clearUndo();
    setChanged(false);
    createName();
  }
  
  /**
   * returns a list with the attributes
   * 
   * @return		a list of the attributes
   */
  public Vector getAttributes() {
    Vector               result;
    int                  i;
    
    result = new Vector();
    for (i = 0; i < getInstances().numAttributes(); i++)
      result.add(getInstances().attribute(i).name());
    Collections.sort(result);
    
    return result;
  }
  
  /**
   * can only reset the changed state to FALSE
   * 
   * @param changed		if false, resets the changed state
   */
  public void setChanged(boolean changed) {
    if (!changed) {
      this.m_Changed = changed;
      createTitle();
    }
  }
  
  /**
   * returns whether the content of the panel was changed
   * 
   * @return		true if the content was changed
   */
  public boolean isChanged() {
    return m_Changed;
  }

  /**
   * returns whether the model is read-only
   * 
   * @return 		true if model is read-only
   */
  public boolean isReadOnly() {
    if (m_TableArff == null)
      return true;
    else
      return ((ArffSortedTableModel) m_TableArff.getModel()).isReadOnly();
  }
  
  /**
   * sets whether the model is read-only
   * 
   * @param value	if true the model is set to read-only
   */
  public void setReadOnly(boolean value) {
    if (m_TableArff != null)
      ((ArffSortedTableModel) m_TableArff.getModel()).setReadOnly(value);
  }

  /**
   * Sets whether to display the attribute index in the header.
   * 
   * @param value	if true then the attribute indices are displayed in the
   * 			table header
   */
  public void setShowAttributeIndex(boolean value) {
    m_ShowAttributeIndex = value;
    if (m_TableArff != null)
      ((ArffSortedTableModel) m_TableArff.getModel()).setShowAttributeIndex(value);
  }
  
  /**
   * Returns whether to display the attribute index in the header.
   * 
   * @return		true if the attribute indices are displayed in the
   * 			table header
   */
  public boolean getShowAttributeIndex() {
    return m_ShowAttributeIndex;
  }
  
  /**
   * returns whether undo support is enabled
   * 
   * @return 		true if undo is enabled
   */
  public boolean isUndoEnabled() {
    return ((ArffSortedTableModel) m_TableArff.getModel()).isUndoEnabled();
  }
  
  /**
   * sets whether undo support is enabled
   * 
   * @param enabled		whether to enable/disable undo support
   */
  public void setUndoEnabled(boolean enabled) {
    ((ArffSortedTableModel) m_TableArff.getModel()).setUndoEnabled(enabled);
  }
  
  /**
   * removes the undo history
   */
  public void clearUndo() {
    ((ArffSortedTableModel) m_TableArff.getModel()).clearUndo();
  }
  
  /**
   * returns whether an undo is possible 
   * 
   * @return		true if undo is possible
   */
  public boolean canUndo() {
    return ((ArffSortedTableModel) m_TableArff.getModel()).canUndo();
  }
  
  /**
   * performs an undo action
   */
  public void undo() {
    if (canUndo()) {
      ((ArffSortedTableModel) m_TableArff.getModel()).undo();
      
      // notify about update
      notifyListener();
    }
  }
  
  /**
   * adds the current state of the instances to the undolist 
   */
  public void addUndoPoint() {
    ((ArffSortedTableModel) m_TableArff.getModel()).addUndoPoint();
        
    // update menu
    setMenu();
  }
  
  /**
   * sets the title (i.e. filename)
   */
  private void createTitle() {
    File              file;
    
    if (m_Filename.equals("")) {
      m_Title = "-none-";
    }
    else if (m_Filename.equals(TAB_INSTANCES)) {
      m_Title = TAB_INSTANCES;
    }
    else {
      try {
        file  = new File(m_Filename);
        m_Title = file.getName();
      }
      catch (Exception e) {
        m_Title = "-none-";
      }
    }
    
    if (isChanged())
      m_Title += " *";
  }
  
  /**
   * sets the relation name
   */
  private void createName() {
    ArffSortedTableModel         model;
    
    model = (ArffSortedTableModel) m_TableArff.getModel();
    if ((model != null) && (model.getInstances() != null))
      m_LabelName.setText("Relation: " + model.getInstances().relationName());
    else
      m_LabelName.setText("");
  }
  
  /**
   * loads the specified file into the table
   * 
   * @param filename		the file to load
   */
  private void loadFile(String filename) {
    ArffSortedTableModel         model;
    
    this.m_Filename = filename;
    
    createTitle();
    
    if (filename.equals("")) {
      model = null;
    }
    else {
      model = new ArffSortedTableModel(filename);
      model.setShowAttributeIndex(getShowAttributeIndex());
    }
    
    m_TableArff.setModel(model);
    setChanged(false);
    createName();
  }
  
  /**
   * calculates the mean of the given numeric column
   */
  private void calcMean() {
    ArffSortedTableModel   model;
    int               i;
    double            mean;
    
    // no column selected?
    if (m_CurrentCol == -1)
      return;
    
    model = (ArffSortedTableModel) m_TableArff.getModel();
    
    // not numeric?
    if (!model.getAttributeAt(m_CurrentCol).isNumeric())
      return;
    
    mean = 0;
    for (i = 0; i < model.getRowCount(); i++)
      mean += model.getInstances().instance(i).value(m_CurrentCol - 1);
    mean = mean / model.getRowCount();
    
    // show result
    ComponentHelper.showMessageBox(
        getParent(), 
        "Mean for attribute...",
        "Mean for attribute '" 
        + m_TableArff.getPlainColumnName(m_CurrentCol) 
        + "':\n\t" + Utils.doubleToString(mean, 3),
        JOptionPane.OK_CANCEL_OPTION,
        JOptionPane.PLAIN_MESSAGE);
  }
  
  /**
   * sets the specified values in a column to a new value
   * 
   * @param o		the menu item
   */
  private void setValues(Object o) {
    String                     msg;
    String                     title;
    String                     value;
    String                     valueNew;
    int                        i;
    ArffSortedTableModel      model;
    
    value    = "";
    valueNew = "";
    
    if (o == menuItemSetMissingValues) {
      title = "Replace missing values..."; 
      msg   = "New value for MISSING values";
    }
    else if (o == menuItemSetAllValues) {
      title = "Set all values..."; 
      msg   = "New value for ALL values";
    }
    else if (o == menuItemReplaceValues) {
      title = "Replace values..."; 
      msg   = "Old value";
    }
    else
      return;
    
    value = ComponentHelper.showInputBox(m_TableArff.getParent(), title, msg, m_LastSearch);
    
    // cancelled?
    if (value == null)
      return;

    m_LastSearch = value;
    
    // replacement
    if (o == menuItemReplaceValues) {
      valueNew = ComponentHelper.showInputBox(m_TableArff.getParent(), title, "New value", m_LastReplace);
      if (valueNew == null)
        return;
      m_LastReplace = valueNew;
    }
    
    model = (ArffSortedTableModel) m_TableArff.getModel();
    model.setNotificationEnabled(false);

    // undo
    addUndoPoint();
    model.setUndoEnabled(false);
    String valueCopy = value;
    // set value
    for (i = 0; i < m_TableArff.getRowCount(); i++) {
      if (o == menuItemSetAllValues) {
        if (valueCopy.equals("NaN") || valueCopy.equals("?")) {
          value = null;
        }
        model.setValueAt(value, i, m_CurrentCol);
      }
      else
        if ( (o == menuItemSetMissingValues) 
            && model.isMissingAt(i, m_CurrentCol) )
          model.setValueAt(value, i, m_CurrentCol);
        else if ( (o == menuItemReplaceValues) 
            && model.getValueAt(i, m_CurrentCol).toString().equals(value) )
          model.setValueAt(valueNew, i, m_CurrentCol);
    }
    model.setUndoEnabled(true);
    model.setNotificationEnabled(true);
    model.notifyListener(new TableModelEvent(model, 0, model.getRowCount(), m_CurrentCol, TableModelEvent.UPDATE));
    
    // refresh
    m_TableArff.repaint();
  }
  
  /**
   * deletes the currently selected attribute
   */
  public void deleteAttribute() {
    ArffSortedTableModel   model;
    
    // no column selected?
    if (m_CurrentCol == -1)
      return;
    
    model = (ArffSortedTableModel) m_TableArff.getModel();

    // really an attribute column?
    if (model.getAttributeAt(m_CurrentCol) == null)
      return;
    
    // really?
    if (ComponentHelper.showMessageBox(
        getParent(), 
        "Confirm...",
        "Do you really want to delete the attribute '" 
        + model.getAttributeAt(m_CurrentCol).name() + "'?",
        JOptionPane.YES_NO_OPTION,
        JOptionPane.QUESTION_MESSAGE) != JOptionPane.YES_OPTION)
      return;
    
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    model.deleteAttributeAt(m_CurrentCol);
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
  }
  
  /**
   * deletes the chosen attributes
   */
  public void deleteAttributes() {
    ListSelectorDialog    dialog;
    ArffSortedTableModel       model;
    Object[]              atts;
    int[]                 indices;
    int                   i;
    JList                 list;
    int                   result;
    
    list   = new JList(getAttributes());
    dialog = new ListSelectorDialog(null, list);
    result = dialog.showDialog();
    
    if (result != ListSelectorDialog.APPROVE_OPTION)
      return;
    
    atts = list.getSelectedValues();
    
    // really?
    if (ComponentHelper.showMessageBox(
        getParent(), 
        "Confirm...",
        "Do you really want to delete these " 
        + atts.length + " attributes?",
        JOptionPane.YES_NO_OPTION,
        JOptionPane.QUESTION_MESSAGE) != JOptionPane.YES_OPTION)
      return;
    
    model   = (ArffSortedTableModel) m_TableArff.getModel();
    indices = new int[atts.length];
    for (i = 0; i < atts.length; i++)
      indices[i] = model.getAttributeColumn(atts[i].toString());
    
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    model.deleteAttributes(indices);
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
  }
  
  /**
   * sets the current attribute as class attribute, i.e. it moves it to the end
   * of the attributes
   */
  public void attributeAsClass() {
    ArffSortedTableModel   model;
    
    // no column selected?
    if (m_CurrentCol == -1)
      return;
    
    model   = (ArffSortedTableModel) m_TableArff.getModel();

    // really an attribute column?
    if (model.getAttributeAt(m_CurrentCol) == null)
      return;
    
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    model.attributeAsClassAt(m_CurrentCol);
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
  }
  
  /**
   * renames the current attribute
   */
  public void renameAttribute() {
    ArffSortedTableModel   model;
    String            newName;
    
    // no column selected?
    if (m_CurrentCol == -1)
      return;
    
    model   = (ArffSortedTableModel) m_TableArff.getModel();

    // really an attribute column?
    if (model.getAttributeAt(m_CurrentCol) == null)
      return;
    
    newName = ComponentHelper.showInputBox(getParent(), "Rename attribute...", "Enter new Attribute name", model.getAttributeAt(m_CurrentCol).name());
    if (newName == null)
      return;
    
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    model.renameAttributeAt(m_CurrentCol, newName);
    setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
  }
  
  /**
   * deletes the currently selected instance
   */
  public void deleteInstance() {
    int               index;
    
    index = m_TableArff.getSelectedRow();
    if (index == -1)
      return;
    
    ((ArffSortedTableModel) m_TableArff.getModel()).deleteInstanceAt(index);
  }
  
  /**
   * deletes all the currently selected instances
   */
  public void deleteInstances() {
    int[]             indices;
    
    if (m_TableArff.getSelectedRow() == -1)
      return;
    
    indices = m_TableArff.getSelectedRows();
    ((ArffSortedTableModel) m_TableArff.getModel()).deleteInstances(indices);
  }
  
  /**
   * sorts the instances via the currently selected column
   */
  public void sortInstances() {
    if (m_CurrentCol == -1)
      return;
    
    ((ArffSortedTableModel) m_TableArff.getModel()).sortInstances(m_CurrentCol);
  }
  
  /**
   * copies the content of the selection to the clipboard
   */
  public void copyContent() {
    StringSelection      selection;
    Clipboard            clipboard;
    
    selection = getTable().getStringSelection();
    if (selection == null)
      return;
    
    clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    clipboard.setContents(selection, selection);
  }
  
  /**
   * searches for a string in the cells
   */
  public void search() {
    String              searchString;
    
    // display dialog
    searchString = ComponentHelper.showInputBox(getParent(), "Search...", "Enter the string to search for", m_LastSearch);
    if (searchString != null)
      m_LastSearch = searchString;
    
    getTable().setSearchString(searchString);
  }
  
  /**
   * clears the search, i.e. resets the found cells
   */
  public void clearSearch() {
    getTable().setSearchString("");
  }
  
  /**
   * calculates the optimal column width for the current column
   */
  public void setOptimalColWidth() {
    // no column selected?
    if (m_CurrentCol == -1)
      return;

    JTableHelper.setOptimalColumnWidth(getTable(), m_CurrentCol);
  }
  
  /**
   * calculates the optimal column widths for all columns
   */
  public void setOptimalColWidths() {
    JTableHelper.setOptimalColumnWidth(getTable());
  }
  
  /**
   * invoked when an action occurs
   * 
   * @param e		the action event
   */
  public void actionPerformed(ActionEvent e) {
    Object          o;
    
    o = e.getSource();
    
    if (o == menuItemMean)
      calcMean();
    else if (o == menuItemSetAllValues)
      setValues(menuItemSetAllValues);
    else if (o == menuItemSetMissingValues)
      setValues(menuItemSetMissingValues);
    else if (o == menuItemReplaceValues)
      setValues(menuItemReplaceValues);
    else if (o == menuItemRenameAttribute)
      renameAttribute();
    else if (o == menuItemAttributeAsClass)
      attributeAsClass();
    else if (o == menuItemDeleteAttribute)
      deleteAttribute();
    else if (o == menuItemDeleteAttributes)
      deleteAttributes();
    else if (o == menuItemDeleteSelectedInstance)
      deleteInstance();
    else if (o == menuItemDeleteAllSelectedInstances)
      deleteInstances();
    else if (o == menuItemSortInstances)
      sortInstances();
    else if (o == menuItemSearch)
      search();
    else if (o == menuItemClearSearch)
      clearSearch();
    else if (o == menuItemUndo)
      undo();
    else if (o == menuItemCopy)
      copyContent();
    else if (o == menuItemOptimalColWidth)
      setOptimalColWidth();
    else if (o == menuItemOptimalColWidths)
      setOptimalColWidths();
  }
  
  /**
   * Invoked when a mouse button has been pressed and released on a component
   * 
   * @param e		the mouse event
   */
  public void mouseClicked(MouseEvent e) {
    int		col;
    boolean	popup;
    
    col   = m_TableArff.columnAtPoint(e.getPoint());
    popup =    ((e.getButton() == MouseEvent.BUTTON3) && (e.getClickCount() == 1))
            || ((e.getButton() == MouseEvent.BUTTON1) && (e.getClickCount() == 1) && e.isAltDown() && !e.isControlDown() && !e.isShiftDown());
    popup = popup && (getInstances() != null);
    
    if (e.getSource() == m_TableArff.getTableHeader()) {
      m_CurrentCol = col;
      
      // Popup-Menu
      if (popup) {
        e.consume();
        setMenu();
        initPopupMenus();
        m_PopupHeader.show(e.getComponent(), e.getX(), e.getY());
      }
    }
    else if (e.getSource() == m_TableArff) {
      // Popup-Menu
      if (popup) {
        e.consume();
        setMenu();
        initPopupMenus();
        m_PopupRows.show(e.getComponent(), e.getX(), e.getY());
      }
    }
    
    // highlihgt column
    if (    (e.getButton() == MouseEvent.BUTTON1)  
         && (e.getClickCount() == 1) 
         && (!e.isAltDown())
         && (col > -1) ) {
      m_TableArff.setSelectedColumn(col);
    }
  }
  
  /**
   * Invoked when the mouse enters a component.
   * 
   * @param e		the mouse event
   */
  public void mouseEntered(MouseEvent e) {
  }
  
  /**
   * Invoked when the mouse exits a component
   * 
   * @param e		the mouse event
   */
  public void mouseExited(MouseEvent e) {
  }
  
  /**
   * Invoked when a mouse button has been pressed on a component
   * 
   * @param e		the mouse event
   */
  public void mousePressed(MouseEvent e) {
  }
  
  /**
   * Invoked when a mouse button has been released on a component.
   * 
   * @param e		the mouse event
   */
  public void mouseReleased(MouseEvent e) {
  }
  
  /**
   * Invoked when the target of the listener has changed its state.
   * 
   * @param e		the change event
   */
  public void stateChanged(ChangeEvent e) {
    m_Changed = true;
    createTitle();
    notifyListener();
  }
  
  /**
   * notfies all listener of the change
   */
  public void notifyListener() {
    Iterator                iter;
    
    iter = m_ChangeListeners.iterator();
    while (iter.hasNext())
      ((ChangeListener) iter.next()).stateChanged(new ChangeEvent(this));
  }
  
  /**
   * Adds a ChangeListener to the panel
   * 
   * @param l		the listener to add
   */
  public void addChangeListener(ChangeListener l) {
    m_ChangeListeners.add(l);
  }
  
  /**
   * Removes a ChangeListener from the panel
   * 
   * @param l		the listener to remove
   */
  public void removeChangeListener(ChangeListener l) {
    m_ChangeListeners.remove(l);
  }
}

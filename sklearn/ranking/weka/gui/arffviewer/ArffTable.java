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
 * ArffTable.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.arffviewer;

import weka.core.Attribute;
import weka.core.Instances;
import weka.core.labelranking.PreferenceAttribute;
import weka.gui.ComponentHelper;
import weka.gui.JTableHelper;
import weka.gui.ViewerDialog;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.datatransfer.StringSelection;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;

import javax.swing.AbstractCellEditor;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TableModelEvent;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableModel;

/**
 * A specialized JTable for the Arff-Viewer.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.8 $ 
 */
public class ArffTable
  extends JTable {
  
  /** for serialization */
  static final long serialVersionUID = -2016200506908637967L;

  /**
   * a special Editor for editing the relation attribute.
   */
  protected class RelationalCellEditor
    extends AbstractCellEditor
    implements TableCellEditor {

    /** for serialization */
    private static final long serialVersionUID = 657969163293205963L;
    
    /** the button for opening the dialog */
    protected JButton m_Button;
    
    /** the current instances */
    protected Instances m_CurrentInst;
    
    /** the row index this editor is for */
    protected int m_RowIndex;
    
    /** the column index this editor is for */
    protected int m_ColumnIndex;
    
    /**
     * initializes the editor
     * 
     * @param rowIndex		the row index
     * @param columnIndex	the column index
     */
    public RelationalCellEditor(int rowIndex, int columnIndex) {
      super();

      m_CurrentInst = getInstancesAt(rowIndex, columnIndex);
      m_RowIndex    = rowIndex;
      m_ColumnIndex = columnIndex;
      
      m_Button = new JButton("...");
      m_Button.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent evt) {
          ViewerDialog        dialog;
          int                 result;
          
          dialog = new ViewerDialog(null);
          dialog.setTitle(
              "Relational attribute Viewer - " 
              + ((ArffSortedTableModel) getModel()).getInstances().attribute(m_ColumnIndex - 1).name());
          result = dialog.showDialog(m_CurrentInst);
          if (result == ViewerDialog.APPROVE_OPTION) {
            m_CurrentInst = dialog.getInstances();
            fireEditingStopped();
          }
          else {
            fireEditingCanceled();
          }
        }
      });
    }

    /**
     * returns the underlying instances at the given position
     * 
     * @param rowIndex		the row index
     * @param columnIndex	the column index
     * @return 			the corresponding instances
     */
    protected Instances getInstancesAt(int rowIndex, int columnIndex) {
      Instances			result;
      ArffSortedTableModel	model;
      double			value;
      
      model = (ArffSortedTableModel) getModel();
      value = model.getInstancesValueAt(rowIndex, columnIndex);
      result = model.getInstances().attribute(columnIndex - 1).relation((int) value);
      
      return result;
    }
    
    /**
     * Sets an initial value for the editor. This will cause the editor to 
     * stopEditing and lose any partially edited value if the editor is 
     * editing when this method is called.
     * 
     * @param table		the table this editor belongs to
     * @param value		the value to edit
     * @param isSelected	whether the cell is selected
     * @param row		the row index
     * @param column		the column index
     * @return			the 
     */
    public Component getTableCellEditorComponent(JTable table,
                                                 Object value,
                                                 boolean isSelected,
                                                 int row,
                                                 int column) {
      return m_Button;
    }

    /**
     * Returns the value contained in the editor.
     * 
     * @return		the value contained in the editor
     */
    public Object getCellEditorValue() {
      return m_CurrentInst;
    }
  }
  
  /** the search string */
  private String m_SearchString;
  /** the listeners for changes */
  private HashSet m_ChangeListeners;
  
  /**
   * initializes with no model
   */
  public ArffTable() {
    this(new ArffSortedTableModel(""));
  }
  
  /**
   * initializes with the given model
   * 
   * @param model		the model to use
   */
  public ArffTable(TableModel model) {
    super(model);
    
    setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
  }
  
  /**
   * sets the new model
   * 
   * @param model		the model to use
   */
  public void setModel(TableModel model) {
    ArffSortedTableModel      arffModel;
    
    // initialize the search
    m_SearchString = null;
    
    // init the listeners
    if (m_ChangeListeners == null)
      m_ChangeListeners = new HashSet();
    
    super.setModel(model);
    
    if (model == null)
      return;
    
    if (!(model instanceof ArffSortedTableModel))
      return;
    
    arffModel = (ArffSortedTableModel) model;
    arffModel.addMouseListenerToHeader(this);
    arffModel.addTableModelListener(this);
    arffModel.sort(0);
    setLayout();
    setSelectedColumn(0);
    
    // disable column moving
    if (getTableHeader() != null)
      getTableHeader().setReorderingAllowed(false);
  }

  /**
   * returns the cell editor for the given cell
   * 
   * @param row		the row index
   * @param column	the column index
   * @return		the cell editor
   */
  public TableCellEditor getCellEditor(int row, int column) {
    TableCellEditor		result;
    
    // relational attribute?
    if (    (getModel() instanceof ArffSortedTableModel) 
	 && (((ArffSortedTableModel) getModel()).getType(column) == Attribute.RELATIONAL) )
      result = new RelationalCellEditor(row, column);
    // default
    else
      result = super.getCellEditor(row, column);
    
    return result;
  }

  /**
   * returns whether the model is read-only
   * 
   * @return 		true if model is read-only
   */
  public boolean isReadOnly() {
    return ((ArffSortedTableModel) getModel()).isReadOnly();
  }
  
  /**
   * sets whether the model is read-only
   * 
   * @param value	if true the model is set to read-only
   */
  public void setReadOnly(boolean value) {
    ((ArffSortedTableModel) getModel()).setReadOnly(value);
  }
  
  /**
   * sets the cell renderer and calcs the optimal column width
   */
  private void setLayout() {
    ArffSortedTableModel      arffModel;
    int                  i;
    JComboBox            combo;
    Enumeration          enm;
    
    arffModel = (ArffSortedTableModel) getModel();
    
    for (i = 0; i < getColumnCount(); i++) {
      // optimal colwidths (only according to header!)
      JTableHelper.setOptimalHeaderWidth(this, i);
      
      // CellRenderer
      getColumnModel().getColumn(i).setCellRenderer(
          new ArffTableCellRenderer());
      
      // CellEditor
      if (i > 0) {
        if (arffModel.getType(i) == Attribute.NOMINAL) {
          combo = new JComboBox();
          combo.addItem(null);
          enm  = arffModel.getInstances().attribute(i - 1).enumerateValues();
          while (enm.hasMoreElements())
            combo.addItem(enm.nextElement());
          getColumnModel().getColumn(i).setCellEditor(new DefaultCellEditor(combo));
        }
        else {
          getColumnModel().getColumn(i).setCellEditor(null);
        }
      }
    }
  }
  
  /**
   * returns the basically the attribute name of the column and not the
   * HTML column name via getColumnName(int)
   * 
   * @param columnIndex		the column index
   * @return 			the plain name
   */
  public String getPlainColumnName(int columnIndex) {
    ArffSortedTableModel      arffModel;
    String               result;
    
    result = "";
    
    if (getModel() == null)
      return result;
    if (!(getModel() instanceof ArffSortedTableModel))  
      return result;
    
    arffModel = (ArffSortedTableModel) getModel();
    
    if ( (columnIndex >= 0) && (columnIndex < getColumnCount()) ) {
      if (columnIndex == 0)
        result = "No.";
      else
        result = arffModel.getAttributeAt(columnIndex).name();
    }
    
    return result;
  }
  
  /**
   * returns the selected content in a StringSelection that can be copied to
   * the clipboard and used in Excel, if nothing is selected the whole table
   * is copied to the clipboard
   * 
   * @return			the current selection
   */
  public StringSelection getStringSelection() {
    StringSelection         result;
    int[]                   indices;
    int                     i;
    int                     n;
    StringBuffer            tmp;
    
    result = null;
    
    // nothing selected? -> all
    if (getSelectedRow() == -1) {
      // really?
      if (ComponentHelper.showMessageBox(
            getParent(),
            "Question...",
            "Do you really want to copy the whole table?",
            JOptionPane.YES_NO_OPTION,
            JOptionPane.QUESTION_MESSAGE ) != JOptionPane.YES_OPTION)
        return result;
      
      indices = new int[getRowCount()];
      for (i = 0; i < indices.length; i++)
        indices[i] = i;
    }
    else {
      indices = getSelectedRows();
    }
    
    // get header
    tmp = new StringBuffer();
    for (i = 0; i < getColumnCount(); i++) {
      if (i > 0)
        tmp.append("\t");
      tmp.append(getPlainColumnName(i));
    }
    tmp.append("\n");
    
    // get content
    for (i = 0; i < indices.length; i++) {
      for (n = 0; n < getColumnCount(); n++) {
        if (n > 0)
          tmp.append("\t");
        tmp.append(getValueAt(indices[i], n).toString());
      }
      tmp.append("\n");
    }
    
    result = new StringSelection(tmp.toString());
    
    return result;
  }
  
  /**
   * sets the search string to look for in the table, NULL or "" disables
   * the search
   * 
   * @param searchString	the search string to use
   */
  public void setSearchString(String searchString) {
    this.m_SearchString = searchString;
    repaint();
  }
  
  /**
   * returns the search string, can be NULL if no search string is set
   * 
   * @return			the current search string
   */
  public String getSearchString() {
    return m_SearchString;
  }
  
  /**
   * sets the selected column
   * 
   * @param index		the column to select
   */
  public void setSelectedColumn(int index) {
    getColumnModel().getSelectionModel().clearSelection();
    getColumnModel().getSelectionModel().setSelectionInterval(index, index);
    resizeAndRepaint();
    if (getTableHeader() != null)
      getTableHeader().resizeAndRepaint();
  }
  
  /**
   * This fine grain notification tells listeners the exact range of cells, 
   * rows, or columns that changed.
   * 
   * @param e		the table event
   */
  public void tableChanged(TableModelEvent e) {
    super.tableChanged(e);
    
    setLayout();
    notifyListener();
  }
  
  /**
   * notfies all listener of the change
   */
  private void notifyListener() {
    Iterator                iter;
    
    iter = m_ChangeListeners.iterator();
    while (iter.hasNext())
      ((ChangeListener) iter.next()).stateChanged(new ChangeEvent(this));
  }
  
  /**
   * Adds a ChangeListener to the panel
   * 
   * @param l			the listener to add
   */
  public void addChangeListener(ChangeListener l) {
    m_ChangeListeners.add(l);
  }
  
  /**
   * Removes a ChangeListener from the panel
   * 
   * @param l			the listener to remove
   */
  public void removeChangeListener(ChangeListener l) {
    m_ChangeListeners.remove(l);
  }
}

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
 * JTableHelper.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.Component;
import java.awt.Point;
import java.awt.Rectangle;
import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

/**
 * A helper class for JTable, e.g. calculating the optimal colwidth.
 *
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.2 $ 
 */

public class JTableHelper {
  // the table to work with
  private JTable          jtable;
  
  /**
   * initializes the object
   */
  public JTableHelper(JTable jtable) {
    this.jtable = jtable;
  }
  
  /**
   * returns the JTable
   */
  public JTable getJTable() {
    return jtable;
  }
  
  /**
   * calcs the optimal column width of the given column
   */
  public int calcColumnWidth(int col) {
    return calcColumnWidth(getJTable(), col);
  }
  
  /**
   *  Calculates the optimal width for the column of the given table. The
   *  calculation is based on the preferred width of the header and cell
   *  renderer.
   *  <br>
   *  Taken from the newsgoup de.comp.lang.java with some modifications.<br>
   *  Taken from FOPPS/EnhancedTable - http://fopps.sourceforge.net/<br>
   *
   *  @param table    the table to calculate the column width
   *  @param col      the column to calculate the widths
   *  @return         the width, -1 if error
   */
  public static int calcColumnWidth(JTable table, int col) {
    int width = calcHeaderWidth(table, col);
    if (width == -1)
      return width;
    
    TableColumnModel columns = table.getColumnModel();
    TableModel data = table.getModel();
    int rowCount = data.getRowCount();
    TableColumn column = columns.getColumn(col);
    try {
      for (int row = rowCount - 1; row >= 0; --row) {
        Component c = table.prepareRenderer(
            table.getCellRenderer(row, col),
            row, col);
        width = Math.max(width, c.getPreferredSize().width + 10);
      }
    } 
    catch (Exception e) {
      e.printStackTrace();
    }
    
    return width;
  }
  
  /**
   * calcs the optimal header width of the given column
   */
  public int calcHeaderWidth(int col) {
    return calcHeaderWidth(getJTable(), col);
  }
  
  /**
   *  Calculates the optimal width for the header of the given table. The
   *  calculation is based on the preferred width of the header renderer.
   *
   *  @param table    the table to calculate the column width
   *  @param col      the column to calculate the widths
   *  @return         the width, -1 if error
   */
  public static int calcHeaderWidth(JTable table, int col) {
    if (table == null) 
      return -1;
    
    if (col < 0 || col > table.getColumnCount()) {
      System.out.println("invalid col " + col);
      return -1;
    }
    
    JTableHeader header = table.getTableHeader();
    TableCellRenderer defaultHeaderRenderer = null;
    if (header != null) defaultHeaderRenderer = header.getDefaultRenderer();
    TableColumnModel columns = table.getColumnModel();
    TableModel data = table.getModel();
    TableColumn column = columns.getColumn(col);
    int width = -1;
    TableCellRenderer h = column.getHeaderRenderer();
    if (h == null) h = defaultHeaderRenderer;
    if (h != null) {
      // Not explicitly impossible
      Component c = h.getTableCellRendererComponent(
          table,
          column.getHeaderValue(),
          false, false, -1, col);
      width = c.getPreferredSize().width + 5;
    }
    
    return width;
  }
  
  /**
   * sets the optimal column width for the given column
   */
  public void setOptimalColumnWidth(int col) {
    setOptimalColumnWidth(getJTable(), col);
  }
  
  /**
   * sets the optimal column width for the given column
   
   */
  public static void setOptimalColumnWidth(JTable jtable, int col) {
    int            width;
    TableColumn    column;
    JTableHeader   header;
    
    if ( (col >= 0) && (col < jtable.getColumnModel().getColumnCount()) ) {
      width = calcColumnWidth(jtable, col);
      
      if (width >= 0) {
        header = jtable.getTableHeader();
        column = jtable.getColumnModel().getColumn(col);
        column.setPreferredWidth(width);
        jtable.sizeColumnsToFit(-1);
        header.repaint();
      }
    }
  }
  
  /**
   * sets the optimal column width for all columns
   */
  public void setOptimalColumnWidth() {
    setOptimalColumnWidth(getJTable());
  }
  
  /**
   * sets the optimal column width for alls column if the given table
   */
  public static void setOptimalColumnWidth(JTable jtable) {
    int            i;
    
    for (i = 0; i < jtable.getColumnModel().getColumnCount(); i++)
      setOptimalColumnWidth(jtable, i);
  }
  
  /**
   * sets the optimal header width for the given column
   */
  public void setOptimalHeaderWidth(int col) {
    setOptimalHeaderWidth(getJTable(), col);
  }
  
  /**
   * sets the optimal header width for the given column
   
   */
  public static void setOptimalHeaderWidth(JTable jtable, int col) {
    int            width;
    TableColumn    column;
    JTableHeader   header;
    
    if ( (col >= 0) && (col < jtable.getColumnModel().getColumnCount()) ) {
      width = calcHeaderWidth(jtable, col);
      
      if (width >= 0) {
        header = jtable.getTableHeader();
        column = jtable.getColumnModel().getColumn(col);
        column.setPreferredWidth(width);
        jtable.sizeColumnsToFit(-1);
        header.repaint();
      }
    }
  }
  
  /**
   * sets the optimal header width for all columns
   */
  public void setOptimalHeaderWidth() {
    setOptimalHeaderWidth(getJTable());
  }
  
  /**
   * sets the optimal header width for alls column if the given table
   */
  public static void setOptimalHeaderWidth(JTable jtable) {
    int            i;
    
    for (i = 0; i < jtable.getColumnModel().getColumnCount(); i++)
      setOptimalHeaderWidth(jtable, i);
  }
  
  /** 
   * Assumes table is contained in a JScrollPane.
   * Scrolls the cell (rowIndex, vColIndex) so that it is visible
   * within the viewport.
   */
  public void scrollToVisible(int row, int col) {
    scrollToVisible(getJTable(), row, col);
  }
  
  /** 
   * Assumes table is contained in a JScrollPane.
   * Scrolls the cell (rowIndex, vColIndex) so that it is visible
   * within the viewport.
   */
  public static void scrollToVisible(JTable table, int row, int col) {
    if (!(table.getParent() instanceof JViewport)) 
      return;
    
    JViewport viewport = (JViewport) table.getParent();
    
    // This rectangle is relative to the table where the
    // northwest corner of cell (0,0) is always (0,0).
    Rectangle rect = table.getCellRect(row, col, true);
    
    // The location of the viewport relative to the table
    Point pt = viewport.getViewPosition();
    
    // Translate the cell location so that it is relative
    // to the view, assuming the northwest corner of the
    // view is (0,0)
    rect.setLocation(rect.x - pt.x, rect.y - pt.y);
    
    // Scroll the area into view
    viewport.scrollRectToVisible(rect);
  }
}

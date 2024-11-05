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
 * ResultSetTable.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql;

import weka.gui.JTableHelper;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JTable;
import javax.swing.table.TableColumnModel;

/**
 * Represents an extended JTable, containing a table model based on a ResultSet
 * and the corresponding query.
 *
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.2 $
 */
public class ResultSetTable
  extends JTable {

  /** for serialization */
  private static final long serialVersionUID = -3391076671854464137L;

  /** the query the table model is based on */
  protected String m_Query;

  /** the connect string with which the query was run */
  protected String m_URL;

  /** the user that was used to connect to the DB */
  protected String m_User;

  /** the password that was used to connect to the DB */
  protected String m_Password;

  /**
   * initializes the table
   * @param url         the database URL
   * @param user        the database user
   * @param pw          the database password
   * @param query       the query
   */
  public ResultSetTable(String url, String user, String pw, String query, 
                        ResultSetTableModel model) {
    super(model);

    m_URL      = url;
    m_User     = user;
    m_Password = pw;
    m_Query    = query;
    
    setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

    // optimal colwidths (only according to header!)/cell renderer
    for (int i = 0; i < getColumnCount(); i++) {
      JTableHelper.setOptimalHeaderWidth(this, i);
      getColumnModel().getColumn(i).setCellRenderer(
          new ResultSetTableCellRenderer());
    }
    
    // double click on column displays optimal colwidth
    final JTable table = this;
    getTableHeader().addMouseListener(new MouseAdapter() {
      public void mouseClicked(MouseEvent e) {
        TableColumnModel columnModel = getColumnModel();
        int viewColumn = columnModel.getColumnIndexAtX(e.getX());
        int column = convertColumnIndexToModel(viewColumn);
        if (    (e.getButton() == MouseEvent.BUTTON1)
             && (e.getClickCount() == 2)
             && (column != -1) )
          JTableHelper.setOptimalColumnWidth(table, column);
      }
    });
    getTableHeader().setToolTipText("double left click on column displays the column with optimal width");
  }

  /**
   * returns the database URL that produced the table model
   */
  public String getURL() {
    return m_URL;
  }

  /**
   * returns the user that produced the table model
   */
  public String getUser() {
    return m_User;
  }

  /**
   * returns the password that produced the table model
   */
  public String getPassword() {
    return m_Password;
  }

  /**
   * returns the query that produced the table model
   */
  public String getQuery() {
    return m_Query;
  }

  /**
   * frees up the memory
   */
  public void finalize() throws Throwable {
    if (getModel() != null)
      ((ResultSetTableModel) getModel()).finalize();
    
    super.finalize();

    System.gc();
  }
}

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
 * QueryPanel.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql;

import weka.gui.ListSelectorDialog;
import weka.gui.sql.event.ConnectionEvent;
import weka.gui.sql.event.ConnectionListener;
import weka.gui.sql.event.HistoryChangedEvent;
import weka.gui.sql.event.HistoryChangedListener;
import weka.gui.sql.event.QueryExecuteEvent;
import weka.gui.sql.event.QueryExecuteListener;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.ResultSet;
import java.util.HashSet;
import java.util.Iterator;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

/**
 * Represents a panel for entering an SQL query.
 *
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.4 $
 */
public class QueryPanel 
  extends JPanel 
  implements ConnectionListener, CaretListener {

  /** for serialization. */
  private static final long serialVersionUID = 4348967824619706636L;

  /** the name of the history. */
  public final static String HISTORY_NAME = "query";

  /** the name for the max rows in the history. */
  public final static String MAX_ROWS = "max_rows";
  
  /** the parent of this panel. */
  protected JFrame m_Parent;

  /** the textarea for the query. */
  protected JTextArea m_TextQuery;

  /** the execute button. */
  protected JButton m_ButtonExecute = new JButton("Execute");

  /** the clear button. */
  protected JButton m_ButtonClear = new JButton("Clear");

  /** the history button. */
  protected JButton m_ButtonHistory = new JButton("History...");

  /** the spinner for the maximum number of rows. */
  protected JSpinner m_SpinnerMaxRows = new JSpinner();

  /** the connection listeners. */
  protected HashSet m_QueryExecuteListeners;

  /** the history listeners. */
  protected HashSet m_HistoryChangedListeners;

  /** for working on the database. */
  protected DbUtils m_DbUtils;

  /** whether we have a connection to a database or not. */
  protected boolean m_Connected;

  /** the query history. */
  protected DefaultListModel m_History = new DefaultListModel();
  
  /**
   * initializes the panel.
   * 
   * @param parent        the parent of this panel
   */
  public QueryPanel(JFrame parent) {
    super();
    
    m_Parent                  = parent;
    m_QueryExecuteListeners   = new HashSet();
    m_HistoryChangedListeners = new HashSet();
    m_DbUtils                 = null;
    m_Connected               = false;

    createPanel();
  }

  /**
   * creates the panel with all its components.
   */
  protected void createPanel() {
    JPanel              panel;
    JPanel              panel2;
    JPanel              panel3;
    SpinnerNumberModel  model;
    
    setLayout(new BorderLayout());
    
    // textarea
    m_TextQuery = new JTextArea();
    m_TextQuery.addCaretListener(this);
    m_TextQuery.setFont(
        new Font("Monospaced", Font.PLAIN, m_TextQuery.getFont().getSize()));
    add(new JScrollPane(m_TextQuery), BorderLayout.CENTER);

    // buttons
    panel = new JPanel(new BorderLayout());
    add(panel, BorderLayout.EAST);
    m_ButtonExecute.setMnemonic('E');
    m_ButtonExecute.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  execute();
	}
      });
    panel.add(m_ButtonExecute, BorderLayout.NORTH);
    panel2 = new JPanel(new BorderLayout());
    panel.add(panel2, BorderLayout.CENTER);
    m_ButtonClear.setMnemonic('r');
    m_ButtonClear.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  clear();
	}
      });
    panel2.add(m_ButtonClear, BorderLayout.NORTH);
    panel3 = new JPanel(new BorderLayout());
    panel2.add(panel3, BorderLayout.CENTER);
    m_ButtonHistory.addActionListener(new ActionListener() {
	public void actionPerformed(ActionEvent e) {
	  showHistory();
	}
      });
    panel3.add(m_ButtonHistory, BorderLayout.NORTH);

    // limit
    panel3 = new JPanel(new FlowLayout());
    panel3.add(new JLabel("max. rows"));
    panel3.add(m_SpinnerMaxRows);
    panel2.add(panel3, BorderLayout.SOUTH);
    model = (SpinnerNumberModel) m_SpinnerMaxRows.getModel();
    model.setMaximum(new Integer(Integer.MAX_VALUE));
    model.setMinimum(new Integer(0));
    model.setValue(new Integer(100));
    model.setStepSize(new Integer(100));
    m_SpinnerMaxRows.setMinimumSize(
        new Dimension(50, m_SpinnerMaxRows.getHeight()));
    m_SpinnerMaxRows.setToolTipText("with 0 all rows are retrieved");
      
    // set initial state
    setButtons();
  }

  /**
   * sets the focus in a designated control.
   */
  public void setFocus() {
    m_TextQuery.requestFocus();
  }

  /**
   * sets the buttons according to the connected-state.
   */
  protected void setButtons() {
    boolean isEmpty;
    
    isEmpty = m_TextQuery.getText().trim().equals("");
    
    m_ButtonExecute.setEnabled((m_Connected) && (!isEmpty));
    m_ButtonClear.setEnabled(!isEmpty);
    m_ButtonHistory.setEnabled(m_History.size() > 0);
  }
  
  /**
   * This method gets called when the connection is either established
   * or disconnected.
   * 
   * @param evt		the event
   */
  public void connectionChange(ConnectionEvent evt) {
    m_Connected = evt.isConnected();
    m_DbUtils   = evt.getDbUtils();
    setButtons();
  }

  /**
   * executes the current query.
   */
  public void execute() {
    Exception     ex;
    ResultSet     rs;
    
    // not connected?
    if (!m_ButtonExecute.isEnabled())
      return;

    // no query?
    if (m_TextQuery.getText().trim().equals(""))
      return;

    // close old resultset
    try {
      if (m_DbUtils.getResultSet() != null)
        m_DbUtils.close();
    }
    catch (Exception e) {
      // ignore (if no resultset present we get an unncessary NullPointerEx.)
    }

    ex = null;
    rs = null;
    
    try {
      if (m_DbUtils.execute(getQuery())) {
        rs = m_DbUtils.getResultSet();
        // add to history
        addHistory(getQuery());
      }
    }
    catch (Exception e) {
      ex = new Exception(e.getMessage());
    }

    notifyQueryExecuteListeners(rs, ex);

    setButtons();
  }

  /**
   * clears the textarea.
   */
  public void clear() {
    m_TextQuery.setText("");
    m_SpinnerMaxRows.setValue(new Integer(100));
  }

  /**
   * adds the given string to the history (removes duplicates).
   * 
   * @param s           the string to add
   */
  protected void addHistory(String s) {
    if (s.equals(""))
      return;
    
    // no duplicates!
    if (m_History.contains(s))
      m_History.removeElement(s);

    m_History.add(0, s);
    
    // send notification
    notifyHistoryChangedListeners();
  }

  /**
   * sets the local history to the given one.
   * 
   * @param history     the history to use
   */
  public void setHistory(DefaultListModel history) {
    int           i;
    
    m_History.clear();
    for (i = 0; i < history.size(); i++)
      m_History.addElement(history.get(i));

    setButtons();
  }

  /**
   * returns the history.
   * 
   * @return        the current history
   */
  public DefaultListModel getHistory() {
    return m_History;
  }

  /**
   * displays the query history.
   */
  public void showHistory() {
    JList                 list;
    ListSelectorDialog    dialog;

    list   = new JList(m_History);
    dialog = new ListSelectorDialog(m_Parent, list);
    
    if (dialog.showDialog() == ListSelectorDialog.APPROVE_OPTION) {
      if (list.getSelectedValue() != null)
        setQuery(list.getSelectedValue().toString());
    }

    setButtons();
  }

  /**
   * sets the query in the textarea.
   * 
   * @param query         the query to display
   */
  public void setQuery(String query) {
    m_TextQuery.setText(query);
  }

  /**
   * returns the currently displayed query.
   * 
   * @return		the query
   */
  public String getQuery() {
    return m_TextQuery.getText();
  }

  /**
   * sets the maximum number of rows to display. 0 means unlimited.
   * 
   * @param rows	the maximum number of rows
   */
  public void setMaxRows(int rows) {
    if (rows >= 0)
      m_SpinnerMaxRows.setValue(new Integer(rows));
  }

  /**
   * returns the current value for the maximum number of rows. 0 means
   * unlimited.
   * 
   * @return		the maximum number of rows
   */
  public int getMaxRows() {
    return ((Integer) m_SpinnerMaxRows.getValue()).intValue();
  }

  /**
   * adds the given listener to the list of listeners.
   * 
   * @param l       the listener to add to the list
   */
  public void addQueryExecuteListener(QueryExecuteListener l) {
    m_QueryExecuteListeners.add(l);
  }

  /**
   * removes the given listener from the list of listeners.
   * 
   * @param l       the listener to remove
   */
  public void removeQueryExecuteListener(QueryExecuteListener l) {
    m_QueryExecuteListeners.remove(l);
  }

  /**
   * notifies the listeners of the event.
   * 
   * @param rs		the resultset
   * @param ex		the exception
   */
  protected void notifyQueryExecuteListeners(ResultSet rs, Exception ex) {
    Iterator              iter;
    QueryExecuteListener  l;

    iter = m_QueryExecuteListeners.iterator();
    while (iter.hasNext()) {
      l = (QueryExecuteListener) iter.next();
      l.queryExecuted(
          new QueryExecuteEvent(
            this, m_DbUtils, getQuery(), getMaxRows(), rs, ex));
    }
  }

  /**
   * adds the given listener to the list of listeners.
   * 
   * @param l       the listener to add to the list
   */
  public void addHistoryChangedListener(HistoryChangedListener l) {
    m_HistoryChangedListeners.add(l);
  }

  /**
   * removes the given listener from the list of listeners.
   * 
   * @param l       the listener to remove
   */
  public void removeHistoryChangedListener(HistoryChangedListener l) {
    m_HistoryChangedListeners.remove(l);
  }

  /**
   * notifies the history listeners of the event.
   */
  protected void notifyHistoryChangedListeners() {
    Iterator                iter;
    HistoryChangedListener  l;

    iter = m_HistoryChangedListeners.iterator();
    while (iter.hasNext()) {
      l = (HistoryChangedListener) iter.next();
      l.historyChanged(
          new HistoryChangedEvent(this, HISTORY_NAME, getHistory()));
    }
  }

  /**
   * Called when the caret position is updated.
   * 
   * @param event	the event
   */
  public void caretUpdate(CaretEvent event) {
    setButtons();
  }
}


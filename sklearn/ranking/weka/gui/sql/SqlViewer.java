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
 * SqlViewer.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql;

import weka.core.Memory;
import weka.gui.LookAndFeel;
import weka.gui.sql.event.ConnectionEvent;
import weka.gui.sql.event.ConnectionListener;
import weka.gui.sql.event.HistoryChangedEvent;
import weka.gui.sql.event.HistoryChangedListener;
import weka.gui.sql.event.QueryExecuteEvent;
import weka.gui.sql.event.QueryExecuteListener;
import weka.gui.sql.event.ResultChangedEvent;
import weka.gui.sql.event.ResultChangedListener;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Properties;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Represents a little tool for querying SQL databases.
 *
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 4724 $
 */
public class SqlViewer 
  extends    JPanel 
  implements ConnectionListener, 
             HistoryChangedListener,
             QueryExecuteListener, 
             ResultChangedListener {

  /** for serialization. */
  private static final long serialVersionUID = -4395028775566514329L;

  /** the name of the history file (in the home directory). */
  protected final static String HISTORY_FILE = "SqlViewerHistory.props";

  /** the width property in the history file. */
  public final static String WIDTH = "width";

  /** the height property in the history file. */
  public final static String HEIGHT = "height";
  
  /** the parent of this panel. */
  protected JFrame m_Parent;

  /** the connection panel. */
  protected ConnectionPanel m_ConnectionPanel;

  /** the query panel. */
  protected QueryPanel m_QueryPanel;

  /** the result panel. */
  protected ResultPanel m_ResultPanel;

  /** the info panel. */
  protected InfoPanel m_InfoPanel;

  /** the connect string with which the query was run. */
  protected String m_URL;

  /** the user that was used to connect to the DB. */
  protected String m_User;

  /** the password that was used to connect to the DB. */
  protected String m_Password;

  /** the currently selected query. */
  protected String m_Query;

  /** stores the history. */
  protected Properties m_History;
  
  /**
   * initializes the SqlViewer.
   * 
   * @param parent        the parent of this panel
   */
  public SqlViewer(JFrame parent) {
    super();
   
    m_Parent   = parent;
    m_URL      = "";
    m_User     = "";
    m_Password = "";
    m_Query    = "";
    m_History  = new Properties();
   
    createPanel();
  }

  /**
   * builds the interface.
   */
  protected void createPanel() {
    JPanel          panel;
    JPanel          panel2;
    
    setLayout(new BorderLayout());

    // connection
    m_ConnectionPanel = new ConnectionPanel(m_Parent);
    panel = new JPanel(new BorderLayout());
    add(panel, BorderLayout.NORTH);
    panel.setBorder(BorderFactory.createCompoundBorder(
                    BorderFactory.createTitledBorder("Connection"),
                    BorderFactory.createEmptyBorder(0, 5, 5, 5)));
    panel.add(m_ConnectionPanel, BorderLayout.CENTER);

    // query
    m_QueryPanel = new QueryPanel(m_Parent);
    panel = new JPanel(new BorderLayout());
    add(panel, BorderLayout.CENTER);
    panel2 = new JPanel(new BorderLayout());
    panel2.setBorder(BorderFactory.createCompoundBorder(
                     BorderFactory.createTitledBorder("Query"),
                     BorderFactory.createEmptyBorder(0, 5, 5, 5)));
    panel2.add(m_QueryPanel, BorderLayout.NORTH);
    panel.add(panel2, BorderLayout.NORTH);

    // result
    m_ResultPanel = new ResultPanel(m_Parent);
    m_ResultPanel.setQueryPanel(m_QueryPanel);
    panel2 = new JPanel(new BorderLayout());
    panel2.setBorder(BorderFactory.createCompoundBorder(
                     BorderFactory.createTitledBorder("Result"),
                     BorderFactory.createEmptyBorder(0, 5, 5, 5)));
    panel2.add(m_ResultPanel, BorderLayout.CENTER);
    panel.add(panel2, BorderLayout.CENTER);

    // info
    m_InfoPanel = new InfoPanel(m_Parent);
    panel = new JPanel(new BorderLayout());
    add(panel, BorderLayout.SOUTH);
    panel.setBorder(BorderFactory.createCompoundBorder(
                    BorderFactory.createTitledBorder("Info"),
                    BorderFactory.createEmptyBorder(0, 5, 5, 5)));
    panel.add(m_InfoPanel, BorderLayout.CENTER);

    // listeners
    addConnectionListener(this);
    addConnectionListener(m_QueryPanel);
    addQueryExecuteListener(this);
    addQueryExecuteListener(m_ResultPanel);
    addResultChangedListener(this);
    addHistoryChangedListener(this);

    // history
    loadHistory(true);
  }
  
  /**
   * This method gets called when the connection is either established
   * or disconnected.
   * 
   * @param evt		the event
   */
  public void connectionChange(ConnectionEvent evt) {
    if (evt.getType() == ConnectionEvent.DISCONNECT) {
      m_InfoPanel.append(   "disconnect from: " 
                          + evt.getDbUtils().getDatabaseURL(),
                          "information_small.gif" );
    }
    else {
      m_InfoPanel.append(   "connecting to: " 
                          + evt.getDbUtils().getDatabaseURL() 
                          + " = " + evt.isConnected(),
                          "information_small.gif" );
    }

    // did an exception happen?
    if (evt.getException() != null)
      m_InfoPanel.append("exception: " + evt.getException(), "error_small.gif");

    // set focus
    if (evt.isConnected())
      m_QueryPanel.setFocus();
    else
      m_ConnectionPanel.setFocus();
  }
  
  /**
   * This method gets called when a query has been executed.
   * 
   * @param evt		the event
   */
  public void queryExecuted(QueryExecuteEvent evt) {
    ResultSetHelper   helper;
    
    if (evt.failed()) {
      m_InfoPanel.append("Query:" + evt.getQuery(), "error_small.gif");
      m_InfoPanel.append("exception: " + evt.getException(), "error_small.gif");
    }
    else {
      m_InfoPanel.append("Query: " + evt.getQuery(), "information_small.gif");
      try {
        if (evt.hasResult()) {
          helper = new ResultSetHelper(evt.getResultSet());
          if ((evt.getMaxRows() > 0) && (helper.getRowCount() >= evt.getMaxRows()))
            m_InfoPanel.append(helper.getRowCount() + " rows selected (" 
                + evt.getMaxRows() + " displayed).", 
                "information_small.gif");
          else if (helper.getRowCount() == -1)
            m_InfoPanel.append("Unknown number of rows selected (due to JDBC driver restrictions).", 
              "information_small.gif");
          else
            m_InfoPanel.append(helper.getRowCount() + " rows selected.", 
                "information_small.gif");
        }

        // save max rows
        loadHistory(false);
        m_History.setProperty(
            QueryPanel.MAX_ROWS, Integer.toString(evt.getMaxRows()));
        saveHistory();
      }
      catch (Exception e) {
        e.printStackTrace();
      }
    }
  }
  
  /**
   * This method gets called when a query has been executed.
   * 
   * @param evt		the event
   */
  public void resultChanged(ResultChangedEvent evt) {
    m_URL      = evt.getURL();
    m_User     = evt.getUser();
    m_Password = evt.getPassword();
    m_Query    = evt.getQuery();
  }
  
  /**
   * This method gets called when a history is modified.
   * It saves the history immediately to the users home directory.
   * 
   * @param evt		the event
   */
  public void historyChanged(HistoryChangedEvent evt) {
    // load history, in case some other process changed it!
    loadHistory(false);
    
    m_History.setProperty( 
        evt.getHistoryName(), modelToString(evt.getHistory()));

    // save it
    saveHistory();
  }

  /**
   * returns the filename of the history file.
   * 
   * @return		the history file
   */
  protected String getHistoryFilename() {
    return   System.getProperties().getProperty("user.home")
           + File.separatorChar
           + HISTORY_FILE;
  }

  /**
   * transforms the given, comma-separated string into a DefaultListModel.
   * 
   * @param s     the string to break up and transform into a list model
   * @return      the generated DefaultListModel
   */
  protected DefaultListModel stringToModel(String s) {
    DefaultListModel    result;
    String              tmpStr;
    int                 i;
    boolean             quote;
    String[]            find;
    String[]            replace;
    int                 index;

    result = new DefaultListModel();
    
    // get rid of doubled quotes, \\n, etc.
    find    = new String[]{"\"\"", "\\n", "\\r", "\\t"};
    replace = new String[]{"\"",   "\n",  "\r",  "\t"};
    for (i = 0; i < find.length; i++) {
      tmpStr = "";
      while (s.length() > 0) {
        index = s.indexOf(find[i]);
        if (index > -1) {
          tmpStr += s.substring(0, index) + replace[i];
          s       = s.substring(index + 2);
        }
        else {
          tmpStr += s;
          s       = "";
        }
      }
      s = tmpStr;
    }

    quote  = false;
    tmpStr = "";
    for (i = 0; i < s.length(); i++) {
      if (s.charAt(i) == '"') {
        quote = !quote;
        tmpStr += "" + s.charAt(i);
      }
      else if (s.charAt(i) == ',') {
        if (quote) {
          tmpStr += "" + s.charAt(i);
        }
        else {
          if (tmpStr.startsWith("\""))
            tmpStr = tmpStr.substring(1, tmpStr.length() - 1);
          result.addElement(tmpStr);
          tmpStr = "";
        }
      }
      else {
        tmpStr += "" + s.charAt(i);
      }
    }
    
    // add last element
    if (!tmpStr.equals("")) {
      if (tmpStr.startsWith("\""))
        tmpStr = tmpStr.substring(1, tmpStr.length() - 1);
      result.addElement(tmpStr);
    }

    return result;
  }

  /**
   * converts the given model into a comma-separated string.
   * 
   * @param m       the model to convert
   * @return        the string representation of the model
   */
  protected String modelToString(DefaultListModel m) {
    String      result;
    String      tmpStr;
    int         i;
    int         n;
    boolean     quote;

    result = "";

    for (i = 0; i < m.size(); i++) {
      if (i > 0)
        result += ",";
      
      tmpStr = m.get(i).toString();
      quote  = (tmpStr.indexOf(",") > -1) || (tmpStr.indexOf(" ") > -1);
      
      if (quote)
        result += "\"";
      
      for (n = 0; n < tmpStr.length(); n++) {
        // double quotes
        if (tmpStr.charAt(n) == '"')
          result += "" + "\"\"";
        // normal character
        else
          result += "" + tmpStr.charAt(n);
      }
      
      if (quote)
        result += "\"";
    }

    return result;
  }

  /**
   * loads the history properties of the SqlViewer in the user's home directory.
   * 
   * @param set       whether to set the read properties in the panels or not
   * @see #HISTORY_FILE
   */
  protected void loadHistory(boolean set) {
    BufferedInputStream		str;
    File			file;
    int				width;
    int				height;

    try {
      file = new File(getHistoryFilename());
      if (file.exists()) {
        str = new BufferedInputStream(
            new FileInputStream(getHistoryFilename()));
        m_History.load(str);
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }

    // set the histories
    if (set) {
      m_ConnectionPanel.setHistory(
          stringToModel(
            m_History.getProperty(ConnectionPanel.HISTORY_NAME, "")));
      m_QueryPanel.setHistory(
          stringToModel(
            m_History.getProperty(QueryPanel.HISTORY_NAME, "")));
      m_QueryPanel.setMaxRows(
          Integer.parseInt(m_History.getProperty(QueryPanel.MAX_ROWS, "100")));

      width  = Integer.parseInt(m_History.getProperty(WIDTH, "0"));
      height = Integer.parseInt(m_History.getProperty(HEIGHT, "0"));
      if ((width != 0) && (height != 0))
	setPreferredSize(new Dimension(width, height));
    }
  }

  /**
   * saves the history properties of the SqlViewer in the user's home directory.
   * 
   * @see #HISTORY_FILE
   */
  protected void saveHistory() {
    BufferedOutputStream    str;
    
    try {
      str = new BufferedOutputStream(
          new FileOutputStream(getHistoryFilename()));
      m_History.store(str, "SQL-Viewer-History");
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * obtains the size of the panel and saves it in the history.
   * 
   * @see #saveHistory()
   */
  public void saveSize() {
    m_History.setProperty(WIDTH, "" + getSize().width);
    m_History.setProperty(HEIGHT, "" + getSize().height);

    saveHistory();
  }
  
  /**
   * calls the clear method of all sub-panels to set back to default values
   * and free up memory.
   */
  public void clear() {
    m_ConnectionPanel.clear();
    m_QueryPanel.clear();
    m_ResultPanel.clear();
    m_InfoPanel.clear();
  }

  /**
   * returns the database URL from the currently active tab in the ResultPanel,
   * otherwise an empty string.
   * 
   * @see 		ResultPanel
   * @return		the currently selected tab's URL
   */
  public String getURL() {
    return m_URL;
  }

  /**
   * returns the user from the currently active tab in the ResultPanel,
   * otherwise an empty string.
   * 
   * @see		ResultPanel
   * @return		the currently selected tab's user
   */
  public String getUser() {
    return m_User;
  }

  /**
   * returns the password from the currently active tab in the ResultPanel,
   * otherwise an empty string.
   * 
   * @see 		ResultPanel
   * @return		the currently selected tab's password
   */
  public String getPassword() {
    return m_Password;
  }

  /**
   * returns the query from the currently active tab in the ResultPanel,
   * otherwise an empty string.
   * 
   * @see		ResultPanel
   * @return		the currently selected tab's query
   */
  public String getQuery() {
    return m_Query;
  }

  /**
   * adds the given listener to the list of listeners.
   * 
   * @param l		the listener to add to the list
   */
  public void addConnectionListener(ConnectionListener l) {
    m_ConnectionPanel.addConnectionListener(l);
  }

  /**
   * removes the given listener from the list of listeners.
   * 
   * @param l		the listener to remove
   */
  public void removeConnectionListener(ConnectionListener l) {
    m_ConnectionPanel.removeConnectionListener(l);
  }

  /**
   * adds the given listener to the list of listeners.
   * 
   * @param l		the listener to add to the list
   */
  public void addQueryExecuteListener(QueryExecuteListener l) {
    m_QueryPanel.addQueryExecuteListener(l);
  }

  /**
   * removes the given listener from the list of listeners.
   * 
   * @param l		the listener to remove
   */
  public void removeQueryExecuteListener(QueryExecuteListener l) {
    m_QueryPanel.removeQueryExecuteListener(l);
  }

  /**
   * adds the given listener to the list of listeners.
   * 
   * @param l		the listener to add to the list
   */
  public void addResultChangedListener(ResultChangedListener l) {
    m_ResultPanel.addResultChangedListener(l);
  }

  /**
   * removes the given listener from the list of listeners.
   * 
   * @param l		the listener to remove
   */
  public void removeResultChangedListener(ResultChangedListener l) {
    m_ResultPanel.removeResultChangedListener(l);
  }

  /**
   * adds the given listener to the list of listeners.
   * 
   * @param l		the listener to add to the list
   */
  public void addHistoryChangedListener(HistoryChangedListener l) {
    m_ConnectionPanel.addHistoryChangedListener(l);
    m_QueryPanel.addHistoryChangedListener(l);
  }

  /**
   * removes the given listener from the list of listeners.
   * 
   * @param l		the listener to remove
   */
  public void removeHistoryChangedListener(HistoryChangedListener l) {
    m_ConnectionPanel.removeHistoryChangedListener(l);
    m_QueryPanel.removeHistoryChangedListener(l);
  }

  /** for monitoring the Memory consumption. */
  private static Memory m_Memory = new Memory(true);
  
  /** the sql viewer. */
  private static SqlViewer m_Viewer;
  
  /**
   * starts the SQL-Viewer interface.
   * 
   * @param args	the commandline arguments - ignored
   */
  public static void main(String[] args) {
    weka.core.logging.Logger.log(weka.core.logging.Logger.Level.INFO, "Logging started");
    LookAndFeel.setLookAndFeel();
    
    try {
      // uncomment to disable the memory management:
      //m_Memory.setEnabled(false);

      final JFrame jf = new JFrame("Weka SQL-Viewer");
      m_Viewer = new SqlViewer(jf);
      jf.getContentPane().setLayout(new BorderLayout());
      jf.getContentPane().add(m_Viewer, BorderLayout.CENTER);
      jf.addWindowListener(new WindowAdapter() {
        public void windowClosing(WindowEvent e) {
          m_Viewer.saveSize();
          jf.dispose();
          System.exit(0);
        }
      });
      jf.pack();
      jf.setSize(800, 600);
      jf.setVisible(true);

      Thread memMonitor = new Thread() {
        public void run() {
          while (true) {
            try {
              this.sleep(4000);

              System.gc();

              if (m_Memory.isOutOfMemory()) {
                // clean up
                jf.dispose();
                m_Viewer = null;
                System.gc();

                // stop threads
                m_Memory.stopThreads();

                // display error
                System.err.println("\ndisplayed message:");
                m_Memory.showOutOfMemory();
                System.err.println("\nexiting");
                System.exit(-1);
              }

            } 
            catch (InterruptedException ex) { 
              ex.printStackTrace(); 
            }
          }
        }
      };

      memMonitor.setPriority(Thread.MAX_PRIORITY);
      memMonitor.start();
    } 
    catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

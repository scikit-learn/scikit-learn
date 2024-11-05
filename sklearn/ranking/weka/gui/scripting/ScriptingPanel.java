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
 * ScriptingPanel.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import weka.core.Tee;
import weka.gui.PropertyDialog;
import weka.gui.ReaderToTextPane;
import weka.gui.scripting.event.TitleUpdatedEvent;
import weka.gui.scripting.event.TitleUpdatedListener;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.util.HashSet;
import java.util.Iterator;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JTextPane;

/**
 * Abstract ancestor for scripting panels.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public abstract class ScriptingPanel
  extends JPanel
  implements TitleUpdatedListener {

  /** for serialization. */
  private static final long serialVersionUID = 7593091442691911406L;

  /** The new output stream for System.out. */
  protected PipedOutputStream m_POO;

  /** The new output stream for System.err. */
  protected PipedOutputStream m_POE;

  /** The thread that sends output from m_POO to the output box. */
  protected ReaderToTextPane m_OutRedirector;

  /** The thread that sends output from m_POE to the output box. */
  protected ReaderToTextPane m_ErrRedirector;
  
  /** whether debug mode is on. */
  protected boolean m_Debug;
  
  /** the listeners for the changes in the title. */
  protected HashSet<TitleUpdatedListener> m_TitleUpdatedListeners;

  /**
   * Default constructor.
   */
  public ScriptingPanel() {
    super();
    
    initialize();
    initGUI();
    initFinish();
  }
  
  /**
   * For initializing member variables.
   */
  protected void initialize() {
    m_POO            = new PipedOutputStream();
    m_POE            = new PipedOutputStream();
    m_Debug          = false;
    m_TitleUpdatedListeners = new HashSet<TitleUpdatedListener>();
  }
  
  /**
   * Sets up the GUI after initializing the members.
   * The JTextArea returned via <code>getOutputArea()</code> must be setup here.
   * 
   * @see	#initialize()
   * @see	#getOutput()
   */
  protected void initGUI() {
  }
  
  /**
   * Finishes up after initializing members and setting up the GUI.
   * Redirects stdout and stderr using <code>getOutputArea()</code>.
   * 
   * @see	#initialize()
   * @see	#initGUI()
   * @see	#getOutput()
   */
  protected void initFinish() {
    // Redirect System.out to the text area
    try {
      PipedInputStream pio = new PipedInputStream(m_POO);
      Tee teeOut = new Tee(System.out);
      System.setOut(teeOut);
      teeOut.add(new PrintStream(m_POO));
      Reader reader = new InputStreamReader(pio);
      m_OutRedirector = new ReaderToTextPane(reader, getOutput(), Color.BLACK);
      m_OutRedirector.start();
    }
    catch (Exception e) {
      System.err.println("Error redirecting stdout");
      e.printStackTrace();
      m_OutRedirector = null;
    }

    // Redirect System.err to the text area
    try {
      PipedInputStream pie = new PipedInputStream(m_POE);
      Tee teeErr = new Tee(System.err);
      System.setErr(teeErr);
      teeErr.add(new PrintStream(m_POE));
      Reader reader = new InputStreamReader(pie);
      m_ErrRedirector = new ReaderToTextPane(reader, getOutput(), Color.RED);
      m_ErrRedirector.start();
    }
    catch (Exception e) {
      System.err.println("Error redirecting stderr");
      e.printStackTrace();
      m_ErrRedirector = null;
    }

    addTitleUpdatedListener(this);
  }
  
  /**
   * Returns an icon to be used in a frame.
   * 
   * @return		the icon
   */
  public abstract ImageIcon getIcon();
  
  /**
   * Returns the current title for the frame/dialog.
   * 
   * @return		the title
   */
  public abstract String getTitle();
  
  /**
   * Returns the text area that is used for displaying output on stdout
   * and stderr.
   * 
   * @return		the JTextArea
   */
  public abstract JTextPane getOutput();
  
  /**
   * Returns the menu bar to to be displayed in the frame.
   * 
   * @return		the menu bar, null if not applicable
   */
  public abstract JMenuBar getMenuBar();
  
  /**
   * Turns on/off debugging mode.
   * 
   * @param value	if true, debug mode is turned on
   */
  public void setDebug(boolean value) {
    m_Debug = value;
  }
  
  /**
   * Returns whether debugging mode is on.
   * 
   * @return		true if debug mode is turned on
   */
  public boolean getDebug() {
    return m_Debug;
  }
  
  /**
   * Adds the listener to the internal list.
   * 
   * @param l		the listener to add
   */
  public void addTitleUpdatedListener(TitleUpdatedListener l) {
    m_TitleUpdatedListeners.add(l);
  }
  
  /**
   * Removes the listener from the internal list.
   * 
   * @param l		the listener to remove
   */
  public void removeTitleUpdatedListener(TitleUpdatedListener l) {
    m_TitleUpdatedListeners.remove(l);
  }
  
  /**
   * Sends the event to all listeners for title updates.
   * 
   * @param e		the event to send
   */
  protected void notifyTitleUpdatedListeners(TitleUpdatedEvent e) {
    Iterator<TitleUpdatedListener>	iter;
    
    iter = m_TitleUpdatedListeners.iterator();
    while (iter.hasNext())
      iter.next().titleUpdated(e);
  }
  
  /**
   * Gets called when the title of the frame/dialog needs updating.
   * 
   * @param event	the event that got sent
   */
  public void titleUpdated(TitleUpdatedEvent event) {
    if (PropertyDialog.getParentDialog(ScriptingPanel.this) != null)
      PropertyDialog.getParentDialog(ScriptingPanel.this).setTitle(getTitle());
    else if (PropertyDialog.getParentFrame(ScriptingPanel.this) != null)
      PropertyDialog.getParentFrame(ScriptingPanel.this).setTitle(getTitle());
  }
  
  /**
   * Displays the panel in a frame.
   * 
   * @param panel	the panel to display
   * @param args	currently ignored commandline parameters
   */
  public static void showPanel(ScriptingPanel panel, String[] args) {
    showPanel(panel, args, 800, 600);
  }
  
  /**
   * Displays the panel in a frame.
   * 
   * @param panel	the panel to display
   * @param args	currently ignored commandline parameters
   * @param width	the width of the frame
   * @param height	the height of the frame
   */
  public static void showPanel(ScriptingPanel panel, String[] args, int width, int height) {
    try {
      JFrame frame = new JFrame();
      frame.getContentPane().setLayout(new BorderLayout());
      frame.getContentPane().add(panel, BorderLayout.CENTER);
      frame.setJMenuBar(panel.getMenuBar());
      frame.setSize(new Dimension(width, height));
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.setTitle(panel.getTitle());
      frame.setIconImage(panel.getIcon().getImage());
      frame.setLocationRelativeTo(null);
      if ((args.length > 0) && (panel instanceof FileScriptingPanel))
	((FileScriptingPanel) panel).open(new File(args[0]));
      frame.setVisible(true);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}

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
 * Script.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.SerializedObject;
import weka.core.Utils;
import weka.core.WekaException;
import weka.gui.ExtensionFileFilter;
import weka.gui.scripting.event.ScriptExecutionEvent;
import weka.gui.scripting.event.ScriptExecutionListener;
import weka.gui.scripting.event.ScriptExecutionEvent.Type;

import java.io.File;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.Document;

/**
 * A simple helper class for loading, saving scripts.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public abstract class Script
  implements OptionHandler, Serializable {

  /** for serialization. */
  private static final long serialVersionUID = 5053328052680586401L;

  /**
   * The Thread for running a script.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5142 $
   */
  public abstract static class ScriptThread
    extends Thread {
    
    /** the owning script. */
    protected Script m_Owner;
    
    /** commandline arguments. */
    protected String[] m_Args;
    
    /** whether the thread was stopped. */
    protected boolean m_Stopped;
    
    /**
     * Initializes the thread.
     * 
     * @param owner	the owning script
     * @param args	the commandline arguments
     */
    public ScriptThread(Script owner, String[] args) {
      super();
      
      m_Owner = owner;
      m_Args  = args.clone();
    }

    /**
     * Returns the owner.
     * 
     * @return		the owning script
     */
    public Script getOwner() {
      return m_Owner;
    }
    
    /**
     * Returns the commandline args.
     * 
     * @return		the arguments
     */
    public String[] getArgs() {
      return m_Args;
    }
    
    /**
     * Performs the actual run.
     */
    protected abstract void doRun();
    
    /**
     * Executes the script.
     */
    public void run() {
      m_Stopped = false;
      
      getOwner().notifyScriptFinishedListeners(new ScriptExecutionEvent(m_Owner, Type.STARTED));
      try {
	doRun();
	if (!m_Stopped)
	  getOwner().notifyScriptFinishedListeners(new ScriptExecutionEvent(m_Owner, Type.FINISHED));
      }
      catch (Exception e) {
	e.printStackTrace();
	getOwner().notifyScriptFinishedListeners(new ScriptExecutionEvent(m_Owner, Type.ERROR, e));
      }
      getOwner().m_ScriptThread = null;
    }
    
    /**
     * Stops the script execution.
     */
    public void stopScript() {
      if (isAlive()) {
	m_Stopped = true;
	try {
	  stop();
	}
	catch (Exception e) {
	  // ignored
	}
      }
    }
  }

  /** the backup extension. */
  public final static String BACKUP_EXTENSION = ".bak";
  
  /** the document this script is a wrapper around. */
  protected Document m_Document;
  
  /** the filename of the script. */
  protected File m_Filename;
  
  /** the newline used on this platform. */
  protected String m_NewLine;
  
  /** whether the script is modified. */
  protected boolean m_Modified;
  
  /** the current script thread. */
  protected transient ScriptThread m_ScriptThread;
  
  /** optional listeners when the script finishes. */
  protected HashSet<ScriptExecutionListener> m_FinishedListeners;
  
  /**
   * Initializes the script.
   */
  public Script() {
    this(null);
  }
  
  /**
   * Initializes the script.
   * 
   * @param doc		the document to use as basis
   */
  public Script(Document doc) {
    this(doc, null);
  }
  
  /**
   * Initializes the script. Automatically loads the specified file, if not
   * null.
   * 
   * @param doc		the document to use as basis
   * @param file	the file to load (if not null)
   */
  public Script(Document doc, File file) {
    initialize();
    
    m_Document = doc;
    
    if (m_Document != null) {
      m_Document.addDocumentListener(new DocumentListener() {
	public void changedUpdate(DocumentEvent e) {
	  m_Modified = true;
	}
	public void insertUpdate(DocumentEvent e) {
	  m_Modified = true;
	}
	public void removeUpdate(DocumentEvent e) {
	  m_Modified = true;
	}
      });
    }
    
    if (file != null)
      open(file);
  }
  
  /**
   * Initializes the script.
   */
  protected void initialize() {
    m_Filename          = null;
    m_NewLine           = System.getProperty("line.separator");
    m_Modified          = false;
    m_ScriptThread      = null;
    m_FinishedListeners = new HashSet<ScriptExecutionListener>();
  }

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    return new Vector().elements();
  }

  /**
   * Parses a given list of options. 
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
  }

  /**
   * Gets the current settings of the script.
   *
   * @return 		an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    return new String[0];
  }
  
  /**
   * Returns the extension filters for this type of script.
   * 
   * @return		the filters
   */
  public abstract ExtensionFileFilter[] getFilters();
  
  /**
   * Returns the default extension. Gets automatically added to files if
   * their name doesn't end with this.
   * 
   * @return		the default extension (incl. the dot)
   * @see		#saveAs(File)
   */
  public abstract String getDefaultExtension();
  
  /**
   * Returns the current filename.
   * 
   * @return		the filename, null if no file loaded/saved
   */
  public File getFilename() {
    return m_Filename;
  }
  
  /**
   * Returns the new line string in use.
   * 
   * @return		the new line string
   */
  public String getNewLine() {
    return m_NewLine;
  }
  
  /**
   * Returns whether the script is modified.
   * 
   * @return		true if the script is modified
   */
  public boolean isModified() {
    return m_Modified;
  }
  
  /**
   * Returns the content.
   * 
   * @return		the content or null in case of an error
   */
  public String getContent() {
    String	result;
    
    if (m_Document == null)
      return "";
    
    try {
      synchronized(m_Document) {
	result = m_Document.getText(0, m_Document.getLength());
      }
    }
    catch (Exception e) {
      e.printStackTrace();
      result = null;
    }
    
    return result;
  }
  
  /**
   * Sets the content.
   * 
   * @param value	the new content
   */
  public void setContent(String value) {
    if (m_Document == null)
      return;
    
    try {
      m_Document.insertString(0, value, null);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Checks whether the extension of the file is a known one.
   * 
   * @param file	the file to check
   * @return		true if the exetnsion is known
   */
  protected boolean checkExtension(File file) {
    boolean			result;
    int				i;
    int				n;
    ExtensionFileFilter[]	filters;
    String[]			exts;

    result = false;
    filters  = getFilters();
    for (i = 0; i < filters.length; i++) {
      exts = filters[i].getExtensions();
      for (n = 0; n < exts.length; n++) {
	if (file.getName().endsWith(exts[n])) {
	  result = true;
	  break;
	}
      }
      if (result)
	break;
    }
    
    return result;
  }
  
  /**
   * Empties the document.
   */
  public void empty() {
    if (m_Document != null) {
      try {
	m_Document.remove(0, m_Document.getLength());
      }
      catch (Exception e) {
	// ignored
      }
    }
    
    m_Modified = false;
    m_Filename = null;
  }
  
  /**
   * Tries to open the file.
   * 
   * @param file	the file to open
   * @return		true if successfully read
   */
  public boolean open(File file) {
    boolean	result;
    String	content;

    if (m_Document == null)
      return true;
    
    // Warn if extension unwknown
    if (!checkExtension(file))
      System.err.println("Extension of file '" + file + "' is unknown!");
    
    try {
      // clear old content
      m_Document.remove(0, m_Document.getLength());
      
      // add new content
      content = ScriptUtils.load(file);
      if (content == null)
	throw new WekaException("Error reading content of file '" + file + "'!");
      m_Document.insertString(0, content, null);
      
      m_Modified = false;
      m_Filename = file;
      result     = true;
    }
    catch (Exception e) {
      e.printStackTrace();
      try {
	m_Document.remove(0, m_Document.getLength());
      }
      catch (Exception ex) {
	// ignored
      }
      result     = false;
      m_Filename = null;
    }
    
    return result;
  }
  
  /**
   * Saves the file under with the current filename.
   * 
   * @return		true if successfully written
   */
  public boolean save() {
    if (m_Filename == null)
      return false;
    else
      return saveAs(m_Filename);
  }
  
  /**
   * Saves the file under with the given filename (and updates the internal
   * filename).
   * 
   * @param file	the filename to write the content to
   * @return		true if successfully written
   */
  public boolean saveAs(File file) {
    boolean	result;
    File	backupFile;
    
    if (m_Document == null)
      return true;
    
    // correct extension?
    if (!checkExtension(file))
      file = new File(file.getPath() + getDefaultExtension());

    // backup previous file
    if (file.exists()) {
      backupFile = new File(file.getPath() + BACKUP_EXTENSION);
      try {
	ScriptUtils.copy(file, backupFile);
      }
      catch (Exception e) {
	e.printStackTrace();
      }
    }
    
    // save current content
    try {
      result     = ScriptUtils.save(file, m_Document.getText(0, m_Document.getLength()));
      m_Filename = file;
      m_Modified = false;
    }
    catch (Exception e) {
      e.printStackTrace();
      result = false;
    }
    
    return result;
  }

  /**
   * Returns whether scripts can be executed.
   * 
   * @return		true if scripts can be executed
   */
  protected abstract boolean canExecuteScripts();

  /**
   * Returns a new thread to execute.
   * 
   * @param args	optional commandline arguments
   * @return		the new thread object
   */
  public abstract ScriptThread newThread(String[] args);
  
  /**
   * Performs pre-execution checks:
   * <ul>
   * 	<li>whether a script is currently running.</li>
   * 	<li>whether script has changed and needs saving</li>
   * 	<li>whether a filename is set (= empty content)</li>
   * </ul>
   * Throws exceptions if checks not met.
   * 
   * @param args	optional commandline arguments
   * @throws Exception	if checks fail
   */
  protected void preCheck(String[] args) throws Exception {
    if (m_ScriptThread != null)
      throw new Exception("A script is currently running!");
    if (m_Modified)
      throw new Exception("The Script has been modified!");
    if (m_Filename == null)
      throw new Exception("The Script contains no content?");
  }
  
  /**
   * Executes the script.
   * 
   * @param args	optional commandline arguments
   */
  protected void execute(String[] args) {
    m_ScriptThread = newThread(args);
    try {
      m_ScriptThread.start();
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /**
   * Executes the script.
   * 
   * @param args	optional commandline arguments, can be null
   * @throws Exception	if checks or execution fail
   */
  public void start(String[] args) throws Exception {
    if (args == null)
      args = new String[0];
    
    preCheck(args);
    
    execute(args);
  }
  
  /**
   * Stops the execution of the script.
   */
  public void stop() {
    if (isRunning()) {
      m_ScriptThread.stopScript();
      m_ScriptThread = null;
      notifyScriptFinishedListeners(new ScriptExecutionEvent(this, Type.STOPPED));
    }
  }
  
  /**
   * Executes the script without loading it first.
   * 
   * @param file	the script to execute
   * @param args	the commandline parameters for the script
   */
  public void run(File file, String[] args) {
    Script	script;
    
    try {
      script = (Script) new SerializedObject(this).getObject();
      script.m_Filename = file;
      script.m_Modified = false;
      script.start(args);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  /**
   * Returns whether the script is still running.
   * 
   * @return		true if the script is still running
   */
  public boolean isRunning() {
    return (m_ScriptThread != null);
  }
  
  /**
   * Adds the given listener to its internal list.
   * 
   * @param l		the listener to add
   */
  public void addScriptFinishedListener(ScriptExecutionListener l) {
    m_FinishedListeners.add(l);
  }
  
  /**
   * Removes the given listener from its internal list.
   * 
   * @param l		the listener to remove
   */
  public void removeScriptFinishedListener(ScriptExecutionListener l) {
    m_FinishedListeners.remove(l);
  }
  
  /**
   * Notifies all listeners.
   * 
   * @param e		the event to send to all listeners
   */
  protected void notifyScriptFinishedListeners(ScriptExecutionEvent e) {
    Iterator<ScriptExecutionListener>	iter;
    
    iter = m_FinishedListeners.iterator();
    while (iter.hasNext())
      iter.next().scriptFinished(e);
  }
  
  /**
   * Returns the content as string.
   * 
   * @return		the current content
   */
  public String toString() {
    String	result;
    
    try {
      if (m_Document == null)
	result = "";
      else
	result = m_Document.getText(0, m_Document.getLength());
    }
    catch (Exception e) {
      result = "";
    }
    
    return result.toString();
  }

  /**
   * Make up the help string giving all the command line options.
   *
   * @param script 	the script to include options for
   * @return 		a string detailing the valid command line options
   */
  protected static String makeOptionString(Script script) {
    StringBuffer 	result;
    Enumeration 	enm;
    Option 		option;
    
    result = new StringBuffer("");

    result.append("\nHelp requested:\n\n");
    result.append("-h or -help\n");
    result.append("\tDisplays this help screen.\n");
    result.append("-s <file>\n");
    result.append("\tThe script to execute.\n");

    enm = script.listOptions();
    while (enm.hasMoreElements()) {
      option = (Option) enm.nextElement();
      result.append(option.synopsis() + '\n');
      result.append(option.description() + "\n");
    }

    result.append("\n");
    result.append("Any additional options are passed on to the script as\n");
    result.append("command-line parameters.\n");
    result.append("\n");
    
    return result.toString();
  }
  
  /**
   * Runs the specified script. All options that weren't "consumed" (like 
   * "-s" for the script filename), will be used as commandline arguments for 
   * the actual script.
   * 
   * @param script	the script object to use
   * @param args	the commandline arguments
   * @throws Exception	if execution fails
   */
  public static void runScript(Script script, String[] args) throws Exception {
    String		tmpStr;
    File		scriptFile;
    Vector<String>	options;
    int			i;
    
    if (Utils.getFlag('h', args) || Utils.getFlag("help", args)) {
      System.out.println(makeOptionString(script));
    }
    else {
      // process options
      tmpStr = Utils.getOption('s', args);
      if (tmpStr.length() == 0)
        throw new WekaException("No script supplied!");
      else
	scriptFile = new File(tmpStr);
      script.setOptions(args);
      
      // remove empty elements from array
      options = new Vector<String>();
      for (i = 0; i < args.length; i++) {
	if (args[i].length() > 0)
	  options.add(args[i]);
      }
      
      // run script
      script.run(scriptFile, options.toArray(new String[options.size()]));
    }
  }
}

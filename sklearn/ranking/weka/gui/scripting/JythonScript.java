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
 * JythonScript.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import weka.core.Utils;
import weka.core.scripting.Jython;
import weka.gui.ExtensionFileFilter;

import java.io.File;

import javax.swing.text.Document;

/**
 * Represents a <a href="http://www.jython.org/" target="_blank">Jython</a> script.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public class JythonScript
  extends Script {
  
  /** for serialization. */
  private static final long serialVersionUID = 3469648507172973169L;

  /**
   * Executes a Jython script in a thread.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5142 $
   */
  public static class JythonThread
    extends ScriptThread {
    
    /**
     * Initializes the thread.
     * 
     * @param owner	the owning script
     * @param args	the commandline arguments
     */
    public JythonThread(Script owner, String[] args) {
      super(owner, args);
    }

    /**
     * Performs the actual run.
     */
    protected void doRun() {
      Jython	jython;
      Class[]	classes;
      Object[]	params;
      String	argv;
      String	arg;
      int	i;
      
      classes = new Class[]{String.class};
      params  = new Object[]{m_Owner.getFilename().getPath()};
      argv    = "sys.argv = ['" + Utils.backQuoteChars(m_Owner.getFilename().getPath()) + "'";
      for (i = 0; i < getArgs().length; i++) {
	arg   = Utils.backQuoteChars(getArgs()[i]);
	argv += ", '" + arg + "'";
      }
      argv   += "]";
      
      jython  = new Jython();
      
      // set commandline parameters
      jython.invoke("exec", new Class[]{String.class}, new Object[]{"import sys"});
      jython.invoke("exec", new Class[]{String.class}, new Object[]{argv});
      
      jython.invoke("execfile", classes, params);
    }
  }
  
  /**
   * Initializes the script.
   */
  public JythonScript() {
    super();
  }
  
  /**
   * Initializes the script.
   * 
   * @param doc		the document to use as basis
   */
  public JythonScript(Document doc) {
    super(doc);
  }
  
  /**
   * Initializes the script. Automatically loads the specified file, if not
   * null.
   * 
   * @param doc		the document to use as basis
   * @param file	the file to load (if not null)
   */
  public JythonScript(Document doc, File file) {
    super(doc, file);
  }
  
  /**
   * Returns the extension filters for this type of script.
   * 
   * @return		the filters
   */
  public ExtensionFileFilter[] getFilters() {
    ExtensionFileFilter[]	result;
    
    result = new ExtensionFileFilter[1];
    result[0] = new ExtensionFileFilter(getDefaultExtension(), "Jython script (*" + getDefaultExtension() + ")");
    
    return result;
  }
  
  /**
   * Returns the default extension. Gets automatically added to files if
   * their name doesn't end with this.
   * 
   * @return		the default extension (incl. the dot)
   */
  public String getDefaultExtension() {
    return ".py";
  }

  /**
   * Returns whether scripts can be executed, i.e., Jython is present.
   * 
   * @return		true if scripts can be executed
   */
  protected boolean canExecuteScripts() {
    return Jython.isPresent();
  }
  
  /**
   * Performs pre-execution checks.
   * <p/>
   * This method checks whether Jython is available (throws an exception if not).
   * 
   * @param args	optional commandline arguments
   * @throws Exception	if checks fail
   */
  protected void preCheck(String[] args) throws Exception {
    super.preCheck(args);
    
    if (!Jython.isPresent())
      throw new Exception("Jython classes are not present in CLASSPATH!");
  }

  /**
   * Returns a new thread to execute.
   * 
   * @param args	optional commandline arguments
   * @return		the new thread object
   */
  public ScriptThread newThread(String[] args) {
    return new JythonThread(this, args);
  }
  
  /**
   * Runs the script from commandline. Use "-h" to list all options.
   * 
   * @param args	the commandline arguments
   * @throws Exception	if execution fails
   */
  public static void main(String[] args) throws Exception {
    runScript(new JythonScript(), args);
  }
}

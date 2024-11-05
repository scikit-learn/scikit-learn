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
 * GroovyScript.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.scripting;

import weka.core.scripting.Groovy;
import weka.gui.ExtensionFileFilter;

import java.io.File;

import javax.swing.text.Document;

/**
 * Represents a <a href="http://groovy.codehaus.org/" target="_blank">Groovy</a> script.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5142 $
 */
public class GroovyScript
  extends Script {
  
  /** for serialization. */
  private static final long serialVersionUID = -3708517162415549420L;

  /**
   * Executes a Groovy script in a thread.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 5142 $
   */
  public static class GroovyThread
    extends ScriptThread {
    
    /**
     * Initializes the thread.
     * 
     * @param owner	the owning script
     * @param args	the commandline arguments
     */
    public GroovyThread(Script owner, String[] args) {
      super(owner, args);
    }

    /**
     * Tests whether the groovy object has a certain method.
     * 
     * @param groovy	the Groovy object to inspect
     * @param name	the method to look for
     * @return		true if the object has the given method
     */
    protected boolean hasMethod(Object groovy, String name) {
      boolean	result;
      try {
	groovy.getClass().getMethod(name, new Class[]{String[].class});
	result = true;
      }
      catch (Exception e) {
	result = false;
      }
      
      return result;
    }
    
    /**
     * Performs the actual run.
     */
    protected void doRun() {
      Object	groovy;
      
      groovy = Groovy.newInstance(m_Owner.getFilename(), Object.class);
      if (hasMethod(groovy, "run"))
	Groovy.invoke(groovy, "run", new Class[]{String[].class}, new Object[]{getArgs()});
      else if (hasMethod(groovy, "main"))
	Groovy.invoke(groovy, "main", new Class[]{String[].class}, new Object[]{getArgs()});
      else
	throw new IllegalStateException("Neither 'run' nor 'main' method found!");
    }
  }
  
  /**
   * Initializes the script.
   */
  public GroovyScript() {
    super();
  }
  
  /**
   * Initializes the script.
   * 
   * @param doc		the document to use as basis
   */
  public GroovyScript(Document doc) {
    super(doc);
  }
  
  /**
   * Initializes the script. Automatically loads the specified file, if not
   * null.
   * 
   * @param doc		the document to use as basis
   * @param file	the file to load (if not null)
   */
  public GroovyScript(Document doc, File file) {
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
    result[0] = new ExtensionFileFilter(getDefaultExtension(), "Groovy script (*" + getDefaultExtension() + ")");
    
    return result;
  }
  
  /**
   * Returns the default extension. Gets automatically added to files if
   * their name doesn't end with this.
   * 
   * @return		the default extension (incl. the dot)
   */
  public String getDefaultExtension() {
    return ".groovy";
  }

  /**
   * Returns whether scripts can be executed, i.e., Groovy is present.
   * 
   * @return		true if scripts can be executed
   */
  protected boolean canExecuteScripts() {
    return Groovy.isPresent();
  }
  
  /**
   * Performs pre-execution checks.
   * <p/>
   * This method checks whether Groovy is available (throws an exception if not).
   * 
   * @param args	optional commandline arguments
   * @throws Exception	if checks fail
   */
  protected void preCheck(String[] args) throws Exception {
    super.preCheck(args);
    
    if (!Groovy.isPresent())
      throw new Exception("Groovy classes are not present in CLASSPATH!");
  }

  /**
   * Returns a new thread to execute.
   * 
   * @param args	optional commandline arguments
   * @return		the new thread object
   */
  public ScriptThread newThread(String[] args) {
    return new GroovyThread(this, args);
  }
  
  /**
   * Runs the script from commandline. Use "-h" to list all options.
   * 
   * @param args	the commandline arguments
   * @throws Exception	if execution fails
   */
  public static void main(String[] args) throws Exception {
    runScript(new GroovyScript(), args);
  }
}

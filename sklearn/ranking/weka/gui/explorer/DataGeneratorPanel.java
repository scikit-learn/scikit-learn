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
 * DataGeneratorPanel.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.explorer;

import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.datagenerators.DataGenerator;
import weka.gui.GenericObjectEditor;
import weka.gui.Logger;
import weka.gui.PropertyPanel;
import weka.gui.SysErrLog;

import java.awt.BorderLayout;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;

import javax.swing.JOptionPane;
import javax.swing.JPanel;

/** 
 * A panel for generating artificial data via DataGenerators.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 1.4 $
 */
public class DataGeneratorPanel
  extends JPanel {

  /** for serialization */
  private static final long serialVersionUID = -2520408165350629380L;

  /** the GOE for the generators */
  protected GenericObjectEditor m_GeneratorEditor = new GenericObjectEditor();

  /** the generated Instances */
  protected Instances m_Instances = null;

  /** the generated output (as text) */
  protected StringWriter m_Output = new StringWriter();

  /** The destination for log/status messages */
  protected Logger m_Log = new SysErrLog();

  /** register the classes */
  static {
     GenericObjectEditor.registerEditors();
  }
  
  /**
   * creates the panel
   */
  public DataGeneratorPanel() {
    setLayout(new BorderLayout());
   
    add(new PropertyPanel(m_GeneratorEditor), BorderLayout.CENTER);

    // editor
    m_GeneratorEditor.setClassType(DataGenerator.class);
    m_GeneratorEditor.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent e) {
	repaint();
      }
    });
    
    // set default generator
    setGenerator(null);
  }

  /**
   * Sets the Logger to receive informational messages
   *
   * @param value 	the Logger that will now get info messages
   */
  public void setLog(Logger value) {
    m_Log = value;
  }

  /**
   * returns the generated instances, null if the process was cancelled.
   *
   * @return the generated Instances
   */
  public Instances getInstances() {
    return m_Instances;
  }

  /**
   * returns the generated output as text
   * 
   * @return		the generated output
   */
  public String getOutput() {
    return m_Output.toString();
  }

  /**
   * sets the generator to use initially
   * 
   * @param value	the data generator to use
   */
  public void setGenerator(DataGenerator value) {
    if (value != null)
      m_GeneratorEditor.setValue(value);
    else
      m_GeneratorEditor.setValue(
          new weka.datagenerators.classifiers.classification.RDG1());
  }

  /**
   * returns the currently selected DataGenerator
   * 
   * @return		the current data generator
   */
  public DataGenerator getGenerator() {
    return (DataGenerator) m_GeneratorEditor.getValue();
  }

  /**
   * generates the instances, returns TRUE if successful
   * 
   * @return		TRUE if successful
   * @see #getInstances()
   */
  public boolean execute() {
    DataGenerator     generator;
    boolean           result;
    String            relName;
    String            cname;
    String            cmd;
    
    result    = true;
    generator = (DataGenerator) m_GeneratorEditor.getValue();
    relName   = generator.getRelationName();

    cname = generator.getClass().getName().replaceAll(".*\\.", "");
    cmd = generator.getClass().getName();
    if (generator instanceof OptionHandler)
      cmd += " " + Utils.joinOptions(((OptionHandler) generator).getOptions());
    
    try {
      m_Log.logMessage("Started " + cname);
      m_Log.logMessage("Command: " + cmd);
      m_Output = new StringWriter();
      generator.setOutput(new PrintWriter(m_Output));
      DataGenerator.makeData(generator, generator.getOptions());
      m_Instances = new Instances(new StringReader(getOutput()));
      m_Log.logMessage("Finished " + cname);
    }
    catch (Exception e) {
      e.printStackTrace();
      JOptionPane.showMessageDialog(
          this, "Error generating data:\n" + e.getMessage(), 
          "Error", JOptionPane.ERROR_MESSAGE);
      m_Instances = null;
      m_Output    = new StringWriter();
      result      = false;
    }

    generator.setRelationName(relName);

    return result;
  }
}

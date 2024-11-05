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
 * AbstractPlotInstances.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.explorer;

import weka.core.Instances;
import weka.core.OptionHandler;
import weka.gui.visualize.PlotData2D;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

/**
 * Abstract superclass for generating plottable instances.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5714 $
 */
public abstract class AbstractPlotInstances
  implements Serializable, OptionHandler {

  /** for serialization. */
  private static final long serialVersionUID = 2375155184845160908L;

  /** the full dataset. */
  protected Instances m_Instances;
  
  /** the plotable instances. */
  protected Instances m_PlotInstances;

  /** whether processing has been finished up already. */
  protected boolean m_FinishUpCalled;
  
  /**
   * Initializes the container.
   */
  public AbstractPlotInstances() {
    initialize();
  }
  
  /**
   * Initializes the member variables.
   */
  protected void initialize() {
    m_Instances      = null;
    m_PlotInstances  = null;
    m_FinishUpCalled = false;
  }

  /**
   * Returns an enumeration of all the available options.
   *
   * @return 		an enumeration of all available options.
   */
  public Enumeration listOptions() {
    return new Vector().elements();
  }

  /**
   * Sets the OptionHandler's options using the given list. All options
   * will be set (or reset) during this call (i.e. incremental setting
   * of options is not possible).
   *
   * @param options 	the list of options as an array of strings
   * @throws Exception 	if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {
  }

  /**
   * Gets the current option settings for the OptionHandler.
   *
   * @return the list of current option settings as an array of strings
   */
  public String[] getOptions() {
    return new String[0];
  }
  
  /**
   * Sets up the structure for the plot instances.
   */
  protected abstract void determineFormat();  
  
  /**
   * Sets the instances that are the basis for the plot instances.
   * 
   * @param value	the training data to initialize with
   */
  public void setInstances(Instances value) {
    m_Instances = value;
  }
  
  /**
   * Returns the training data.
   * 
   * @return		the training data
   */
  public Instances getInstances() {
    return m_Instances;
  }
  
  /**
   * Default implementation only ensures that a dataset has been set.
   */
  protected void check() {
    if (m_Instances == null)
      throw new IllegalStateException("No instances set!");
  }
  
  /**
   * Performs checks, sets up the structure for the plot instances.
   */
  public void setUp() {
    m_FinishUpCalled = false;
    
    check();
    determineFormat();
  }
  
  /**
   * Performs optional post-processing.
   */
  protected void finishUp() {
    m_FinishUpCalled = true;
  }
  
  /**
   * Returns the generated plot instances.
   * 
   * @return		the instances to plot
   */
  public Instances getPlotInstances() {
    if (!m_FinishUpCalled)
      finishUp();
    
    return m_PlotInstances;
  }
  /**
   * Assembles and returns the plot. The relation name of the dataset gets
   * added automatically.
   * 
   * @param name	the name of the plot
   * @return		the plot
   * @throws Exception	if plot generation fails
   */
  protected abstract PlotData2D createPlotData(String name) throws Exception;
  
  /**
   * Assembles and returns the plot. The relation name of the dataset gets
   * added automatically.
   * 
   * @param name	the name of the plot
   * @return		the plot
   * @throws Exception	if plot generation fails
   */
  public PlotData2D getPlotData(String name) throws Exception {
    if (!m_FinishUpCalled)
      finishUp();
    
    return createPlotData(name);
  }
  
  /**
   * For freeing up memory. Plot data cannot be generated after this call!
   */
  public void cleanUp() {
    m_Instances      = null;
    m_PlotInstances  = null;
    m_FinishUpCalled = false;
  }
}

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
 * ExplorerDefaults.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.gui.explorer;

import weka.core.Utils;

import java.io.Serializable;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Properties;
import java.util.Vector;

/**
 * This class offers get methods for the default Explorer settings in 
 * the props file <code>weka/gui/explorer/Explorer.props</code>.
 *
 * @author  FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6624 $
 * @see #PROPERTY_FILE
 */
public class ExplorerDefaults
  implements Serializable {
  
  /** for serialization. */
  private static final long serialVersionUID = 4954795757927524225L;

  /** The name of the properties file. */
  public final static String PROPERTY_FILE = "weka/gui/explorer/Explorer.props";

  /** Properties associated with the explorer options. */
  protected static Properties PROPERTIES;
  static {
    try {
      PROPERTIES = Utils.readProperties(PROPERTY_FILE);
    }
    catch (Exception e) {
      System.err.println("Problem reading properties. Fix before continuing.");
      e.printStackTrace();
      PROPERTIES = new Properties();
    }
  }

  /**
   * returns the value for the specified property, if non-existent then the
   * default value.
   *
   * @param property      the property to retrieve the value for
   * @param defaultValue  the default value for the property
   * @return              the value of the specified property
   */
  public static String get(String property, String defaultValue) {
    return PROPERTIES.getProperty(property, defaultValue);
  }
  
  public static void set(String property, String value) {
    PROPERTIES.setProperty(property, value);
  }
  
  /**
   * returns the associated properties file.
   * 
   * @return              the props file
   */
  public final static Properties getProperties() {
    return PROPERTIES;
  }
  
  /**
   * Tries to instantiate the class stored for this property, optional 
   * options will be set as well. Returns null if unsuccessful.
   * 
   * @param property      the property to get the object for
   * @param defaultValue  the default object spec string
   * @return              if successful the fully configured object, null
   *                      otherwise
   */
  protected static Object getObject(String property, String defaultValue) {
    return getObject(property, defaultValue, Object.class);
  }
  
  /**
   * Tries to instantiate the class stored for this property, optional 
   * options will be set as well. Returns null if unsuccessful.
   * 
   * @param property      the property to get the object for
   * @param defaultValue  the default object spec string
   * @param cls           the class the object must be derived from
   * @return              if successful the fully configured object, null
   *                      otherwise
   */
  protected static Object getObject(String property, String defaultValue, Class cls) {
    Object	result;
    String	tmpStr;
    String[]	tmpOptions;
    
    result = null;
    
    try {
      tmpStr     = get(property, defaultValue);
      tmpOptions = Utils.splitOptions(tmpStr);
      if (tmpOptions.length != 0) {
	tmpStr        = tmpOptions[0];
	tmpOptions[0] = "";
	result        = Utils.forName(cls, tmpStr, tmpOptions);
      }
    }
    catch (Exception e) {
      e.printStackTrace();
      result = null;
    }
    
    return result;
  }

  /**
   * returns if the GOEs in the Explorer will be initialized based on the
   * data that is loaded into the Explorer.
   * 
   * @return		true if the GOEs get initialized
   */
  public static boolean getInitGenericObjectEditorFilter() {
    return Boolean.parseBoolean(get("InitGenericObjectEditorFilter", "false"));
  }
  
  /**
   * returns an array with the classnames of all the additional panels to 
   * display as tabs in the Explorer.
   * 
   * @return		the classnames
   */
  public static String[] getTabs() {
    String[]	result;
    String	tabs;
    
    // read and split on comma
    tabs   = get("Tabs", "weka.gui.explorer.ClassifierPanel,weka.gui.explorer.ClustererPanel,weka.gui.explorer.AssociationsPanel,weka.gui.explorer.AttributeSelectionPanel,weka.gui.explorer.VisualizePanel");
    result = tabs.split(",");
    
    return result;
  }
  
  /**
   * Returns the initial directory for the file chooser used for opening
   * datasets. 
   * <p/>
   * The following placeholders are recognized:
   * <pre>
   *   %t - the temp directory
   *   %h - the user's home directory
   *   %c - the current directory
   *   %% - gets replaced by a single percentage sign
   * </pre>
   * 
   * @return		the default directory
   */
  public static String getInitialDirectory() {
    String	result;
    
    result = get("InitialDirectory", "%c");
    result = result.replaceAll("%t", System.getProperty("java.io.tmpdir"));
    result = result.replaceAll("%h", System.getProperty("user.home"));
    result = result.replaceAll("%c", System.getProperty("user.dir"));
    result = result.replaceAll("%%", System.getProperty("%"));
    
    return result;
  }
  
  /**
   * returns the default filter (fully configured) for the preprocess panel.
   * 
   * @return		the default filter, null if none
   */
  public static Object getFilter() {
    return getObject("Filter", "", weka.filters.Filter.class);
  }
  
  /**
   * returns the default classifier (fully configured) for the classify panel.
   * 
   * @return		the default classifier, ZeroR by default
   */
  public static Object getClassifier() {
    Object	result;
    
    result = getObject(
		"Classifier", 
		weka.classifiers.rules.ZeroR.class.getName(), 
		weka.classifiers.Classifier.class);
    if (result == null)
      result = new weka.classifiers.rules.ZeroR();
    
    return result;
  }
  
  /**
   * returns the default classifier test mode for the classify panel.
   * 
   * @return		the default test mode
   */
  public static int getClassifierTestMode() {
    return Integer.parseInt(get("ClassifierTestMode", "1"));
  }
  
  /**
   * returns the default number of folds of the CV in the classify panel.
   * 
   * @return		the default number of folds
   */
  public static int getClassifierCrossvalidationFolds() {
    return Integer.parseInt(get("ClassifierCrossvalidationFolds", "10"));
  }
  
  /**
   * returns the default classifier test mode for the classify panel (0-99).
   * 
   * @return		the default precentage split
   */
  public static int getClassifierPercentageSplit() {
    return Integer.parseInt(get("ClassifierPercentageSplit", "66"));
  }

  /**
   * returns whether the built model is output.
   * 
   * @return		true if the built model is output
   */
  public static boolean getClassifierOutputModel() {
    return Boolean.parseBoolean(get("ClassifierOutputModel", "true"));
  }

  /**
   * returns whether additional per-class stats of the classifier are output.
   * 
   * @return		true if stats are output
   */
  public static boolean getClassifierOutputPerClassStats() {
    return Boolean.parseBoolean(get("ClassifierOutputPerClassStats", "true"));
  }

  /**
   * returns whether entropy-based evaluation meastures of the classifier 
   * are output.
   * 
   * @return		true if output
   */
  public static boolean getClassifierOutputEntropyEvalMeasures() {
    return Boolean.parseBoolean(get("ClassifierOutputEntropyEvalMeasures", "false"));
  }

  /**
   * returns whether the confusion matrix for the classifier is output.
   * 
   * @return		true if matrix is output
   */
  public static boolean getClassifierOutputConfusionMatrix() {
    return Boolean.parseBoolean(get("ClassifierOutputConfusionMatrix", "true"));
  }

  /**
   * returns whether the predictions of the classifier are output as well.
   * 
   * @return		true if predictions are output as well
   */
  public static boolean getClassifierOutputPredictions() {
    return Boolean.parseBoolean(get("ClassifierOutputPredictions", "false"));
  }

  /**
   * returns the string with the additional indices to output alongside the
   * predictions. 
   * 
   * @return		the indices, 0 if none are output
   */
  public static String getClassifierOutputAdditionalAttributes() {
    return get("ClassifierOutputAdditionalAttributes", "");
  }

  /**
   * returns whether the predictions of the classifier are stored for 
   * visualization.
   * 
   * @return		true if predictions are stored
   */
  public static boolean getClassifierStorePredictionsForVis() {
    return Boolean.parseBoolean(get("ClassifierStorePredictionsForVis", "true"));
  }

  /**
   * returns whether the evaluation of the classifier is done cost-sensitively.
   * 
   * @return		true if cost-sensitively done
   */
  public static boolean getClassifierCostSensitiveEval() {
    return Boolean.parseBoolean(get("ClassifierCostSensitiveEval", "false"));
  }
  
  /**
   * returns the default random seed value for the classifier for the 
   * classify panel.
   * 
   * @return		the default random seed
   */
  public static int getClassifierRandomSeed() {
    return Integer.parseInt(get("ClassifierRandomSeed", "1"));
  }

  /**
   * returns whether the order is preserved in case of the percentage split 
   * in the classify tab.
   * 
   * @return		true if order is preserved
   */
  public static boolean getClassifierPreserveOrder() {
    return Boolean.parseBoolean(get("ClassifierPreserveOrder", "false"));
  }

  /**
   * returns whether the source of a sourcable Classifier is output
   * in the classify tab.
   * 
   * @return		true if the source code is output
   */
  public static boolean getClassifierOutputSourceCode() {
    return Boolean.parseBoolean(get("ClassifierOutputSourceCode", "false"));
  }

  /**
   * returns the default classname for a sourcable Classifier in the classify tab.
   * 
   * @return		the default classname
   */
  public static String getClassifierSourceCodeClass() {
    return get("ClassifierSourceCodeClass", "Foobar");
  }
  
  /**
   * Returns an instance of the class used for generating plot instances
   * for displaying the classifier errors.
   * 
   * @return		an instance of the class
   */
  public static ClassifierErrorsPlotInstances getClassifierErrorsPlotInstances() {
    ClassifierErrorsPlotInstances	result;
    String			classname;
    String[]			options;
    
    try {
      options = Utils.splitOptions(get("ClassifierErrorsPlotInstances", "weka.gui.explorer.ClassifierErrorsPlotInstances"));
      classname  = options[0];
      options[0] = "";
      result     = (ClassifierErrorsPlotInstances) Utils.forName(ClassifierErrorsPlotInstances.class, classname, options);
    }
    catch (Exception e) {
      e.printStackTrace();
      result = new ClassifierErrorsPlotInstances();
    }
    
    return result;
  }
  
  /**
   * Returns the minimum size in pixels for plots of plotting classifier errors 
   * of numeric attributes.
   * 
   * @return		the size
   */
  public static int getClassifierErrorsMinimumPlotSizeNumeric() {
    return Integer.parseInt(get("ClassifierErrorsMinimumPlotSizeNumeric", "1"));
  }
  
  /**
   * Returns the maximum size in pixels for plots of plotting classifier errors 
   * of numeric attributes.
   * 
   * @return		the size
   */
  public static int getClassifierErrorsMaximumPlotSizeNumeric() {
    return Integer.parseInt(get("ClassifierErrorsMaximumPlotSizeNumeric", "20"));
  }
  
  /**
   * returns the default clusterer (fully configured) for the clusterer panel.
   * 
   * @return		the default clusterer, EM by default
   */
  public static Object getClusterer() {
    Object	result;
    
    result = getObject(
		"Clusterer", 
		weka.clusterers.EM.class.getName(), 
		weka.clusterers.Clusterer.class);
    if (result == null)
      result = new weka.clusterers.EM();
    
    return result;
  }
  
  /**
   * returns the default cluster test mode for the cluster panel.
   * 
   * @return		the default test mode
   */
  public static int getClustererTestMode() {
    return Integer.parseInt(get("ClustererTestMode", "3"));
  }

  /**
   * returns whether the clusters are storeed for visualization purposes
   * in the cluster panel.
   * 
   * @return		true if clusters are stored
   */
  public static boolean getClustererStoreClustersForVis() {
    return Boolean.parseBoolean(get("ClustererStoreClustersForVis", "true"));
  }
  
  /**
   * Returns an instance of the class used for generating plot instances
   * for displaying the cluster assignments.
   * 
   * @return		an instance of the class
   */
  public static ClustererAssignmentsPlotInstances getClustererAssignmentsPlotInstances() {
    ClustererAssignmentsPlotInstances	result;
    String			classname;
    String[]			options;
    
    try {
      options = Utils.splitOptions(get("ClustererAssignmentsPlotInstances", "weka.gui.explorer.ClustererAssignmentsPlotInstances"));
      classname  = options[0];
      options[0] = "";
      result     = (ClustererAssignmentsPlotInstances) Utils.forName(ClustererAssignmentsPlotInstances.class, classname, options);
    }
    catch (Exception e) {
      e.printStackTrace();
      result = new ClustererAssignmentsPlotInstances();
    }
    
    return result;
  }
  
  /**
   * returns the default associator (fully configured) for the associations panel.
   * 
   * @return		the default associator, Apriori by default
   */
  public static Object getAssociator() {
    Object	result;
    
    result = getObject(
		"Associator", 
		weka.associations.Apriori.class.getName(), 
		weka.associations.Associator.class);
    if (result == null)
      result = new weka.associations.Apriori();
    
    return result;
  }
  
  /**
   * returns the default attribute evalautor (fully configured) for the 
   * attribute selection panel.
   * 
   * @return		the default attribute evaluator, CfsSubsetEval by default
   */
  public static Object getASEvaluator() {
    Object	result;
    
    result = getObject(
		"ASEvaluation", 
		weka.attributeSelection.CfsSubsetEval.class.getName(), 
		weka.attributeSelection.ASEvaluation.class);
    if (result == null)
      result = new weka.attributeSelection.CfsSubsetEval();
    
    return result;
  }
  
  /**
   * returns the default attribute selection search scheme (fully configured) 
   * for the attribute selection panel.
   * 
   * @return		the default search scheme, BestFirst by default
   */
  public static Object getASSearch() {
    Object	result;
    
    result = getObject(
		"ASSearch", 
		weka.attributeSelection.BestFirst.class.getName(), 
		weka.attributeSelection.ASSearch.class);
    if (result == null)
      result = new weka.attributeSelection.BestFirst();
    
    return result;
  }
  
  /**
   * returns the default attribute selection test mode for the attribute
   * selection panel.
   * 
   * @return		the default test mode
   */
  public static int getASTestMode() {
    return Integer.parseInt(get("ASTestMode", "0"));
  }
  
  /**
   * returns the default number of folds of the CV in the attribute selection 
   * panel.
   * 
   * @return		the default number of folds
   */
  public static int getASCrossvalidationFolds() {
    return Integer.parseInt(get("ASCrossvalidationFolds", "10"));
  }
  
  /**
   * returns the default random seed value in the attribute selection panel.
   * 
   * @return		the default random seed
   */
  public static int getASRandomSeed() {
    return Integer.parseInt(get("ASRandomSeed", "1"));
  }
  
  /**
   * only for testing - prints the content of the props file.
   * 
   * @param args	commandline parameters - ignored
   */
  public static void main(String[] args) {
    Enumeration		names;
    String		name;
    Vector		sorted;
    
    System.out.println("\nExplorer defaults:");
    names = PROPERTIES.propertyNames();

    // sort names
    sorted = new Vector();
    while (names.hasMoreElements())
      sorted.add(names.nextElement());
    Collections.sort(sorted);
    names = sorted.elements();
    
    // output
    while (names.hasMoreElements()) {
      name = names.nextElement().toString();
      System.out.println("- " + name + ": " + PROPERTIES.getProperty(name, ""));
    }
    System.out.println();
  }
}

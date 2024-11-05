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
 * AssociatorEvaluation.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.associations;

import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Drawable;
import weka.core.converters.ConverterUtils.DataSource;

import java.util.Enumeration;
import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * Class for evaluating Associaters.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6496 $
 */
public class AssociatorEvaluation
  implements RevisionHandler {

  /** the result string */
  protected StringBuffer m_Result;
    
  /**
   * default constructor
   */
  public AssociatorEvaluation() {
    super();
    
    m_Result = new StringBuffer();
  }
  
  /**
   * Generates an option string to output on the commandline.
   * 
   * @param associator	the associator to generate the string for
   * @return		the option string
   */
  protected static String makeOptionString(Associator associator) {
    StringBuffer	text;
    
    text = new StringBuffer();   
    
    // general options
    text.append("\nGeneral options:\n\n");
    text.append("-t <training file>\n");
    text.append("\tThe name of the training file.\n");
    text.append("-g <name of graph file>\n");
    text.append("\tOutputs the graph representation (if supported) of the associator to a file.\n");
    
    // associator specific options, if any
    if (associator instanceof OptionHandler) {
      text.append(
	  "\nOptions specific to " 
	  + associator.getClass().getName().replaceAll(".*\\.", "") + ":\n\n");
      
      Enumeration enm = ((OptionHandler) associator).listOptions();
      while (enm.hasMoreElements()) {
	Option option = (Option) enm.nextElement();
	text.append(option.synopsis() + "\n");
	text.append(option.description() + "\n");
      }
    }
    
    return text.toString();
  }

  /**
   * Evaluates an associator with the options given in an array of strings.
   *
   * @param associatorString 	class of associator as a string
   * @param options 		the array of string containing the options
   * @throws Exception 		if model could not be evaluated successfully
   * @return 			a string describing the results 
   */
  public static String evaluate(String associatorString, String[] options) throws Exception {
    Associator associator;	 

    // Create associator
    try {
      associator = (Associator) Class.forName(associatorString).newInstance();
    }
    catch (Exception e) {
      throw new Exception("Can't find class with name " + associatorString + '.');
    }
    
    return evaluate(associator, options);
  }
  
  /**
   * Evaluates the associator with the given commandline options and returns
   * the evaluation string.
   * 
   * @param associator	the Associator to evaluate
   * @param options	the commandline options
   * @return		the generated output string
   * @throws Exception	if evaluation fails
   */
  public static String evaluate(Associator associator, String[] options) 
    throws Exception {

    String trainFileString = "";
    String graphFileName = "";
    AssociatorEvaluation eval;
    DataSource loader;

    // help?
    if (Utils.getFlag('h', options))
      throw new Exception("\nHelp requested.\n" + makeOptionString(associator));
    
    try {
      // general options
      trainFileString = Utils.getOption('t', options);
      if (trainFileString.length() == 0) 
	throw new Exception("No training file given!");
      loader = new DataSource(trainFileString);

      graphFileName = Utils.getOption('g', options);

      // associator specific options
      if (associator instanceof OptionHandler) {
        ((OptionHandler) associator).setOptions(options);
      }
      
      // left-over options?
      Utils.checkForRemainingOptions(options);
    }
    catch (Exception e) {
      throw new Exception(
	  "\nWeka exception: " 
	  + e.getMessage() + "\n" 
	  + makeOptionString(associator));
    }
    
    // load file and build associations
    eval = new AssociatorEvaluation();
    String results = eval.evaluate(associator, new Instances(loader.getDataSet()));

    // If associator is drawable output string describing graph
    if ((associator instanceof Drawable) && (graphFileName.length() != 0)) {
      BufferedWriter writer = new BufferedWriter(new FileWriter(graphFileName));
      writer.write(((Drawable) associator).graph());
      writer.newLine();
      writer.flush();
      writer.close();
    }
    
    return results;
  }
  
  /**
   * Evaluates the associator with the given commandline options and returns
   * the evaluation string.
   * 
   * @param associator	the Associator to evaluate
   * @param data	the data to run the associator with
   * @return		the generated output string
   * @throws Exception	if evaluation fails
   */
  public String evaluate(Associator associator, Instances data) 
    throws Exception {
    
    long startTime;
    long endTime;
    
    // build associations
    startTime = System.currentTimeMillis();
    associator.buildAssociations(data);
    endTime = System.currentTimeMillis();

    m_Result = new StringBuffer(associator.toString());
    m_Result.append("\n=== Evaluation ===\n\n");
    m_Result.append("Elapsed time: " + (((double) (endTime - startTime)) / 1000) + "s");
    m_Result.append("\n");
    
    return m_Result.toString();
  }

  /**
   * Tests whether the current evaluation object is equal to another
   * evaluation object
   *
   * @param obj the object to compare against
   * @return true if the two objects are equal
   */
  public boolean equals(Object obj) {
    if ((obj == null) || !(obj.getClass().equals(this.getClass())))
      return false;
    
    AssociatorEvaluation cmp = (AssociatorEvaluation) obj;
    
    // TODO: better comparison???
    String associatingResults1 = m_Result.toString().replaceAll("Elapsed time.*", "");
    String associatingResults2 = cmp.m_Result.toString().replaceAll("Elapsed time.*", "");
    if (!associatingResults1.equals(associatingResults2)) 
      return false;
    
    return true;
  }
  
  /**
   * returns a summary string of the evaluation with a no title
   * 
   * @return		the summary string
   */
  public String toSummaryString() {
    return toSummaryString("");
  }
  
  /**
   * returns a summary string of the evaluation with a default title
   * 
   * @param title	the title to print before the result
   * @return		the summary string
   */
  public String toSummaryString(String title) {
    StringBuffer	result;
    
    result = new StringBuffer(title);
    if (title.length() != 0)
      result.append("\n");
    result.append(m_Result);
    
    return result.toString();
  }
  
  /**
   * returns the current result
   * 
   * @return		the currently stored result
   * @see		#toSummaryString()
   */
  public String toString() {
    return toSummaryString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6496 $");
  }

  /**
   * A test method for this class. Just extracts the first command line
   * argument as an associator class name and calls evaluate.
   * 
   * @param args 	an array of command line arguments, the first of which
   * 			must be the class name of an associator.
   */
  public static void main(String[] args) {
    try {
      if (args.length == 0) {
	throw new Exception(
	    "The first argument must be the class name of a kernel");
      }
      String associator = args[0];
      args[0] = "";
      System.out.println(evaluate(associator, args));
    }
    catch (Exception ex) {
      ex.printStackTrace();
      System.err.println(ex.getMessage());
    }
  }
}

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
 *    Run.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka;

import java.util.ArrayList;

/**
 * Helper class that executes Weka schemes from the command line. Performs
 * Suffix matching on the scheme name entered by the user - e.g.<br><br>
 * 
 * java weka.Run NaiveBayes <br><br>
 * 
 * will prompt the user to choose among weka.classifiers.bayes.ComplementNaiveBayes,
 * weka.classifiers.bayes.NaiveBayes, weka.classifiers.bayes.NaiveBayesMultinomial,
 * weka.classifiers.bayes.NaiveBayesMultinomialUpdateable, weka.classifiers.bayes.NaiveBayesSimple,
 * weka.classifiers.bayes.NaiveBayesUpdateable
 * 
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: 6676 $
 *
 */
public class Run {
  
  public enum SchemeType {
    CLASSIFIER("classifier"),
    CLUSTERER("clusterer"),
    ASSOCIATOR("association rules"),
    ATTRIBUTE_SELECTION("attribute selection"),
    FILTER("filter");
    
    private final String m_stringVal;
    
    SchemeType(String name) {
      m_stringVal = name;
    }
    
    public String toString() {
      return m_stringVal;
    }
  }
  
  /**
   * Main method for this class. -help or -h prints usage info.
   * 
   * @param args
   */
  public static void main(String[] args) {
    try {
      if (args.length == 0 || args[0].equalsIgnoreCase("-h") ||
          args[0].equalsIgnoreCase("-help")) {
        System.err.println("Usage:\n\tweka.Run [-no-scan] [-no-load] <scheme name [scheme options]>");
        System.exit(1);
      }
      boolean noScan = false;
      boolean noLoad = false;
      if (args[0].equals("-list-packages")) {
        weka.core.WekaPackageManager.loadPackages(true);
        System.exit(0);
      } else if (args[0].equals("-no-load")) {
        noLoad = true;
        if (args.length > 1) {
          if (args[1].equals("-no-scan")) {
            noScan = true;
          }
        }
      } else if (args[0].equals("-no-scan")) {
        noScan = true;
        if (args.length > 1) {
          if (args[1].equals("-no-load")) {
            noLoad = true;
          }
        }
      }
      
      if (!noLoad) {
        weka.core.WekaPackageManager.loadPackages(false);
      }
      
      int schemeIndex = 0;
      if (noLoad && noScan) {
        schemeIndex = 2;
      } else if (noLoad || noScan) {
        schemeIndex = 1;
      }
      
      String schemeToRun = null;
      String[] options = null;
      
      if (schemeIndex >= args.length) {
        System.err.println("No scheme name given.");
        System.exit(1);
      }
      schemeToRun = args[schemeIndex];
      options = new String[args.length - schemeIndex - 1];
      if (options.length > 0) {
        System.arraycopy(args, schemeIndex + 1, options, 0, options.length);
      }
      
           
      if (!noScan) {     
        ArrayList<String> matches = weka.core.ClassDiscovery.find(schemeToRun);
        ArrayList<String> prunedMatches = new ArrayList<String>();
        // prune list for anything that isn't a runnable scheme      
        for (int i = 0; i < matches.size(); i++) {
          try {
            Object scheme = java.beans.Beans.instantiate((new Run()).getClass().getClassLoader(),
                matches.get(i));          
            if (scheme instanceof weka.classifiers.Classifier ||
                scheme instanceof weka.clusterers.Clusterer ||
                scheme instanceof weka.associations.Associator ||
                scheme instanceof weka.attributeSelection.ASEvaluation ||
                scheme instanceof weka.filters.Filter) {
              prunedMatches.add(matches.get(i));
            }
          } catch (Exception ex) {
            // ignore any classes that we can't instantiate due to no no-arg constructor
          }
        }

        if (prunedMatches.size() == 0) {
          System.err.println("Can't find scheme " + schemeToRun + ", or it is not runnable.");
          System.exit(1);
        } else if (prunedMatches.size() > 1) {
          java.io.BufferedReader br = 
            new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
          boolean done = false;
          while (!done) {
            System.out.println("Select a scheme to run, or <return> to exit:");
            for (int i = 0; i < prunedMatches.size(); i++) {
              System.out.println("\t" + (i+1) + ") " + prunedMatches.get(i));
            }
            System.out.print("\nEnter a number > ");
            String choice = null;
            int schemeNumber = 0;
            try {
              choice = br.readLine();
              if (choice.equals("")) {
                System.exit(0);
              } else {
                schemeNumber = Integer.parseInt(choice);
                schemeNumber--;
                if (schemeNumber >= 0 && schemeNumber < prunedMatches.size()) {
                  schemeToRun = prunedMatches.get(schemeNumber);
                  done = true;
                }
              }
            } catch (java.io.IOException ex) {
              // ignore
            }
          }
        } else {
          schemeToRun = prunedMatches.get(0);
        }
      }

      Object scheme = null;
      try {
        scheme = java.beans.Beans.instantiate((new Run()).getClass().getClassLoader(),
            schemeToRun);
      } catch (Exception ex) {
        System.err.println(schemeToRun + " is not runnable!");
        System.exit(1);
      }
      // now see which interfaces/classes this scheme implements/extends
      ArrayList<SchemeType> types = new ArrayList<SchemeType>();      
      if (scheme instanceof weka.classifiers.Classifier) {
        types.add(SchemeType.CLASSIFIER);
      }
      if (scheme instanceof weka.clusterers.Clusterer) {
        types.add(SchemeType.CLUSTERER);
      }
      if (scheme instanceof weka.associations.Associator) {
        types.add(SchemeType.ASSOCIATOR);
      }
      if (scheme instanceof weka.attributeSelection.ASEvaluation) {
        types.add(SchemeType.ATTRIBUTE_SELECTION);
      }
      if (scheme instanceof weka.filters.Filter) {
        types.add(SchemeType.FILTER);
      }
      
      SchemeType selectedType = null;
      if (types.size() == 0) {
        System.err.println("" + schemeToRun + " is not runnable!");
        System.exit(1);
      }
      if (types.size() == 1) {
        selectedType = types.get(0);
      } else {
        java.io.BufferedReader br = 
          new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
        boolean done = false;
        while (!done) {
          System.out.println("" + schemeToRun + " can be executed as any of the following:");
          for (int i = 0; i < types.size(); i++) {
            System.out.println("\t" + (i+1) + ") " + types.get(i));
          }
          System.out.print("\nEnter a number > ");
          String choice = null;
          int typeNumber = 0;
          try {
            choice = br.readLine();
            if (choice.equals("")) {
              System.exit(0);
            } else {
              typeNumber = Integer.parseInt(choice);
              typeNumber--;
              if (typeNumber >= 0 && typeNumber < types.size()) {
                selectedType = types.get(typeNumber);
                done = true;
              }
            }
          } catch (java.io.IOException ex) {
            // ignore
          }
        }
      }
            
      if (selectedType == SchemeType.CLASSIFIER) {
        weka.classifiers.AbstractClassifier.runClassifier((weka.classifiers.Classifier)scheme, options);
      } else if (selectedType == SchemeType.CLUSTERER) {
        weka.clusterers.AbstractClusterer.runClusterer((weka.clusterers.Clusterer)scheme, options);
      } else if (selectedType == SchemeType.ATTRIBUTE_SELECTION) {
        weka.attributeSelection.ASEvaluation.runEvaluator((weka.attributeSelection.ASEvaluation)scheme, options);
      } else if (selectedType == SchemeType.ASSOCIATOR) {
        weka.associations.AbstractAssociator.runAssociator((weka.associations.Associator)scheme, options);
      } else if (selectedType == SchemeType.FILTER) {
        weka.filters.Filter.runFilter((weka.filters.Filter)scheme, options);
      }
    } 
    catch (Exception e) {
      if (    ((e.getMessage() != null) && (e.getMessage().indexOf("General options") == -1))
	   || (e.getMessage() == null) )
	e.printStackTrace();
      else
	System.err.println(e.getMessage());
    }
  }
}
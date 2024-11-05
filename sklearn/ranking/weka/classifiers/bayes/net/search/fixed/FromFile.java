/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * FromFile.java
 * Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 * 
 */
package weka.classifiers.bayes.net.search.fixed;

import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.net.BIFReader;
import weka.classifiers.bayes.net.ParentSet;
import weka.classifiers.bayes.net.search.SearchAlgorithm;
import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;

import java.util.Enumeration;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * The FromFile reads the structure of a Bayes net from a file in BIFF format.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -B &lt;BIF File&gt;
 *  Name of file containing network structure in BIF format
 * </pre>
 * 
 <!-- options-end -->
 * 
 * @author Remco Bouckaert
 * @version $Revision: 1.8 $
 */
public class FromFile 
	extends SearchAlgorithm {
  
  	/** for serialization */
  	static final long serialVersionUID = 7334358169507619525L;
  
	/** name of file to read structure from **/
	String m_sBIFFile = "";
	    
	/**
	 * Returns a string describing this object
	 * @return a description of the classifier suitable for
	 * displaying in the explorer/experimenter gui
	 */
	public String globalInfo() {
	  return 
	      "The FromFile reads the structure of a Bayes net from a file "
	    + "in BIFF format.";
	}

	/**
	 * 
	 * @param bayesNet
	 * @param instances the instances to work with
	 * @throws Exception if attribute from BIF file could not be found
	 */
	public void buildStructure (BayesNet bayesNet, Instances instances) throws Exception {
		// read network structure in BIF format
		BIFReader bifReader = new BIFReader();
		bifReader.processFile(m_sBIFFile);
		// copy parent sets
        for (int iAttribute = 0; iAttribute < instances.numAttributes(); iAttribute++) {
            int iBIFAttribute = bifReader.getNode(bayesNet.getNodeName(iAttribute));
            ParentSet bifParentSet = bifReader.getParentSet(iBIFAttribute);
        	for (int iBIFParent = 0; iBIFParent < bifParentSet.getNrOfParents(); iBIFParent++) {
        	    String sParent = bifReader.getNodeName(bifParentSet.getParent(iBIFParent));
        	    int iParent = 0;
        	    while (iParent < instances.numAttributes() && !bayesNet.getNodeName(iParent).equals(sParent)) {
        	        iParent++;
        	    }
        	    if (iParent >= instances.numAttributes()) {
        	        throw new Exception("Could not find attribute " + sParent + " from BIF file in data");
        	    }
        		bayesNet.getParentSet(iAttribute).addParent(iParent, instances);
        	}
        }
	} // buildStructure

    /**
     * Set name of network in BIF file to read structure from
     * 
     * @param sBIFFile the name of the BIF file
     */
    public void setBIFFile(String sBIFFile) {
    	m_sBIFFile = sBIFFile;
    }

    /**
     * Get name of network in BIF file to read structure from
     * @return BIF file name
     */
    public String getBIFFile() {
        return m_sBIFFile;
    }

	/**
	 * Returns an enumeration describing the available options.
	 *
	 * @return an enumeration of all the available options.
	 */
	public Enumeration listOptions() {
	  Vector newVector = new Vector();

	  newVector.addElement(new Option("\tName of file containing network structure in BIF format\n", 
					 "B", 1, "-B <BIF File>"));

          Enumeration en = super.listOptions();
          while (en.hasMoreElements())
            newVector.addElement(en.nextElement());
          
	  return newVector.elements();
	}

	/**
	 * Parses a given list of options. <p/>
	 *
	 <!-- options-start -->
	 * Valid options are: <p/>
	 * 
	 * <pre> -B &lt;BIF File&gt;
	 *  Name of file containing network structure in BIF format
	 * </pre>
	 * 
	 <!-- options-end -->
	 *
	 * @param options the list of options as an array of strings
	 * @throws Exception if an option is not supported
	 */
	public void setOptions(String[] options) throws Exception {
	  setBIFFile( Utils.getOption('B', options));
          
          super.setOptions(options);
	}

	/**
	 * Gets the current settings of the search algorithm.
	 *
	 * @return an array of strings suitable for passing to setOptions
	 */
	public String [] getOptions() {
          String[] superOptions = super.getOptions();
	  String [] options  = new String [2 + superOptions.length];
	  int current = 0;

          options[current++] = "-B";
	  options[current++] = "" + getBIFFile();

          // insert options from parent class
          for (int iOption = 0; iOption < superOptions.length; iOption++) {
                  options[current++] = superOptions[iOption];
          }

	  // Fill up rest with empty strings, not nulls!
	  while (current < options.length) {
		options[current++] = "";
	  }
	  return options;
	}

	/**
	 * Returns the revision string.
	 * 
	 * @return		the revision
	 */
	public String getRevision() {
	  return RevisionUtils.extract("$Revision: 1.8 $");
	}

} // class FromFile

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
 *    ClusterEvaluation.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package  weka.clusterers;

import weka.core.Drawable;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Range;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.MethodDescriptor;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/**
 * Class for evaluating clustering models.<p/>
 *
 * Valid options are: <p/>
 *
 * -t name of the training file <br/>
 * Specify the training file. <p/>
 *
 * -T name of the test file <br/>
 * Specify the test file to apply clusterer to. <p/>
 *
 * -d name of file to save clustering model to <br/>
 * Specify output file. <p/>
 *
 * -l name of file to load clustering model from <br/>
 * Specifiy input file. <p/>
 *
 * -p attribute range <br/>
 * Output predictions. Predictions are for the training file if only the
 * training file is specified, otherwise they are for the test file. The range
 * specifies attribute values to be output with the predictions.
 * Use '-p 0' for none. <p/>
 *
 * -x num folds <br/>
 * Set the number of folds for a cross validation of the training data.
 * Cross validation can only be done for distribution clusterers and will
 * be performed if the test file is missing. <p/>
 * 
 * -s num <br/>
 * Sets the seed for randomizing the data for cross-validation. <p/>
 *
 * -c class <br/>
 * Set the class attribute. If set, then class based evaluation of clustering
 * is performed. <p/>
 * 
 * -g name of graph file <br/>
 * Outputs the graph representation of the clusterer to the file. Only for
 * clusterer that implemented the <code>weka.core.Drawable</code> interface.
 * <p/>
 *
 * @author   Mark Hall (mhall@cs.waikato.ac.nz)
 * @version  $Revision: 6021 $
 * @see	     weka.core.Drawable
 */
public class ClusterEvaluation 
  implements Serializable, RevisionHandler {

  /** for serialization */
  static final long serialVersionUID = -830188327319128005L;
  
  /** the clusterer */
  private Clusterer m_Clusterer;

  /** holds a string describing the results of clustering the training data */
  private StringBuffer m_clusteringResults;

  /** holds the number of clusters found by the clusterer */
  private int m_numClusters;

  /** holds the assigments of instances to clusters for a particular testing
      dataset */
  private double[] m_clusterAssignments;

  /** holds the average log likelihood for a particular testing dataset
     if the clusterer is a DensityBasedClusterer */
  private double m_logL;

  /** will hold the mapping of classes to clusters (for class based 
      evaluation) */
  private int[] m_classToCluster = null;

  /**
   * set the clusterer
   * @param clusterer the clusterer to use
   */
  public void setClusterer(Clusterer clusterer) {
    m_Clusterer = clusterer;
  }

  /**
   * return the results of clustering.
   * @return a string detailing the results of clustering a data set
   */
  public String clusterResultsToString() {
    return m_clusteringResults.toString();
  }

  /**
   * Return the number of clusters found for the most recent call to
   * evaluateClusterer
   * @return the number of clusters found
   */
  public int getNumClusters() {
    return m_numClusters;
  }

  /**
   * Return an array of cluster assignments corresponding to the most
   * recent set of instances clustered.
   * @return an array of cluster assignments
   */
  public double[] getClusterAssignments() {
    return m_clusterAssignments;
  }

  /**
   * Return the array (ordered by cluster number) of minimum error class to
   * cluster mappings
   * @return an array of class to cluster mappings
   */
  public int[] getClassesToClusters() {
    return m_classToCluster;
  }

  /**
   * Return the log likelihood corresponding to the most recent
   * set of instances clustered.
   *
   * @return a <code>double</code> value
   */
  public double getLogLikelihood() {
    return m_logL;
  }

  /**
   * Constructor. Sets defaults for each member variable. Default Clusterer
   * is EM.
   */
  public ClusterEvaluation () {
    setClusterer(new SimpleKMeans());
    m_clusteringResults = new StringBuffer();
    m_clusterAssignments = null;
  }

  /**
   * Evaluate the clusterer on a set of instances. Calculates clustering
   * statistics and stores cluster assigments for the instances in
   * m_clusterAssignments
   * 
   * @param test the set of instances to cluster
   * @throws Exception if something goes wrong
   */
  public void evaluateClusterer(Instances test) throws Exception {
    evaluateClusterer(test, "");
  }

  /**
   * Evaluate the clusterer on a set of instances. Calculates clustering
   * statistics and stores cluster assigments for the instances in
   * m_clusterAssignments
   * 
   * @param test the set of instances to cluster
   * @param testFileName the name of the test file for incremental testing, 
   * if "" or null then not used
   * @throws Exception if something goes wrong
   */
  public void evaluateClusterer(Instances test, String testFileName) throws Exception {
    int i = 0;
    int cnum;
    double loglk = 0.0;
    int cc = m_Clusterer.numberOfClusters();
    m_numClusters = cc;
    double[] instanceStats = new double[cc];
    Instances testRaw = null;
    boolean hasClass = (test.classIndex() >= 0);
    int unclusteredInstances = 0;
    Vector<Double> clusterAssignments = new Vector<Double>();
    Filter filter = null;
    DataSource source = null;
    Instance inst;

    if (testFileName == null)
      testFileName = "";
    
    // load data
    if (testFileName.length() != 0)
      source = new DataSource(testFileName);
    else
      source = new DataSource(test);
    testRaw = source.getStructure(test.classIndex());
    
    // If class is set then do class based evaluation as well
    if (hasClass) {
      if (testRaw.classAttribute().isNumeric())
	throw new Exception("ClusterEvaluation: Class must be nominal!");

      filter = new Remove();
      ((Remove) filter).setAttributeIndices("" + (testRaw.classIndex() + 1));
      ((Remove) filter).setInvertSelection(false);
      filter.setInputFormat(testRaw);
    }
    
    i = 0;
    while (source.hasMoreElements(testRaw)) {
      // next instance
      inst = source.nextElement(testRaw);
      if (filter != null) {
	filter.input(inst);
	filter.batchFinished();
	inst = filter.output();
      }
      
      cnum = -1;
      try {
	if (m_Clusterer instanceof DensityBasedClusterer) {
	  loglk += ((DensityBasedClusterer)m_Clusterer).
	    logDensityForInstance(inst);
	  cnum = m_Clusterer.clusterInstance(inst); 
	  clusterAssignments.add((double) cnum);
	}
	else {
	  cnum = m_Clusterer.clusterInstance(inst);
	  clusterAssignments.add((double) cnum);
	}
      }
      catch (Exception e) {
	clusterAssignments.add(-1.0);
	unclusteredInstances++;
      }
      
      if (cnum != -1) {
	instanceStats[cnum]++;
      }
    }
    
    double sum = Utils.sum(instanceStats);
    loglk /= sum;
    m_logL = loglk;
    m_clusterAssignments = new double [clusterAssignments.size()];
    for (i = 0; i < clusterAssignments.size(); i++) {
      m_clusterAssignments[i] = clusterAssignments.get(i);
    }
    int numInstFieldWidth = (int)((Math.log(clusterAssignments.size())/Math.log(10))+1);
    
    m_clusteringResults.append(m_Clusterer.toString());
    m_clusteringResults.append("Clustered Instances\n\n");
    int clustFieldWidth = (int)((Math.log(cc)/Math.log(10))+1);
    for (i = 0; i < cc; i++) {
      if (instanceStats[i] > 0)
	m_clusteringResults.append(Utils.doubleToString((double)i, 
							clustFieldWidth, 0) 
				   + "      " 
				   + Utils.doubleToString(instanceStats[i],
							  numInstFieldWidth, 0) 
				   + " (" 
				   + Utils.doubleToString((instanceStats[i] / 
							   sum * 100.0)
							  , 3, 0) + "%)\n");
    }
    
    if (unclusteredInstances > 0)
      m_clusteringResults.append("\nUnclustered instances : "
				 +unclusteredInstances);

    if (m_Clusterer instanceof DensityBasedClusterer)
      m_clusteringResults.append("\n\nLog likelihood: " 
				 + Utils.doubleToString(loglk, 1, 5) 
				 + "\n");       
    
    if (hasClass) {
      evaluateClustersWithRespectToClass(test, testFileName);
    }
  }

  /**
   * Evaluates cluster assignments with respect to actual class labels.
   * Assumes that m_Clusterer has been trained and tested on 
   * inst (minus the class).
   * 
   * @param inst the instances (including class) to evaluate with respect to
   * @param fileName the name of the test file for incremental testing, 
   * if "" or null then not used
   * @throws Exception if something goes wrong
   */
  private void evaluateClustersWithRespectToClass(Instances inst, String fileName)
    throws Exception {
    
    
    
    int numClasses = inst.classAttribute().numValues();
    int[][] counts = new int [m_numClusters][numClasses];
    int[] clusterTotals = new int[m_numClusters];
    double[] best = new double[m_numClusters+1];
    double[] current = new double[m_numClusters+1];
    DataSource source = null;
    Instances instances = null;
    Instance instance = null;
    int i;
    int numInstances;
        

    if (fileName == null)
      fileName = "";
    
    if (fileName.length() != 0) {
      source = new DataSource(fileName);
    }
    else
      source = new DataSource(inst);
    instances = source.getStructure(inst.classIndex());

    i = 0;
    while (source.hasMoreElements(instances)) {
      instance = source.nextElement(instances);
      if (m_clusterAssignments[i] >= 0) {
        counts[(int)m_clusterAssignments[i]][(int)instance.classValue()]++;
        clusterTotals[(int)m_clusterAssignments[i]]++;        
      }
      i++;
    }
    numInstances = i;
   
    best[m_numClusters] = Double.MAX_VALUE;
    mapClasses(m_numClusters, 0, counts, clusterTotals, current, best, 0);

    m_clusteringResults.append("\n\nClass attribute: "
			+inst.classAttribute().name()
			+"\n");
    m_clusteringResults.append("Classes to Clusters:\n");
    String matrixString = toMatrixString(counts, clusterTotals, new Instances(inst, 0));
    m_clusteringResults.append(matrixString).append("\n");

    int Cwidth = 1 + (int)(Math.log(m_numClusters) / Math.log(10));
    // add the minimum error assignment
    for (i = 0; i < m_numClusters; i++) {
      if (clusterTotals[i] > 0) {
	m_clusteringResults.append("Cluster "
				   +Utils.doubleToString((double)i,Cwidth,0));
	m_clusteringResults.append(" <-- ");
	
	if (best[i] < 0) {
	  m_clusteringResults.append("No class\n");
	} else {
	  m_clusteringResults.
	    append(inst.classAttribute().value((int)best[i])).append("\n");
	}
      }
    }
    m_clusteringResults.append("\nIncorrectly clustered instances :\t"
			       +best[m_numClusters]+"\t"
			       +(Utils.doubleToString((best[m_numClusters] / 
						       numInstances * 
						       100.0), 8, 4))
			       +" %\n");

    // copy the class assignments
    m_classToCluster = new int [m_numClusters];
    for (i = 0; i < m_numClusters; i++) {
      m_classToCluster[i] = (int)best[i];
    }
  }

  /**
   * Returns a "confusion" style matrix of classes to clusters assignments
   * @param counts the counts of classes for each cluster
   * @param clusterTotals total number of examples in each cluster
   * @param inst the training instances (with class)
   * @return the "confusion" style matrix as string
   * @throws Exception if matrix can't be generated
   */
  private String toMatrixString(int[][] counts, int[] clusterTotals,
				Instances inst) 
    throws Exception {
    StringBuffer ms = new StringBuffer();

    int maxval = 0;
    for (int i = 0; i < m_numClusters; i++) {
      for (int j = 0; j < counts[i].length; j++) {
	if (counts[i][j] > maxval) {
	  maxval = counts[i][j];
	}
      }
    }

    int Cwidth = 1 + Math.max((int)(Math.log(maxval) / Math.log(10)),
			      (int)(Math.log(m_numClusters) / Math.log(10)));

    ms.append("\n");
    
    for (int i = 0; i < m_numClusters; i++) {
      if (clusterTotals[i] > 0) {
	ms.append(" ").append(Utils.doubleToString((double)i, Cwidth, 0));
      }
    }
    ms.append("  <-- assigned to cluster\n");
    
    for (int i = 0; i< counts[0].length; i++) {

      for (int j = 0; j < m_numClusters; j++) {
	if (clusterTotals[j] > 0) {
	  ms.append(" ").append(Utils.doubleToString((double)counts[j][i], 
						     Cwidth, 0));
	}
      }
      ms.append(" | ").append(inst.classAttribute().value(i)).append("\n");
    }

    return ms.toString();
  }

  /**
   * Finds the minimum error mapping of classes to clusters. Recursively
   * considers all possible class to cluster assignments.
   * 
   * @param numClusters the number of clusters
   * @param lev the cluster being processed
   * @param counts the counts of classes in clusters
   * @param clusterTotals the total number of examples in each cluster
   * @param current the current path through the class to cluster assignment
   * tree
   * @param best the best assignment path seen
   * @param error accumulates the error for a particular path
   */
  public static void mapClasses(int numClusters, int lev, int[][] counts, int[] clusterTotals,
			  double[] current, double[] best, int error) {
    // leaf
    if (lev == numClusters) {
      if (error < best[numClusters]) {
	best[numClusters] = error;
	for (int i = 0; i < numClusters; i++) {
	  best[i] = current[i];
	}
      }
    } else {
      // empty cluster -- ignore
      if (clusterTotals[lev] == 0) {
	current[lev] = -1; // cluster ignored
	mapClasses(numClusters, lev+1, counts, clusterTotals, current, best,
		   error);
      } else {
	// first try no class assignment to this cluster
	current[lev] = -1; // cluster assigned no class (ie all errors)
	mapClasses(numClusters, lev+1, counts, clusterTotals, current, best,
		   error+clusterTotals[lev]);
	// now loop through the classes in this cluster
	for (int i = 0; i < counts[0].length; i++) {
	  if (counts[lev][i] > 0) {
	    boolean ok = true;
	    // check to see if this class has already been assigned
	    for (int j = 0; j < lev; j++) {
	      if ((int)current[j] == i) {
		ok = false;
		break;
	      }
	    }
	    if (ok) {
	      current[lev] = i;
	      mapClasses(numClusters, lev+1, counts, clusterTotals, current, best, 
			 (error + (clusterTotals[lev] - counts[lev][i])));
	    }
	  }
	}
      }
    }
  }

  /**
   * Evaluates a clusterer with the options given in an array of
   * strings. It takes the string indicated by "-t" as training file, the
   * string indicated by "-T" as test file.
   * If the test file is missing, a stratified ten-fold
   * cross-validation is performed (distribution clusterers only).
   * Using "-x" you can change the number of
   * folds to be used, and using "-s" the random seed.
   * If the "-p" option is present it outputs the classification for
   * each test instance. If you provide the name of an object file using
   * "-l", a clusterer will be loaded from the given file. If you provide the
   * name of an object file using "-d", the clusterer built from the
   * training data will be saved to the given file.
   *
   * @param clusterer machine learning clusterer
   * @param options the array of string containing the options
   * @throws Exception if model could not be evaluated successfully
   * @return a string describing the results 
   */
  public static String evaluateClusterer(Clusterer clusterer, String[] options)
    throws Exception {
    
    int seed = 1, folds = 10;
    boolean doXval = false;
    Instances train = null;
    Random random;
    String trainFileName, testFileName, seedString, foldsString;
    String objectInputFileName, objectOutputFileName, attributeRangeString;
    String graphFileName;
    String[] savedOptions = null;
    boolean printClusterAssignments = false;
    Range attributesToOutput = null;
    StringBuffer text = new StringBuffer();
    int theClass = -1; // class based evaluation of clustering
    boolean updateable = (clusterer instanceof UpdateableClusterer);
    DataSource source = null;
    Instance inst;

    if (Utils.getFlag('h', options) || Utils.getFlag("help", options)) {
      
      // global info requested as well?
      boolean globalInfo = Utils.getFlag("synopsis", options) ||
        Utils.getFlag("info", options);
      
      throw  new Exception("Help requested." 
          + makeOptionString(clusterer, globalInfo));
    }
    
    try {
      // Get basic options (options the same for all clusterers
      //printClusterAssignments = Utils.getFlag('p', options);
      objectInputFileName = Utils.getOption('l', options);
      objectOutputFileName = Utils.getOption('d', options);
      trainFileName = Utils.getOption('t', options);
      testFileName = Utils.getOption('T', options);
      graphFileName = Utils.getOption('g', options);

      // Check -p option
      try {
	attributeRangeString = Utils.getOption('p', options);
      }
      catch (Exception e) {
	throw new Exception(e.getMessage() + "\nNOTE: the -p option has changed. " +
			    "It now expects a parameter specifying a range of attributes " +
			    "to list with the predictions. Use '-p 0' for none.");
      }
      if (attributeRangeString.length() != 0) {
	printClusterAssignments = true;
	if (!attributeRangeString.equals("0")) 
	  attributesToOutput = new Range(attributeRangeString);
      }

      if (trainFileName.length() == 0) {
        if (objectInputFileName.length() == 0) {
          throw  new Exception("No training file and no object " 
			       + "input file given.");
        }

        if (testFileName.length() == 0) {
          throw  new Exception("No training file and no test file given.");
        }
      }
      else {
	if ((objectInputFileName.length() != 0) 
	    && (printClusterAssignments == false)) {
	  throw  new Exception("Can't use both train and model file " 
			       + "unless -p specified.");
	}
      }

      seedString = Utils.getOption('s', options);

      if (seedString.length() != 0) {
	seed = Integer.parseInt(seedString);
      }

      foldsString = Utils.getOption('x', options);

      if (foldsString.length() != 0) {
	folds = Integer.parseInt(foldsString);
	doXval = true;
      }
    }
    catch (Exception e) {
      throw  new Exception('\n' + e.getMessage() 
			   + makeOptionString(clusterer, false));
    }

    try {
      if (trainFileName.length() != 0) {
	source = new DataSource(trainFileName);
	train  = source.getStructure();

	String classString = Utils.getOption('c',options);
	if (classString.length() != 0) {
	  if (classString.compareTo("last") == 0)
	    theClass = train.numAttributes();
	  else if (classString.compareTo("first") == 0)
	    theClass = 1;
	  else
	    theClass = Integer.parseInt(classString);

	  if (theClass != -1) {
	    if (doXval || testFileName.length() != 0)
	      throw new Exception("Can only do class based evaluation on the "
		  +"training data");

	    if (objectInputFileName.length() != 0)
	      throw new Exception("Can't load a clusterer and do class based "
		  +"evaluation");

	    if (objectOutputFileName.length() != 0)
	      throw new Exception(
		  "Can't do class based evaluation and save clusterer");
	  }
	}
	else {
	  // if the dataset defines a class attribute, use it
	  if (train.classIndex() != -1) {
	    theClass = train.classIndex() + 1;
	    System.err.println(
		"Note: using class attribute from dataset, i.e., attribute #" 
		+ theClass);
	  }
	}

	if (theClass != -1) {
	  if (theClass < 1 || theClass > train.numAttributes())
	    throw new Exception("Class is out of range!");

	  if (!train.attribute(theClass - 1).isNominal())
	    throw new Exception("Class must be nominal!");
	  
	  train.setClassIndex(theClass - 1);
	}
      }
    }
    catch (Exception e) {
      throw  new Exception("ClusterEvaluation: " + e.getMessage() + '.');
    }

    // Save options
    if (options != null) {
      savedOptions = new String[options.length];
      System.arraycopy(options, 0, savedOptions, 0, options.length);
    }

    if (objectInputFileName.length() != 0)
      Utils.checkForRemainingOptions(options);

    // Set options for clusterer
    if (clusterer instanceof OptionHandler)
      ((OptionHandler)clusterer).setOptions(options);

    Utils.checkForRemainingOptions(options);

    Instances trainHeader = train;
    if (objectInputFileName.length() != 0) {
      // Load the clusterer from file
      //      clusterer = (Clusterer) SerializationHelper.read(objectInputFileName);
      java.io.ObjectInputStream ois = 
        new java.io.ObjectInputStream(
        new java.io.BufferedInputStream(
        new java.io.FileInputStream(objectInputFileName)));
      clusterer = (Clusterer) ois.readObject();
      // try and get the training header
      try {
        trainHeader = (Instances) ois.readObject();
      } catch (Exception ex) {
        // don't moan if we cant
      }
    }
    else {
      // Build the clusterer if no object file provided
      if (theClass == -1) {
	if (updateable) {
	  clusterer.buildClusterer(source.getStructure());
	  while (source.hasMoreElements(train)) {
	    inst = source.nextElement(train);
	    ((UpdateableClusterer) clusterer).updateClusterer(inst);
	  }
	  ((UpdateableClusterer) clusterer).updateFinished();
	}
	else {
	  clusterer.buildClusterer(source.getDataSet());
	}
      }
      else {
	Remove removeClass = new Remove();
	removeClass.setAttributeIndices("" + theClass);
	removeClass.setInvertSelection(false);
	removeClass.setInputFormat(train);
	if (updateable) {
	  Instances clusterTrain = Filter.useFilter(train, removeClass);
	  clusterer.buildClusterer(clusterTrain);
          trainHeader = clusterTrain;
	  while (source.hasMoreElements(train)) {
	    inst = source.nextElement(train);
	    removeClass.input(inst);
	    removeClass.batchFinished();
	    Instance clusterTrainInst = removeClass.output();
	    ((UpdateableClusterer) clusterer).updateClusterer(clusterTrainInst);
	  }
	  ((UpdateableClusterer) clusterer).updateFinished();
	}
	else {
	  Instances clusterTrain = Filter.useFilter(source.getDataSet(), removeClass);
	  clusterer.buildClusterer(clusterTrain);
          trainHeader = clusterTrain;
	}
	ClusterEvaluation ce = new ClusterEvaluation();
	ce.setClusterer(clusterer);
	ce.evaluateClusterer(train, trainFileName);
	
	return "\n\n=== Clustering stats for training data ===\n\n" +
	  ce.clusterResultsToString();
      }
    }

    /* Output cluster predictions only (for the test data if specified,
       otherwise for the training data */
    if (printClusterAssignments) {
      return printClusterings(clusterer, trainFileName, testFileName, attributesToOutput);
    }

    text.append(clusterer.toString());
    text.append("\n\n=== Clustering stats for training data ===\n\n" 
		+ printClusterStats(clusterer, trainFileName));

    if (testFileName.length() != 0) {
      // check header compatibility
      DataSource test = new DataSource(testFileName);
      Instances testStructure = test.getStructure();
      if (!trainHeader.equalHeaders(testStructure)) {
        throw new Exception("Training and testing data are not compatible\n" + trainHeader.equalHeadersMsg(testStructure));
      }

      text.append("\n\n=== Clustering stats for testing data ===\n\n" 
		  + printClusterStats(clusterer, testFileName));
    }

    if ((clusterer instanceof DensityBasedClusterer) && 
	(doXval == true) && 
	(testFileName.length() == 0) && 
	(objectInputFileName.length() == 0)) {
      // cross validate the log likelihood on the training data
      random = new Random(seed);
      random.setSeed(seed);
      train = source.getDataSet();
      train.randomize(random);
      text.append(
	  crossValidateModel(
	      clusterer.getClass().getName(), train, folds, savedOptions, random));
    }

    // Save the clusterer if an object output file is provided
    if (objectOutputFileName.length() != 0) {
      //SerializationHelper.write(objectOutputFileName, clusterer);
      saveClusterer(objectOutputFileName, clusterer, trainHeader);
    }

    // If classifier is drawable output string describing graph
    if ((clusterer instanceof Drawable) && (graphFileName.length() != 0)) {
      BufferedWriter writer = new BufferedWriter(new FileWriter(graphFileName));
      writer.write(((Drawable) clusterer).graph());
      writer.newLine();
      writer.flush();
      writer.close();
    }
    
    return  text.toString();
  }

  private static void saveClusterer(String fileName, 
                             Clusterer clusterer, 
                             Instances header) throws Exception {
    java.io.ObjectOutputStream oos = 
      new java.io.ObjectOutputStream(
      new java.io.BufferedOutputStream(
      new java.io.FileOutputStream(fileName)));

    oos.writeObject(clusterer);
    if (header != null) {
      oos.writeObject(header);
    }
    oos.flush();
    oos.close();
  }

  /**
   * Perform a cross-validation for DensityBasedClusterer on a set of instances.
   *
   * @param clusterer the clusterer to use
   * @param data the training data
   * @param numFolds number of folds of cross validation to perform
   * @param random random number seed for cross-validation
   * @return the cross-validated log-likelihood
   * @throws Exception if an error occurs
   */
  public static double crossValidateModel(DensityBasedClusterer clusterer,
					  Instances data,
					  int numFolds,
					  Random random) throws Exception {
    Instances train, test;
    double foldAv = 0;;
    data = new Instances(data);
    data.randomize(random);
    //    double sumOW = 0;
    for (int i = 0; i < numFolds; i++) {
      // Build and test clusterer
      train = data.trainCV(numFolds, i, random);

      clusterer.buildClusterer(train);

      test = data.testCV(numFolds, i);
      
      for (int j = 0; j < test.numInstances(); j++) {
	try {
	  foldAv += ((DensityBasedClusterer)clusterer).
	    logDensityForInstance(test.instance(j));
	  //	  sumOW += test.instance(j).weight();
	  //	double temp = Utils.sum(tempDist);
	} catch (Exception ex) {
	  // unclustered instances
	}
      }
    }
   
    //    return foldAv / sumOW;
    return foldAv / data.numInstances();
  }

  /**
   * Performs a cross-validation 
   * for a DensityBasedClusterer clusterer on a set of instances.
   *
   * @param clustererString a string naming the class of the clusterer
   * @param data the data on which the cross-validation is to be 
   * performed 
   * @param numFolds the number of folds for the cross-validation
   * @param options the options to the clusterer
   * @param random a random number generator
   * @return a string containing the cross validated log likelihood
   * @throws Exception if a clusterer could not be generated 
   */
  public static String crossValidateModel (String clustererString, 
					   Instances data, 
					   int numFolds, 
					   String[] options,
					   Random random)
    throws Exception {
    Clusterer clusterer = null;
    String[] savedOptions = null;
    double CvAv = 0.0;
    StringBuffer CvString = new StringBuffer();

    if (options != null) {
      savedOptions = new String[options.length];
    }

    data = new Instances(data);

    // create clusterer
    try {
      clusterer = (Clusterer)Class.forName(clustererString).newInstance();
    }
    catch (Exception e) {
      throw  new Exception("Can't find class with name " 
			   + clustererString + '.');
    }

    if (!(clusterer instanceof DensityBasedClusterer)) {
      throw  new Exception(clustererString 
			   + " must be a distrinbution " 
			   + "clusterer.");
    }

    // Save options
    if (options != null) {
      System.arraycopy(options, 0, savedOptions, 0, options.length);
    }

    // Parse options
    if (clusterer instanceof OptionHandler) {
      try {
	((OptionHandler)clusterer).setOptions(savedOptions);
	Utils.checkForRemainingOptions(savedOptions);
      }
      catch (Exception e) {
	throw  new Exception("Can't parse given options in " 
			     + "cross-validation!");
      }
    }
    CvAv = crossValidateModel((DensityBasedClusterer)clusterer, data, numFolds, random);

    CvString.append("\n" + numFolds 
		    + " fold CV Log Likelihood: " 
		    + Utils.doubleToString(CvAv, 6, 4) 
		    + "\n");
    return  CvString.toString();
  }


  // ===============
  // Private methods
  // ===============
  /**
   * Print the cluster statistics for either the training
   * or the testing data.
   *
   * @param clusterer the clusterer to use for generating statistics.
   * @param fileName the file to load
   * @return a string containing cluster statistics.
   * @throws Exception if statistics can't be generated.
   */
  private static String printClusterStats (Clusterer clusterer, 
					   String fileName)
    throws Exception {
    StringBuffer text = new StringBuffer();
    int i = 0;
    int cnum;
    double loglk = 0.0;
    int cc = clusterer.numberOfClusters();
    double[] instanceStats = new double[cc];
    int unclusteredInstances = 0;

    if (fileName.length() != 0) {
      DataSource source = new DataSource(fileName);
      Instances structure = source.getStructure();
      Instance inst;
      while (source.hasMoreElements(structure)) {
	inst = source.nextElement(structure);
	try {
	  cnum = clusterer.clusterInstance(inst);

	  if (clusterer instanceof DensityBasedClusterer) {
	    loglk += ((DensityBasedClusterer)clusterer).
	      logDensityForInstance(inst);
	    //	    temp = Utils.sum(dist);
	  }
	  instanceStats[cnum]++;
	}
	catch (Exception e) {
	  unclusteredInstances++;
	}
	i++;
      }

      /*
      // count the actual number of used clusters
      int count = 0;
      for (i = 0; i < cc; i++) {
	if (instanceStats[i] > 0) {
	  count++;
	}
      }
      if (count > 0) {
	double[] tempStats = new double [count];
	count=0;
	for (i=0;i<cc;i++) {
	  if (instanceStats[i] > 0) {
	    tempStats[count++] = instanceStats[i];
	}
	}
	instanceStats = tempStats;
	cc = instanceStats.length;
	} */

      int clustFieldWidth = (int)((Math.log(cc)/Math.log(10))+1);
      int numInstFieldWidth = (int)((Math.log(i)/Math.log(10))+1);
      double sum = Utils.sum(instanceStats);
      loglk /= sum;
      text.append("Clustered Instances\n");

      for (i = 0; i < cc; i++) {
	if (instanceStats[i] > 0) {
	  text.append(Utils.doubleToString((double)i, 
					   clustFieldWidth, 0) 
		      + "      " 
		      + Utils.doubleToString(instanceStats[i], 
					     numInstFieldWidth, 0) 
		      + " (" 
		    + Utils.doubleToString((instanceStats[i]/sum*100.0)
					   , 3, 0) + "%)\n");
	}
      }
      if (unclusteredInstances > 0) {
	text.append("\nUnclustered Instances : "+unclusteredInstances);
      }

      if (clusterer instanceof DensityBasedClusterer) {
	text.append("\n\nLog likelihood: " 
		    + Utils.doubleToString(loglk, 1, 5) 
		    + "\n");
      }
    }

    return text.toString();
  }


  /**
   * Print the cluster assignments for either the training
   * or the testing data.
   *
   * @param clusterer the clusterer to use for cluster assignments
   * @param trainFileName the train file
   * @param testFileName an optional test file
   * @param attributesToOutput the attributes to print
   * @return a string containing the instance indexes and cluster assigns.
   * @throws Exception if cluster assignments can't be printed
   */
  private static String printClusterings (Clusterer clusterer, String trainFileName,
					  String testFileName, Range attributesToOutput)
    throws Exception {

    StringBuffer text = new StringBuffer();
    int i = 0;
    int cnum;
    DataSource source = null;
    Instance inst;
    Instances structure;
    
    if (testFileName.length() != 0)
      source = new DataSource(testFileName);
    else
      source = new DataSource(trainFileName);
    
    structure = source.getStructure();
    while (source.hasMoreElements(structure)) {
      inst = source.nextElement(structure);
      try {
	cnum = clusterer.clusterInstance(inst);
	
	text.append(i + " " + cnum + " "
	    + attributeValuesString(inst, attributesToOutput) + "\n");
      }
      catch (Exception e) {
	/*	  throw  new Exception('\n' + "Unable to cluster instance\n" 
	 + e.getMessage()); */
	text.append(i + " Unclustered "
	    + attributeValuesString(inst, attributesToOutput) + "\n");
      }
      i++;
    }
    
    return text.toString();
  }

  /**
   * Builds a string listing the attribute values in a specified range of indices,
   * separated by commas and enclosed in brackets.
   *
   * @param instance the instance to print the values from
   * @param attRange the range of the attributes to list
   * @return a string listing values of the attributes in the range
   */
  private static String attributeValuesString(Instance instance, Range attRange) {
    StringBuffer text = new StringBuffer();
    if (attRange != null) {
      boolean firstOutput = true;
      attRange.setUpper(instance.numAttributes() - 1);
      for (int i=0; i<instance.numAttributes(); i++)
	if (attRange.isInRange(i)) {
	  if (firstOutput) text.append("(");
	  else text.append(",");
	  text.append(instance.toString(i));
	  firstOutput = false;
	}
      if (!firstOutput) text.append(")");
    }
    return text.toString();
  }

  /**
   * Make up the help string giving all the command line options
   *
   * @param clusterer the clusterer to include options for
   * @return a string detailing the valid command line options
   */
  private static String makeOptionString (Clusterer clusterer,
                                          boolean globalInfo) {
    StringBuffer optionsText = new StringBuffer("");
    // General options
    optionsText.append("\n\nGeneral options:\n\n");
    optionsText.append("-h or -help\n");
    optionsText.append("\tOutput help information.\n");
    optionsText.append("-synopsis or -info\n");
    optionsText.append("\tOutput synopsis for clusterer (use in conjunction "
        + " with -h)\n");
    optionsText.append("-t <name of training file>\n");
    optionsText.append("\tSets training file.\n");
    optionsText.append("-T <name of test file>\n");
    optionsText.append("\tSets test file.\n");
    optionsText.append("-l <name of input file>\n");
    optionsText.append("\tSets model input file.\n");
    optionsText.append("-d <name of output file>\n");
    optionsText.append("\tSets model output file.\n");
    optionsText.append("-p <attribute range>\n");
    optionsText.append("\tOutput predictions. Predictions are for " 
		       + "training file" 
		       + "\n\tif only training file is specified," 
		       + "\n\totherwise predictions are for the test file."
		       + "\n\tThe range specifies attribute values to be output"
		       + "\n\twith the predictions. Use '-p 0' for none.\n");
    optionsText.append("-x <number of folds>\n");
    optionsText.append("\tOnly Distribution Clusterers can be cross validated.\n");
    optionsText.append("-s <random number seed>\n");
    optionsText.append("\tSets the seed for randomizing the data in cross-validation\n");
    optionsText.append("-c <class index>\n");
    optionsText.append("\tSet class attribute. If supplied, class is ignored");
    optionsText.append("\n\tduring clustering but is used in a classes to");
    optionsText.append("\n\tclusters evaluation.\n");
    if (clusterer instanceof Drawable) {
      optionsText.append("-g <name of graph file>\n");
      optionsText.append("\tOutputs the graph representation of the clusterer to the file.\n");
    }

    // Get scheme-specific options
    if (clusterer instanceof OptionHandler) {
      optionsText.append("\nOptions specific to " 
			 + clusterer.getClass().getName() + ":\n\n");
      Enumeration enu = ((OptionHandler)clusterer).listOptions();

      while (enu.hasMoreElements()) {
	Option option = (Option)enu.nextElement();
	optionsText.append(option.synopsis() + '\n');
	optionsText.append(option.description() + "\n");
      }
    }
    
    // Get global information (if available)
    if (globalInfo) {
      try {
        String gi = getGlobalInfo(clusterer);
        optionsText.append(gi);
      } catch (Exception ex) {
        // quietly ignore
      }
    }

    return  optionsText.toString();
  }
  
  /**
   * Return the global info (if it exists) for the supplied clusterer
   * 
   * @param clusterer the clusterer to get the global info for
   * @return the global info (synopsis) for the clusterer
   * @throws Exception if there is a problem reflecting on the clusterer
   */
  protected static String getGlobalInfo(Clusterer clusterer) throws Exception {
    BeanInfo bi = Introspector.getBeanInfo(clusterer.getClass());
    MethodDescriptor[] methods;
    methods = bi.getMethodDescriptors();
    Object[] args = {};
    String result = "\nSynopsis for " + clusterer.getClass().getName()
      + ":\n\n";
    
    for (int i = 0; i < methods.length; i++) {
      String name = methods[i].getDisplayName();
      Method meth = methods[i].getMethod();
      if (name.equals("globalInfo")) {
        String globalInfo = (String)(meth.invoke(clusterer, args));
        result += globalInfo;
        break;
      }
    }
    
    return result;
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
    
    ClusterEvaluation cmp = (ClusterEvaluation) obj;
    
    if ((m_classToCluster != null) != (cmp.m_classToCluster != null)) return false;
    if (m_classToCluster != null) {
      for (int i = 0; i < m_classToCluster.length; i++) {
        if (m_classToCluster[i] != cmp.m_classToCluster[i])
  	return false;
      }
    }
    
    if ((m_clusterAssignments != null) != (cmp.m_clusterAssignments != null)) return false;
    if (m_clusterAssignments != null) {
      for (int i = 0; i < m_clusterAssignments.length; i++) {
        if (m_clusterAssignments[i] != cmp.m_clusterAssignments[i])
  	return false;
      }
    }

    if (Double.isNaN(m_logL) != Double.isNaN(cmp.m_logL)) return false;
    if (!Double.isNaN(m_logL)) {
      if (m_logL != cmp.m_logL) return false;
    }
    
    if (m_numClusters != cmp.m_numClusters) return false;
    
    // TODO: better comparison? via members?
    String clusteringResults1 = m_clusteringResults.toString().replaceAll("Elapsed time.*", "");
    String clusteringResults2 = cmp.m_clusteringResults.toString().replaceAll("Elapsed time.*", "");
    if (!clusteringResults1.equals(clusteringResults2)) return false;
    
    return true;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6021 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param args the options
   */
  public static void main (String[] args) {
    try {
      if (args.length == 0) {
	throw  new Exception("The first argument must be the name of a " 
			     + "clusterer");
      }

      String ClustererString = args[0];
      args[0] = "";
      Clusterer newClusterer = AbstractClusterer.forName(ClustererString, null);
      System.out.println(evaluateClusterer(newClusterer, args));
    }
    catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
}

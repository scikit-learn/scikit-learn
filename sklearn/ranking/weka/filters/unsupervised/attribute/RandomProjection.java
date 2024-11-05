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
 *    RandomProjection.java
 *    Copyright (C) 2003 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.attribute;

import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance; 
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Reduces the dimensionality of the data by projecting it onto a lower dimensional subspace using a random matrix with columns of unit length (i.e. It will reduce the number of attributes in the data while preserving much of its variation like PCA, but at a much less computational cost).<br/>
 * It first applies the  NominalToBinary filter to convert all attributes to numeric before reducing the dimension. It preserves the class attribute.<br/>
 * <br/>
 * For more information, see:<br/>
 * <br/>
 * Dmitriy Fradkin, David Madigan: Experiments with random projections for machine learning. In: KDD '03: Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining, New York, NY, USA, 517-522, 003.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;inproceedings{Fradkin003,
 *    address = {New York, NY, USA},
 *    author = {Dmitriy Fradkin and David Madigan},
 *    booktitle = {KDD '03: Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining},
 *    pages = {517-522},
 *    publisher = {ACM Press},
 *    title = {Experiments with random projections for machine learning},
 *    year = {003}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -N &lt;number&gt;
 *  The number of dimensions (attributes) the data should be reduced to
 *  (default 10; exclusive of the class attribute, if it is set).</pre>
 * 
 * <pre> -D [SPARSE1|SPARSE2|GAUSSIAN]
 *  The distribution to use for calculating the random matrix.
 *  Sparse1 is:
 *    sqrt(3)*{-1 with prob(1/6), 0 with prob(2/3), +1 with prob(1/6)}
 *  Sparse2 is:
 *    {-1 with prob(1/2), +1 with prob(1/2)}
 * </pre>
 * 
 * <pre> -P &lt;percent&gt;
 *  The percentage of dimensions (attributes) the data should
 *  be reduced to (exclusive of the class attribute, if it is set). This -N
 *  option is ignored if this option is present or is greater
 *  than zero.</pre>
 * 
 * <pre> -M
 *  Replace missing values using the ReplaceMissingValues filter</pre>
 * 
 * <pre> -R &lt;num&gt;
 *  The random seed for the random number generator used for
 *  calculating the random matrix (default 42).</pre>
 * 
 <!-- options-end -->
 *
 * @author Ashraf M. Kibriya (amk14@cs.waikato.ac.nz) 
 * @version $Revision: 6749 $ [1.0 - 22 July 2003 - Initial version (Ashraf M. Kibriya)]
 */
public class RandomProjection 
  extends Filter 
  implements UnsupervisedFilter, OptionHandler, TechnicalInformationHandler {

  /** for serialization */
  static final long serialVersionUID = 4428905532728645880L;

  /** Stores the number of dimensions to reduce the data to */
  protected int m_k = 10;

  /** Stores the dimensionality the data should be reduced to as percentage of the original dimension */
  protected double m_percent = 0.0;

  /** Is the random matrix will be computed using 
      Gaussian distribution or not */
  protected boolean m_useGaussian = false;

  /** distribution type: sparse 1 */
  public static final int SPARSE1 = 1;
  /** distribution type: sparse 2 */
  public static final int SPARSE2 = 2;
  /** distribution type: gaussian */
  public static final int GAUSSIAN = 3;

  /** The types of distributions that can be used for 
  calculating the random matrix */
  public static final Tag [] TAGS_DSTRS_TYPE = {
    new Tag(SPARSE1, "Sparse1"),
    new Tag(SPARSE2, "Sparse2"),
    new Tag(GAUSSIAN, "Gaussian"),
  };

  /** Stores the distribution to use for calculating the
      random matrix */
  protected int m_distribution = SPARSE1;
 
  /** Should the missing values be replaced using 
      unsupervised.ReplaceMissingValues filter */
  protected boolean m_useReplaceMissing = false;

  /** Keeps track of output format if it is defined or not */
  protected boolean m_OutputFormatDefined = false;

  /** The NominalToBinary filter applied to the data before this filter */
  protected Filter m_ntob; // = new weka.filters.unsupervised.attribute.NominalToBinary();

  /** The ReplaceMissingValues filter */
  protected Filter m_replaceMissing;
    
  /** Stores the random seed used to generate the random matrix */
  protected long m_rndmSeed = 42;

  /** The random matrix */
  protected double m_rmatrix[][];

  /** The random number generator used for generating the random matrix */
  protected Random m_random;

  /**
   * Returns an enumeration describing the available options.
   *
   * @return an enumeration of all the available options.
   */
  public Enumeration listOptions() {

    Vector newVector = new Vector(2);

    newVector.addElement(new Option(
	      "\tThe number of dimensions (attributes) the data should be reduced to\n"
             +"\t(default 10; exclusive of the class attribute, if it is set).",
	      "N", 1, "-N <number>"));

    newVector.addElement(new Option(
	      "\tThe distribution to use for calculating the random matrix.\n"
	     +"\tSparse1 is:\n"
	     +"\t  sqrt(3)*{-1 with prob(1/6), 0 with prob(2/3), +1 with prob(1/6)}\n"
	     +"\tSparse2 is:\n"
	     +"\t  {-1 with prob(1/2), +1 with prob(1/2)}\n",
	      "D", 1, "-D [SPARSE1|SPARSE2|GAUSSIAN]"));

    //newVector.addElement(new Option(
    //	      "\tUse Gaussian distribution for calculating the random matrix.",
    //	      "G", 0, "-G"));

    newVector.addElement(new Option(
	      "\tThe percentage of dimensions (attributes) the data should\n"
	     +"\tbe reduced to (exclusive of the class attribute, if it is set). This -N\n"
	     +"\toption is ignored if this option is present or is greater\n"
	     +"\tthan zero.",
	      "P", 1, "-P <percent>"));

    newVector.addElement(new Option(
	      "\tReplace missing values using the ReplaceMissingValues filter",
	      "M", 0, "-M"));

    newVector.addElement(new Option(
	      "\tThe random seed for the random number generator used for\n"
	     +"\tcalculating the random matrix (default 42).",
	      "R", 0, "-R <num>"));
 
    return newVector.elements();
  }

  /**
   * Parses a given list of options. <p/>
   * 
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -N &lt;number&gt;
   *  The number of dimensions (attributes) the data should be reduced to
   *  (default 10; exclusive of the class attribute, if it is set).</pre>
   * 
   * <pre> -D [SPARSE1|SPARSE2|GAUSSIAN]
   *  The distribution to use for calculating the random matrix.
   *  Sparse1 is:
   *    sqrt(3)*{-1 with prob(1/6), 0 with prob(2/3), +1 with prob(1/6)}
   *  Sparse2 is:
   *    {-1 with prob(1/2), +1 with prob(1/2)}
   * </pre>
   * 
   * <pre> -P &lt;percent&gt;
   *  The percentage of dimensions (attributes) the data should
   *  be reduced to (exclusive of the class attribute, if it is set). This -N
   *  option is ignored if this option is present or is greater
   *  than zero.</pre>
   * 
   * <pre> -M
   *  Replace missing values using the ReplaceMissingValues filter</pre>
   * 
   * <pre> -R &lt;num&gt;
   *  The random seed for the random number generator used for
   *  calculating the random matrix (default 42).</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported
   */
  public void setOptions(String[] options) throws Exception {


    String mString = Utils.getOption('P', options);
    if (mString.length() != 0) {
	setPercent((double) Double.parseDouble(mString)); //setNumberOfAttributes((int) Integer.parseInt(mString));
    } else {
        setPercent(0);
	mString = Utils.getOption('N', options);
	if (mString.length() != 0) 
	    setNumberOfAttributes(Integer.parseInt(mString));	    
	else	    
	    setNumberOfAttributes(10);
    }    
    
    mString = Utils.getOption('R', options);
    if(mString.length()!=0) {
	setRandomSeed( Long.parseLong(mString) );
    }

    mString = Utils.getOption('D', options);
    if(mString.length()!=0) {
	if(mString.equalsIgnoreCase("sparse1"))
	   setDistribution( new SelectedTag(SPARSE1, TAGS_DSTRS_TYPE) );
	else if(mString.equalsIgnoreCase("sparse2"))
	   setDistribution( new SelectedTag(SPARSE2, TAGS_DSTRS_TYPE) );
	else if(mString.equalsIgnoreCase("gaussian"))
	   setDistribution( new SelectedTag(GAUSSIAN, TAGS_DSTRS_TYPE) );	   
    }

    if(Utils.getFlag('M', options))
	setReplaceMissingValues(true);
    else
	setReplaceMissingValues(false);


   //if(Utils.getFlag('G', options))
   //    setUseGaussian(true);
   //else
   //    setUseGaussian(false);
    
  }

  /**
   * Gets the current settings of the filter.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String [] getOptions() {

    String [] options = new String [10];
    int current = 0;

    //if (getUseGaussian()) {
    //  options[current++] = "-G";
    //}

    if (getReplaceMissingValues()) {
      options[current++] = "-M";
    }

    if (getPercent() == 0) {
      options[current++] = "-N";
      options[current++] = "" + getNumberOfAttributes();
    }
    else {
      options[current++] = "-P";
      options[current++] = "" + getPercent();
    }
    
    options[current++] = "-R";
    options[current++] = "" + getRandomSeed();
    
    SelectedTag t = getDistribution();
    options[current++] = "-D";
    options[current++] = ""+t.getSelectedTag().getReadable();


    while (current < options.length) {
      options[current++] = "";
    }

    return options;
  }
    
   
  /**
   * Returns a string describing this filter
   *
   * @return a description of the filter suitable for
   * displaying in the explorer/experimenter gui
   */
  public String globalInfo() {

    return "Reduces the dimensionality of the data by projecting"
	 + " it onto a lower dimensional subspace using a random"
	 + " matrix with columns of unit length (i.e. It will reduce"
	 + " the number of attributes in the data while preserving"
	 + " much of its variation like PCA, but at a much less"
	 + " computational cost).\n"
	 + "It first applies the  NominalToBinary filter to" 
	 + " convert all attributes to numeric before reducing the"
	 + " dimension. It preserves the class attribute.\n\n"
	 + "For more information, see:\n\n"
	 + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.INPROCEEDINGS);
    result.setValue(Field.AUTHOR, "Dmitriy Fradkin and David Madigan");
    result.setValue(Field.TITLE, "Experiments with random projections for machine learning");
    result.setValue(Field.BOOKTITLE, "KDD '03: Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining");
    result.setValue(Field.YEAR, "003");
    result.setValue(Field.PAGES, "517-522");
    result.setValue(Field.PUBLISHER, "ACM Press");
    result.setValue(Field.ADDRESS, "New York, NY, USA");
    
    return result;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String numberOfAttributesTipText() {

    return "The number of dimensions (attributes) the data should"
         + " be reduced to.";
  }

  /** 
   * Sets the number of attributes (dimensions) the data should be reduced to
   * 
   * @param newAttNum the goal for the dimensions
   */
  public void setNumberOfAttributes(int newAttNum) {
      m_k = newAttNum;
  }
  
  /** 
   * Gets the current number of attributes (dimensionality) to which the data 
   * will be reduced to.
   *  
   * @return the number of dimensions
   */
  public int getNumberOfAttributes() {
      return m_k;
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String percentTipText() {

      return  " The percentage of dimensions (attributes) the data should"
            + " be reduced to  (inclusive of the class attribute). This "
	    + " NumberOfAttributes option is ignored if this option is"
	    + " present or is greater than zero.";
  }

  /** 
   * Sets the percent the attributes (dimensions) of the data should be reduced to
   * 
   * @param newPercent the percentage of attributes
   */
  public void setPercent(double newPercent) {
      if(newPercent > 0)
	  newPercent /= 100;
      m_percent = newPercent;
  }

  /** 
   * Gets the percent the attributes (dimensions) of the data will be reduced to
   * 
   * @return the percentage of attributes
   */
  public double getPercent() {
      return m_percent * 100;
  }


  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String randomSeedTipText() {
      return  "The random seed used by the random"
	     +" number generator used for generating"
	     +" the random matrix ";
  }

  /** 
   * Sets the random seed of the random number generator
   * 
   * @param seed the random seed value
   */
  public void setRandomSeed(long seed) {
      m_rndmSeed = seed;
  }

  /** 
   * Gets the random seed of the random number generator
   * 
   * @return the random seed value
   */
  public long getRandomSeed() {
      return m_rndmSeed;
  }


  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String  distributionTipText() {
      return "The distribution to use for calculating the random matrix.\n"
	    +"Sparse1 is:\n"
	    +" sqrt(3) * { -1 with prob(1/6), \n"
	    +"               0 with prob(2/3),  \n"
            +"              +1 with prob(1/6) } \n"
	    +"Sparse2 is:\n"
	    +" { -1 with prob(1/2), \n"
	    +"   +1 with prob(1/2) } ";
      
  }
  /** 
   * Sets the distribution to use for calculating the random matrix
   * 
   * @param newDstr the distribution to use
   */
  public void setDistribution(SelectedTag newDstr) {

      if (newDstr.getTags() == TAGS_DSTRS_TYPE) {
	  m_distribution = newDstr.getSelectedTag().getID();
      }
  }

  /** 
   * Returns the current distribution that'll be used for calculating the 
   * random matrix
   * 
   * @return the current distribution
   */
  public SelectedTag getDistribution() {
      return new SelectedTag(m_distribution, TAGS_DSTRS_TYPE);
  }

  /**
   * Returns the tip text for this property
   *
   * @return tip text for this property suitable for
   * displaying in the explorer/experimenter gui
   */
  public String replaceMissingValuesTipText() {

    return "If set the filter uses weka.filters.unsupervised.attribute.ReplaceMissingValues"
	 + " to replace the missing values";
  }

  /** 
   * Sets either to use replace missing values filter or not
   * 
   * @param t if true then the replace missing values is used
   */
  public void setReplaceMissingValues(boolean t) {
      m_useReplaceMissing = t;
  }

  /** 
   * Gets the current setting for using ReplaceMissingValues filter
   * 
   * @return true if the replace missing values filter is used
   */
  public boolean getReplaceMissingValues() {
      return m_useReplaceMissing;
  }

  /** 
   * Returns the Capabilities of this filter.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();

    // attributes
    result.enableAllAttributes();
    result.enable(Capability.MISSING_VALUES);
    
    // class
    result.enableAllClasses();
    result.enable(Capability.MISSING_CLASS_VALUES);
    result.enable(Capability.NO_CLASS);

    //RANKING BEGIN
    result.disable(Capability.RANKING);
    result.disable(Capability.PREFERENCE_ATTRIBUTE);
    //RANKING END
    
    return result;
  }

  /**
   * Sets the format of the input instances.
   *
   * @param instanceInfo an Instances object containing the input 
   * instance structure (any instances contained in the object are 
   * ignored - only the structure is required).
   * @return true if the outputFormat may be collected immediately
   * @throws Exception if the input format can't be set 
   * successfully
   */
  public boolean setInputFormat(Instances instanceInfo) throws Exception {      
    super.setInputFormat(instanceInfo);
    /*
    if (instanceInfo.classIndex() < 0) {
      throw new UnassignedClassException("No class has been assigned to the instances");
    }
    */
    
    for(int i=0; i<instanceInfo.numAttributes(); i++) {        
	if( i!=instanceInfo.classIndex() && (instanceInfo.attribute(i).isNominal() || instanceInfo.attribute(i).isRanking() ) ) {
            if(instanceInfo.classIndex()>=0)
                m_ntob = new weka.filters.supervised.attribute.NominalToBinary();
            else
                m_ntob = new weka.filters.unsupervised.attribute.NominalToBinary();
            
            break;
	}
    }

    //r.setSeed(m_rndmSeed); //in case the setRandomSeed() is not
                           //called we better set the seed to its 
                           //default value of 42.
    boolean temp=true;
    if(m_replaceMissing!=null) {
	m_replaceMissing = new weka.filters.unsupervised.attribute.ReplaceMissingValues();
	if(m_replaceMissing.setInputFormat(instanceInfo))
	    temp=true;
	else
	    temp=false;
    }
    
    if(m_ntob!=null) {
	if(m_ntob.setInputFormat(instanceInfo)) {
	    setOutputFormat();
	    return temp && true;
	}
	else { 
	    return false;
	}
    }
    else {
	setOutputFormat();
	return temp && true;
    }
  }

   
  /**
   * Input an instance for filtering.
   *
   * @param instance the input instance
   * @return true if the filtered instance may now be
   * collected with output().
   * @throws IllegalStateException if no input format has been set
   */
  public boolean input(Instance instance) throws Exception {

    Instance newInstance=null;

    if (getInputFormat()==null) {
	throw new IllegalStateException("No input instance format defined");
    }
    if(m_NewBatch) {
      resetQueue();
      //if(ntob!=null) 
      //	  ntob.m_NewBatch=true;
      m_NewBatch = false;
    }
    
    boolean replaceDone=false;
    if(m_replaceMissing!=null) {
	if(m_replaceMissing.input(instance)) {
	    if(m_OutputFormatDefined == false)
		setOutputFormat();
	    newInstance = m_replaceMissing.output();
	    replaceDone = true;
	}
	else
	    return false;;
    }

    if(m_ntob!=null) {
	if(replaceDone==false)
	    newInstance = instance;
	if(m_ntob.input(newInstance)) {
	    if(m_OutputFormatDefined == false) 
		setOutputFormat();
	    newInstance = m_ntob.output();
	    newInstance = convertInstance(newInstance);
	    push(newInstance);
	    return true;	
	}
	else {
	    return false;
	}
    }
    else {
	if(replaceDone==false)
	    newInstance = instance;
	newInstance = convertInstance(newInstance);
	push(newInstance);
	return true;
    }
  }


  /**
   * Signify that this batch of input to the filter is finished.
   *
   * @return true if there are instances pending output
   * @throws NullPointerException if no input structure has been defined,
   * @throws Exception if there was a problem finishing the batch.
   */
  public boolean batchFinished() throws Exception {
      if (getInputFormat() == null) {
	  throw new NullPointerException("No input instance format defined");
      }
      
      boolean conversionDone=false;
      if(m_replaceMissing!=null) {
	  if(m_replaceMissing.batchFinished()) {
	      Instance newInstance, instance;
	      
	      while((instance=m_replaceMissing.output())!=null) {
		  if(!m_OutputFormatDefined)
		      setOutputFormat();
		  if(m_ntob!=null) {
		      m_ntob.input(instance);
		  }
		  else {
		      newInstance = convertInstance(instance);
		      push(newInstance);
		  }
	      }

	      if(m_ntob!=null) {
		  if(m_ntob.batchFinished()) {
		      //Instance newInstance, instance;
		      while((instance=m_ntob.output())!=null) {
			  if(!m_OutputFormatDefined)
			      setOutputFormat();
			  newInstance = convertInstance(instance);
			  push(newInstance);
		      }
		      m_ntob = null;		      
		  }
	      }
	      m_replaceMissing = null;
	      conversionDone=true;
	  }
      }

      if(conversionDone==false && m_ntob!=null) {
	  if(m_ntob.batchFinished()) {
	      Instance newInstance, instance;
	      while((instance=m_ntob.output())!=null) {
		  if(!m_OutputFormatDefined)
		      setOutputFormat();
		  newInstance = convertInstance(instance);
		  push(newInstance);
	      }
	      m_ntob = null;
	  }
      }
      m_OutputFormatDefined=false;
      return super.batchFinished();
  }
    

  /** Sets the output format */  
  protected void setOutputFormat() {
      Instances currentFormat;
      if(m_ntob!=null) {
	  currentFormat = m_ntob.getOutputFormat();
      }
      else 
	  currentFormat = getInputFormat();
      
      if(m_percent>0)
	  { m_k = (int) ((getInputFormat().numAttributes()-1)*m_percent); 
	  // System.out.print("numAtts: "+currentFormat.numAttributes());
	  // System.out.print("percent: "+m_percent);
	  // System.out.print("percent*numAtts: "+(currentFormat.numAttributes()*m_percent));
	  // System.out.println("m_k: "+m_k);
	  }

      Instances newFormat;
      int newClassIndex=-1;
      FastVector attributes = new FastVector();
      for(int i=0; i<m_k; i++) {
	  attributes.addElement( new Attribute("K"+(i+1)) );
      }
      if(currentFormat.classIndex()!=-1)  {  //if classindex is set
	  //attributes.removeElementAt(attributes.size()-1);
	  attributes.addElement(currentFormat.attribute(currentFormat.classIndex()));
	  newClassIndex = attributes.size()-1;
      }

      newFormat = new Instances(currentFormat.relationName(), attributes, 0);
      if(newClassIndex!=-1)
	  newFormat.setClassIndex(newClassIndex);
      m_OutputFormatDefined=true;

      m_random = new Random();
      m_random.setSeed(m_rndmSeed);

      m_rmatrix = new double[m_k][currentFormat.numAttributes()];
      if(m_distribution==GAUSSIAN) {
	  for(int i=0; i<m_rmatrix.length; i++) 
	      for(int j=0; j<m_rmatrix[i].length; j++) 
		  m_rmatrix[i][j] = m_random.nextGaussian();
      }
      else {
	  boolean useDstrWithZero = (m_distribution==SPARSE1);
	  for(int i=0; i<m_rmatrix.length; i++) 
	      for(int j=0; j<m_rmatrix[i].length; j++) 
		  m_rmatrix[i][j] = rndmNum(useDstrWithZero);
      }

      setOutputFormat(newFormat);
  }

  /**
   * converts a single instance to the required format
   *
   * @param currentInstance     the instance to convert
   * @return                    the converted instance
   */
  protected Instance convertInstance(Instance currentInstance) {

      Instance newInstance;
      double vals[] = new double[getOutputFormat().numAttributes()];
      int classIndex = (m_ntob==null) ? getInputFormat().classIndex():m_ntob.getOutputFormat().classIndex();

      for(int i = 0; i < m_k; i++) {
        vals[i] = computeRandomProjection(i,classIndex,currentInstance);
      }
      if (classIndex != -1) {
        vals[m_k] = currentInstance.value(classIndex);
      }

      newInstance = new DenseInstance(currentInstance.weight(), vals);
      newInstance.setDataset(getOutputFormat());

      return newInstance;
  }


  /**
   * computes one random projection for a given instance (skip missing values)
   *
   * @param rpIndex     offset the new random projection attribute
   * @param classIndex  classIndex of the input instance
   * @param instance    the instance to convert
   * @return    the random sum
   */

  protected double computeRandomProjection(int rpIndex, int classIndex, Instance instance) {

    double sum = 0.0;
    for(int i = 0; i < instance.numValues(); i++) {
      int index = instance.index(i);
      if (index != classIndex) {
        double value = instance.valueSparse(i);
        if (!Utils.isMissingValue(value)) {
          sum += m_rmatrix[rpIndex][index] * value;
        }
      }
    }
    return sum;
  }

  private static final int weights[] = {1, 1, 4};
  private static final int vals[] = {-1, 1, 0};
  private static final int weights2[] = {1, 1};
  private static final int vals2[] = {-1, 1};
  private static final double sqrt3 = Math.sqrt(3);

  /**
   * returns a double x such that <br/>
   *      x = sqrt(3) * { -1 with prob. 1/6, 0 with prob. 2/3, 1 with prob. 1/6 }
   *      
   * @param useDstrWithZero
   * @return the generated number
   */
  protected double rndmNum(boolean useDstrWithZero) {
      if(useDstrWithZero)
	  return sqrt3 * vals[weightedDistribution(weights)];
      else
	  return vals2[weightedDistribution(weights2)];
  }

  /** 
   * Calculates a weighted distribution
   * 
   * @param weights the weights to use
   * @return
   */
  protected int weightedDistribution(int [] weights) {
      int sum=0; 
      
      for(int i=0; i<weights.length; i++) 
	  sum += weights[i];
      
      int val = (int)Math.floor(m_random.nextDouble()*sum);
      
      for(int i=0; i<weights.length; i++) {
	  val -= weights[i];
	  if(val<0)
	      return i;
      }
      return -1;
  }  
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6749 $");
  }

  /**
   * Main method for testing this class.
   *
   * @param argv should contain arguments to the filter: 
   * use -h for help
   */
  public static void main(String [] argv) {
    runFilter(new RandomProjection(), argv);
  }
}

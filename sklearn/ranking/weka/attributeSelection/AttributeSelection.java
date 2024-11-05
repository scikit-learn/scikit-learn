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
 *    AttributeSelection.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.MethodDescriptor;
import java.beans.PropertyDescriptor;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.util.Enumeration;
import java.util.Random;

/** 
 * Attribute selection class. Takes the name of a search class and
 * an evaluation class on the command line. <p/>
 *
 * Valid options are: <p/>
 *
 * -h <br/>
 * Display help. <p/>
 *
 * -i &lt;name of input file&gt; <br/>
 * Specify the training data file. <p/>
 * 
 * -c &lt;class index&gt; <br/>
 * The index of the attribute to use as the class. <p/>
 * 
 * -s &lt;search method&gt; <br/>
 * The full class name of the search method followed by search method options
 * (if any).<br/>
 * Eg. -s "weka.attributeSelection.BestFirst -N 10" <p/>
 *
 * -x &lt;number of folds&gt; <br/>
 * Perform a cross validation. <p/>
 *
 * -n &lt;random number seed&gt; <br/>
 * Specify a random number seed. Use in conjuction with -x. (Default = 1). <p/>
 * 
 * ------------------------------------------------------------------------ <p/>
 * 
 * Example usage as the main of an attribute evaluator (called FunkyEvaluator):
 * <pre>
 * public static void main(String [] args) {
 *   runEvaluator(new FunkyEvaluator(), args);
 * }
 * </pre>
 * <p/>
 *
 * ------------------------------------------------------------------------ <p/>
 *
 * @author   Mark Hall (mhall@cs.waikato.ac.nz)
 * @version  $Revision: 1.47 $
 */
public class AttributeSelection 
  implements Serializable, RevisionHandler {
  
  /** for serialization */
  static final long serialVersionUID = 4170171824147584330L;

  /** the instances to select attributes from */
  private Instances m_trainInstances;

  /** the attribute/subset evaluator */
  private ASEvaluation m_ASEvaluator;

  /** the search method */
  private ASSearch m_searchMethod;

  /** the number of folds to use for cross validation */
  private int m_numFolds;

  /** holds a string describing the results of the attribute selection */
  private StringBuffer m_selectionResults;

  /** rank features (if allowed by the search method) */
  private boolean m_doRank;
  
  /** do cross validation */
  private boolean m_doXval;

  /** seed used to randomly shuffle instances for cross validation */
  private int m_seed;

  /** number of attributes requested from ranked results */
  private int m_numToSelect;
  
  /** the selected attributes */
  private int [] m_selectedAttributeSet;

  /** the attribute indexes and associated merits if a ranking is produced */
  private double [][] m_attributeRanking;

  /** if a feature selection run involves an attribute transformer */
  private AttributeTransformer m_transformer = null;
  
  /** the attribute filter for processing instances with respect to
      the most recent feature selection run */
  private Remove m_attributeFilter = null;

  /** hold statistics for repeated feature selection, such as
      under cross validation */
  private double [][] m_rankResults = null;
  private double [] m_subsetResults = null;
  private int m_trials = 0;

  /**
   * Return the number of attributes selected from the most recent
   * run of attribute selection
   * @return the number of attributes selected
   */
  public int numberAttributesSelected() throws Exception {
    int [] att = selectedAttributes();
    return att.length-1;
  }

  /**
   * get the final selected set of attributes.
   * @return an array of attribute indexes
   * @exception Exception if attribute selection has not been performed yet
   */
  public int [] selectedAttributes () throws Exception {
    if (m_selectedAttributeSet == null) {
      throw new Exception("Attribute selection has not been performed yet!");
    }
    return m_selectedAttributeSet;
  }

  /**
   * get the final ranking of the attributes.
   * @return a two dimensional array of ranked attribute indexes and their
   * associated merit scores as doubles.
   * @exception Exception if a ranking has not been produced
   */
  public double [][] rankedAttributes () throws Exception {
    if (m_attributeRanking == null) {
      throw new Exception("Ranking has not been performed");
    }
    return m_attributeRanking;
  }

  /**
   * set the attribute/subset evaluator
   * @param evaluator the evaluator to use
   */
  public void setEvaluator (ASEvaluation evaluator) {
    m_ASEvaluator = evaluator;
  }

  /**
   * set the search method
   * @param search the search method to use
   */
  public void setSearch (ASSearch search) {
    m_searchMethod = search;

    if (m_searchMethod instanceof RankedOutputSearch) {
      setRanking(((RankedOutputSearch)m_searchMethod).getGenerateRanking());
    }
  }

  /**
   * set the number of folds for cross validation
   * @param folds the number of folds
   */
  public void setFolds (int folds) {
    m_numFolds = folds;
  }

  /**
   * produce a ranking (if possible with the set search and evaluator)
   * @param r true if a ranking is to be produced
   */
  public void setRanking (boolean r) {
    m_doRank = r;
  }

  /**
   * do a cross validation
   * @param x true if a cross validation is to be performed
   */
  public void setXval (boolean x) {
    m_doXval = x;
  }

  /**
   * set the seed for use in cross validation
   * @param s the seed
   */
  public void setSeed (int s) {
    m_seed = s;
  }

  /**
   * get a description of the attribute selection
   * @return a String describing the results of attribute selection
   */
  public String toResultsString() {
    return m_selectionResults.toString();
  }

  /**
   * reduce the dimensionality of a set of instances to include only those 
   * attributes chosen by the last run of attribute selection.
   * @param in the instances to be reduced
   * @return a dimensionality reduced set of instances
   * @exception Exception if the instances can't be reduced
   */
  public Instances reduceDimensionality(Instances in) throws Exception {
    if (m_attributeFilter == null) {
      throw new Exception("No feature selection has been performed yet!");
    }

    if (m_transformer != null) {
      Instances transformed = new Instances(m_transformer.transformedHeader(),
					    in.numInstances());
      for (int i=0;i<in.numInstances();i++) {
	transformed.add(m_transformer.convertInstance(in.instance(i)));
      }
      return Filter.useFilter(transformed, m_attributeFilter);
    }

    return Filter.useFilter(in, m_attributeFilter);
  }

  /**
   * reduce the dimensionality of a single instance to include only those 
   * attributes chosen by the last run of attribute selection.
   * @param in the instance to be reduced
   * @return a dimensionality reduced instance
   * @exception Exception if the instance can't be reduced
   */
  public Instance reduceDimensionality(Instance in) throws Exception {
    if (m_attributeFilter == null) {
      throw new Exception("No feature selection has been performed yet!");
    }
    if (m_transformer != null) {
      in = m_transformer.convertInstance(in);
    }
    m_attributeFilter.input(in);
    m_attributeFilter.batchFinished();
    Instance result = m_attributeFilter.output();
    return result;
  }

  /**
   * constructor. Sets defaults for each member varaible. Default
   * attribute evaluator is CfsSubsetEval; default search method is
   * BestFirst.
   */
  public AttributeSelection () {
    setFolds(10);
    setRanking(false);
    setXval(false);
    setSeed(1);
    setEvaluator(new CfsSubsetEval());
    setSearch(new GreedyStepwise());
    m_selectionResults = new StringBuffer();
    m_selectedAttributeSet = null;
    m_attributeRanking = null;
  }

  /**
   * Perform attribute selection with a particular evaluator and
   * a set of options specifying search method and input file etc.
   *
   * @param ASEvaluator an evaluator object
   * @param options an array of options, not only for the evaluator
   * but also the search method (if any) and an input data file
   * @return the results of attribute selection as a String
   * @exception Exception if no training file is set
   */
  public static String SelectAttributes (ASEvaluation ASEvaluator, 
					 String[] options)
    throws Exception {
    String trainFileName, searchName;
    Instances train = null;
    ASSearch searchMethod = null;
    String[] optionsTmp = (String[]) options.clone();
    boolean helpRequested = false;

    try {
      // get basic options (options the same for all attribute selectors
      trainFileName = Utils.getOption('i', options);
      helpRequested = Utils.getFlag('h', optionsTmp);

      if (helpRequested || (trainFileName.length() == 0)) {
        searchName = Utils.getOption('s', optionsTmp);
        if (searchName.length() != 0) {
          String[] searchOptions = Utils.splitOptions(searchName);
          searchMethod = (ASSearch)Class.forName(searchOptions[0]).newInstance();
        }

        if (helpRequested)
          throw new Exception("Help requested.");
        else
          throw new Exception("No training file given.");
      }
    }
    catch (Exception e) {
      throw  new Exception('\n' + e.getMessage() 
			   + makeOptionString(ASEvaluator, searchMethod));
    }

    DataSource source = new DataSource(trainFileName);
    train = source.getDataSet();
    return SelectAttributes(ASEvaluator, options, train);
  }

  /**
   * returns a string summarizing the results of repeated attribute
   * selection runs on splits of a dataset.
   * @return a summary of attribute selection results
   * @exception Exception if no attribute selection has been performed.
   */
  public String CVResultsString () throws Exception {
    StringBuffer CvString = new StringBuffer();
    
    if ((m_subsetResults == null && m_rankResults == null) ||
	( m_trainInstances == null)) {
      throw new Exception("Attribute selection has not been performed yet!");
    }

    int fieldWidth = (int)(Math.log(m_trainInstances.numAttributes()) +1.0);

    CvString.append("\n\n=== Attribute selection " + m_numFolds 
		    + " fold cross-validation ");

    if (!(m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) && 
	!(m_ASEvaluator instanceof UnsupervisedAttributeEvaluator) &&
	//RANKING BEGIN
	(m_trainInstances.classAttribute().isNominal()||m_trainInstances.classAttribute().isRanking())) {
    	//RANKING END
	CvString.append("(stratified), seed: ");
	CvString.append(m_seed+" ===\n\n");
    }
    else {
      CvString.append("seed: "+m_seed+" ===\n\n");
    }

    if ((m_searchMethod instanceof RankedOutputSearch) && (m_doRank == true)) {
      CvString.append("average merit      average rank  attribute\n");

      // calcualte means and std devs
      for (int i = 0; i < m_rankResults[0].length; i++) {
	m_rankResults[0][i] /= m_numFolds; // mean merit
	double var = m_rankResults[0][i]*m_rankResults[0][i]*m_numFolds;
	var = (m_rankResults[2][i] - var);
	var /= m_numFolds;

	if (var <= 0.0) {
	  var = 0.0;
	  m_rankResults[2][i] = 0;
	}
	else {
	  m_rankResults[2][i] = Math.sqrt(var);
	}

	m_rankResults[1][i] /= m_numFolds; // mean rank
	var = m_rankResults[1][i]*m_rankResults[1][i]*m_numFolds;
	var = (m_rankResults[3][i] - var);
	var /= m_numFolds;

	if (var <= 0.0) {
	  var = 0.0;
	  m_rankResults[3][i] = 0;
	}
	else {
	  m_rankResults[3][i] = Math.sqrt(var);
	}
      }

      // now sort them by mean rank
      int[] s = Utils.sort(m_rankResults[1]);
      for (int i=0; i<s.length; i++) {
	if (m_rankResults[1][s[i]] > 0) {
	  CvString.append(Utils.doubleToString(Math.
					       abs(m_rankResults[0][s[i]]),
					       6, 3) 
			  + " +-" 
			  + Utils.doubleToString(m_rankResults[2][s[i]], 6, 3) 
			  + "   " 
			  + Utils.doubleToString(m_rankResults[1][s[i]],
						 fieldWidth+2, 1) 
			  + " +-" 
			  + Utils.doubleToString(m_rankResults[3][s[i]], 5, 2) 
			+"  "
			  + Utils.doubleToString(((double)(s[i] + 1)), 
						 fieldWidth, 0)
			  + " " 
			  + m_trainInstances.attribute(s[i]).name() 
			  + "\n");
	}
      }
    }
    else {
      CvString.append("number of folds (%)  attribute\n");

      for (int i = 0; i < m_subsetResults.length; i++) {
	if ((m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) ||
	    (i != m_trainInstances.classIndex())) {
	  CvString.append(Utils.doubleToString(m_subsetResults[i], 12, 0) 
			  + "(" 
			  + Utils.doubleToString((m_subsetResults[i] / 
						  m_numFolds * 100.0)
						 , 3, 0) 
			  + " %)  " 
			  + Utils.doubleToString(((double)(i + 1)),
						 fieldWidth, 0)
			  + " " 
			  + m_trainInstances.attribute(i).name() 
			  + "\n");
	}
      }
    }

    return CvString.toString();
  }

  /**
   * Select attributes for a split of the data. Calling this function
   * updates the statistics on attribute selection. CVResultsString()
   * returns a string summarizing the results of repeated calls to
   * this function. Assumes that splits are from the same dataset---
   * ie. have the same number and types of attributes as previous
   * splits.
   *
   * @param split the instances to select attributes from
   * @exception Exception if an error occurs
   */
  public void selectAttributesCVSplit(Instances split) throws Exception {
    double[][] attributeRanking = null;

    // if the train instances are null then set equal to this split.
    // If this is the case then this function is more than likely being
    // called from outside this class in order to obtain CV statistics
    // and all we need m_trainIstances for is to get at attribute names
    // and types etc.
    if (m_trainInstances == null) {
      m_trainInstances = split;
    }

    // create space to hold statistics
    if (m_rankResults == null && m_subsetResults == null) {
      m_subsetResults = new double[split.numAttributes()];
      m_rankResults = new double[4][split.numAttributes()];
    }

    m_ASEvaluator.buildEvaluator(split);
    // Do the search
    int[] attributeSet = m_searchMethod.search(m_ASEvaluator, 
					       split);
    // Do any postprocessing that a attribute selection method might 
    // require
    attributeSet = m_ASEvaluator.postProcess(attributeSet);
    
    if ((m_searchMethod instanceof RankedOutputSearch) && 
	(m_doRank == true)) {
      attributeRanking = ((RankedOutputSearch)m_searchMethod).
	rankedAttributes();
      
      // System.out.println(attributeRanking[0][1]);
      for (int j = 0; j < attributeRanking.length; j++) {
	// merit
	m_rankResults[0][(int)attributeRanking[j][0]] += 
	  attributeRanking[j][1];
	// squared merit
	m_rankResults[2][(int)attributeRanking[j][0]] += 
	  (attributeRanking[j][1]*attributeRanking[j][1]);
	// rank
	m_rankResults[1][(int)attributeRanking[j][0]] += (j + 1);
	// squared rank
	m_rankResults[3][(int)attributeRanking[j][0]] += (j + 1)*(j + 1);
	// += (attributeRanking[j][0] * attributeRanking[j][0]);
      }
    } else {
      for (int j = 0; j < attributeSet.length; j++) {
	m_subsetResults[attributeSet[j]]++;
      }
    }

    m_trials++;
  }

  /**
   * Perform a cross validation for attribute selection. With subset
   * evaluators the number of times each attribute is selected over
   * the cross validation is reported. For attribute evaluators, the
   * average merit and average ranking + std deviation is reported for
   * each attribute.
   *
   * @return the results of cross validation as a String
   * @exception Exception if an error occurs during cross validation
   */
  public String CrossValidateAttributes () throws Exception {
    Instances cvData = new Instances(m_trainInstances);
    Instances train;

    Random random = new Random(m_seed);
    cvData.randomize(random);

    if (!(m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) && 
	!(m_ASEvaluator instanceof UnsupervisedAttributeEvaluator)) {
    	//RANKING BEGIN
      if (cvData.classAttribute().isNominal() || cvData.classAttribute().isRanking()) {
    	//RANKING END
	cvData.stratify(m_numFolds);
      }

    }

    for (int i = 0; i < m_numFolds; i++) {
      // Perform attribute selection
      train = cvData.trainCV(m_numFolds, i, random);
      selectAttributesCVSplit(train);
    }

    return  CVResultsString();
  }

  /**
   * Perform attribute selection on the supplied training instances.
   *
   * @param data the instances to select attributes from
   * @exception Exception if there is a problem during selection
   */
  public void SelectAttributes (Instances data) throws Exception {
    int [] attributeSet;
    
    m_transformer = null;
    m_attributeFilter = null;
    m_trainInstances = data;
    
    if (m_doXval == true && (m_ASEvaluator instanceof AttributeTransformer)) {
      throw new Exception("Can't cross validate an attribute transformer.");
    }

    if (m_ASEvaluator instanceof SubsetEvaluator &&
	m_searchMethod instanceof Ranker) {
      throw new Exception(m_ASEvaluator.getClass().getName()
			  +" must use a search method other than Ranker");
    }

    if (m_ASEvaluator instanceof AttributeEvaluator &&
	!(m_searchMethod instanceof Ranker)) {
      //      System.err.println("AttributeEvaluators must use a Ranker search "
      //			 +"method. Switching to Ranker...");
      //      m_searchMethod = new Ranker();
      throw new Exception("AttributeEvaluators must use the Ranker search "
			  + "method");
    }

    if (m_searchMethod instanceof RankedOutputSearch) {
      m_doRank = ((RankedOutputSearch)m_searchMethod).getGenerateRanking();
    }

    if (m_ASEvaluator instanceof UnsupervisedAttributeEvaluator ||
	m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) {
      // unset the class index
      //      m_trainInstances.setClassIndex(-1);
    } else {
      // check that a class index has been set
      if (m_trainInstances.classIndex() < 0) {
	m_trainInstances.setClassIndex(m_trainInstances.numAttributes()-1);
      }
    }
    
    // Initialize the attribute evaluator
    m_ASEvaluator.buildEvaluator(m_trainInstances);
    if (m_ASEvaluator instanceof AttributeTransformer) {
      m_trainInstances = 
	((AttributeTransformer)m_ASEvaluator).transformedHeader();
      m_transformer = (AttributeTransformer)m_ASEvaluator;
    }
    int fieldWidth = (int)(Math.log(m_trainInstances.numAttributes()) +1.0);

    // Do the search
    attributeSet = m_searchMethod.search(m_ASEvaluator, 
					 m_trainInstances);

    // try and determine if the search method uses an attribute transformer---
    // this is a bit of a hack to make things work properly with RankSearch
    // using PrincipalComponents as its attribute ranker
     try {
       BeanInfo bi = Introspector.getBeanInfo(m_searchMethod.getClass());
       PropertyDescriptor properties[];
       MethodDescriptor methods[];
       //       methods = bi.getMethodDescriptors();
       properties = bi.getPropertyDescriptors();
       for (int i=0;i<properties.length;i++) {
	 String name = properties[i].getDisplayName();
	 Method meth = properties[i].getReadMethod();
	 Object retType = meth.getReturnType();
	 if (retType.equals(ASEvaluation.class)) {
	   Class args [] = { };
	   ASEvaluation tempEval = (ASEvaluation)(meth.invoke(m_searchMethod,
							      (Object[])args));
	   if (tempEval instanceof AttributeTransformer) {
	     // grab the transformed data header
	     m_trainInstances = 
	       ((AttributeTransformer)tempEval).transformedHeader();
	     m_transformer = (AttributeTransformer)tempEval;
	   }
	 }
       }
     } catch (IntrospectionException ex) {
       System.err.println("AttributeSelection: Couldn't "
			  +"introspect");
     }
     
     
     // Do any postprocessing that a attribute selection method might require
     attributeSet = m_ASEvaluator.postProcess(attributeSet);
     if (!m_doRank) {
       m_selectionResults.append(printSelectionResults());
     }
     
    if ((m_searchMethod instanceof RankedOutputSearch) && m_doRank == true) {
      m_attributeRanking = 
	((RankedOutputSearch)m_searchMethod).rankedAttributes();
      m_selectionResults.append(printSelectionResults());
      m_selectionResults.append("Ranked attributes:\n");

      // retrieve the number of attributes to retain
      m_numToSelect = 
	((RankedOutputSearch)m_searchMethod).getCalculatedNumToSelect();

      // determine fieldwidth for merit
      int f_p=0;
      int w_p=0;
     
      for (int i = 0; i < m_numToSelect; i++) {
	double precision = (Math.abs(m_attributeRanking[i][1]) - 
			    (int)(Math.abs(m_attributeRanking[i][1])));
        double intPart = (int)(Math.abs(m_attributeRanking[i][1]));

	if (precision > 0) {
	  precision = Math.abs((Math.log(Math.abs(precision)) / 
				Math.log(10)))+3;
	}
	if (precision > f_p) {
	  f_p = (int)precision;
	}

        if (intPart == 0) {
          if (w_p < 2) {
            w_p = 2;
          }
        } else if ((Math.abs((Math.log(Math.abs(m_attributeRanking[i][1])) 
		       / Math.log(10)))+1) > w_p) {
	  if (m_attributeRanking[i][1] > 0) {
	    w_p = (int)Math.abs((Math.log(Math.abs(m_attributeRanking[i][1]))
				 / Math.log(10)))+1;
	  }
	}
      }

      for (int i = 0; i < m_numToSelect; i++) {
	m_selectionResults.
	  append(Utils.doubleToString(m_attributeRanking[i][1],
				      f_p+w_p+1,f_p) 
		 + Utils.doubleToString((m_attributeRanking[i][0] + 1),
					fieldWidth+1,0) 
		 + " " 
		 + m_trainInstances.
		 attribute((int)m_attributeRanking[i][0]).name() 
		 + "\n");
      }
       
      // set up the selected attributes array - usable by a filter or
      // whatever
      if (m_trainInstances.classIndex() >= 0) {
        if ((!(m_ASEvaluator instanceof UnsupervisedSubsetEvaluator)        
             && !(m_ASEvaluator instanceof UnsupervisedAttributeEvaluator)) ||
            m_ASEvaluator instanceof AttributeTransformer) {
          // one more for the class
          m_selectedAttributeSet = new int[m_numToSelect + 1];
          m_selectedAttributeSet[m_numToSelect] = 
              m_trainInstances.classIndex();
        } else {
          m_selectedAttributeSet = new int[m_numToSelect];
        }
      } else {
	m_selectedAttributeSet = new int[m_numToSelect];
      }

      m_selectionResults.append("\nSelected attributes: ");

      for (int i = 0; i < m_numToSelect; i++) {
	m_selectedAttributeSet[i] = (int)m_attributeRanking[i][0];
	 
	if (i == m_numToSelect - 1) {
	  m_selectionResults.append(((int)m_attributeRanking[i][0] + 1) 
				    + " : " 
				    + (i + 1) 
				    + "\n");
	}
	else {
	  m_selectionResults.append(((int)m_attributeRanking[i][0] + 1));
	  m_selectionResults.append(",");
	}
      }
    } else {
      // set up the selected attributes array - usable by a filter or
      // whatever
      if ((!(m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) 
	   && !(m_ASEvaluator instanceof UnsupervisedAttributeEvaluator)) || 
	  m_trainInstances.classIndex() >= 0) 
	// one more for the class
	{
	  m_selectedAttributeSet = new int[attributeSet.length + 1];
	  m_selectedAttributeSet[attributeSet.length] = 
	    m_trainInstances.classIndex();
	}
      else {
	m_selectedAttributeSet = new int[attributeSet.length];
      }
      
      for (int i = 0; i < attributeSet.length; i++) {
	m_selectedAttributeSet[i] = attributeSet[i];
      }

      m_selectionResults.append("Selected attributes: ");

      for (int i = 0; i < attributeSet.length; i++) {
	if (i == (attributeSet.length - 1)) {
	  m_selectionResults.append((attributeSet[i] + 1) 
				    + " : " 
				    + attributeSet.length 
				    + "\n");
	}
	else {
	  m_selectionResults.append((attributeSet[i] + 1) + ",");
	}
      }

      for (int i=0;i<attributeSet.length;i++) {
	m_selectionResults.append("                     "
				  +m_trainInstances
				  .attribute(attributeSet[i]).name()
				  +"\n");
      }
    }

    // Cross validation should be called from here
    if (m_doXval == true) {
      m_selectionResults.append(CrossValidateAttributes()); 
    }

    // set up the attribute filter with the selected attributes
    if (m_selectedAttributeSet != null && !m_doXval) {
      m_attributeFilter = new Remove();
      m_attributeFilter.setAttributeIndicesArray(m_selectedAttributeSet);
      m_attributeFilter.setInvertSelection(true);
      m_attributeFilter.setInputFormat(m_trainInstances);  
    }

    // Save space
    m_trainInstances = new Instances(m_trainInstances, 0);
  }

  /**
   * Perform attribute selection with a particular evaluator and
   * a set of options specifying search method and options for the
   * search method and evaluator.
   *
   * @param ASEvaluator an evaluator object
   * @param options an array of options, not only for the evaluator
   * but also the search method (if any) and an input data file
   * @param train the input instances
   * @return the results of attribute selection as a String
   * @exception Exception if incorrect options are supplied
   */
  public static String SelectAttributes (ASEvaluation ASEvaluator, 
					 String[] options, 
					 Instances train)
    throws Exception {
    int seed = 1, folds = 10;
    String foldsString, seedString, searchName;
    String classString;
    String searchClassName;
    String[] searchOptions = null; //new String [1];
    ASSearch searchMethod = null;
    boolean doCrossVal = false;
    int classIndex = -1;
    boolean helpRequested = false;
    AttributeSelection trainSelector = new AttributeSelection();

    try {
      if (Utils.getFlag('h', options)) {
	helpRequested = true;
      }

      // does data already have a class attribute set?
      if (train.classIndex() != -1)
	classIndex = train.classIndex() + 1;

      // get basic options (options the same for all attribute selectors
      classString = Utils.getOption('c', options);

      if (classString.length() != 0) {
	if (classString.equals("first")) {
	  classIndex = 1;
	} else if (classString.equals("last")) {
	  classIndex = train.numAttributes();
	} else {
	  classIndex = Integer.parseInt(classString);
	}
      }

      if ((classIndex != -1) && 
	  ((classIndex == 0) || (classIndex > train.numAttributes()))) {
	throw  new Exception("Class index out of range.");
      }

      if (classIndex != -1) {
	train.setClassIndex(classIndex - 1);
      }
      else {
	//	classIndex = train.numAttributes();
	//	train.setClassIndex(classIndex - 1);
      }
      
      foldsString = Utils.getOption('x', options);

      if (foldsString.length() != 0) {
	folds = Integer.parseInt(foldsString);
	doCrossVal = true;
      }
      
      trainSelector.setFolds(folds);
      trainSelector.setXval(doCrossVal);

      seedString = Utils.getOption('n', options);

      if (seedString.length() != 0) {
	seed = Integer.parseInt(seedString);
      }

      trainSelector.setSeed(seed);

      searchName = Utils.getOption('s', options);

      if ((searchName.length() == 0) && 
	  (!(ASEvaluator instanceof AttributeEvaluator))) {
	throw  new Exception("No search method given.");
      }

      if (searchName.length() != 0) {
	searchName = searchName.trim();
	// split off any search options
	int breakLoc = searchName.indexOf(' ');
	searchClassName = searchName;
	String searchOptionsString = "";

	if (breakLoc != -1) {
	  searchClassName = searchName.substring(0, breakLoc);
	  searchOptionsString = searchName.substring(breakLoc).trim();
	  searchOptions = Utils.splitOptions(searchOptionsString);
	}
      }
      else {
	try {
	  searchClassName = new String("weka.attributeSelection.Ranker");
	  searchMethod = (ASSearch)Class.
	    forName(searchClassName).newInstance();
	}
	catch (Exception e) {
	  throw  new Exception("Can't create Ranker object");
	}
      }

      // if evaluator is a subset evaluator
      // create search method and set its options (if any)
      if (searchMethod == null) {
	searchMethod = ASSearch.forName(searchClassName, searchOptions);
      }

      // set the search method
      trainSelector.setSearch(searchMethod);
    }
    catch (Exception e) {
      throw  new Exception('\n' + e.getMessage() 
			   + makeOptionString(ASEvaluator, searchMethod));
    }

    try {
      // Set options for ASEvaluator
      if (ASEvaluator instanceof OptionHandler) {
	((OptionHandler)ASEvaluator).setOptions(options);
      }

      /* // Set options for Search method
	 if (searchMethod instanceof OptionHandler)
	 {
	 if (searchOptions != null)
	 {
	 ((OptionHandler)searchMethod).setOptions(searchOptions);
	 }
	 }
	 Utils.checkForRemainingOptions(searchOptions); */
    }
    catch (Exception e) {
      throw  new Exception("\n" + e.getMessage() 
			   + makeOptionString(ASEvaluator, searchMethod));
    }

    try {
      Utils.checkForRemainingOptions(options);
    }
    catch (Exception e) {
      throw  new Exception('\n' + e.getMessage() 
			   + makeOptionString(ASEvaluator, searchMethod));
    }

    if (helpRequested) {
      System.out.println(makeOptionString(ASEvaluator, searchMethod));
      System.exit(0);
    }

    // set the attribute evaluator
    trainSelector.setEvaluator(ASEvaluator);

    // do the attribute selection
    trainSelector.SelectAttributes(train);

    // return the results string
    return trainSelector.toResultsString();
  }


  /**
   * Assembles a text description of the attribute selection results.
   *
   * @return a string describing the results of attribute selection.
   */
  private String printSelectionResults () {
    StringBuffer text = new StringBuffer();
    text.append("\n\n=== Attribute Selection on all input data ===\n\n" 
		+ "Search Method:\n");
    text.append(m_searchMethod.toString());
    text.append("\nAttribute ");

    if (m_ASEvaluator instanceof SubsetEvaluator) {
      text.append("Subset Evaluator (");
    }
    else {
      text.append("Evaluator (");
    }

    if (!(m_ASEvaluator instanceof UnsupervisedSubsetEvaluator) 
	&& !(m_ASEvaluator instanceof UnsupervisedAttributeEvaluator)) {
      text.append("supervised, ");
      text.append("Class (");

      if (m_trainInstances.attribute(m_trainInstances.classIndex())
	  .isNumeric()) {
	text.append("numeric): ");
      }
      else {
	text.append("nominal): ");
      }

      text.append((m_trainInstances.classIndex() + 1) 
		  + " " 
		  + m_trainInstances.attribute(m_trainInstances
					       .classIndex()).name() 
		  + "):\n");
    }
    else {
      text.append("unsupervised):\n");
    }

    text.append(m_ASEvaluator.toString() + "\n");
    return  text.toString();
  }


  /**
   * Make up the help string giving all the command line options
   *
   * @param ASEvaluator the attribute evaluator to include options for
   * @param searchMethod the search method to include options for
   * @return a string detailing the valid command line options
   * @throws Exception if something goes wrong
   */
  private static String makeOptionString (ASEvaluation ASEvaluator, 
					  ASSearch searchMethod)
    throws Exception {
    
    StringBuffer optionsText = new StringBuffer("");
    // General options
    optionsText.append("\n\nGeneral options:\n\n");
    optionsText.append("-h\n\tdisplay this help\n");
    optionsText.append("-i <name of input file>\n");
    optionsText.append("\tSets training file.\n");
    optionsText.append("-c <class index>\n");
    optionsText.append("\tSets the class index for supervised attribute\n");
    optionsText.append("\tselection. Default=last column.\n");
    optionsText.append("-s <class name>\n");
    optionsText.append("\tSets search method for subset evaluators.\n");
    optionsText.append("-x <number of folds>\n");
    optionsText.append("\tPerform a cross validation.\n");
    optionsText.append("-n <random number seed>\n");
    optionsText.append("\tUse in conjunction with -x.\n");

    // Get attribute evaluator-specific options
    if (ASEvaluator instanceof OptionHandler) {
      optionsText.append("\nOptions specific to " 
			 + ASEvaluator.getClass().getName() 
			 + ":\n\n");
      Enumeration enu = ((OptionHandler)ASEvaluator).listOptions();

      while (enu.hasMoreElements()) {
	Option option = (Option)enu.nextElement();
	optionsText.append(option.synopsis() + '\n');
	optionsText.append(option.description() + "\n");
      }
    }

    if (searchMethod != null) {
      if (searchMethod instanceof OptionHandler) {
	optionsText.append("\nOptions specific to " 
			   + searchMethod.getClass().getName() 
			   + ":\n\n");
	Enumeration enu = ((OptionHandler)searchMethod).listOptions();

	while (enu.hasMoreElements()) {
	  Option option = (Option)enu.nextElement();
	  optionsText.append(option.synopsis() + '\n');
	  optionsText.append(option.description() + "\n");
	}
      }
    }
    else {
      if (ASEvaluator instanceof SubsetEvaluator) {
	System.out.println("No search method given.");
      }
    }

    return  optionsText.toString();
  }


  /**
   * Main method for testing this class.
   *
   * @param args the options
   */
  public static void main (String[] args) {
    try {
      if (args.length == 0) {
	throw  new Exception("The first argument must be the name of an " 
			     + "attribute/subset evaluator");
      }

      String EvaluatorName = args[0];
      args[0] = "";
      ASEvaluation newEval = ASEvaluation.forName(EvaluatorName, null);
      System.out.println(SelectAttributes(newEval, args));
    }
    catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.47 $");
  }
}

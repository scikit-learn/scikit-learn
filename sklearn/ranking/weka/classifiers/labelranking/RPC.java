package weka.classifiers.labelranking;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Random;
import java.util.Vector;
import weka.classifiers.AbstractClassifier;
import weka.classifiers.Evaluation;
import weka.classifiers.SingleClassifierEnhancer;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.WeightedInstancesHandler;
import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Utils;
import weka.core.converters.XArffLoader;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;

/**
 * The RPC class will learn Label Rankings by using the 
 * "Ranking by Pairwise Comparison" (RPC) algorithm. 
 * For more information, see:
 * "Label Ranking by Learning Pairwise Preferences"
 * By Eyke H�llermeier, Johannes F�rnkranz, Weiwei Cheng and Klaus Brinker
 * @author Weiwei Cheng
 * @author Alexander Balz
 */


public class RPC extends SingleClassifierEnhancer implements OptionHandler, WeightedInstancesHandler{

	Classifier[][] cl;
    Classifier[][] clCopy;
    protected int voting=0;
    protected String classifierString = defaultClassifierString();
 //   protected int labelSequence=0;
    
    /** Voting mode selection: soft voting */
    public static final int SELECTION_SOFT = 0;
    /** Voting mode selection: binary voting */
    public static final int SELECTION_BIN = 1;
    
    /** Choosing between ordering and ranking: Ordering */
 //   public static final int SELECTION_ORDERING = 1;
    /** Choosing between ordering and ranking: Ordering */
 //   public static final int SELECTION_RANKING = 0;

    /** Attribute selection methods */
    public static final Tag [] TAGS_SELECTION = {
      new Tag(SELECTION_SOFT, "Soft voting"),
      new Tag(SELECTION_BIN, "Binary voting")
    };
    /*
    public static final Tag [] TAGS_SELECTION_ORDERING = {
        new Tag(SELECTION_ORDERING, "Ordering"),
        new Tag(SELECTION_RANKING, "Ranking")
    };*/
    
	/**
	 * 
	 */
	private static final long serialVersionUID = 7652692558884463154L;

	
	protected void setClassifierString(String cls){
		classifierString = cls;
	}
	
	protected String getClassifierString(){
		return classifierString;
	}
	
	public RPC() throws Exception{
		String[] options = new String[]{"-W",getClassifierString(),"-V",String.valueOf(voting)};
	    			
		super.setOptions(options);
		setOptions(options);
	}
    
	/**
	 * Capabilities necessary for the ranker.
	 */
	 public Capabilities getCapabilities() {
		    Capabilities result = super.getCapabilities();   // returns the object from weka.classifiers.Classifier
		    
		    result.disableAll();
		    // attributes
		    result.enable(Capability.NOMINAL_ATTRIBUTES);
		    result.enable(Capability.NUMERIC_ATTRIBUTES);
		    result.enable(Capability.DATE_ATTRIBUTES);
		    result.enable(Capability.MISSING_VALUES);
		    result.enable(Capability.PREFERENCE_ATTRIBUTE);
		    result.enableAllAttributes();
		 
		    // class
		    //result.enable(Capability.NOMINAL_CLASS);
		    result.enable(Capability.MISSING_CLASS_VALUES);
		    result.enable(Capability.RANKING);
		 
		    return result;
	}

	 /**
	  * Displays a short description of the classifier.
	  * @return
	  */
	 public String globalInfo(){
		 return "Label Ranking by pairwise comparison.";
	 }
	 
	 /**
	  * Builds the classifier with given dataset.
	  * @throws Exception 
	  */
	 @Override
	 public void buildClassifier(Instances data) throws Exception{
		 cl = getClassifier(data,data.numAttributes()-1);
	 }
	  
	 /**
	  * returns information about the built classifier
	  */
	 public String toString(){
		 String text="Base classifier: "+getClassifierString()+ "Voting type: ";
		 if(voting==1)text+="binary";
		 else text+="soft";
		 text+="Instances include ";
		 //if(labelSequence==0)text+= "rankings.";
		 //else text+="orderings.";
		 
		 return text;
	 }
	  
	 /** Only for numeric predictions, so not really what we were looking for.
     * Used for prediction: returns a double value containing
	 * the first label in a label ranking.
	 * @throws Exception 
	 *
	 @Override
	 public double classifyInstance(Instance inst) throws Exception {
		 double[] feature = new double[inst.numAttributes()-1]; 
		 for(int i=0; i<feature.length; i++){
			 feature[i]=inst.value(i);
		 }
		 
		 double[] labelRank;
		 if(getVoting().getSelectedTag().getID()==0)
			 labelRank = getLabelRanking(cl.clone(),feature,"soft");			 
		 else
			 labelRank = getLabelRanking(cl.clone(),feature,"bin");	
		 
		 int[] rnk = new int[labelRank.length];
		 for(int k=0; k<labelRank.length; k++){
			 rnk[k] = (int)labelRank[k];
		 }
		 
		 PreferenceDenseInstance pdi = (PreferenceDenseInstance)inst;
		 int key=0;
		 for(int j=0; ; j++){
			 if(pdi.getHashMap().get(j)!=null){
				 key = j;
				 break;
			 }
		 }
		 pdi.insertPrediction(key, rnk);
		 
		 //returning label with highest order in ranking.
		 for(int k=0; k<labelRank.length; k++){
			 if(labelRank[k]==0.0)return (double)k;
		 }
		 return labelRank[0];
	 }*/
	 	 
	 /**
	  * Inside of the output of distributionForInstance, we store the Label Rankings
	  * predicted by the Ranking algorithm, so we will not have distributions over
	  * the separate labels like in nominal case.
	  * @param a single Instance of a data set.
	  * @return a predicted Label Ranking of that Instance.
	  */
	 @Override
	 public double[] distributionForInstance(Instance inst) throws Exception{
	//	 RankUtilities.ordering = labelSequence;
		 double[] feature = new double[inst.numAttributes()-1]; 
		 for(int i=0; i<feature.length; i++){
			 feature[i]=inst.value(i);
		 }
		 
		 double[] labelRank;
		 if(getVoting().getSelectedTag().getID()==0)
			 labelRank = getLabelRanking(cl.clone(),feature,"soft");			 
		 else
			 labelRank = getLabelRanking(cl.clone(),feature,"bin");	
		 
		 int[] rnk = new int[labelRank.length];
		 for(int k=0; k<labelRank.length; k++){
			 rnk[k] = (int)labelRank[k];
		 }
		 
		 //Prediction Hashmap and Preferences Hashmap must possess the same key
		 //for being able to relate them.
		 PreferenceDenseInstance pdi = (PreferenceDenseInstance)inst;
		 int key=0;
		 for(int j=0; ; j++){
			 if(pdi.getHashMap().get(j)!=null){
				 key = j;
				 break;
			 }
		 }
		 pdi.insertPrediction(key, rnk);
//		 inst = pdi;
		 
		 return labelRank;
	 }
 	 
	    /**
	     * Create c(c-1)/2 base classifiers (c is number of labels) with predefined methods.
	     * @param tuples The training data
	     * @param len_f The length of feature vector
	     * @return c(c-1)/2 base classifiers with predefined methods (c is number of labels) 
	     * @throws Exception
	     */
	    public Classifier[][] getClassifier(Instances tuples, int len_f) throws Exception{
	    	//The rankings have to be read out again and put inside of the label(s) array.
	        double[] feature = new double[len_f];
	        double[] label = new double[tuples.getLabels().size()];
	        
	        ArrayList<double[]> features = new ArrayList<double[]>();
	        ArrayList<double[]> labels = new ArrayList<double[]>();
       
	        for(int i=0; i<tuples.size()-1; i++){
	            for(int j=0; j<len_f; j++){
	                Arrays.fill(feature, j, j+1, tuples.get(i).value(j));
	            }   
	            features.add(feature);
	            feature = new double[len_f];//Necessary!
	            
	            //pairwise rankings are extracted.
	            for(int z=len_f; z<tuples.get(0).numAttributes();z++){
	            	PreferenceDenseInstance pdi;
	            	int[][] rankIndices = new int[1][1];
	            	if(tuples.get(i) instanceof PreferenceDenseInstance){
	            		pdi = (PreferenceDenseInstance) tuples.get(i);
	            		for(int ii=0; ; ii++){
	            			if(pdi.getHashMap().get(ii)!=null){
	            				rankIndices = pdi.getHashMap().get(ii);
	            				break;
	            			}        			
	            		}

	            		//generating an (in)complete ranking if possible.
	            		if(RankUtilities.noPartialRanking(rankIndices)){
	            			label = RankUtilities.triangleLabels;	            		
	            		}
	            		//when having a partial ranking, we give all included pairwise rankings to RPC.
	            		else{
	            			for(int n=0; n<rankIndices.length; n++){
	            				for(int o=0; o<rankIndices.length; o++){
	            					if(rankIndices[n][o]==1){
	            						labels.add(new double[]{n,o});
	            					}
	            				}
	            			}
	            			label[0]=Double.POSITIVE_INFINITY;
	            		}
	            	}
	            }	
	    		           	
	                       	            
	            /** transform between order and rank  
	             *  (depends on the data set. if not needed, simply add array label)
	             *  This is already done in the RankUtilities class, because it is needed
	             *  for the ranking measures, too.
	             *  Start here*/
/*	            double[] label_trans = new double[label.length];
	            for(int k=0; k<label.length; k++){
	                label_trans[(int)label[k]] = k;
	            }
	            label = label_trans;*/
	            /** End here*/
	            //excluding missing class values. Now, the altered array bounds have to be taken into account.
	            if(label[0]!=Double.POSITIVE_INFINITY){
	            	labels.add(label);
	            }
//	            label = new double[tuples.get(0).numAttributes()-len_f];//Necessary!
	            if(tuples.getLabels().size()!=0)
	            	label = new double[tuples.getLabels().size()];	   
	            else
	            	label = new double[RankUtilities.labels.size()];
	        }

	        //Label ranks are needed here.(NOT order of labels)
	        //Logistic Regression for L1L2 L1L3 ..... (index start from 1)
	        Classifier[][] clfer = new AbstractClassifier[label.length][label.length+1];
	        
	        FastVector attributes = new FastVector();
	        FastVector classAttValues = new FastVector();
	        
	        for(int i=0; i<feature.length; i++){
	            attributes.addElement(new Attribute("A" + (int)(i+1)));
	        }
	        
	        //list all possible values for class attribute 
	        for(int i=0; i<2; i++){
	            classAttValues.addElement(new Double(i).toString());
	        }
	        	        
	        attributes.addElement(new Attribute("class", classAttValues));
	        
	        for(int m=1; m<label.length; m++){
	            for(int n=m+1; n<label.length+1; n++){
	                
	                Instances data = new Instances("data", attributes, 0);
	                data.setClassIndex(data.numAttributes()-1);
	                //data.setClassIndex(data.);
	                
//	                for(int i=0; i<tuples.size(); i++){
	                  for(int i=0; i<features.size(); i++){
	                    Instance inst = new DenseInstance(data.numAttributes());
	                    for(int j=0; j<feature.length; j++){
	                        inst.setValue(j, features.get(i)[j]);
	                    }
	                    	                                        
	                    //compare the ranks between two labels, determine the class value 	                    
	                      try{
	                    	  if(labels.get(i).length>1 && m-1<labels.get(i).length && n-1<labels.get(i).length && labels.get(i)[m-1]<labels.get(i)[n-1]){
	                    		  inst.setValue(data.classIndex(), 1);
	                    	  }
	                    	  else{
	                    		  inst.setValue(data.classIndex(), 0);
	                    	  }
	                      }
	                      catch(Exception e){
	                    	  //e.printStackTrace();
	                    	  //If we land here, it is tried to classify a ranking that has been thrown
	                    	  //out of the adjacence matrix for receiving the upper triangular matrix.
	                    	  //Because they don't exist anymore, we do nothing in here.
	                      }
	                    	data.add(inst);	                    
	                }
	                clfer[m][n] 
	                		= AbstractClassifier.makeCopy(getClassifier());
	                clfer[m][n].buildClassifier(data);
	            }
	        }
	        
	        return clfer;	        
	    }
	    
	    /**
	     * 
	     * @param clfer The predefined binary base classifiers 
	     * @param feature The feature vector of input instance
	     * @param option Binary voting or soft voting
	     * @return A label ranking of given instance (Not an order of labels!)
	     * @throws Exception
	     */
	    public double[] getLabelRanking(Classifier[][] clfer, double[] feature, String option) throws Exception{
	    	
	        double[] label = new double[clfer.length];
	        double[] labelranking = new double[clfer.length];
	        double[] sum_vote = new double[clfer.length+1];
	        double[] temp = new double[clfer.length];
	        
	        //transform feature to instance
	        Instance instance = new DenseInstance(feature.length);
	        for(int i=0; i<feature.length; i++){
	            instance.setValue(i, feature[i]);
	        }
	        
	        /*create dummy instances (necessary for certain classifier, e.g. J48)*/
	        FastVector attributes = new FastVector();
	        FastVector classAttValues = new FastVector();
	        for(int i=0; i<feature.length; i++){
	            attributes.addElement(new Attribute("A" + (int)(i+1)));
	        }
	        
	        //list all possible values for class attribute 
	        
	        for(int i=0; i<2; i++){
	            classAttValues.addElement(new Double(i).toString());
	        }
	         
	        
	        //add them into attribute list
	        attributes.addElement(new Attribute("class", classAttValues));
	        
	        Instances dummyInstances = new Instances("data", attributes, 0);
	        dummyInstances.setClassIndex(dummyInstances.numAttributes()-1);
	        
	        instance.setDataset(dummyInstances);
	        
	        /**
	         * Voting
	         */       
	        try{
	        	if(clfer.length==0)
	        			clfer = cl.clone();
	        for(int i=1; i<clfer[0].length; i++){
	            for(int j=1; j<clfer[0].length; j++){
	                if(i==j){
	                    continue;
	                }
	                else if(i<j){
	                    try{
	                        if(option.equals("bin")){
	                            sum_vote[i] = sum_vote[i] + clfer[i][j].classifyInstance(instance);
	                        }
	                        else{
	                        	sum_vote[i] = sum_vote[i] + clfer[i][j].distributionForInstance(instance)[1];
	                        }
	                    }catch(Exception e){//if there is not data for current model, vote zero for both.
	                        e.printStackTrace();
	                    }
	                    
	                }    
	                else{
	                    try{
	                        if(option.equals("bin")){
	                            sum_vote[i] = sum_vote[i] - clfer[j][i].classifyInstance(instance) + 1;
	                        }
	                        else{
	                            sum_vote[i] = sum_vote[i] - clfer[j][i].distributionForInstance(instance)[1] + 1;
	                        }                        
	                    }catch(Exception e){
	                        e.printStackTrace();
	                    }                    
	                }
	            }
	        }
	        }catch(Exception e){
	        	e.printStackTrace();
	        	e.getCause();	        	
	        }
	        
	        for(int i=0; i<temp.length; i++){
	            temp[i] = sum_vote[i+1]; //+  Double.MIN_VALUE * (i+1);
	        }
	        
	        int[] a = Utils.sort(temp);//sort in ascending order, return origin index 
	        for(int i=0; i<temp.length; i++){
	            label[temp.length-i-1] = a[i]+1;
	        }
	          
	        //from order to rank
	        for(int i=0; i<labelranking.length; i++){
	            labelranking[(int)label[i]-1] = i;
	        }
	        
	        return labelranking;
	    }
	    
	    /**
	     * String describing default classifier.
	     * 
	     * @return the default classifier classname
	     */
	    protected String defaultClassifierString() {	      
	      return "weka.classifiers.functions.Logistic";
	    }
	    
	    /**
	     * Method to set the  base Classifier.
	     */
	    public void setClassifier(Classifier cls){
	    	m_Classifier = cls;
	    }
	    
	    /**
	     * Returning the base classifier.
	     */
	    public Classifier getClassifier(){
	    	return m_Classifier;
	    }
	    
	    /**
	     * Returns an enumeration describing the available options.
	     * 
	     * <!-- options-start -->
	     *  <pre>
	     *  -W 		The base classifier to be used.
	     *  -V 		The voting strategy in RPC. "bin" is for binary voting. 
	     *  		Everything else: soft-voting.
	     *  </pre>
	     * <!-- options-end -->
	     *
	     * @return an enumeration of all the available options.
	     */
	    public Enumeration<Option> listOptions() {
	      Vector<Option> newVector = new Vector<Option>(4);
	      
//	      newVector.addElement(new Option("\tSet the base classifier", "W", 1, "-W"));
	      newVector.addElement(new Option("\tSet voting mode", "V",1,"-V"));
	      newVector.addElement(new Option(
	              "",
	              "", 0, "\nOptions specific to classifier "
	              + getVoting() + ":"));
	     /* newVector.addElement(new Option("\tSet ordering/ranking mode", "O",1,"-O"));
	      newVector.addElement(new Option(
	              "",
	              "", 0, "\nOptions specific to classifier "
	              + getLabelSequence() + ":"));*/
	      
	      Enumeration<?> enu = super.listOptions();
	      while (enu.hasMoreElements()) {
	        newVector.addElement((Option) enu.nextElement());
	      }
	      
	      return newVector.elements();
	    }

	    /**
	     * Parses a given list of options. <p/>
	     *
	     <!-- options-start -->
	     * Valid options are: <p/>
	     * 
	     * <pre> -W
	     *  Full name of base classifier.
	     *  (default: weka.classifiers.trees.DecisionStump)</pre>
	     *  <pre> -V
	     *  The voting strategy.
	     *  (default: soft-voting)</pre>
	     * 
	     <!-- options-end -->
	     *
	     * Options after -- are passed to the designated classifier.<p>
	     *
	     * @param options the list of options as an array of strings
	     * @throws Exception if an option is not supported
	     */
	    public void setOptions(String[] options) throws Exception {
	      super.setOptions(options);
	      setVoting(new SelectedTag(Integer.parseInt(Utils.getOption("V", options)), TAGS_SELECTION));
	      //setLabelSequence(new SelectedTag(Integer.parseInt(Utils.getOption("O", options)), TAGS_SELECTION_ORDERING) );
	    }

	    /**
	     * Gets the current settings of the Classifier.
	     *
	     * @return an array of strings suitable for passing to setOptions
	     */
	    public String[] getOptions() {
	    	   String[] superOptions = super.getOptions();
	    	   String[] options = new String[superOptions.length + 3];
	    	    
	    	    int current = 0;
	    	    	    	    
	    	    	options[current++]="-V";
	    	    	options[current++]=String.valueOf(getVoting().getSelectedTag().getID());
	    	    	//options[current++]="-O";
	    	    	//options[current++]=String.valueOf(getLabelSequence().getSelectedTag().getID());
	    	    	
	    	    System.arraycopy(superOptions, 0, options, current, 
	    	        superOptions.length);

	    	    current += superOptions.length;
	    	    while (current < options.length) {
	    	      options[current++] = "";
	    	    }
	    	    return options;
	    }

	    
	    /**
	     * Returns the tip text for this property
	     * @return tip text for this property suitable for
	     * displaying in the explorer/experimenter gui
	     */
	    public String votingTipText() {
	      return "Set the voting strategy for the base classifier. In case of a " +
	      		 "binary voting, the classifyInstance method of the base classifier" +
	      		 "will be used for voting. Otherwise, voting will be done by base classifier's" +
	      		 "DistributionForInstance method.";
	    }


	    /**
	     * Gets the selected voting method.
	     * @return the method to use.
	     */
	    public SelectedTag getVoting() {	      
	      return new SelectedTag(voting, TAGS_SELECTION);
	    }
	    
	    /**
	     * Sets the voting mode.
	     * @param method
	     */
	    public void setVoting(SelectedTag method){
	    	if (method.getTags() == TAGS_SELECTION) {
		        voting = method.getSelectedTag().getID();
		      }
	    }
	    
    
	    public static void main (String[] args)throws Exception{
	    	//Making an evaluation like described on the paper.
	    	Evaluation eval;
	    	Instances inst;
	    	String path;
	    	RPC rpc;

	    		path = args[0];	
	    		rpc = new RPC();	
	    		
	    		System.out.println(rpc.getCapabilities());
	    	
	    		inst = rpc.init(path, rpc, "iris_dense.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using iris_dense data set:");	
	    		double[] res = rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path,  rpc,"wine_dense.xarff");	
	    		eval = new Evaluation(inst);
	    		System.out.println("Using wine data set:");	
	    		res = rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path,  rpc,"glass_dense.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using glass data set:");	
	    		res = rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path, rpc,"vowel_dense.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using vowel data set:");	
	    		res = rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path, rpc,"vehicle_dense.xarff");	
	    		eval = new Evaluation(inst);
	    		System.out.println("Using vehicle data set:");	
	    		res =rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		//Ordering data sets...
	    		/*
	    		inst = rpc.init(path, rpc,"spo.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using spo dataset:");	
	    		res = rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path, rpc,"heat.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using heat dataset:");	
	    		res =rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path, rpc,"dtt.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using dtt dataset:");	
	    		res =rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path, rpc,"cold.xarff");	
	    		eval = new Evaluation(inst);
	    		System.out.println("Using cold dataset:");	
	    		res = rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    		
	    		inst = rpc.init(path, rpc,"diau.xarff");		
	    		eval = new Evaluation(inst);
	    		System.out.println("Using diau dataset:");	
	    		res =rpc.createStats(rpc,eval, inst);
	    		System.out.println("kendall's tau: "+res[0]+"\t spearman correlation: "+res[1]+"\n\n");
	    			*/
	    	}
	    	
	    	/**
	    	 * Computing kendall's tau and spearman correlation with 5 repeats.
	    	 * @param d
	    	 * @param inst
	    	 * @return
	    	 * @throws Exception
	    	 */
	    	public double[] createStats(RPC rpc, Evaluation eval, Instances inst) throws Exception{
	    		int seed=5;
	    		for(int i=1; i<seed+1; i++){
	    				eval.crossValidateModel("weka.classifiers.labelranking.RPC", inst, 10, new String[]{"-V","0"}, new Random(i));
	    		}
	    		return new double[]{eval.tau(), eval.spearman()[1]};
	    	}
	    	
	    	/**
	    	 * Setting up the instances and our classifier.
	    	 * @param d
	    	 * @param filename
	    	 * @throws Exception
	    	 */
	    	public Instances init(String path, RPC rpc, String filename) throws Exception{
	    		File file;
	    		XArffLoader xarff = new XArffLoader();
	    		file = new File(path+filename);
	    		xarff.setFile(file);
	    		Instances inst = xarff.getDataSet();
	    		inst.setClassIndex(inst.numAttributes()-1);
	    		rpc.buildClassifier(inst);
	    		return inst;
	    	}
	    
	    
	    	 
}

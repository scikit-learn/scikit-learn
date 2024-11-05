package weka.classifiers.labelranking;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import weka.classifiers.AbstractClassifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.converters.XArffLoader;
import weka.core.labelranking.PreferenceDenseInstance;
import weka.core.labelranking.RankUtilities;
import weka.core.neighboursearch.LinearNNSearch;

/**
 * In this class, a label ranker based on the Plackett-Luce Model
 * by means of minorization-maximization technique is implemented. 
 * For more details, see: "Label Ranking Methods based on the Plackett-Luce Model"
 * by Weiwei Cheng, Krzysztof Demczyński and Eyke Hüllermeier.
 * This file contains a transferring to java of the plackmm algorithm written
 * in MATLAB. The original version can be found on 
 * http://www.stat.psu.edu/~dhunter/code/btmatlab/
 * @author Alexander Balz
 */
public class PLMM extends AbstractClassifier implements OptionHandler, WeightedInstancesHandler {
	protected int KNN=5;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -7304847934573818076L;
	
	Instances instances = null;
	public static ArrayList<int[]> matlabDataset = new ArrayList<int[]>();
	double[] model;
	ArrayList<int[]> instData;
	int M;
	int N;
	int P;
	
	/**
	 * Constructor - taking the options given by the super class.
	 * @throws Exception
	 */
	public PLMM() throws Exception{
		String[] options = new String[]{"-K", String.valueOf(KNN)};
		super.setOptions(options);
		setOptions(options);
	}
	
	
	/**
	 * Converting the Instances' input to
	 * a new input form like in the original MATLAB version.
	 * Then, the model is built with the Plackett-Luce-Model.
	 */
	@Override
	public void buildClassifier(Instances data) throws Exception {
		instances = data;
		instData = new ArrayList<int[]>();
		PreferenceDenseInstance pdi;
		int[][] hashMap = null;
		double[] ranking = null;
		
		//First, the conversion to another input format is done.
		for(int i=0; i<data.numInstances(); i++){
			if(data.get(i) instanceof PreferenceDenseInstance)
				pdi = (PreferenceDenseInstance)data.get(i);
			else
				throw new Exception("The Instances do not contain rankings.");
			
			
			//Obtain the hashMap of pairwise rankings.
			int key=0;
			while(true){
				if(pdi.getHashMap().get(key)!=null){
					hashMap = pdi.getHashMap().get(key);
					break;
				}
				key++;
			}
			
			if(RankUtilities.cycles(hashMap)){
				throw new Exception("The input data set contains cycles, they need to be cleared before continuing.");
			}
			
			if(RankUtilities.noPartialRanking(hashMap)){
				ranking = RankUtilities.triangleLabels;
			}
			
			for(int j=0; j<ranking.length; j++){
				int[] newData = new int[3];
				newData[0]=j;
				newData[1]=i;
				for(int k=0; k<ranking.length; k++){
					if(ranking[k]==j){
						newData[2]=k;
					}
				}
				instData.add(newData);
			}
			
		}
		
//		PLMMDataset.fillList("/home/alex/Desktop/PLMValues.txt");
		if(matlabDataset.size()!=0){
			instData = matlabDataset;
		}
		
		//Now, we begin to build the model by means of Plackett-Luce.
		M = max(instData,0)+1;		//The amount of "individuals"
		N = max(instData,1)+1;		//The amount of "contests"
		P = max(instData,2)+1;		//largest amount of "individuals per contest"
		
		
		double[][] f = new double[P][N];	//f[i][j] = one who placed i in contest j
		double[][] r = new double[M][N];	//r[i][j] = place of i in contest j.
		
		for(int i=0; i<instData.size(); i++){
			if(instData.get(i)[0]==0)
				f[instData.get(i)[2]][instData.get(i)[1]]=instData.get(i)[0];
			else
				f[instData.get(i)[2]][instData.get(i)[1]]=instData.get(i)[0]+1;
			if(instData.get(i)[2]+P*(instData.get(i)[1])==0)
				r[instData.get(i)[0]][instData.get(i)[1]]=instData.get(i)[2]+P*(instData.get(i)[1]);
			else
				r[instData.get(i)[0]][instData.get(i)[1]]=instData.get(i)[2]+P*(instData.get(i)[1])+1;

		}
		
		
		double[][] r2 = new double[r.length][r[0].length];
		
		for(int rx=0; rx<r.length; rx++){
			for(int ry=0; ry<r[0].length; ry++){
				r2[rx][ry] = r[rx][ry];
			}
		}
		
		
		double[] w = new double[M];  // w[i] = # times i placed higher than last.
		double[] pp = sumThreshold(f,0);  //pp[j] = #individuals in contest.
		
		for(int i=0; i<N; i++){
			double[] tmp = matrixFromTo(f,(pp[i]-1),i);
			//w = vectorQuery(w,tmp);
			for(int y=0; y<tmp.length; y++){
				if((int) tmp[y]==0)
					w[(int) tmp[y]]+=1;
				else
					w[(int) tmp[y]-1]+=1;
			}
		}
		
		int[] temp = new int[pp.length];
		int current = P;
		for(int p=1; p<=(N-1)*P; p++){
			temp[p]=current;
			if(current>=(N-1)*P) break;
			current+=P;
		}
		for(int q=0; q<pp.length; q++){
			pp[q]=pp[q]+temp[q];
		}
		
		double[] gamma = ones(M);			//(unscaled) initial gamma vector.
		
		
		double[] dgamma= new double[]{1};
		int iterations=0;
		
		while(norm(dgamma,2)>1e-09){
			iterations++;
			double[][]g = new double[f.length][f[0].length];
			
			for(int i=0; i<g.length; i++){
				for(int j=0; j<g[0].length; j++){
					g[i][j]=f[i][j];
				}
			}
			
			for(int y=0; y<f[0].length; y++){
				for(int x=0; x<f.length; x++){
					if(f[x][y]>0){
						g[x][y]= gamma[(int)f[x][y]-1];
					}
				}
			}
			double[][] tmpMatrix = new double[P][g[0].length];
			
			int newX=0;
			for(int x=P-1; x>=0; x--){
				for(int y=0; y<g[0].length; y++){
					tmpMatrix[newX][y]=g[x][y];
				}
				newX++;
			}
			
			g=cumSum(tmpMatrix);
			
			tmpMatrix = new double[P][g[0].length];
			
		    newX=0;
			for(int x=P-1; x>=0; x--){
				for(int y=0; y<g[0].length; y++){
					tmpMatrix[newX][y]=g[x][y];
				}
				newX++;
			}
			
			g=tmpMatrix;
			
				int ii=1;	
				int position=0;
				for(int k=0; k<g[0].length; k++){
					for(int j=0; j<g.length; j++){
						if(pp[position]==ii){
							g[j][k]=0;
							position++;
						}
						ii++;
					}
				}
				
			//}
			
			for(int y=0; y<g[0].length; y++){
				for(int x=0; x<g.length; x++){
					if(g[x][y]>0){
						g[x][y] = 1/g[x][y];
					}
				}
			}
			
			g=cumSum(g);
			
		
			for(int y=0; y<r[0].length; y++){
				for(int x=0; x<r.length; x++){
					if(r[x][y]>0){
						r2[x][y]= getPosition((int)r[x][y]-1,g);
						//	g[(int)r[x][y]-1];
					}
				}
			}
			
			
			
			double[] sr2 = rowSum(r2);
			
			double[] newgamma = new double[sr2.length];
			
			for(int i=0; i<newgamma.length; i++){
				newgamma[i]=w[i]/sr2[i];
			}
			
			dgamma = new double[newgamma.length];
			
			for(int a=0; a<dgamma.length; a++){
				dgamma[a]=newgamma[a]-gamma[a];
			}
	
			gamma = newgamma;
		}
		
		double sum =0;
		for(int b=0; b<gamma.length; b++){
			sum+=gamma[b];
		}
		
		model = new double[gamma.length];
		for(int c=0; c<model.length; c++){
			model[c]=gamma[c]/sum;
		}
		
		System.out.println(printModel()+"\n\n");
//		System.out.println("Model length: "+model.length);
	}
	
	
	/**
	 * Predicts a label ranking for instance inst.
	 */
	 @Override
	 public double[] distributionForInstance(Instance inst) throws Exception{
		 int actualElement=0;
		 double[] prediction = new double[instances.getLabels().size()];
		 LinearNNSearch lns = new LinearNNSearch(instances);
		 Instances neighbours = lns.kNearestNeighbours(inst, KNN);
		 //The result of buildClassifier can be found in the 'model' array.
		 buildClassifier(neighbours);
		 	 
		 while(true){
			 boolean allNegative=true;
			 //1. Element of max: maximum value in array, 2. element: its index inside of the array.
			 double[] max= new double[]{0,0};
			 //Sorting the elements, the labels with highest probability will be preferred.
			 for(int i=0; i<model.length; i++){
				 if(model[i]>0) allNegative = false;
				 if(model[i]>max[0]){
					 max[0]=model[i];
					 max[1]=i;
				 }
			 }
			 if(allNegative)break;
			 prediction[actualElement] = max[1];
			 actualElement++;
			 model[(int)max[1]]=-1;
		 }
		 actualElement=0;
		 return prediction;
	 }

	
	/**
	 * Evaluating on multiple data sets.
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception{
		PLMM pl = new PLMM();
		
		insertMatlabData("/home/alex/Desktop/PLMValues.txt");
		XArffLoader xarff = new XArffLoader();
		xarff.setFile(new File("/home/alex/Desktop/RankingDatasets/iris_dense.xarff"));
		Instances dat = xarff.getDataSet();
		
		pl.buildClassifier(dat);
//		double[] out = pl.model;
		
//		for(int i=0; i<out.length; i++){
//			System.out.println(out[i]+"\n");
//		}
	}

	@Override
	public Enumeration<?> listOptions() {
		Enumeration<?> out = super.listOptions();
		return out;
	}

	@Override
	public void setOptions(String[] options) throws Exception {
		super.setOptions(options);
		Integer.parseInt(Utils.getOption("K", options));	
	}

	@Override
	public String[] getOptions() {
		int current=0;
		String[] opt = new String[super.getOptions().length+2];
		 for(int i=0; i<super.getOptions().length; i++){
			 opt[i]=super.getOptions()[i];
			 current++;
		 }
		
		opt[current++]="-K";
    	opt[current++]=String.valueOf(KNN);
		return opt;
	}
		
	@Override
	public String toString(){
		return getClass().getName() + '@' + Integer.toHexString(hashCode());
	}
	
	public String globalInfo(){
		return "PLMM is an instance-based classifier using the Plackett-Luce model" +
				"in order to build a classifier model. Predictions are made by applying a " +
				"k-nearest neighbour algorithm based on minorization-maximization technique" +
				"to classify a query instance.";
	}
		
	@Override
	public Capabilities getCapabilities(){
		   Capabilities result = super.getCapabilities();   
		   // returns the object from weka.classifiers.Classifier
		    
		    result.disableAll();
		    // attributes
		    result.enable(Capability.NOMINAL_ATTRIBUTES);
		    result.enable(Capability.NUMERIC_ATTRIBUTES);
		    result.enable(Capability.DATE_ATTRIBUTES);
		    result.enable(Capability.MISSING_VALUES);
		    result.enable(Capability.PREFERENCE_ATTRIBUTE);
		    result.enableAllAttributes();
		 
		    // class
		    result.enable(Capability.MISSING_CLASS_VALUES);
		    result.enable(Capability.RANKING);
		 
		    return result;
	}
	
	/**
	 * Returns the maximum value in the column 
	 * of a given list.
	 * @param list
	 * @param column
	 * @return
	 */
	public int max(ArrayList<int[]> list, int column){
		int max=0;
		for(int i=0; i<list.size(); i++){
			int tmp = list.get(i)[column];
			if(tmp>max){
				max = tmp;
			}
		}
		return max;
	}

	/**
	 * Counting the elements per column which are bigger
	 * than a given value an putting in into a new array.
	 * (Like in MATLAB for sum(f>0) with threshold 0.
	 * @param list
	 * @return
	 */	
	public double[] sumThreshold(double[][] input, int threshold){
		double[] output = new double[input[0].length];
		for(int j=0; j<input[0].length; j++){
			for(int i=0; i<input.length; i++){
				if(input[i][j]>=threshold){
					output[j]++;
				}
			}
		}
		return output;
	}
	
	/**
	 * determining the maximum length of 
	 * the arrays contained in the list.
	 * @param list
	 * @return
	 */
	public int longestArray(ArrayList<int[]> list){
		int longest=0;
		for(int i=0; i<list.size(); i++){
			if(list.get(i).length>longest){
				longest = list.get(i).length;
			}
		}
		return longest;
	}
	
	/**
	 * Returns a vector filled with ones.
	 * @param M
	 * @return
	 */
	public double[] ones(int M){
		double[] vector = new double[M];
		for(int i=0; i<M; i++){
			vector[i]=1;
		}
		return vector;
	}
	
	
	/**
	 * Those two norm-functions are equivalent to their MATLAB
	 * versions.
	 * @param dgamma
	 * @return
	 *
	public double norm(int dgamma){
		return dgamma;
	}*/
	
	public double norm(double[] dgamm, int p){
		double sum=0;
		for(int i=0; i<dgamm.length; i++){
	//		for(int j=0; j<dgamma[0].length; j++){
				sum+= Math.pow(Math.abs(dgamm[i])	,p);
		//	}
		}
		sum = Math.pow(sum, (1/(double)p));
		return sum;
	}
	
	/**
	 * Returns the cumulative sum of a matrix.
	 * @param input
	 * @return
	 */
	public double[][] cumSum(double[][] input){
		double[] sums = new double[input[0].length];
		
		for(int i=0; i<input.length; i++){
			for(int j=0; j<input[0].length; j++){
				sums[j]+=input[i][j];
				input[i][j] = sums[j];
			}
		}
		
		return input;
	}
	
	/**
	 * Returns a visualization of a matrix.
	 * @param input
	 * @return
	 */
	public String matrixToString(double[][] input){
		String matrix ="";
		for(int i=0; i<input.length; i++){
			for(int j=0; j<input[0].length; j++){
				matrix+=input[i][j]+",";
			}
			matrix = matrix.substring(0, matrix.length()-1);
			matrix+="\n";
		}
		return matrix;
	}
	
	/**
	 * Equivalent to a sum(input,2) call in MATLAB. 
	 * adding the elements of each row in a matrix and
	 * store it into a vector.
	 * @param input
	 * @return
	 */
	public double[] rowSum(double[][] input){
		double[] output = new double[input.length];
		for(int i=0; i<input.length; i++){
			for(int j=0; j<input[0].length; j++){
				output[i]+=input[i][j];
			}
		}
		return output;
	}
	
	/**
	 * Returns "to" elements of a specified
	 * matrix column.
	 * @param matrix
	 * @param to
	 * @param column
	 * @return
	 */
	public double[] matrixFromTo(double[][] matrix, double to, int column){
		double[] output = new double[(int)to];
		for(int x=0; x<to; x++){
			try{
				output[x] = matrix[x][column];
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		
		return output;
	}
		   
    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String KNNTipText() {
      return "Sets the k nearest neighbours for prediction of new rankings.";
    }

    /**
     * Gets the selected voting method.
     * @return the method to use.
     */
    public int getKNN() {	      
      return KNN;
    }
    
    /**
     * Sets the voting mode.
     * @param method
     */
    public void setKNN(int input){
    	KNN = input;
    }
    
    

	/**
	 * If desired, it is possible to insert a data set  in *.txt format
	 * which is structured like a matlab-dataset, like the nascar2002 
	 * dataset from http://www.stat.psu.edu/~dhunter/code/btmatlab/nascar2002.txt
	 * @param fileName
	 * @throws IOException
	 */
	public static void insertMatlabData(String fileName) throws IOException{
		String tmp="";
		File file = new File(fileName);
		
		BufferedReader fileIn = new BufferedReader(new FileReader(file));
		
		while((tmp = fileIn.readLine())!=null){
			String elem = "";
			int[] data = new int[3];
			int datCount=0;
			for(int i=0; i<tmp.length(); i++){
				if(tmp.charAt(i)==' '){
					data[datCount]=Integer.parseInt(elem)-1;
					datCount++;
					elem="";
				}
				else if(tmp.charAt(i)!='\n'){					
					elem+=tmp.charAt(i);
				}
				if(tmp.charAt(i)=='\n' || i==tmp.length()-1){
					data[datCount]=Integer.parseInt(elem)-1;
				}
			}
			matlabDataset.add(data);
		}
	}	 
	
	
	/**
	 * Returns the pos. element of a matrix, iterating column-wise
	 * and then row-wise.
	 * @param pos
	 * @param matrix
	 * @return
	 */
	public double getPosition(int pos, double[][] matrix){
		
		int actual=0;
		for(int j=0; j<matrix[0].length; j++){
			for(int i=0; i<matrix.length; i++){
				actual++;
				if(actual==pos){
					return matrix[i][j];
				}
			}
		}
		
		return 0;
	}


	/**
	 * Returns a string containing the generated model.
	 * @return
	 */
	public String printModel(){
		String out="";
		for(int i=0; i<model.length; i++){
			out+=model[i]+", ";
		}
		return out;
	}


}

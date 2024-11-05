package weka.core.labelranking;

import java.util.List;

/**
 * The RankUtilities class provides some useful functions
 * for handling label rankings like cycle detection
 * and a function to find partial rankings.
 * @author Alexander Balz
 * 
 */
public class RankUtilities{
	//Variable to check if we got partial rankings in our dataset.
	public static boolean noPartialRank=true;
	public static double[] triangleLabels;
	public static int[][] triangle;
	//Are rankings or orderings inside of the instances? (0 for rankings, 1 for orderings).
	//public static int ordering=0;
	//a safety copy if the instances class lost labels, like in filtering case.
	public static List<String[]> labels;
	/**
	 * a value of 0 means that the value with equal coordinates in pairLabels has not yet been processed.
	 * 1 means that a node is being processed and 2 means that this node has already been processed.
	 */	
	private static int[][] visited;
	
	private static int[][] pairs;
	
	/**
	 * Cycle detection is based on a depth-first search, marking all visited values in a separate matrix.
	 * It consists of 2 methods cycles and findCycles.
	 * cycles determines the first node to be processed and calls the 2. method to do the depth-first search.
	 * Afterwards, cycles tests, if our rankings contain several graphs that have no connection to each other
	 * and tries to find cycles in each one of them.
	 * @param pairLabels
	 * @return
	 * @throws Exception 
	 */
	public static boolean cycles(int[][] pairLabels){
		visited = new int[pairLabels.length][pairLabels.length];
		pairs = pairLabels;

		int[] firstNode = new int[2];
		boolean first=false;
		//find first node in pairLabels matrix.
		for(int x=0; x<pairs.length; x++ ){			
			for(int y=0; y<pairs.length; y++){
				if(pairs[x][y]==1){
					firstNode[0]=x;
					firstNode[1]=y;
					first=true;
					break;
				}
			}
			if(first)break;
		}
		
		try{
			findCycles(firstNode);
		}
		catch(Exception e){return true;}
		
		//If a graph consists of branches without connection to each other, we process these branches here.		
		while(true){
			boolean allNodesProcessed=true;
			for(int i=0; i<pairs.length; i++){
				for(int j=0; j<pairs.length; j++){
					if(i==j&&pairs[i][j]==1)return true;
					if(pairs[i][j]==1 && visited[i][j]!=2){
						
						try{
							findCycles(new int[]{i,j});
						}
						catch(Exception e){
							return true;
						}
						allNodesProcessed = false;
					}
				}
			}
			if(allNodesProcessed)break;
		}
		return false;
	}
	
	/**
	 * Here, our depth-first search and marking processed nodes are done.
	 * @param pairLabels
	 * @param node
	 * @return
	 * @throws Exception 
	 */
	private static void findCycles(int[] node) throws Exception{
		if(visited[node[0]][node[1]]==1)throw new Exception("cycle detected!");
		else if(visited[node[0]][node[1]]==0){
			visited[node[0]][node[1]]=1;
			//Each neighbor-node will be visited.
			for(int i=0; i<visited.length; i++){
				if(pairs[node[1]][i]==1){
					findCycles(new int[]{node[1],i});						
				}
				visited[node[0]][node[1]]=2;
			}
		}
	}
	
	
	/**
	 * This function tests, if a ranking is complete or not
	 * by counting the ones in the input matrix per line.
	 * These numbers have to be unique if we have a complete ranking. 
	 * @param the pairwise label rankings in matrix form
	 * @return a boolean value if ranking is complete or not.
	 */
	public static boolean noPartialRanking(int[][] pairLabels){
		boolean unique = true;
		int[] countsLine = new int[pairLabels.length];
		int[] countsCol = new int[pairLabels.length];
		double[] pLabels = new double[pairLabels.length];
		int[][] newPairs;
		int newLength=pairLabels.length;
		int[][]triangular = new int[1][1];
		
		//Calculating transitive closure so we can see circles in the pairwise error matrix later.	
		//Problem about this was that not only cycles will appear in the matrix, but other ranking pairs, too.
//		if(RankUtilities.cycles(pairLabels)){
//			RankUtilities.transitiveClosure(pairLabels);
//		}
		
		//defining the original order of labels in the matrix.
		for(int i=0; i<pLabels.length; i++){
			pLabels[i] = i;
		}
		
		//counting ones per line
		for(int i=0; i<pairLabels.length; i++){
			for(int j=0; j<pairLabels.length; j++){
				if(pairLabels[i][j]==1)countsLine[i]++;
			}
		}
		//counting ones per column.
		for(int i=0; i<pairLabels.length; i++){
			for(int j=0; j<pairLabels.length; j++){
				if(pairLabels[j][i]==1)countsCol[i]++;
			}
		}
		
		//finding empty lines and columns, calculating length of new array.
		for(int k=0; k<countsLine.length; k++){
			if(countsCol[k]==0 && countsLine[k]==0){				
				for(int n=0; n<pairLabels.length; n++){
					if(pLabels[n]!=Double.POSITIVE_INFINITY && pLabels[n]>pLabels[k]){
					  pLabels[n]--;
					}
				}
				pLabels[k]=Double.POSITIVE_INFINITY;
				newLength--;
			}
		}	
		
		//Deleting empty lines and columns, generating a new matrix without them.
		newPairs = new int[newLength][newLength];
		
		int infty=0; 
		int inftyrow=0;
		
		for(int a=0; a<pairLabels.length; a++){
			if(pLabels[a]==Double.POSITIVE_INFINITY){
				inftyrow++;
				continue;
			}
			for(int b=0; b<pairLabels.length; b++){
				if(pLabels[b]!=Double.POSITIVE_INFINITY /*&& a<newPairs.length && (b-infty)<newPairs.length*/){
					newPairs[a-inftyrow][b-infty]=pairLabels[a][b];
				}
				else{
					infty++;
				}
			}
			infty=0;
		}
		
		countsLine = new int[newPairs.length];
		countsCol = new int[newPairs.length];
		
		//counting ones per line for our new array.
		for(int i=0; i<newPairs.length; i++){
			for(int j=0; j<newPairs.length; j++){
				if(newPairs[i][j]==1)countsLine[i]++;
			}
		}
		//counting ones per column for our new array.
		for(int i=0; i<newPairs.length; i++){
			for(int j=0; j<newPairs.length; j++){
				if(newPairs[j][i]==1)countsCol[i]++;
			}
		}
		
		
		//checking if number of ones per line is unique or not. If it is unique, we have to
		//go on trying to calculate an upper triangular matrix. If it is not, this ranking has to be incomplete.
		for(int x=0; x<countsLine.length; x++){
			for(int y=0; y<countsLine.length; y++){
				if(x!=y){
					if(countsLine[x]==countsLine[y])unique=false;
				}
			}
		}
		
		
				
		//here, we change the order of our labels. Ptr. shows to the first position and is replaced with the label including the highest amount of ones
		//in its line. Then ptr is increased to the next label and so on.
		if(unique){
						
			int ptr=0;
			for(int z=countsLine.length-1; z>=0; z--){
				infty=0;
				for(int zz=0; zz<pLabels.length; zz++){
					if(pLabels[zz]==Double.POSITIVE_INFINITY){
						infty++;
						continue;
					}
					if(z==countsLine[zz-infty]){
						pLabels[zz/*+infty*/]=ptr;
						ptr++;
					}
				}
			}
		
			//After this, in the triangular[][] field, a triangular upper matrix is held in case of an (in)complete ranking; otherwise it isn't.
			//To do that, we have to change label positions in our matrix. Instead, we generate a new one - triangular, so the original one
			//will not be manipulated.
			triangular = new int[newPairs.length][newPairs.length];			
					
			
			
			boolean containsInfinity=false;
			for(int inf=0; inf<pLabels.length; inf++){
				if(pLabels[inf]==Double.POSITIVE_INFINITY)containsInfinity=true;
			}
			
			if(containsInfinity==true){
				int iPlus=0;
				int jPlus=0;
				for(int i=0; i<newPairs.length; i++){
				
					while(pLabels[iPlus]==Double.POSITIVE_INFINITY){
						iPlus++;
					}
				
					for(int j=0; j<newPairs.length; j++){
						while(pLabels[jPlus]==Double.POSITIVE_INFINITY){
							jPlus++;
						}
						if((int)pLabels[iPlus]<triangular.length && (int)pLabels[jPlus]<triangular.length){
							triangular[(int)pLabels[iPlus]][(int)pLabels[jPlus]] = newPairs[i][j];
						}
						jPlus++;
					}
					jPlus=0;
					iPlus++;
				}
			}
			else{
				for(int i=0; i<newPairs.length; i++){
					for(int j=0; j<newPairs.length; j++){
						triangular[(int)pLabels[i]][(int)pLabels[j]]=newPairs[i][j];
					}
				}
			}
			
			int zeroes=0;
			
			for(int i=0; i<triangular.length; i++){
				for(int j=0; j<triangular.length; j++){
					if(j<=zeroes){
						if(triangular[i][j]!=0)return false;
					}
					else{
						if(triangular[i][j]!=1)return false;
					}
				}
				zeroes++;
			}
			
		}
		
		triangle = triangular;
		triangleLabels = pLabels;
		
		return unique;
	}
	
	/**
	 * Calculating the transitive Closure to derive all pairwise rankings
	 * if a given ranking was incomplete.
	 * @param the pairwise label rankings in matrix form
	 *
	public static void transitiveClosure(int[][] pairLabels){
		for(int k=0; k<pairLabels.length; k++){
			for(int i=0; i<pairLabels.length; i++){
				if(pairLabels[i][k]==1){
					for(int j=0; j<pairLabels.length; j++){
						if(pairLabels[k][j]==1)pairLabels[i][j]=1;
					}
				}
			}
		}
	}*/
	
	/**
	 * Converts a double[] with ranking information into a String representation.
	 */
	public String doubleRankingToString(double[] readLabels){		
		String ranking="";
		for(int rn=0; rn<labels.size(); rn++){
			for(int in=0; in<readLabels.length; in++){				
				if(readLabels[in]==rn){
					ranking+=labels.get(in)+">";
				}
			}
		}
		ranking = ranking.substring(0,ranking.length()-1);
		return ranking;
	}
	
	     
}

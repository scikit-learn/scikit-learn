/*
 * Created on 12.02.2007
 */
package weka.classifiers.labelranking;

import java.util.ArrayList;
import java.util.Arrays;
import weka.core.Utils;

/**
 * @author Weiwei Cheng
 */
public class TreeLabelRankingEvaluation {
	
	/**
	 * @param a The first ranking
	 * @param b The second ranking
	 * @return Spearman's rank correlation coefficient. 
	 * A rank correlation coefficient is in the interval [-1,1] where:
	 * 		If the agreement between the two rankings is perfect (i.e., the two rankings are the same) the coefficient has value 1.
	 * 		If the disagreement between the two rankings is perfect (i.e., one ranking is the reverse of the other) the coefficient has value -1.
	 * 		For all other arrangements the value lies between -1 and 1, and increasing values imply increasing agreement between the rankings.
	 * 		If the rankings are completely independent, the coefficient has value 0.
	 * @throws Exception 
	 */
	public static double spearman(ArrayList<Double> a, ArrayList<Double> b) throws Exception{
		
		if(a.size()!=b.size()){
			throw new Exception("The sizes of given rankings are not same.");
		}
		
		else{
			
			double dif_squ_sum = 0;
			for(int i=0; i<a.size(); i++)
				dif_squ_sum = dif_squ_sum + (a.get(i) - b.get(i)) * (a.get(i) - b.get(i));
			//Careful! a.size()*a.size()*a.size() is huge.
			return 1-((6*dif_squ_sum/a.size())/(a.size()+1))/(a.size()-1);//rho
		}

	}
	
	/**
	 * Note: the index of the ranking is started from 0.
	 * 
	 * @param a The first ranking
	 * @param b The second ranking
	 * @return Kendall tau rank correlation coefficient.
  	 * A rank correlation coefficient is in the interval [-1,1] where:
  	 * 		If the agreement between the two rankings is perfect (i.e., the two rankings are the same) the coefficient has value 1.
  	 * 		If the disagreement between the two rankings is perfect (i.e., one ranking is the reverse of the other) the coefficient has value -1.
  	 * 		For all other arrangements the value lies between -1 and 1, and increasing values imply increasing agreement between the rankings.
  	 * 		If the rankings are completely independent, the coefficient has value 0.
	 * @throws Exception 
	 */
	public static double kendall(ArrayList<Double> a, ArrayList<Double> b) throws Exception{
		
		double tau = -1;
		
		if(a.size()!=b.size()){
			throw new Exception("The sizes of given rankings are not same.");
		}
		
		else{
			/**
			 * 1. Sort ArrayList a and get the new vision of b.
			 */			
			ArrayList<Double> b_adjust = new ArrayList<Double>(a.size());
			for(int i=0; i<a.size(); i++){
				b_adjust.add(b.get(a.indexOf((double)i)));
			}
			
			/**
			 * 2. Compute the times of pairwise adjacent transpositions.
			 */
			int transposition = 0;
			for(int i=0; i<a.size()-1; i++)
				for(int j=i+1; j<a.size(); j++){
    				if(b_adjust.get(i)>b_adjust.get(j))
    					transposition++;
			}
			
			/**
			 * 3. Compute tau.
			 */
			tau = 1-(4D*transposition/(a.size())/(a.size()-1));
			return tau;
		}
	
	}
	
	/**
	 * This function is the same with kendall(ArrayList<Double> a, ArrayList<Double> b), 
	 * only takes different type of arguments. Note: the index of the ranking is started from 0. 
	 * @param a
	 * @param b
	 * @return Kendall tau rank correlation coefficient.
	 * @throws Exception 
	 */
    public static double kendall(double[] a, double[] b) throws Exception {

        double tau = -1;
        
        if(a.length!=b.length){
            throw new Exception("The sizes of given rankings are not same.");
        }
        
        else{
            /**
             * 1. Sort array a and get the new vision of b.
             */         
            double[] b_adjust = new double[a.length];
            for(int i=0; i<a.length; i++){
                int[] sort = Utils.sort(a);
                b_adjust[i] = b[sort[i]];
            }
            
            /**
             * 2. Compute the times of pairwise adjacent transpositions.
             */
            int transposition = 0;
            for(int i=0; i<a.length-1; i++)
                for(int j=i+1; j<a.length; j++){
                    if(b_adjust[i]>b_adjust[j])
                        transposition++;
            }
            
            /**
             * 3. Compute tau.
             */
            tau = 1-(4D*transposition/(a.length)/(a.length-1));
            return tau;
        }
        
    }
	
	/**
	 * 
	 * @param x The first ranking
	 * @param y The second ranking
	 * @return Pearson correlation coefficient
	 * @throws Exception 
	 */
	public static double pearson(ArrayList<Double> x, ArrayList<Double> y) throws Exception{
		
		double r = -1;
		
		if(x.size()!=y.size()){
			throw new Exception("The sizes of given rankings are not same.");
		}
		
		else{
			double sumX = 0;
			double sumY = 0;
			double sumXY = 0;
			double sumXX = 0;
			double sumYY = 0;
			int len = x.size();
			
			for(int i=0; i<len; i++){
				sumX = sumX + x.get(i);
				sumY = sumY + y.get(i);
				sumXY = sumXY + x.get(i)*y.get(i);
				sumXX = sumXX + x.get(i)*x.get(i);
				sumYY = sumYY + y.get(i)*y.get(i);
			}
				
			r = (sumXY - sumX*sumY/len)/(Math.sqrt(sumXX-sumX*sumX/len)*Math.sqrt(sumYY-sumY*sumY/len));
			return r;
		}
	}
	
	/**
	 * Get the size of the symmetric difference of two sets, i.e., |Delta(X,Y)|. 
	 * X,Y are sets of different pairs within ranking x and y, respectively. 
	 * Delta(X,Y) is equal to the union of both relative complements OR
	 * the union of X and Y, minus their intersection.
	 * 
	 * A simple test shows that Option 2 is faster.
	 * 
	 * @param x the first ranking.
	 * @param y the second ranking.
	 * @return
	 */
	public static int symDiffSize(double[] x, double[] y){
	    
	    //get the sets of pairs of x and y.
        double[][] X = new double[x.length*(x.length-1)/2][2];
        double[][] Y = new double[y.length*(y.length-1)/2][2];
        
        int count=0;
        for(int j=0; j<x.length-1; j++){
            for(int k=j+1; k<x.length; k++){
                X[count][0] = x[j];
                X[count][1] = x[k];
                count++;
            }
        }
        
        count=0;
        for(int j=0; j<y.length-1; j++){
            for(int k=j+1; k<y.length; k++){
                Y[count][0] = y[j];
                Y[count][1] = y[k];
                count++;
            }
        }
        
        //Option 1: 
        //get the union and intersection of pairs. 
//	    ArrayList<double[]> union = new ArrayList<double[]>(); 
//	    ArrayList<double[]> intersection = new ArrayList<double[]>();
//	    
//	    loopX:
//	    for(int i=0; i<X.length; i++){
//	        for(int j=0; j<Y.length; j++){
//	            if(Arrays.equals(X[i], Y[j])){
//	                union.add(X[i]);
//	                intersection.add(X[i]);
//	                continue loopX;
//	            }
//	        }
//	        union.add(X[i]);
//	    }
//        
//	    loopY:
//	    for(int i=0; i<Y.length; i++){
//	        for(int j=0; j<union.size(); j++){
//	            if(Arrays.equals(Y[i], union.get(j))){
//	                continue loopY;
//	            }
//	        }
//	        union.add(Y[i]);
//	    }
//
//        for(int i=0; i<intersection.size(); i++){
//            union.remove(intersection.get(i));
//        }
//
//        return union.size();
        
        //Option 2: 
        //get X-Y and Y-X.
        ArrayList<double[]> complementXY = new ArrayList<double[]>();
        ArrayList<double[]> complementYX = new ArrayList<double[]>();
        
        loopXY:
        for(int i=0; i<X.length; i++){
            for(int j=0; j<Y.length; j++){
                if(Arrays.equals(X[i], Y[j])){
                    continue loopXY;
                }
            }
            complementXY.add(X[i]);
        }
        
        loopYX:
        for(int i=0; i<Y.length; i++){
            for(int j=0; j<X.length; j++){
                if(Arrays.equals(X[j], Y[i])){
                    continue loopYX;
                }
            }
            complementYX.add(Y[i]);
        }
        
        return complementXY.size()+complementYX.size();
        
	} 
	
	/**
	 * 
	 * @param x The first top K (predicted top K)
	 * @param y The second top K (true top K)
	 * @return kofK
	 * kofK is a measurement between two top K ranking. 
	 * 		It is equal to the number of the same entries of 
	 * 		both top K ranking divided by K.
	 * @throws Exception 
	 */
	public static double kofK(int[] x, int[] y) throws Exception{
		
		double k = 0;
		
		int K = x.length;
		
		if(K!=y.length){
			throw new Exception("The lengths of given rankings are not the same.");
		}
		
		else{
			for(int i=0; i<K; i++)
				for(int j=0; j<K; j++){
					if(x[i]==y[j]){
						k++;
						break;
					}
			}
			return k/K;
		}
	}
	
    /**
     * 
     * @param k  k-combination
     * @param n  with n elements 
     * @return The number of k-combinations from a set with n elements.
     * @throws Exception
     */
    public static long combination(int k, int n) throws Exception{
        
        if(k>n){
            throw new Exception("The first parameter is bigger than the second.");
        }
        
        return factorial(n)/(factorial(n-k)*factorial(k));
    }
    
    /**
     * 
     * @param n
     * @return the factorial of n.
     * @throws Exception 
     */
    public static long factorial(int n) throws Exception{
        long f = 1;
        if(n<0){
            throw new Exception("n is negative.");
        }
        else if(n>20){
            throw new Exception("Overflow. n is too big.");
        }
        for(int i=1; i<=n; i++){
            f =  f*i;
        }
        return f;
    }
    
    /**
     * 
     * @param a Matrix a
     * @param b Matrix b
     * @return Multiplication of a and b.
     * @throws Exception If the two matrices don't match.
     */
    public static double[][] matrixMultiplication(double[][] a, double[][] b) throws Exception{
        
        if(a[0].length!=b.length)
            throw new Exception("The first matrix's column and the second matrix's row don't match.");
        
        double[][] mx = new double[a.length][b[0].length];
        
        for(int m=0; m<mx.length; m++){
            for(int n=0; n<mx[0].length; n++){
               for(int i=0; i<b.length; i++){
            	   mx[m][n] = mx[m][n] + a[m][i]*b[i][n];
               }
            }
    	}
        
        return mx;
    }
    
    /**
     * 
     * @param a The first array.
     * @param b The second array.
     * @return The inner product between two arrays a and b.
     * @throws Exception If two arrays don't match.
     */
    public static double dotProduct(double[] a, double[] b) throws Exception{
        
        if(a.length!=b.length)
            throw new Exception("The lengths of arrays don't match.");

        double result = 0;
        for(int i=0; i<a.length; i++)
            result = result + a[i]*b[i];
        
        return result;
    }
    
    /**
     * cosine measure 
     * 
     * @param a The first vector.
     * @param b The second vector.
     * @return Cosine of the angle between two vectors. 
     */
    public static double cosine(double[] a, double[] b) throws Exception{
        
        if(a.length!=b.length)
            throw new Exception("The lengths of these two vectors don't match.");
        
        double dotproduct = 0;
        double a_norm_sqr = 0;
        double b_norm_sqr = 0;
        
        for(int i=0; i<a.length; i++){
            dotproduct = dotproduct + a[i]*b[i];
            a_norm_sqr = a_norm_sqr + a[i]*a[i];
            b_norm_sqr = b_norm_sqr + b[i]*b[i];
        }
        
        return dotproduct/Math.sqrt(a_norm_sqr)/Math.sqrt(b_norm_sqr);
    }
    
    /**
     * 
     * @param a double array
     * @return standard variance of array a 
     */
    public static double standarddeviation(double[] a){
//    	double sum = 0 ;
//    	double sqrsum = 0;
//    	for(int i=0; i<a.length; i++){
//    		sum = sum + a[i]; 
//    		sqrsum = sqrsum + a[i]*a[i];
//    	}
//    	
//    	return Math.sqrt(sqrsum/(a.length-1) - sum/a.length * sum/a.length - sum/a.length * sum/a.length/(a.length-1));
        return Math.sqrt(Utils.variance(a));
    }
    
    /**
     * 
     * @param a double array
     * @param b double array
     * @return t-statistic from a paired t-test
     * @throws Exception 
     */
    public static double pairedttest(double[] a, double[] b) throws Exception{
        
        if(a.length!=b.length)
            throw new Exception("The lengths of these two vectors don't match.");  

        double[] diff = new double[a.length];
        double diffsum = 0;
        
        for(int i=0; i<a.length; i++){
            diff[i] = b[i] - a[i];
            diffsum = diffsum + diff[i];
        }
        
        return diffsum/a.length * Math.sqrt(a.length) / standarddeviation(diff);
    }

    /**
     * 
     * @param a 
     * @param b
     * @return the euclidean distance between vector a and vector b.
     * @throws Exception
     */
    public static double euclideanDistance(double[] a, double[] b) throws Exception{
        
        if(a.length!=b.length){
            throw new Exception("The lengths of given vectors are not same.");
        }
        
        double squsum = 0;
        for(int i=0; i<a.length; i++){
            squsum = squsum + (a[i]-b[i])*(a[i]-b[i]);
        }
        
        return Math.sqrt(squsum);
    }
    
    /**
     * 
     * @param a
     * @param b
     * @param p
     * @return Minkowski distance between a and b of order p.
     * @throws Exception
     */
    public static double minkowskiDistance(double[] a, double[] b, double p) throws Exception{
        
        if(a.length!=b.length){
            throw new Exception("The lengths of given vectors are not same.");
        }
        if(p<0){
            throw new Exception("The degree p must be positive.");
        }
        
        double sum = 0;
        for(int i=0; i<a.length; i++){
            sum = sum + Math.pow(Math.abs(a[i]-b[i]), p);
        }
        
        return Math.pow(sum, 1/p);
    }
    
    /**
     * 
     * @param a object 1
     * @param b object 2
     * @param w weight vector
     * @return the weighted Manhattan distance
     * @throws Exception If the lengths don't match.
     */
    public static double weightedManhattanDistance(double[] a, double[] b, double[] w) throws Exception{
        
        if(w.length!=a.length||w.length!=b.length)
            throw new Exception("The lengths of vectors don't match.");        
        
        double dis = 0;
        for(int i=0; i<w.length; i++){
            dis = dis + w[i]*Math.abs(a[i]-b[i]);
        }
        return dis;        
    }
  
    /**
     * Calculate the difference (in terms of position) between the true nearest case and predicted nearest case.
     * 
     * @param order_true The true order.
     * @param order_pred The predicted order.
     * @return pos(Case_pred)-pos(Case_true)
     * @throws Exception If the true nearest case cannot be found.
     */
    public static int positionErrorNN(int[] order_true, int[] order_pred) throws Exception{
        for(int i=1; i<order_pred.length; i++){
            if(order_pred[i]==order_true[1]){
                return i-1; 
            }
        }       
        throw new Exception("The desired object cannot be found in predicted ordering.");
    }
    
    /**
     * Calculate the difference (in terms of position) between the true topmost case and predicted topmost case.
     * 
     * @param order_true The true order.
     * @param order_pred The predicted order.
     * @return pos(Case_pred)-pos(Case_true)
     * @throws Exception If the true nearest case cannot be found.
     */
    public static int positionError(int[] order_true, int[] order_pred) throws Exception{
        for(int i=0; i<order_pred.length; i++){
            if(order_pred[i]==order_true[0]){
                return i-1; 
            }
        }       
        throw new Exception("The desired object cannot be found in predicted ordering.");
    }
    
	/**
	 * The main method is for test purpose.
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
//		ArrayList<Double> al1 = new ArrayList<Double>();
//		ArrayList<Double> al2 = new ArrayList<Double>();
//		al1.add(0.0);al2.add(2.0);
//		al1.add(1.0);al2.add(3.0);
//		al1.add(2.0);al2.add(0.0);
//		al1.add(3.0);al2.add(1.0);
//		al1.add(4.0);al2.add(4.0);
//		al1.add(5.0);al2.add(6.0);
//		al1.add(6.0);al2.add(7.0);
//		al1.add(7.0);al2.add(5.0);
//
//		ArrayList<Double> a = new ArrayList<Double>();
//		ArrayList<Double> b = new ArrayList<Double>();
//		ArrayList<Double> c = new ArrayList<Double>();
//		ArrayList<Double> d = new ArrayList<Double>();	
//		ArrayList<Double> e = new ArrayList<Double>();
//		ArrayList<Double> f = new ArrayList<Double>();
//		ArrayList<Double> g = new ArrayList<Double>();
//		ArrayList<Double> h = new ArrayList<Double>();
//		
//		a.add(0.0); d.add(0.0);
//		a.add(1.0); d.add(2.0);
//		a.add(2.0); d.add(3.0);
//		a.add(3.0); d.add(1.0);
//
//		b.add(0.0); e.add(0.0); g.add(2.0);
//		b.add(1.0); e.add(3.0); g.add(3.0);
//		b.add(3.0); e.add(2.0); g.add(0.0);
//		b.add(2.0); e.add(1.0); g.add(1.0);
//		
//		c.add(0.0); f.add(1.0); h.add(3.0);
//		c.add(2.0); f.add(3.0); h.add(2.0);
//		c.add(1.0); f.add(2.0); h.add(1.0);
//		c.add(3.0); f.add(0.0); h.add(0.0);
//		
//		System.out.println(kendall(al1, al2));
//		System.out.println(spearman(al1, al2));
//		System.out.println(pearson(al1, al2));
//		System.out.println(kendall(a, a));
//		System.out.println(spearman(a, a));
//		System.out.println(pearson(a, a));
//		System.out.println(kendall(a, b));
//		System.out.println(spearman(a, b));
//		System.out.println(pearson(a, b));
//		System.out.println(kendall(a, c));
//		System.out.println(spearman(a, c));
//		System.out.println(pearson(a, c));
//		System.out.println(kendall(a, d));
//		System.out.println(spearman(a, d));
//		System.out.println(pearson(a, d));
//		System.out.println(kendall(a, e));
//		System.out.println(spearman(a, e));
//		System.out.println(pearson(a, e));
//		System.out.println(kendall(a, f));
//		System.out.println(spearman(a, f));
//		System.out.println(pearson(a, f));
//		System.out.println(kendall(a, g));
//		System.out.println(spearman(a, g));
//		System.out.println(pearson(a, g));
//		System.out.println(kendall(a, h));
//		System.out.println(spearman(a, h));
//		System.out.println(pearson(a, h));
//		
//		ArrayList<Double> s1 = new ArrayList<Double>();
//		ArrayList<Double> s2 = new ArrayList<Double>();
//		s1.add(0.0);s2.add(4.0);
//		s1.add(1.0);s2.add(1.0);
//		s1.add(2.0);s2.add(2.0);
//		s1.add(3.0);s2.add(3.0);
//		s1.add(4.0);s2.add(0.0);
//
//		System.out.println(spearman(s1, s2));
//		System.out.println(kendall(s1, s2));
//		System.out.println(pearson(s1, s2));
//		
//		ArrayList<Double> x = new ArrayList<Double>();
//		ArrayList<Double> y = new ArrayList<Double>();
//		x.add(3.0);y.add(4.0);
//		x.add(9.0);y.add(7.0);
//		x.add(2.0);y.add(5.0);
//		x.add(0.0);y.add(1.0);
//		x.add(8.0);y.add(9.0);
//		x.add(1.0);y.add(2.0);
//		x.add(5.0);y.add(8.0);
//		x.add(6.0);y.add(3.0);
//		x.add(7.0);y.add(6.0);
//		x.add(4.0);y.add(0.0);
//		System.out.println(spearman(y, x));
//		System.out.println(kendall(x, y));
//		System.out.println(pearson(x, y));
//		
//		x.clear();
//		y.clear();
//		x.add(83.0);y.add(141.0);
//		x.add(86.0);y.add(162.0);
//		x.add(88.0);y.add(161.0);
//		x.add(92.0);y.add(154.0);
//		x.add(94.0);y.add(171.0);
//		x.add(98.0);y.add(174.0);
//		x.add(101.0);y.add(184.0);
//		x.add(114.0);y.add(190.0);
//		x.add(117.0);y.add(187.0);
//		x.add(121.0);y.add(191.0);
//		System.out.println(pearson(x, y));
//		
//		x.clear();
//		y.clear();
//		x.add(0.0);y.add(2.0);
//		x.add(1.0);y.add(3.0);
//		x.add(2.0);y.add(0.0);
//		x.add(3.0);y.add(1.0);
//		x.add(4.0);y.add(4.0);
//		x.add(5.0);y.add(6.0);
//		x.add(6.0);y.add(7.0);
//		x.add(7.0);y.add(5.0);
//		x.add(8.0);y.add(5.0);
//		x.add(9.0);y.add(10.0);
//		x.add(10.0);y.add(11.0);
//		x.add(11.0);y.add(9.0);
//		System.out.println(kendall(x, y));
//		System.out.println(spearman(x, y));
//		
//		int[] x = {1,2,3,4,5};
//		int[] y = {0,0,0,0,0};
//		System.out.println(kofK(x, y));
		
//		double[][] a = {{1,2,3,4},{4,1,2,3},{4,1,2,3}};
//		double[][] b = {{3,1,3},{2,2,2},{1,5,8},{2,2,3}};
//		double[][] c = matrixMultiplication(a, b);

//        double[] a = {1,1};
//        double[] b = {0,1};
//        double[] c = {5d/9,1};
//        double[] d = {1,0};
//        double[] e = {1,1d/3};
//        double[] f = {1,1};
//        System.out.println(cosine(a,b));
//        System.out.println(cosine(c,d));
//        System.out.println(cosine(e,f));
  
//        double[] a = {106,84,110,91,109,91,111,107,121,105,
//					  99,94,119,88,118,97,103,106,95,106,
//					  85,106,101,105,96,105,107,128,111,101};
//        double[] b = {0.035,0.269,0.486,0.589,0.706};
//        double[] c = {0.373,0.56,0.636,0.729,0.768};
//        System.out.println(pairedttest(b,c));
//        System.out.println(standarddeviation(a));
//        System.out.println(standarddeviation(b));
//        System.out.println(standarddeviation(c));   
//	    
//        double[] a = {1,2,3,4};
//        double[] b = {4,3,2,1};
//        double[] c = {1,2,3,4,5};
//        System.out.println(symDiffSize(a, a));
//        System.out.println(symDiffSize(a, b));
//        System.out.println(symDiffSize(a, c));
//        System.out.println(symDiffSize(b, c));
//        
//        a = new double[100];
//        b = new double[100]; 
//        for(int i=0; i<100; i++){
//            a[i]=i;
//            b[i]=100-i;
//        }
//        long x = System.currentTimeMillis();
//        System.out.println(symDiffSize(a, a));
//        System.out.println(symDiffSize(a, b));
//        System.out.println(symDiffSize(a, c));
//	    System.out.println(System.currentTimeMillis()-x);
	    double[] q = {1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0};
	    double[] a = {1,1,0,0,0,0,0,0,2,0,0,0,1,0,2,0,1,1,1,0,0,1};
	    double[] b = {0,0,1,1,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,1,1,0};
	    double[] c = {4,0,0,0,1,1,0,0,4,1,1,1,0,0,4,1,0,0,0,0,0,0};
	    System.out.println(cosine(q, a));
	    System.out.println(cosine(q, b));
	    System.out.println(cosine(q, c));
    }

}
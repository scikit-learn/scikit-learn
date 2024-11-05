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
 *    EstimatorUtils.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.Vector;
 
/** 
 * Contains static utility functions for Estimators.<p>
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class EstimatorUtils
  implements RevisionHandler {
  
  /** 
   * Find the minimum distance between values
   * @param inst sorted instances, sorted
   * @param attrIndex index of the attribute, they are sorted after
   * @return the minimal distance
   */
  public static double findMinDistance(Instances inst, int attrIndex) {
    double min = Double.MAX_VALUE;
    int numInst = inst.numInstances();
    double diff;
    if (numInst < 2) return min;
    int begin = -1;
    Instance instance = null;
    do { 
      begin++;
      if (begin < numInst) 
	{ instance = inst.instance(begin); }
    } while (begin < numInst && instance.isMissing(attrIndex)); 

    double secondValue = inst.instance(begin).value(attrIndex);
    for (int i = begin; i < numInst && !inst.instance(i).isMissing(attrIndex);  i++) {
      double firstValue = secondValue; 
      secondValue = inst.instance(i).value(attrIndex);
      if (secondValue != firstValue) {
	diff = secondValue - firstValue;
	if (diff < min && diff > 0.0) {
	  min = diff;
	}
      }
    }
    return min;
  }

  /** 
   * Find the minimum and the maximum of the attribute and return it in 
   * the last parameter..
   * @param inst instances used to build the estimator
   * @param attrIndex index of the attribute
   * @param minMax the array to return minimum and maximum in
   * @return number of not missing values
   * @exception Exception if parameter minMax wasn't initialized properly
   */
  public static int getMinMax(Instances inst, int attrIndex, double [] minMax) 
    throws Exception {
    double min = Double.NaN;
    double max = Double.NaN;
    Instance instance = null;
    int numNotMissing = 0;
    if ((minMax == null) || (minMax.length < 2)) {
      throw new Exception("Error in Program, privat method getMinMax");
    }
    
    Enumeration enumInst = inst.enumerateInstances();
    if (enumInst.hasMoreElements()) {
      do {
	instance = (Instance) enumInst.nextElement();
      } while (instance.isMissing(attrIndex) && (enumInst.hasMoreElements()));
      
      // add values if not  missing
      if (!instance.isMissing(attrIndex)) {
	numNotMissing++;
	min = instance.value(attrIndex);
	max = instance.value(attrIndex);
      }
      while (enumInst.hasMoreElements()) {
	instance = (Instance) enumInst.nextElement();
	if (!instance.isMissing(attrIndex)) {
	  numNotMissing++;
	  if (instance.value(attrIndex) < min) {
	    min = (instance.value(attrIndex));
	  } else {
	    if (instance.value(attrIndex) > max) {	      
	      max = (instance.value(attrIndex));
	    }
	  }
	}
      }
    }
    minMax[0] = min;
    minMax[1] = max;
    return numNotMissing;
  }

  /**
   * Returns a dataset that contains all instances of a certain class value.
   *
   * @param data dataset to select the instances from
   * @param attrIndex index of the relevant attribute
   * @param classIndex index of the class attribute
   * @param classValue the relevant class value 
   * @return a dataset with only 
   */
  public static Vector getInstancesFromClass(Instances data, int attrIndex,
					     int classIndex,
					     double classValue, Instances workData) {
    //Oops.pln("getInstancesFromClass classValue"+classValue+" workData"+data.numInstances());
    Vector dataPlusInfo = new Vector(0);
    int num = 0;
    int numClassValue = 0;
    //workData = new Instances(data, 0);
    for (int i = 0; i < data.numInstances(); i++) {
      if (!data.instance(i).isMissing(attrIndex)) {
	num++;
	if (data.instance(i).value(classIndex) == classValue) {
	  workData.add(data.instance(i));
	  numClassValue++;
	}
      }
    } 

    Double alphaFactor = new Double((double)numClassValue/(double)num);
    dataPlusInfo.add(workData);
    dataPlusInfo.add(alphaFactor);
    return dataPlusInfo;
  }


  /**
   * Returns a dataset that contains of all instances of a certain class value.
   * @param data dataset to select the instances from
   * @param classIndex index of the class attribute
   * @param classValue the class value 
   * @return a dataset with only instances of one class value
   */
  public static Instances getInstancesFromClass(Instances data, int classIndex,
						double classValue) {
     Instances workData = new Instances(data, 0);
    for (int i = 0; i < data.numInstances(); i++) {
      if (data.instance(i).value(classIndex) == classValue) {
	workData.add(data.instance(i));
      }
     
    }
    return workData;
  }
  
    
   
  /**
   * Output of an n points of a density curve.
   * Filename is parameter f + ".curv".
   *
   * @param f string to build filename
   * @param est
   * @param min
   * @param max
   * @param numPoints
   * @throws Exception if something goes wrong
   */
  public static void writeCurve(String f, Estimator est, 
				double min, double max,
				int numPoints) throws Exception {

    PrintWriter output = null;
    StringBuffer text = new StringBuffer("");
    
    if (f.length() != 0) {
      // add attribute indexnumber to filename and extension .hist
      String name = f + ".curv";
      output = new PrintWriter(new FileOutputStream(name));
    } else {
      return;
    }

    double diff = (max - min) / ((double)numPoints - 1.0);
    try {
      text.append("" + min + " " + est.getProbability(min) + " \n");

      for (double value = min + diff; value < max; value += diff) {
	text.append("" + value + " " + est.getProbability(value) + " \n");
      }
      text.append("" + max + " " + est.getProbability(max) + " \n");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.out.println(ex.getMessage());
    }
    output.println(text.toString());    

    // close output
    if (output != null) {
      output.close();
    }
  }

  /**
   * Output of an n points of a density curve.
   * Filename is parameter f + ".curv".
   *
   * @param f string to build filename
   * @param est
   * @param classEst
   * @param classIndex
   * @param min
   * @param max
   * @param numPoints
   * @throws Exception if something goes wrong
   */
  public static void writeCurve(String f, Estimator est, 
				Estimator classEst,
				double classIndex,
				double min, double max,
				int numPoints) throws Exception {

    PrintWriter output = null;
    StringBuffer text = new StringBuffer("");
    
    if (f.length() != 0) {
      // add attribute indexnumber to filename and extension .hist
      String name = f + ".curv";
      output = new PrintWriter(new FileOutputStream(name));
    } else {
      return;
    }

    double diff = (max - min) / ((double)numPoints - 1.0);
    try {
      text.append("" + min + " " + 
		  est.getProbability(min) * classEst.getProbability(classIndex)
		  + " \n");

      for (double value = min + diff; value < max; value += diff) {
	text.append("" + value + " " + 
		    est.getProbability(value) * classEst.getProbability(classIndex)
		    + " \n");
      }
      text.append("" + max + " " +
		  est.getProbability(max) * classEst.getProbability(classIndex)
		  + " \n");
    } catch (Exception ex) {
      ex.printStackTrace();
      System.out.println(ex.getMessage());
    }
    output.println(text.toString());    

    // close output
    if (output != null) {
      output.close();
    }
  }

  
  /**
   * Returns a dataset that contains of all instances of a certain value
   * for the given attribute.
   * @param data dataset to select the instances from
   * @param index the index of the attribute  
   * @param v the value 
   * @return a subdataset with only instances of one value for the attribute 
   */
  public static Instances getInstancesFromValue(Instances data, int index,
					  double v) {
    Instances workData = new Instances(data, 0);
    for (int i = 0; i < data.numInstances(); i++) {
      if (data.instance(i).value(index) == v) {
	workData.add(data.instance(i));
      }
    } 
    return workData;
  }

   
  /**
   * Returns a string representing the cutpoints
   */
  public static String cutpointsToString(double [] cutPoints, boolean [] cutAndLeft) {
    StringBuffer text = new StringBuffer("");
    if (cutPoints == null) {
      text.append("\n# no cutpoints found - attribute \n"); 
    } else {
      text.append("\n#* "+cutPoints.length+" cutpoint(s) -\n"); 
      for (int i = 0; i < cutPoints.length; i++) {
	text.append("# "+cutPoints[i]+" "); 
	text.append(""+cutAndLeft[i]+"\n");
      }
      text.append("# end\n");
    }
    return text.toString();
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.4 $");
  }
}

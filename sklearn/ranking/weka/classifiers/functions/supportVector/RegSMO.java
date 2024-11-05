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
 *    RegSMO.java
 *    Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.functions.supportVector;

import weka.core.Instances;
import weka.core.Option;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.util.Enumeration;
import java.util.Vector;

/**
 <!-- globalinfo-start -->
 * Implementation of SMO for support vector regression as described in :<br/>
 * <br/>
 * A.J. Smola, B. Schoelkopf (1998). A tutorial on support vector regression.
 * <p/>
 <!-- globalinfo-end -->
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;misc{Smola1998,
 *    author = {A.J. Smola and B. Schoelkopf},
 *    note = {NeuroCOLT2 Technical Report NC2-TR-1998-030},
 *    title = {A tutorial on support vector regression},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -P &lt;double&gt;
 *  The epsilon for round-off error.
 *  (default 1.0e-12)</pre>
 * 
 * <pre> -L &lt;double&gt;
 *  The epsilon parameter in epsilon-insensitive loss function.
 *  (default 1.0e-3)</pre>
 * 
 * <pre> -W &lt;double&gt;
 *  The random number seed.
 *  (default 1)</pre>
 * 
 <!-- options-end -->
 *
 * @author  Remco Bouckaert (remco@cs.waikato.ac.nz,rrb@xm.co.nz)
 * @version $Revision: 1.4 $
 */
public class RegSMO 
  extends RegOptimizer
  implements TechnicalInformationHandler {
  
  /** for serialization */
  private static final long serialVersionUID = -7504070793279598638L;

  /** tolerance parameter, smaller changes on alpha in inner loop will be ignored **/
  protected double m_eps = 1.0e-12;
  
  /** Precision constant for updating sets */
  protected final static double m_Del = 1e-10; //1000 * Double.MIN_VALUE;
  
  /** error cache containing m_error[i] = SVMOutput(i) - m_target[i] - m_b <br/>
   * note, we don't need m_b in the cache, since if we do, we need to maintain 
   * it when m_b is updated */
  double[] m_error;
  
  /** alpha value for first candidate **/  
  protected double m_alpha1;
  
  /** alpha* value for first candidate **/  
  protected double m_alpha1Star;

  /** alpha value for second candidate **/  
  protected double m_alpha2;
  
  /** alpha* value for second candidate **/
  protected double m_alpha2Star;
  
  /**
   * default constructor
   */
  public RegSMO() {
    super();
  }
  
  /**
   * Returns a string describing classifier
   * 
   * @return 		a description suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String globalInfo() {
    return
        "Implementation of SMO for support vector regression as described "
      + "in :\n\n"
      + getTechnicalInformation().toString();
  }

  /**
   * Returns an instance of a TechnicalInformation object, containing 
   * detailed information about the technical background of this class,
   * e.g., paper reference or book this class is based on.
   * 
   * @return 		the technical information about this class
   */
  public TechnicalInformation getTechnicalInformation() {
    TechnicalInformation 	result;
    
    result = new TechnicalInformation(Type.MISC);
    result.setValue(Field.AUTHOR, "A.J. Smola and B. Schoelkopf");
    result.setValue(Field.TITLE, "A tutorial on support vector regression");
    result.setValue(Field.NOTE, "NeuroCOLT2 Technical Report NC2-TR-1998-030");
    result.setValue(Field.YEAR, "1998");
    
    return result;
  }
  
  /**
   * Returns an enumeration describing the available options
   * 
   * @return an enumeration of all the available options
   */
  public Enumeration listOptions() {
    Vector result = new Vector();
    
    result.addElement(new Option(
	"\tThe epsilon for round-off error.\n" 
	+ "\t(default 1.0e-12)", 
	"P", 1, "-P <double>"));
    
    Enumeration enm = super.listOptions();
    while (enm.hasMoreElements()) {
      result.addElement(enm.nextElement());
    }
    
    return result.elements();
  }
  
  /**
   * Parses a given list of options. <p/>
   *
   <!-- options-start -->
   * Valid options are: <p/>
   * 
   * <pre> -P &lt;double&gt;
   *  The epsilon for round-off error.
   *  (default 1.0e-12)</pre>
   * 
   * <pre> -L &lt;double&gt;
   *  The epsilon parameter in epsilon-insensitive loss function.
   *  (default 1.0e-3)</pre>
   * 
   * <pre> -W &lt;double&gt;
   *  The random number seed.
   *  (default 1)</pre>
   * 
   <!-- options-end -->
   *
   * @param options the list of options as an array of strings
   * @throws Exception if an option is not supported 
   */
  public void setOptions(String[] options) throws Exception {
    String	tmpStr;
    
    tmpStr = Utils.getOption('P', options);
    if (tmpStr.length() != 0) {
      setEpsilon(Double.parseDouble(tmpStr));
    } else {
      setEpsilon(1.0e-12);
    }
    
    super.setOptions(options);
  }
  
  /**
   * Gets the current settings of the classifier.
   *
   * @return an array of strings suitable for passing to setOptions
   */
  public String[] getOptions() {
    int       	i;
    Vector    	result;
    String[]  	options;

    result = new Vector();

    options = super.getOptions();
    for (i = 0; i < options.length; i++)
      result.add(options[i]);
    
    result.add("-P");
    result.add("" + getEpsilon());

    return (String[]) result.toArray(new String[result.size()]);	  
  }
  
  /**
   * Returns the tip text for this property
   * 
   * @return 		tip text for this property suitable for
   * 			displaying in the explorer/experimenter gui
   */
  public String epsilonTipText() {
    return "The epsilon for round-off error (shouldn't be changed).";
  }
  
  /**
   * Get the value of epsilon.
   * 
   * @return 		Value of epsilon.
   */
  public double getEpsilon() {
    return m_eps;
  }
  
  /**
   * Set the value of epsilon.
   * 
   * @param v  		Value to assign to epsilon.
   */
  public void setEpsilon(double v) {
    m_eps = v;
  }
  
  /** initialize various variables before starting the actual optimizer 
   * 
   * @param data 	data set used for learning
   * @throws Exception	if something goes wrong
   */
  protected void init(Instances data) throws Exception {
    super.init(data);
    
    //init error cache
    m_error = new double[m_nInstances];
    for (int i = 0; i < m_nInstances; i++) {
      m_error[i] = -m_target[i];
    }
  }
  
  /** 
   * wrap up various variables to save memeory and do some housekeeping after optimization
   * has finished.
   *
   * @throws Exception 	if something goes wrong
   */
  protected void wrapUp() throws Exception {
    m_error = null;
    super.wrapUp();
  }

  /** 
   * Finds optimal point on line constrained by first (i1) and second (i2) 
   * candidate. Parameters correspond to pseudocode (see technicalinformation)
   * 
   * @param i1
   * @param alpha1
   * @param alpha1Star
   * @param C1
   * @param i2
   * @param alpha2
   * @param alpha2Star
   * @param C2
   * @param gamma
   * @param eta
   * @param deltaPhi
   * @return
   */
  protected boolean findOptimalPointOnLine(int i1, double alpha1, double alpha1Star, double C1, 
      int i2, double alpha2, double alpha2Star, double C2, 
      double gamma, double eta, double deltaPhi) {
    if (eta <= 0) {
      // this may happen due to numeric instability
      // due to Mercer's condition, this should not happen, hence we give up
      return false;
    }
    
    boolean case1 = false;
    boolean case2 = false;
    boolean case3 = false;
    boolean case4 = false;
    boolean finished = false;
    
    //		while !finished 
    //		% this loop is passed at most three times 
    //		% case variables needed to avoid attempting small changes twice 
    while (!finished) {
      //			if (case1 == 0) && 
      //				(alpha1 > 0 || (alpha1* == 0 && deltaPhi > 0)) && 
      //				(alpha2 > 0 || (alpha2* == 0 && deltaPhi < 0)) 
      //				compute L, H (wrt. alpha1, alpha2) 
      //				if L < H 
      //					a2 = alpha2 ? - deltaPhi/eta 
      //					a2 = min(a2, H) 
      //					a2 = max(L, a2) 
      //					a1 = alpha1 ? - (a2 ? alpha2) 
      //					update alpha1, alpha2 if change is larger than some eps 
      //				else 
      //					finished = 1 
      //				endif 
      //				case1 = 1; 
      
      if ((case1 == false) && 
	  (alpha1 > 0 || (alpha1Star == 0 && deltaPhi > 0)) && 
	  (alpha2 > 0 || (alpha2Star == 0 && deltaPhi < 0))) {
	// compute L, H (wrt. alpha1, alpha2) 
	double L = Math.max(0, gamma - C1);
	double H = Math.min(C2, gamma);
	if (L < H) {
	  double a2 = alpha2 - deltaPhi / eta;
	  a2 = Math.min(a2, H);
	  a2 = Math.max(L, a2);
	  // To prevent precision problems
	  if (a2 > C2 - m_Del * C2) {
	    a2 = C2;
	  } else if (a2 <= m_Del * C2) {
	    a2 = 0;
	  }
	  double a1 = alpha1 - (a2 - alpha2);
	  if (a1 > C1 - m_Del * C1) {
	    a1 = C1;
	  } else if (a1 <= m_Del * C1) {
	    a1 = 0;
	  }
	  // update alpha1, alpha2 if change is larger than some eps
	  if (Math.abs(alpha1 - a1) > m_eps) {
	    deltaPhi += eta * (a2 - alpha2);
	    alpha1 = a1;
	    alpha2 = a2;
	  }
	} else {
	  finished = true;
	}
	case1 = true;
      }
      
      //			elseif (case2 == 0) && 
      //				(alpha1 > 0 || (alpha1* == 0 && deltaPhi > 2 epsilon)) && 
      //				(alpha2* > 0 || (alpha2 == 0 && deltaPhi > 2 epsilon)) 
      //				compute L, H (wrt. alpha1, alpha2*) 
      //				if L < H 
      //					a2 = alpha2* + (deltaPhi ?- 2 epsilon)/eta 
      //					a2 = min(a2, H) 
      //					a2 = max(L, a2) 
      //					a1 = alpha1 + (a2 ? alpha2*) 
      //					update alpha1, alpha2* if change is larger than some eps 
      //				else 
      //					finished = 1 
      //				endif 
      //				case2 = 1; 
      
      else if (
	  (case2 == false)
	  && (alpha1 > 0 || (alpha1Star == 0 && deltaPhi > 2 * m_epsilon))
	  && (alpha2Star > 0 || (alpha2 == 0 && deltaPhi > 2 * m_epsilon))) {
	// compute L, H (wrt. alpha1, alpha2*) 
	double L = Math.max(0, -gamma);
	double H = Math.min(C2, -gamma + C1);
	if (L < H) {
	  double a2 = alpha2Star + (deltaPhi - 2 * m_epsilon) / eta;
	  a2 = Math.min(a2, H);
	  a2 = Math.max(L, a2);
	  // To prevent precision problems
	  if (a2 > C2 - m_Del * C2) {
	    a2 = C2;
	  } else if (a2 <= m_Del * C2) {
	    a2 = 0;
	  }
	  double a1 = alpha1 + (a2 - alpha2Star);
	  if (a1 > C1 - m_Del * C1) {
	    a1 = C1;
	  } else if (a1 <= m_Del * C1) {
	    a1 = 0;
	  }
	  // update alpha1, alpha2* if change is larger than some eps 
	  if (Math.abs(alpha1 - a1) > m_eps) {
	    deltaPhi += eta * (-a2 + alpha2Star);
	    alpha1 = a1;
	    alpha2Star = a2;
	  }
	} else {
	  finished = true;
	}
	case2 = true;
      }
      
      //			elseif (case3 == 0) && 
      //				(alpha1* > 0 || (alpha1 == 0 && deltaPhi < -2 epsilon)) && 
      //				(alpha2 > 0 || (alpha2* == 0 && deltaPhi < -2 epsilon)) 
      //				compute L, H (wrt. alpha1*, alpha2) 
      //				if L < H 
      //					a2 = alpha2 ?- (deltaPhi ?+ 2 epsilon)/eta 
      //					a2 = min(a2, H) 
      //					a2 = max(L, a2) 
      //					a1 = alpha1* + (a2 ? alpha2) 
      //					update alpha1*, alpha2 if change is larger than some eps 
      //				else 
      //					finished = 1 
      //				endif 
      //				case3 = 1; 
      
      else if (
	  (case3 == false)
	  && (alpha1Star > 0 || (alpha1 == 0 && deltaPhi < - 2 * m_epsilon))
	  && (alpha2 > 0 || (alpha2Star == 0 && deltaPhi < - 2 * m_epsilon))) {
	// compute L, H (wrt. alpha1*, alpha2)
	double L = Math.max(0, gamma);
	double H = Math.min(C2, C1 + gamma);
	if (L < H) {
	  // note Smola's psuedocode has a minus, where there should be a plus in the following line, Keerthi's is correct
	  double a2 = alpha2 - (deltaPhi + 2 * m_epsilon) / eta;
	  a2 = Math.min(a2, H);
	  a2 = Math.max(L, a2);
	  // To prevent precision problems
	  if (a2 > C2 - m_Del * C2) {
	    a2 = C2;
	  } else if (a2 <= m_Del * C2) {
	    a2 = 0;
	  }
	  double a1 = alpha1Star + (a2 - alpha2);
	  if (a1 > C1 - m_Del * C1) {
	    a1 = C1;
	  } else if (a1 <= m_Del * C1) {
	    a1 = 0;
	  }
	  // update alpha1*, alpha2 if change is larger than some eps 
	  if (Math.abs(alpha1Star - a1) > m_eps) {
	    deltaPhi += eta * (a2 - alpha2);
	    alpha1Star = a1;
	    alpha2 = a2;
	  }
	} else {
	  finished = true;
	}
	case3 = true;
      }
      
      //			elseif (case4 == 0) && 
      //				(alpha1* > 0 || (alpha1 == 0 && deltaPhi < 0)) && 
      //				(alpha2* > 0 || (alpha2 == 0 && deltaPhi > 0)) 
      //				compute L, H (wrt. alpha1*, alpha2*) 
      //				if L < H 
      //					a2 = alpha2* + deltaPhi/eta 
      //					a2 = min(a2, H) 
      //					a2 = max(L, a2) 
      //					a1 = alpha1* ? (a2 ? alpha2*) 
      //					update alpha1*, alpha2* if change is larger than some eps 
      //				else 
      //					finished = 1 
      //				endif 
      //				case4 = 1; 
      //			else 
      //				finished = 1 
      //			endif 
      
      else if ((case4 == false) && 
	  (alpha1Star > 0 || (alpha1 == 0 && deltaPhi < 0)) && 
	  (alpha2Star > 0 || (alpha2 == 0 && deltaPhi > 0))) {
	// compute L, H (wrt. alpha1*, alpha2*) 
	double L = Math.max(0, -gamma - C1);
	double H = Math.min(C2, -gamma);
	if (L < H) {
	  double a2 = alpha2Star + deltaPhi / eta;
	  a2 = Math.min(a2, H);
	  a2 = Math.max(L, a2);
	  // To prevent precision problems
	  if (a2 > C2 - m_Del * C2) {
	    a2 = C2;
	  } else if (a2 <= m_Del * C2) {
	    a2 = 0;
	  }
	  double a1 = alpha1Star - (a2 - alpha2Star);
	  if (a1 > C1 - m_Del * C1) {
	    a1 = C1;
	  } else if (a1 <= m_Del * C1) {
	    a1 = 0;
	  }
	  // update alpha1*, alpha2* if change is larger than some eps 
	  if (Math.abs(alpha1Star - a1) > m_eps) {
	    deltaPhi += eta * (-a2 + alpha2Star);
	    
	    alpha1Star = a1;
	    alpha2Star = a2;
	  }
	} else {
	  finished = true;
	}
	case4 = true;
      } else {
	finished = true;
      }
      
      //			update deltaPhi
      // using 4.36 from Smola's thesis:
      // deltaPhi = deltaPhi - eta * ((alpha1New-alpha1StarNew)-(alpha1-alpha1Star));
      // the update is done inside the loop, saving us to remember old values of alpha1(*)
      //deltaPhi += eta * ((alpha2 - alpha2Star) - dAlpha2Old);
      //dAlpha2Old = (alpha2 - alpha2Star);
      
      //		endwhile 
      
    }
    
    if (Math.abs(alpha1 - m_alpha[i1]) > m_eps
	|| Math.abs(alpha1Star - m_alphaStar[i1]) > m_eps
	|| Math.abs(alpha2 - m_alpha[i2]) > m_eps
	|| Math.abs(alpha2Star - m_alphaStar[i2]) > m_eps) {
      
      if (alpha1 > C1 - m_Del * C1) {
	alpha1 = C1;
      } else if (alpha1 <= m_Del * C1) {
	alpha1 = 0;
      }
      if (alpha1Star > C1 - m_Del * C1) {
	alpha1Star = C1;
      } else if (alpha1Star <= m_Del * C1) {
	alpha1Star = 0;
      }
      if (alpha2 > C2 - m_Del * C2) {
	alpha2 = C2;
      } else if (alpha2 <= m_Del * C2) {
	alpha2 = 0;
      }
      if (alpha2Star > C2 - m_Del * C2) {
	alpha2Star = C2;
      } else if (alpha2Star <= m_Del * C2) {
	alpha2Star = 0;
      }
      
      // store new alpha's
      m_alpha[i1] = alpha1;
      m_alphaStar[i1] = alpha1Star;
      m_alpha[i2] = alpha2;
      m_alphaStar[i2] = alpha2Star;
      
      // update supportvector set
      if (alpha1 != 0 || alpha1Star != 0){
	if (!m_supportVectors.contains(i1)) {
	  m_supportVectors.insert(i1);	
	}
      } else {
	m_supportVectors.delete(i1);
      }
      if (alpha2 != 0 || alpha2Star != 0){
	if (!m_supportVectors.contains(i2)) {
	  m_supportVectors.insert(i2);
	}
      } else {
	m_supportVectors.delete(i2);
      }
      return true;
    }
    
    return false;
  }

  /** 
   * takeStep method from pseudocode.
   * Parameters correspond to pseudocode (see technicalinformation)
   * 
   * @param i1
   * @param i2
   * @param alpha2
   * @param alpha2Star
   * @param phi2
   * @return
   * @throws Exception
   */
  protected int takeStep(int i1, int i2, double alpha2, double alpha2Star, double phi2) throws Exception {
    //		if (i1 == i2) return 0 
    if (i1 == i2) {
      return 0;
    }
    double C1 = m_C * m_data.instance(i1).weight();
    double C2 = m_C * m_data.instance(i2).weight();
    //		alpha1, alpha1* = Lagrange multipliers for i1 
    //		y1 = target[i1] 
    //		phi1 = SVM output on point[i1] ? y1 (in error cache) 
    double alpha1 = m_alpha[i1];
    double alpha1Star = m_alphaStar[i1];
    double y1 = m_target[i1];
    double phi1 = m_error[i1];
    
    //		k11 = kernel(point[i1],point[i1]) 
    //		k12 = kernel(point[i1],point[i2]) 
    //		k22 = kernel(point[i2],point[i2]) 
    //		eta = 2*k12? - k11? - k22 
    //		gamma = alpha1 ?- alpha1* + alpha2 ?- alpha2* 
    
    double k11 = m_kernel.eval(i1, i1, m_data.instance(i1));
    double k12 = m_kernel.eval(i1, i2, m_data.instance(i1));
    double k22 = m_kernel.eval(i2, i2, m_data.instance(i2));
    double eta = -2 * k12 + k11 + k22; // note, Smola's psuedocode has signs swapped, Keerthi's doesn't
    if (eta < 0) {
      // this may happen due to numeric instability
      // due to Mercer's condition, this should not happen, hence we give up
      return 0;
    }
    double gamma = alpha1 - alpha1Star + alpha2 - alpha2Star;
    
    //		% we assume eta < 0. otherwise one has to repeat the complete 
    //		% reasoning similarly (compute objective function for L and H 
    //		% and decide which one is largest 
    //		case1 = case2 = case3 = case4 = finished = 0 
    //		alpha1old = alpha1, alpha1old* = alpha1* 
    //		alpha2old = alpha2, alpha2old* = alpha2* 
    //		deltaPhi = phi1 ?- phi2 
    
    double alpha1old = alpha1;
    double alpha1Starold = alpha1Star;
    double alpha2old = alpha2;
    double alpha2Starold = alpha2Star;
    double deltaPhi = phi2 - phi1;
    
    if (findOptimalPointOnLine(i1, alpha1, alpha1Star, C1, i2, alpha2, alpha2Star, C2, gamma, eta, deltaPhi)) {
      alpha1 = m_alpha[i1];
      alpha1Star = m_alphaStar[i1];
      alpha2 = m_alpha[i2];
      alpha2Star = m_alphaStar[i2];
      
      
      //		Update error cache using new Lagrange multipliers 
      double dAlpha1 = alpha1 - alpha1old - (alpha1Star - alpha1Starold);
      double dAlpha2 = alpha2 - alpha2old - (alpha2Star - alpha2Starold);
      for (int j = 0; j < m_nInstances; j++) {
	if ((j != i1) && (j != i2)/* && m_error[j] != MAXERR*/) {
	  m_error[j] += dAlpha1 * m_kernel.eval(i1, j, m_data.instance(i1)) + dAlpha2 * m_kernel.eval(i2, j, m_data.instance(i2));
	}
      }
      m_error[i1] += dAlpha1 * k11 + dAlpha2 * k12;
      m_error[i2] += dAlpha1 * k12 + dAlpha2 * k22;
      
      //		Update threshold to reflect change in Lagrange multipliers
      double b1 = Double.MAX_VALUE;
      double b2 = Double.MAX_VALUE;
      if ((0 < alpha1 && alpha1 < C1) || (0 < alpha1Star && alpha1Star < C1) ||(0 < alpha2 && alpha2 < C2) || (0 < alpha2Star && alpha2Star < C2)) {
	if (0 < alpha1 && alpha1 < C1) {
	  b1 = m_error[i1] - m_epsilon;
	} else if (0 < alpha1Star && alpha1Star < C1) {
	  b1 = m_error[i1] + m_epsilon;
	}
	if (0 < alpha2 && alpha2 < C2) {
	  b2 = m_error[i2] - m_epsilon;
	} else if (0 < alpha2Star && alpha2Star < C2) {
	  b2 = m_error[i2] + m_epsilon;
	}
	if (b1 < Double.MAX_VALUE) {
	  m_b = b1;
	  if (b2 < Double.MAX_VALUE) {
	    m_b = (b1 + b2) / 2.0;
	  }
	} else if (b2 < Double.MAX_VALUE) {
	  m_b = b2;
	}
      } else if (m_b == 0) {
	// both alpha's are on the boundary, and m_b is not initialized
	m_b = (m_error[i1] + m_error[i2])/2.0;
      }
      
      //		if changes in alpha1(*), alpha2(*) are larger than some eps 
      //			return 1 
      //		else 
      //			return 0 
      //		endif
      return 1;
    } else {
      return 0;
    }
    //	endprocedure 
  }
  
  /** 
   * examineExample method from pseudocode.
   * Parameters correspond to pseudocode (see technicalinformation)
   * 
   * @param i2
   * @return
   * @throws Exception
   */
  protected int examineExample(int i2) throws Exception {
    //	procedure examineExample(i2) 
    //		y2 = target[i2] 
    double y2 = m_target[i2];
    //		alpha2, alpha2* = Lagrange multipliers for i2 
    double alpha2 = m_alpha[i2];
    double alpha2Star = m_alphaStar[i2];
    //		C2, C2* = Constraints for i2 
    double C2 = m_C;
    double C2Star = m_C;
    //		phi2 = SVM output on point[i2] ? y2 (in error cache) 
    double phi2 = m_error[i2];
    // phi2b contains the error, taking the offset in account
    double phi2b = phi2 - m_b;
    //		if ((phi2 > epsilon && alpha2* < C2*) ||
    //			(phi2 < epsilon && alpha2* > 0 ) ||
    //			(-?phi2 > epsilon && alpha2 < C2 ) ||
    //			(?-phi2 > epsilon && alpha2 > 0 )) 
    if ((phi2b > m_epsilon && alpha2Star < C2Star)
	|| (phi2b < m_epsilon && alpha2Star > 0)
	|| (-phi2b > m_epsilon && alpha2 < C2)
	|| (-phi2b > m_epsilon && alpha2 > 0)) {
      
      //			if (number of non?zero & non?C alpha > 1) 
      //				i1 = result of second choice heuristic 
      //				if takeStep(i1,i2) return 1 
      //			endif 
      int i1 = secondChoiceHeuristic(i2);
      if (i1 >= 0 && (takeStep(i1, i2, alpha2, alpha2Star, phi2) > 0)) {
	return 1;
      }
      //			loop over all non?zero and non?C alpha, random start 
      //				i1 = identity of current alpha 
      //				if takeStep(i1,i2) return 1 
      //			endloop 
      for (i1 = 0; i1 < m_target.length; i1++) {
	if ((m_alpha[i1] > 0 && m_alpha[i1] < m_C) || (m_alphaStar[i1] > 0 && m_alphaStar[i1] < m_C)) {
	  if (takeStep(i1, i2, alpha2, alpha2Star, phi2) > 0) {
	    return 1;
	  }
	}
      }
      //			loop over all possible i1, with random start 
      //				i1 = loop variable 
      //				if takeStep(i1,i2) return 1 
      //			endloop 
      for (i1 = 0; i1 < m_target.length; i1++) {
	if (takeStep(i1, i2, alpha2, alpha2Star, phi2) > 0) {
	  return 1;
	}
      }
      //		endif 
    }
    //		return 0 
    return 0;
    //	endprocedure 
  }
  
  /** 
   * applies heuristic for finding candidate that is expected to lead to
   * good gain when applying takeStep together with second candidate.
   * 
   * @param i2 index of second candidate
   * @return
   */
  protected int secondChoiceHeuristic(int i2) {
    // randomly select an index i1 (not equal to i2) with non?zero and non?C alpha, if any
    for (int i = 0; i < 59; i++) {
      int i1 = m_random.nextInt(m_nInstances);
      if ((i1 != i2) && (m_alpha[i1] > 0 && m_alpha[i1] < m_C) || (m_alphaStar[i1] > 0 && m_alphaStar[i1] < m_C)) {
	return i1;
      }
    }
    return -1;
  }
  
  /**
   * finds alpha and alpha* parameters that optimize the SVM target function
   * 
   * @throws Exception
   */
  public void optimize() throws Exception {
    
    //	main routine: 
    //		initialize threshold to zero 
    //		numChanged = 0 
    //		examineAll = 1 
    //		SigFig = -100 
    //		LoopCounter = 0 
    int numChanged = 0;
    int examineAll = 1;
    int sigFig = -100;
    int loopCounter = 0;
    //		while ((numChanged > 0 | examineAll) | (SigFig < 3)) 
    while ((numChanged > 0 || (examineAll > 0)) | (sigFig < 3)) {
      //			LoopCounter++ 
      //			numChanged = 0; 
      loopCounter++;
      numChanged = 0;
      //			if (examineAll) 
      //				loop I over all training examples 
      //				numChanged += examineExample(I)
      //			else 
      //				loop I over examples where alpha is not 0 & not C 
      //				numChanged += examineExample(I) 
      //			endif
      int numSamples = 0;
      if (examineAll > 0) {
	for (int i = 0; i < m_nInstances; i++) {
	  numChanged += examineExample(i);
	}
      } else {
	for (int i = 0; i < m_target.length; i++) {
	  if ((m_alpha[i] > 0 && m_alpha[i] < m_C * m_data.instance(i).weight()) || 
	      (m_alphaStar[i] > 0 && m_alphaStar[i] < m_C * m_data.instance(i).weight())) {
	    numSamples++;
	    numChanged += examineExample(i);
	  }
	}
      }
      //	
      //		if (mod(LoopCounter, 2) == 0) 
      //				MinimumNumChanged = max(1, 0.1*NumSamples) 
      //			else 
      //				MinimumNumChanged = 1 
      //			endif 
      int minimumNumChanged = 1;
      if (loopCounter % 2 == 0) {
	minimumNumChanged = (int) Math.max(1, 0.1 * numSamples);
      }
      
      //			if (examineAll == 1) 
      //				examineAll = 0 
      //			elseif (numChanged < MinimumNumChanged) 
      //				examineAll = 1 
      //			endif 
      if (examineAll == 1) {
	examineAll = 0;
      } else if (numChanged < minimumNumChanged) {
	examineAll = 1;
      }
      
      //		endwhile 
      if (loopCounter == 2500) {
	break;
      }
    }
    //	endmain 
  }
  
  /** 
   * learn SVM parameters from data using Smola's SMO algorithm.
   * Subclasses should implement something more interesting.
   * 
   * @param instances	the data to learn from
   * @throws Exception	if something goes wrong
   */
  public void buildClassifier(Instances instances) throws Exception {
    // initialize variables
    init(instances);
    // solve optimization problem
    optimize();
    // clean up
    wrapUp();
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

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
 *    RELEASE INFORMATION (December 27, 2004)
 *    
 *    FCBF algorithm:
 *      Template obtained from Weka
 *      Developped for Weka by Zheng Alan Zhao   
 *      December 27, 2004
 *
 *    FCBF algorithm is a feature selection method based on Symmetrical Uncertainty 
 *    Measurement for relevance redundancy analysis. The details of FCBF algorithm are 
 *    in L. Yu and H. Liu. Feature selection for high-dimensional data: a fast 
 *    correlation-based filter solution. In Proceedings of the twentieth International 
 *    Conference on Machine Learning, pages 856--863, 2003.
 *    
 *    
 *    CONTACT INFORMATION
 *    
 *    For algorithm implementation:
 *    Zheng Zhao: zhaozheng at asu.edu
 *      
 *    For the algorithm:
 *    Lei Yu: leiyu at asu.edu
 *    Huan Liu: hliu at asu.edu
 *     
 *    Data Mining and Machine Learning Lab
 *    Computer Science and Engineering Department
 *    Fulton School of Engineering
 *    Arizona State University
 *    Tempe, AZ 85287
 *
 *    AttributeSetEvaluator.java
 *
 *    Copyright (C) 2004 Data Mining and Machine Learning Lab, 
 *                       Computer Science and Engineering Department, 
 *       		 Fulton School of Engineering, 
 *                       Arizona State University
 *
 */

package weka.attributeSelection;


/**
 * Abstract attribute set evaluator.
 *
 * @author Zheng Zhao: zhaozheng at asu.edu
 * @version $Revision: 1.3 $
 */
public abstract class AttributeSetEvaluator extends ASEvaluation {
  
    /** for serialization */
    private static final long serialVersionUID = -5744881009422257389L;
  
    // ===============
    // Public methods.
    // ===============

    /**
     * evaluates an individual attribute
     *
     * @param attribute the index of the attribute to be evaluated
     * @return the "merit" of the attribute
     * @exception Exception if the attribute could not be evaluated
     */
    public abstract double evaluateAttribute(int attribute) throws Exception;

  /**
   * Evaluates a set of attributes
   *
   * @param attributes an <code>int[]</code> value
   * @param classAttributes an <code>int[]</code> value
   * @return a <code>double</code> value
   * @exception Exception if an error occurs
   */
  public abstract double evaluateAttribute(int[] attributes, int[] classAttributes) 
    throws Exception;
}

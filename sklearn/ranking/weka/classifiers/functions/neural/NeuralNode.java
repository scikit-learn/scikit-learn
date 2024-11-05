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
 *    NeuralNode.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.functions.neural;

import weka.core.RevisionUtils;

import java.util.Random;

/**
 * This class is used to represent a node in the neuralnet.
 * 
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 5402 $
 */
public class NeuralNode
  extends NeuralConnection {

  /** for serialization */
  private static final long serialVersionUID = -1085750607680839163L;
    
  /** The weights for each of the input connections, and the threshold. */
  private double[] m_weights;
  
  /** The best (lowest error) weights. Only used when validation set is used */
  private double[] m_bestWeights;
  
  /** The change in the weights. */
  private double[] m_changeInWeights;
  
  private Random m_random;

  /** Performs the operations for this node. Currently this
   * defines that the node is either a sigmoid or a linear unit. */
  private NeuralMethod m_methods;

  /** 
   * @param id The string name for this node (used to id this node).
   * @param r A random number generator used to generate initial weights.
   * @param m The methods this node should use to update.
   */
  public NeuralNode(String id, Random r, NeuralMethod m) {
    super(id);
    m_weights = new double[1];
    m_bestWeights = new double[1];
    m_changeInWeights = new double[1];
    
    m_random = r;
    
    m_weights[0] = m_random.nextDouble() * .1 - .05;
    m_changeInWeights[0] = 0;

    m_methods = m;
  }
  
  /**
   * Set how this node should operate (note that the neural method has no
   * internal state, so the same object can be used by any number of nodes.
   * @param m The new method.
   */
  public void setMethod(NeuralMethod m) {
    m_methods = m;
  } 

  public NeuralMethod getMethod() {
    return m_methods;
  }

  /**
   * Call this to get the output value of this unit. 
   * @param calculate True if the value should be calculated if it hasn't been
   * already.
   * @return The output value, or NaN, if the value has not been calculated.
   */
  public double outputValue(boolean calculate) {
    
    if (Double.isNaN(m_unitValue) && calculate) {
      //then calculate the output value;
      m_unitValue = m_methods.outputValue(this);
    }
    
    return m_unitValue;
  }

  
  /**
   * Call this to get the error value of this unit.
   * @param calculate True if the value should be calculated if it hasn't been
   * already.
   * @return The error value, or NaN, if the value has not been calculated.
   */
  public double errorValue(boolean calculate) {

    if (!Double.isNaN(m_unitValue) && Double.isNaN(m_unitError) && calculate) {
      //then calculate the error.
      m_unitError = m_methods.errorValue(this);
    }
    return m_unitError;
  }

  /**
   * Call this to reset the value and error for this unit, ready for the next
   * run. This will also call the reset function of all units that are 
   * connected as inputs to this one.
   * This is also the time that the update for the listeners will be performed.
   */
  public void reset() {
    
    if (!Double.isNaN(m_unitValue) || !Double.isNaN(m_unitError)) {
      m_unitValue = Double.NaN;
      m_unitError = Double.NaN;
      m_weightsUpdated = false;
      for (int noa = 0; noa < m_numInputs; noa++) {
	m_inputList[noa].reset();
      }
    }
  }
  
  /**
   * Call this to have the connection save the current
   * weights.
   */
  public void saveWeights() {
    // copy the current weights
    System.arraycopy(m_weights, 0, m_bestWeights, 0, m_weights.length);
    
    // tell inputs to save weights
    for (int i = 0; i < m_numInputs; i++) {
      m_inputList[i].saveWeights();
    }
  }
  
  /**
   * Call this to have the connection restore from the saved
   * weights.
   */
  public void restoreWeights() {
    // copy the saved best weights back into the weights
    System.arraycopy(m_bestWeights, 0, m_weights, 0, m_weights.length);
    
    // tell inputs to restore weights
    for (int i = 0; i < m_numInputs; i++) {
      m_inputList[i].restoreWeights();
    }
  }

  /**
   * Call this to get the weight value on a particular connection.
   * @param n The connection number to get the weight for, -1 if The threshold
   * weight should be returned.
   * @return The value for the specified connection or if -1 then it should 
   * return the threshold value. If no value exists for the specified 
   * connection, NaN will be returned.
   */
  public double weightValue(int n) {
    if (n >= m_numInputs || n < -1) {
      return Double.NaN;
    }
    return m_weights[n + 1];
  }

  /**
   * call this function to get the weights array.
   * This will also allow the weights to be updated.
   * @return The weights array.
   */
  public double[] getWeights() {
    return m_weights;
  }

  /**
   * call this function to get the chnage in weights array.
   * This will also allow the change in weights to be updated.
   * @return The change in weights array.
   */
  public double[] getChangeInWeights() {
    return m_changeInWeights;
  }

  /**
   * Call this function to update the weight values at this unit.
   * After the weights have been updated at this unit, All the
   * input connections will then be called from this to have their
   * weights updated.
   * @param l The learning rate to use.
   * @param m The momentum to use.
   */
  public void updateWeights(double l, double m) {
    
    if (!m_weightsUpdated && !Double.isNaN(m_unitError)) {
      m_methods.updateWeights(this, l, m);
     
      //note that the super call to update the inputs is done here and
      //not in the m_method updateWeights, because it is not deemed to be
      //required to update the weights at this node (while the error and output
      //value ao need to be recursively calculated)
      super.updateWeights(l, m); //to call all of the inputs.
    }
    
  }

  /**
   * This will connect the specified unit to be an input to this unit.
   * @param i The unit.
   * @param n It's connection number for this connection.
   * @return True if the connection was made, false otherwise.
   */
  protected boolean connectInput(NeuralConnection i, int n) {
    
    //the function that this overrides can do most of the work.
    if (!super.connectInput(i, n)) {
      return false;
    }
    
    //note that the weights are shifted 1 forward in the array so
    //it leaves the numinputs aligned on the space the weight needs to go.
    m_weights[m_numInputs] = m_random.nextDouble() * .1 - .05;
    m_changeInWeights[m_numInputs] = 0;
    
    return true;
  }

  /**
   * This will allocate more space for input connection information
   * if the arrays for this have been filled up.
   */
  protected void allocateInputs() {
    
    NeuralConnection[] temp1 = new NeuralConnection[m_inputList.length + 15];
    int[] temp2 = new int[m_inputNums.length + 15];
    double[] temp4 = new double[m_weights.length + 15];
    double[] temp5 = new double[m_changeInWeights.length + 15];
    double[] temp6 = new double[m_bestWeights.length + 15];

    temp4[0] = m_weights[0];
    temp5[0] = m_changeInWeights[0];
    temp6[0] = m_bestWeights[0];
    for (int noa = 0; noa < m_numInputs; noa++) {
      temp1[noa] = m_inputList[noa];
      temp2[noa] = m_inputNums[noa];
      temp4[noa+1] = m_weights[noa+1];
      temp5[noa+1] = m_changeInWeights[noa+1];
      temp6[noa+1] = m_bestWeights[noa+1];
    }
    
    m_inputList = temp1;
    m_inputNums = temp2;
    m_weights = temp4;
    m_changeInWeights = temp5;
    m_bestWeights = temp6;
  }

  
  

  /**
   * This will disconnect the input with the specific connection number
   * From this node (only on this end however).
   * @param i The unit to disconnect.
   * @param n The connection number at the other end, -1 if all the connections
   * to this unit should be severed (not the same as removeAllInputs).
   * @return True if the connection was removed, false if the connection was 
   * not found.
   */
  protected boolean disconnectInput(NeuralConnection i, int n) {
    
    int loc = -1;
    boolean removed = false;
    do {
      loc = -1;
      for (int noa = 0; noa < m_numInputs; noa++) {
	if (i == m_inputList[noa] && (n == -1 || n == m_inputNums[noa])) {
	  loc = noa;
	  break;
	}
      }
      
      if (loc >= 0) {
	for (int noa = loc+1; noa < m_numInputs; noa++) {
	  m_inputList[noa-1] = m_inputList[noa];
	  m_inputNums[noa-1] = m_inputNums[noa];
	  
	  m_weights[noa] = m_weights[noa+1];
	  m_changeInWeights[noa] = m_changeInWeights[noa+1];
	  
	  m_inputList[noa-1].changeOutputNum(m_inputNums[noa-1], noa-1);
	}
	m_numInputs--;
	removed = true;
      }      
    } while (n == -1 && loc != -1);
    return removed;
  }
  
  /**
   * This function will remove all the inputs to this unit.
   * In doing so it will also terminate the connections at the other end.
   */
  public void removeAllInputs() {
    super.removeAllInputs();
    
    double temp1 = m_weights[0];
    double temp2 = m_changeInWeights[0];

    m_weights = new double[1];
    m_changeInWeights = new double[1];

    m_weights[0] = temp1;
    m_changeInWeights[0] = temp2;
    
  }  
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5402 $");
  }
}

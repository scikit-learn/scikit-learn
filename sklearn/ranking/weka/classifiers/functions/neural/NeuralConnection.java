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
 *    NeuralConnection.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 */

package weka.classifiers.functions.neural;

import weka.core.RevisionHandler;

import java.awt.Color;
import java.awt.Graphics;
import java.io.Serializable;

/** 
 * Abstract unit in a NeuralNetwork.
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 5402 $
 */
public abstract class NeuralConnection
  implements Serializable, RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -286208828571059163L;

  //bitwise flags for the types of unit.

  /** This unit is not connected to any others. */
  public static final int UNCONNECTED = 0;
  
  /** This unit is a pure input unit. */
  public static final int PURE_INPUT = 1;
  
  /** This unit is a pure output unit. */
  public static final int PURE_OUTPUT = 2;
  
  /** This unit is an input unit. */
  public static final int INPUT = 4;
  
  /** This unit is an output unit. */
  public static final int OUTPUT = 8;
  
  /** This flag is set once the unit has a connection. */
  public static final int CONNECTED = 16;



  /////The difference between pure and not is that pure is used to feed 
  /////the neural network the attribute values and the errors on the outputs
  /////Beyond that they do no calculations, and have certain restrictions
  /////on the connections they can make.



  /** The list of inputs to this unit. */
  protected NeuralConnection[] m_inputList;

  /** The list of outputs from this unit. */
  protected NeuralConnection[] m_outputList;

  /** The numbering for the connections at the other end of the input lines. */
  protected int[] m_inputNums;
  
  /** The numbering for the connections at the other end of the out lines. */
  protected int[] m_outputNums;

  /** The number of inputs. */
  protected int m_numInputs;

  /** The number of outputs. */
  protected int m_numOutputs;

  /** The output value for this unit, NaN if not calculated. */
  protected double m_unitValue;

  /** The error value for this unit, NaN if not calculated. */
  protected double m_unitError;
  
  /** True if the weights have already been updated. */
  protected boolean m_weightsUpdated;
  
  /** The string that uniquely (provided naming is done properly) identifies
   * this unit. */
  protected String m_id;

  /** The type of unit this is. */
  protected int m_type;

  /** The x coord of this unit purely for displaying purposes. */
  protected double m_x;
  
  /** The y coord of this unit purely for displaying purposes. */
  protected double m_y;
  

  
  
  /**
   * Constructs The unit with the basic connection information prepared for
   * use. 
   * 
   * @param id the unique id of the unit
   */
  public NeuralConnection(String id) {
    
    m_id = id;
    m_inputList = new NeuralConnection[0];
    m_outputList = new NeuralConnection[0];
    m_inputNums = new int[0];
    m_outputNums = new int[0];

    m_numInputs = 0;
    m_numOutputs = 0;

    m_unitValue = Double.NaN;
    m_unitError = Double.NaN;

    m_weightsUpdated = false;
    m_x = 0;
    m_y = 0;
    m_type = UNCONNECTED;
  }
  
  
  /**
   * @return The identity string of this unit.
   */
  public String getId() {
    return m_id;
  }

  /**
   * @return The type of this unit.
   */
  public int getType() {
    return m_type;
  }

  /**
   * @param t The new type of this unit.
   */
  public void setType(int t) {
    m_type = t;
  }

  /**
   * Call this to reset the unit for another run.
   * It is expected by that this unit will call the reset functions of all 
   * input units to it. It is also expected that this will not be done
   * if the unit has already been reset (or atleast appears to be).
   */
  public abstract void reset();

  /**
   * Call this to get the output value of this unit. 
   * @param calculate True if the value should be calculated if it hasn't been
   * already.
   * @return The output value, or NaN, if the value has not been calculated.
   */
  public abstract double outputValue(boolean calculate);

  /**
   * Call this to get the error value of this unit.
   * @param calculate True if the value should be calculated if it hasn't been
   * already.
   * @return The error value, or NaN, if the value has not been calculated.
   */
  public abstract double errorValue(boolean calculate);
  
  /**
   * Call this to have the connection save the current
   * weights.
   */
  public abstract void saveWeights();
  
  /**
   * Call this to have the connection restore from the saved
   * weights.
   */
  public abstract void restoreWeights();

  /**
   * Call this to get the weight value on a particular connection.
   * @param n The connection number to get the weight for, -1 if The threshold
   * weight should be returned.
   * @return This function will default to return 1. If overridden, it should
   * return the value for the specified connection or if -1 then it should 
   * return the threshold value. If no value exists for the specified 
   * connection, NaN will be returned.
   */
  public double weightValue(int n) {
    return 1;
  }

  /**
   * Call this function to update the weight values at this unit.
   * After the weights have been updated at this unit, All the
   * input connections will then be called from this to have their
   * weights updated.
   * @param l The learning Rate to use.
   * @param m The momentum to use.
   */
  public void updateWeights(double l, double m) {
    
    //the action the subclasses should perform is upto them 
    //but if they coverride they should make a call to this to
    //call the method for all their inputs.
    
    if (!m_weightsUpdated) {
      for (int noa = 0; noa < m_numInputs; noa++) {
	m_inputList[noa].updateWeights(l, m);
      }
      m_weightsUpdated = true;
    }
    
  }

  /**
   * Use this to get easy access to the inputs.
   * It is not advised to change the entries in this list
   * (use the connecting and disconnecting functions to do that)
   * @return The inputs list.
   */
  public NeuralConnection[] getInputs() {
    return m_inputList;
  }

  /**
   * Use this to get easy access to the outputs.
   * It is not advised to change the entries in this list
   * (use the connecting and disconnecting functions to do that)
   * @return The outputs list.
   */
  public NeuralConnection[] getOutputs() {
    return m_outputList;
  }

  /**
   * Use this to get easy access to the input numbers.
   * It is not advised to change the entries in this list
   * (use the connecting and disconnecting functions to do that)
   * @return The input nums list.
   */
  public int[] getInputNums() {
    return m_inputNums;
  }

  /**
   * Use this to get easy access to the output numbers.
   * It is not advised to change the entries in this list
   * (use the connecting and disconnecting functions to do that)
   * @return The outputs list.
   */
  public int[] getOutputNums() {
    return m_outputNums;
  }

  /**
   * @return the x coord.
   */
  public double getX() {
    return m_x;
  }
  
  /**
   * @return the y coord.
   */
  public double getY() {
    return m_y;
  }
  
  /**
   * @param x The new value for it's x pos.
   */
  public void setX(double x) {
    m_x = x;
  }
  
  /**
   * @param y The new value for it's y pos.
   */
  public void setY(double y) {
    m_y = y;
  }
  
  
  /**
   * Call this function to determine if the point at x,y is on the unit.
   * @param g The graphics context for font size info.
   * @param x The x coord.
   * @param y The y coord.
   * @param w The width of the display.
   * @param h The height of the display.
   * @return True if the point is on the unit, false otherwise.
   */
  public boolean onUnit(Graphics g, int x, int y, int w, int h) {

    int m = (int)(m_x * w);
    int c = (int)(m_y * h);
    if (x > m + 10 || x < m - 10 || y > c + 10 || y < c - 10) {
      return false;
    }
    return true;

  }
  
  /**
   * Call this function to draw the node.
   * @param g The graphics context.
   * @param w The width of the drawing area.
   * @param h The height of the drawing area.
   */
  public void drawNode(Graphics g, int w, int h) {
    
    if ((m_type & OUTPUT) == OUTPUT) {
      g.setColor(Color.orange);
    }
    else {
      g.setColor(Color.red);
    }
    g.fillOval((int)(m_x * w) - 9, (int)(m_y * h) - 9, 19, 19);
    g.setColor(Color.gray);
    g.fillOval((int)(m_x * w) - 5, (int)(m_y * h) - 5, 11, 11);
  }

  /**
   * Call this function to draw the node highlighted.
   * @param g The graphics context.
   * @param w The width of the drawing area.
   * @param h The height of the drawing area.
   */
  public void drawHighlight(Graphics g, int w, int h) {
   
    drawNode(g, w, h);
    g.setColor(Color.yellow);
    g.fillOval((int)(m_x * w) - 5, (int)(m_y * h) - 5, 11, 11);
  }

  /** 
   * Call this function to draw the nodes input connections.
   * @param g The graphics context.
   * @param w The width of the drawing area.
   * @param h The height of the drawing area.
   */
  public void drawInputLines(Graphics g, int w, int h) {

    g.setColor(Color.black);
    
    int px = (int)(m_x * w);
    int py = (int)(m_y * h);
    for (int noa = 0; noa < m_numInputs; noa++) {
      g.drawLine((int)(m_inputList[noa].getX() * w)
		 , (int)(m_inputList[noa].getY() * h)
		 , px, py);
    }
  }

  /**
   * Call this function to draw the nodes output connections.
   * @param g The graphics context.
   * @param w The width of the drawing area.
   * @param h The height of the drawing area.
   */
  public void drawOutputLines(Graphics g, int w, int h) {
    
    g.setColor(Color.black);
    
    int px = (int)(m_x * w);
    int py = (int)(m_y * h);
    for (int noa = 0; noa < m_numOutputs; noa++) {
      g.drawLine(px, py
		 , (int)(m_outputList[noa].getX() * w)
		 , (int)(m_outputList[noa].getY() * h));
    }
  }


  /**
   * This will connect the specified unit to be an input to this unit.
   * @param i The unit.
   * @param n It's connection number for this connection.
   * @return True if the connection was made, false otherwise.
   */
  protected boolean connectInput(NeuralConnection i, int n) {
    
    for (int noa = 0; noa < m_numInputs; noa++) {
      if (i == m_inputList[noa]) {
	return false;
      }
    }
    if (m_numInputs >= m_inputList.length) {
      //then allocate more space to it.
      allocateInputs();
    }
    m_inputList[m_numInputs] = i;
    m_inputNums[m_numInputs] = n;
    m_numInputs++;
    return true;
  }
  
  /**
   * This will allocate more space for input connection information
   * if the arrays for this have been filled up.
   */
  protected void allocateInputs() {
    
    NeuralConnection[] temp1 = new NeuralConnection[m_inputList.length + 15];
    int[] temp2 = new int[m_inputNums.length + 15];

    for (int noa = 0; noa < m_numInputs; noa++) {
      temp1[noa] = m_inputList[noa];
      temp2[noa] = m_inputNums[noa];
    }
    m_inputList = temp1;
    m_inputNums = temp2;
  }

  /** 
   * This will connect the specified unit to be an output to this unit.
   * @param o The unit.
   * @param n It's connection number for this connection.
   * @return True if the connection was made, false otherwise.
   */
  protected boolean connectOutput(NeuralConnection o, int n) {
    
    for (int noa = 0; noa < m_numOutputs; noa++) {
      if (o == m_outputList[noa]) {
	return false;
      }
    }
    if (m_numOutputs >= m_outputList.length) {
      //then allocate more space to it.
      allocateOutputs();
    }
    m_outputList[m_numOutputs] = o;
    m_outputNums[m_numOutputs] = n;
    m_numOutputs++;
    return true;
  }
  
  /**
   * Allocates more space for output connection information
   * if the arrays have been filled up.
   */
  protected void allocateOutputs() {
    
    NeuralConnection[] temp1 
      = new NeuralConnection[m_outputList.length + 15];
    
    int[] temp2 = new int[m_outputNums.length + 15];
    
    for (int noa = 0; noa < m_numOutputs; noa++) {
      temp1[noa] = m_outputList[noa];
      temp2[noa] = m_outputNums[noa];
    }
    m_outputList = temp1;
    m_outputNums = temp2;
  }
  
  /**
   * This will disconnect the input with the specific connection number
   * From this node (only on this end however).
   * @param i The unit to disconnect.
   * @param n The connection number at the other end, -1 if all the connections
   * to this unit should be severed.
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
	  //set the other end to have the right connection number.
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
    
    for (int noa = 0; noa < m_numInputs; noa++) {
      //this command will simply remove any connections this node has
      //with the other in 1 go, rather than seperately.
      m_inputList[noa].disconnectOutput(this, -1);
    }
    
    //now reset the inputs.
    m_inputList = new NeuralConnection[0];
    setType(getType() & (~INPUT));
    if (getNumOutputs() == 0) {
      setType(getType() & (~CONNECTED));
    }
    m_inputNums = new int[0];
    m_numInputs = 0;
    
  }

 

  /**
   * Changes the connection value information for one of the connections.
   * @param n The connection number to change.
   * @param v The value to change it to.
   */
  protected void changeInputNum(int n, int v) {
    
    if (n >= m_numInputs || n < 0) {
      return;
    }

    m_inputNums[n] = v;
  }
  
  /**
   * This will disconnect the output with the specific connection number
   * From this node (only on this end however).
   * @param o The unit to disconnect.
   * @param n The connection number at the other end, -1 if all the connections
   * to this unit should be severed.
   * @return True if the connection was removed, false if the connection was
   * not found.
   */  
  protected boolean disconnectOutput(NeuralConnection o, int n) {
    
    int loc = -1;
    boolean removed = false;
    do {
      loc = -1;
      for (int noa = 0; noa < m_numOutputs; noa++) {
	if (o == m_outputList[noa] && (n == -1 || n == m_outputNums[noa])) {
	  loc =noa;
	  break;
	}
      }
      
      if (loc >= 0) {
	for (int noa = loc+1; noa < m_numOutputs; noa++) {
	  m_outputList[noa-1] = m_outputList[noa];
	  m_outputNums[noa-1] = m_outputNums[noa];

	  //set the other end to have the right connection number
	  m_outputList[noa-1].changeInputNum(m_outputNums[noa-1], noa-1);
	}
	m_numOutputs--;
	removed = true;
      }
    } while (n == -1 && loc != -1);
    
    return removed;
  }

  /**
   * This function will remove all outputs to this unit.
   * In doing so it will also terminate the connections at the other end.
   */
  public void removeAllOutputs() {
    
    for (int noa = 0; noa < m_numOutputs; noa++) {
      //this command will simply remove any connections this node has
      //with the other in 1 go, rather than seperately.
      m_outputList[noa].disconnectInput(this, -1);
    }
    
    //now reset the inputs.
    m_outputList = new NeuralConnection[0];
    m_outputNums = new int[0];
    setType(getType() & (~OUTPUT));
    if (getNumInputs() == 0) {
      setType(getType() & (~CONNECTED));
    }
    m_numOutputs = 0;
    
  }

  /**
   * Changes the connection value information for one of the connections.
   * @param n The connection number to change.
   * @param v The value to change it to.
   */
  protected void changeOutputNum(int n, int v) {
    
    if (n >= m_numOutputs || n < 0) {
      return;
    }

    m_outputNums[n] = v;
  }
  
  /**
   * @return The number of input connections.
   */
  public int getNumInputs() {
    return m_numInputs;
  }

  /**
   * @return The number of output connections.
   */
  public int getNumOutputs() {
    return m_numOutputs;
  }


  /**
   * Connects two units together.
   * @param s The source unit.
   * @param t The target unit.
   * @return True if the units were connected, false otherwise.
   */
  public static boolean connect(NeuralConnection s, NeuralConnection t) {
    
    if (s == null || t == null) {
      return false;
    }
    //this ensures that there is no existing connection between these 
    //two units already. This will also cause the current weight there to be 
    //lost
 
    disconnect(s, t);
    if (s == t) {
      return false;
    }
    if ((t.getType() & PURE_INPUT) == PURE_INPUT) {
      return false;   //target is an input node.
    }
    if ((s.getType() & PURE_OUTPUT) == PURE_OUTPUT) {
      return false;   //source is an output node
    }
    if ((s.getType() & PURE_INPUT) == PURE_INPUT 
	&& (t.getType() & PURE_OUTPUT) == PURE_OUTPUT) {      
      return false;   //there is no actual working node in use
    }
    if ((t.getType() & PURE_OUTPUT) == PURE_OUTPUT && t.getNumInputs() > 0) {
      return false; //more than 1 node is trying to feed a particular output
    }

    if ((t.getType() & PURE_OUTPUT) == PURE_OUTPUT &&
	(s.getType() & OUTPUT) == OUTPUT) {
      return false; //an output node already feeding out a final answer
    }

    if (!s.connectOutput(t, t.getNumInputs())) {
      return false;
    }
    if (!t.connectInput(s, s.getNumOutputs() - 1)) {
      
      s.disconnectOutput(t, t.getNumInputs());
      return false;

    }

    //now ammend the type.
    if ((s.getType() & PURE_INPUT) == PURE_INPUT) {
      t.setType(t.getType() | INPUT);
    }
    else if ((t.getType() & PURE_OUTPUT) == PURE_OUTPUT) {
      s.setType(s.getType() | OUTPUT);
    }
    t.setType(t.getType() | CONNECTED);
    s.setType(s.getType() | CONNECTED);
    return true;
  }

  /**
   * Disconnects two units.
   * @param s The source unit.
   * @param t The target unit.
   * @return True if the units were disconnected, false if they weren't
   * (probably due to there being no connection).
   */
  public static boolean disconnect(NeuralConnection s, NeuralConnection t) {
    
    if (s == null || t == null) {
      return false;
    }

    boolean stat1 = s.disconnectOutput(t, -1);
    boolean stat2 = t.disconnectInput(s, -1);
    if (stat1 && stat2) {
      if ((s.getType() & PURE_INPUT) == PURE_INPUT) {
	t.setType(t.getType() & (~INPUT));
      }
      else if ((t.getType() & (PURE_OUTPUT)) == PURE_OUTPUT) {
	s.setType(s.getType() & (~OUTPUT));
      }
      if (s.getNumInputs() == 0 && s.getNumOutputs() == 0) {
	s.setType(s.getType() & (~CONNECTED));
      }
      if (t.getNumInputs() == 0 && t.getNumOutputs() == 0) {
	t.setType(t.getType() & (~CONNECTED));
      }
    }
    return stat1 && stat2;
  }
}













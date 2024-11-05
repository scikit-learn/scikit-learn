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
 *    NeuralMethod.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 */


package weka.classifiers.functions.neural;

import java.io.Serializable;

/**
 * This is an interface used to create classes that can be used by the 
 * neuralnode to perform all it's computations.
 *
 * @author Malcolm Ware (mfw4@cs.waikato.ac.nz)
 * @version $Revision: 1.6 $
 */
public interface NeuralMethod extends Serializable {
  
  /**
   * This function calculates what the output value should be.
   * @param node The node to calculate the value for.
   * @return The value.
   */
  double outputValue(NeuralNode node);

  /**
   * This function calculates what the error value should be.
   * @param node The node to calculate the error for.
   * @return The error.
   */
  double errorValue(NeuralNode node);

  /**
   * This function will calculate what the change in weights should be
   * and also update them.
   * @param node The node to update the weights for.
   * @param learn The learning rate to use.
   * @param momentum The momentum to use.
   */
  void updateWeights(NeuralNode node, double learn, double momentum);

}

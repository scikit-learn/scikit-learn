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
 * RelationalLocator.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;


/**
 * This class locates and records the indices of relational attributes, 
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6228 $
 * @see Attribute#RELATIONAL
 */
public class RelationalLocator
  extends AttributeLocator {

  /** for serialization */

  private static final long serialVersionUID = 4646872277151854732L;

  /**
   * Initializes the RelationalLocator with the given data.
   * 
   * @param data	the data to work on
   */
  public RelationalLocator(Instances data) {
    super(data, Attribute.RELATIONAL);
  }
  
  /**
   * Initializes the RelationalLocator with the given data. 
   * Checks only the given range.
   * 
   * @param data	the data to work on
   * @param fromIndex	the first index to inspect (including)
   * @param toIndex	the last index to check (including)
   */
  public RelationalLocator(Instances data, int fromIndex, int toIndex) {
    super(data, Attribute.RELATIONAL, fromIndex, toIndex);
  }
  
  /**
   * Initializes the RelationalLocator with the given data.
   * Checks only the specified attribute indices.
   * 
   * @param data	the data to work on
   * @param indices	the attribute indices to check
   */
  public RelationalLocator(Instances data, int[] indices) {
    super(data, Attribute.RELATIONAL, indices);
  }

  /**
   * Copies relational values contained in the instance copied to a new
   * dataset. The Instance must already be assigned to a dataset. This
   * dataset and the destination dataset must have the same structure.
   *
   * @param inst 		the Instance containing the relational values 
   * 				to copy.
   * @param destDataset 	the destination set of Instances
   * @param strAtts 		an AttributeLocator containing the indices of 
   * 				any relational attributes in the dataset.  
   */
  public static void copyRelationalValues(Instance inst, Instances destDataset, 
                               AttributeLocator strAtts) {

    if (inst.dataset() == null) {
      throw new IllegalArgumentException("Instance has no dataset assigned!!");
    } else if (inst.dataset().numAttributes() != destDataset.numAttributes()) {
      throw new IllegalArgumentException(
	  "Src and Dest differ in # of attributes: " 
	  + inst.dataset().numAttributes() + " != " + destDataset.numAttributes());
    } 
    copyRelationalValues(inst, true, inst.dataset(), strAtts,
                     destDataset, strAtts);
  }

  /**
   * Takes relational values referenced by an Instance and copies them from a
   * source dataset to a destination dataset. The instance references are
   * updated to be valid for the destination dataset. The instance may have the 
   * structure (i.e. number and attribute position) of either dataset (this
   * affects where references are obtained from). Only works if the number
   * of relational attributes is the same in both indices (implicitly these 
   * relational attributes should be semantically same but just with shifted 
   * positions).
   *
   * @param instance 		the instance containing references to relations 
   * 				in the source dataset that will have references 
   * 				updated to be valid for the destination dataset.
   * @param instSrcCompat 	true if the instance structure is the same as 
   * 				the source, or false if it is the same as the 
   * 				destination (i.e. which of the relational 
   * 				attribute indices contains the correct 
   * 				locations for this instance).
   * @param srcDataset 		the dataset for which the current instance 
   * 				relationvalue references are valid (after any 
   * 				position mapping if needed)
   * @param srcLoc 		an AttributeLocator containing the indices of 
   * 				relational attributes in the source datset.
   * @param destDataset 	the dataset for which the current instance 
   * 				relation references need to be inserted (after 
   * 				any position mapping if needed)
   * @param destLoc 	an AttributeLocator containing the indices of 
   * 				relational attributes in the destination datset.
   */
  public static void copyRelationalValues(Instance instance, boolean instSrcCompat,
                                  Instances srcDataset, AttributeLocator srcLoc,
                                  Instances destDataset, AttributeLocator destLoc) {
	  
    if (srcDataset == destDataset)
      return;
    
    if (srcLoc.getAttributeIndices().length != destLoc.getAttributeIndices().length)
      throw new IllegalArgumentException(
	  "Src and Dest relational indices differ in length: "
	  + srcLoc.getAttributeIndices().length + " != " + destLoc.getAttributeIndices().length);

    if (srcLoc.getLocatorIndices().length != destLoc.getLocatorIndices().length)
      throw new IllegalArgumentException(
	  "Src and Dest locator indices differ in length: "
	  + srcLoc.getLocatorIndices().length + " != " + destLoc.getLocatorIndices().length);

    for (int i = 0; i < srcLoc.getAttributeIndices().length; i++) {
      int instIndex  = instSrcCompat 
      		 	  ? srcLoc.getActualIndex(srcLoc.getAttributeIndices()[i]) 
      		 	  : destLoc.getActualIndex(destLoc.getAttributeIndices()[i]);
      Attribute src  = srcDataset.attribute(srcLoc.getActualIndex(srcLoc.getAttributeIndices()[i]));
      Attribute dest = destDataset.attribute(destLoc.getActualIndex(destLoc.getAttributeIndices()[i]));
      if (!instance.isMissing(instIndex)) {
        int valIndex = dest.addRelation(src.relation((int)instance.value(instIndex)));
        instance.setValue(instIndex, (double)valIndex);
      }
    }
    
    // recurse if necessary
    int[] srcIndices  = srcLoc.getLocatorIndices();
    int[] destIndices = destLoc.getLocatorIndices();
    for (int i = 0; i < srcIndices.length; i++) {
      int index = instSrcCompat
	         ? srcLoc.getActualIndex(srcIndices[i])
	         : destLoc.getActualIndex(destIndices[i]);
      if (instance.isMissing(index))
	continue;
      Instances rel = instSrcCompat
      		         ? instance.relationalValue(index)
      		         : instance.relationalValue(index);
      AttributeLocator srcRelAttsNew = srcLoc.getLocator(srcIndices[i]);
      Instances srcDatasetNew = srcRelAttsNew.getData();
      AttributeLocator destRelAttsNew = destLoc.getLocator(destIndices[i]);
      Instances destDatasetNew = destRelAttsNew.getData();
      for (int n = 0; n < rel.numInstances(); n++) {
        copyRelationalValues(rel.instance(n), instSrcCompat, srcDatasetNew, srcRelAttsNew, destDatasetNew, destRelAttsNew);
      }
    }
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6228 $");
  }
}

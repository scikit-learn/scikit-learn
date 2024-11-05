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
 *    RemoveFrequentValues.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.filters.unsupervised.instance;

import weka.core.Attribute;
import weka.core.AttributeStats;
import weka.core.Capabilities;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.RevisionUtils;
import weka.core.SingleIndex;
import weka.core.UnsupportedAttributeTypeException;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
import weka.filters.UnsupervisedFilter;

import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

/** 
 <!-- globalinfo-start -->
 * Determines which values (frequent or infrequent ones) of an (nominal) attribute are retained and filters the instances accordingly. In case of values with the same frequency, they are kept in the way they appear in the original instances object. E.g. if you have the values "1,2,3,4" with the frequencies "10,5,5,3" and you chose to keep the 2 most common values, the values "1,2" would be returned, since the value "2" comes before "3", even though they have the same frequency.
 * <p/>
 <!-- globalinfo-end -->
 * 
 <!-- options-start -->
 * Valid options are: <p/>
 * 
 * <pre> -C &lt;num&gt;
 *  Choose attribute to be used for selection.</pre>
 * 
 * <pre> -N &lt;num&gt;
 *  Number of values to retain for the sepcified attribute, 
 *  i.e. the ones with the most instances (default 2).</pre>
 * 
 * <pre> -L
 *  Instead of values with the most instances the ones with the 
 *  least are retained.
 * </pre>
 * 
 * <pre> -H
 *  When selecting on nominal attributes, removes header
 *  references to excluded values.</pre>
 * 
 * <pre> -V
 *  Invert matching sense.</pre>
 * 
 <!-- options-end -->
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5499 $
 */
public class RemoveFrequentValues 
   extends Filter 
   implements OptionHandler, UnsupervisedFilter {
  
   /** for serialization */
   static final long serialVersionUID = -2447432930070059511L;

   /** The attribute's index setting. */
   private SingleIndex m_AttIndex = new SingleIndex("last"); 

   /** the number of values to retain. */
   protected int m_NumValues = 2; 
   
   /** whether to retain values with least instances instead of most. */
   protected boolean m_LeastValues = false;
   
   /** whether to invert the matching sense. */
   protected boolean m_Invert = false;
   
   /** Modify header for nominal attributes? */
   protected boolean m_ModifyHeader = false;
   
   /** If m_ModifyHeader, stores a mapping from old to new indexes */
   protected int [] m_NominalMapping;
   
   /** contains the values to retain */
   protected HashSet m_Values = null;
   
   /**
    * Returns a string describing this filter
    * @return a description of the classifier suitable for
    * displaying in the explorer/experimenter gui
    */
   public String globalInfo() {
     return 
         "Determines which values (frequent or infrequent ones) of an "
       + "(nominal) attribute are retained and filters the instances "
       + "accordingly. In case of values with the same frequency, they are "
       + "kept in the way they appear in the original instances object. E.g. "
       + "if you have the values \"1,2,3,4\" with the frequencies \"10,5,5,3\" "
       + "and you chose to keep the 2 most common values, the values \"1,2\" "
       + "would be returned, since the value \"2\" comes before \"3\", even "
       + "though they have the same frequency.";
   }

   /**
    * Returns an enumeration describing the available options.
    *
    * @return an enumeration of all the available options.
    */
   public Enumeration listOptions() {
      Vector newVector = new Vector(5);

      newVector.addElement(new Option(
                "\tChoose attribute to be used for selection.",
                "C", 1, "-C <num>"));
      newVector.addElement(new Option(
                  "\tNumber of values to retain for the sepcified attribute, \n"
                + "\ti.e. the ones with the most instances (default 2).",
                "N", 1, "-N <num>"));
      newVector.addElement(new Option(
                  "\tInstead of values with the most instances the ones with the \n"
                + "\tleast are retained.\n",
                "L", 0, "-L"));
      newVector.addElement(new Option(
                  "\tWhen selecting on nominal attributes, removes header\n"
            	 + "\treferences to excluded values.",
            	 "H", 0, "-H"));
      newVector.addElement(new Option(
            	 "\tInvert matching sense.",
                "V", 0, "-V"));

      return newVector.elements();
   }

   /**
    * Parses a given list of options. <p/>
    * 
    <!-- options-start -->
    * Valid options are: <p/>
    * 
    * <pre> -C &lt;num&gt;
    *  Choose attribute to be used for selection.</pre>
    * 
    * <pre> -N &lt;num&gt;
    *  Number of values to retain for the sepcified attribute, 
    *  i.e. the ones with the most instances (default 2).</pre>
    * 
    * <pre> -L
    *  Instead of values with the most instances the ones with the 
    *  least are retained.
    * </pre>
    * 
    * <pre> -H
    *  When selecting on nominal attributes, removes header
    *  references to excluded values.</pre>
    * 
    * <pre> -V
    *  Invert matching sense.</pre>
    * 
    <!-- options-end -->
    *
    * @param options the list of options as an array of strings
    * @throws Exception if an option is not supported
    */
   public void setOptions(String[] options) throws Exception {
      String attIndex = Utils.getOption('C', options);
      if (attIndex.length() != 0) {
         setAttributeIndex(attIndex);
      } else {
         setAttributeIndex("last");
      }
      
      String numValues = Utils.getOption('N', options);
      if (numValues.length() != 0) {
         setNumValues(Integer.parseInt(numValues));
      } else {
         setNumValues(2);
      }
      
      setUseLeastValues(Utils.getFlag('L', options));

      setModifyHeader(Utils.getFlag('H', options));

      setInvertSelection(Utils.getFlag('V', options));
      
      if (getInputFormat() != null) {
         setInputFormat(getInputFormat());
      }
   }

   /**
    * Gets the current settings of the filter.
    *
    * @return an array of strings suitable for passing to setOptions
    */
   public String[] getOptions() {
      String [] options = new String [7];
      int current = 0;

      options[current++] = "-C";
      options[current++] = "" + (getAttributeIndex());
      options[current++] = "-N";
      options[current++] = "" + (getNumValues());
      if (getUseLeastValues()) {
        options[current++] = "-H";
      }
      if (getModifyHeader()) {
         options[current++] = "-H";
      }
      if (getInvertSelection()) {
         options[current++] = "-V";
      }
      while (current < options.length) {
        options[current++] = "";
      }
      return options;
   }

   /**
    * Returns the tip text for this property
    * @return tip text for this property suitable for
    * displaying in the explorer/experimenter gui
    */
   public String attributeIndexTipText() {
     return "Choose attribute to be used for selection (default last).";
   }

   /**
    * Get the index of the attribute used.
    *
    * @return the index of the attribute
    */
   public String getAttributeIndex() {
     return m_AttIndex.getSingleIndex();
   }

   /**
    * Sets index of the attribute used.
    *
    * @param attIndex the index of the attribute
    */
   public void setAttributeIndex(String attIndex) {
     m_AttIndex.setSingleIndex(attIndex);
   }

   /**
    * Returns the tip text for this property
    * 
    * @return tip text for this property suitable for
    * displaying in the explorer/experimenter gui
    */
   public String numValuesTipText() {
      return "The number of values to retain.";
   }

   /**
    * Gets how many values are retained
    *
    * @return how many values are retained
    */
   public int getNumValues() {
      return m_NumValues;
   }

   /**
    * Sets how many values are retained  
    *
    * @param numValues the number of values to retain
    */
   public void setNumValues(int numValues) {
      m_NumValues = numValues;
   }

   /**
    * Returns the tip text for this property
    * 
    * @return tip text for this property suitable for
    * displaying in the explorer/experimenter gui
    */
   public String useLeastValuesTipText() {
      return "Retains values with least instance instead of most.";
   }

   /**
    * Gets whether to use values with least or most instances
    *
    * @return true if values with least instances are retained
    */
   public boolean getUseLeastValues() {
      return m_LeastValues;
   }

   /**
    * Sets whether to use values with least or most instances  
    *
    * @param leastValues whether values with least or most instances are retained
    */
   public void setUseLeastValues(boolean leastValues) {
      m_LeastValues = leastValues;
   }

   /**
    * Returns the tip text for this property
    * @return tip text for this property suitable for
    * displaying in the explorer/experimenter gui
    */
   public String modifyHeaderTipText() {
      return "When selecting on nominal attributes, removes header references to "
      + "excluded values.";
   }
   
   /**
    * Gets whether the header will be modified when selecting on nominal
    * attributes.
    *
    * @return true if so.
    */
   public boolean getModifyHeader() {
      return m_ModifyHeader;
   }
   
   /**
    * Sets whether the header will be modified when selecting on nominal
    * attributes.
    *
    * @param newModifyHeader true if so.
    */
   public void setModifyHeader(boolean newModifyHeader) {
      m_ModifyHeader = newModifyHeader;
   }
   
   /**
    * Returns the tip text for this property
    * 
    * @return tip text for this property suitable for
    * displaying in the explorer/experimenter gui
    */
   public String invertSelectionTipText() {
      return "Invert matching sense.";
   }

   /**
    * Get whether the supplied columns are to be removed or kept
    *
    * @return true if the supplied columns will be kept
    */
   public boolean getInvertSelection() {
      return m_Invert;
   }

   /**
    * Set whether selected values should be removed or kept. If true the 
    * selected values are kept and unselected values are deleted. 
    *
    * @param invert the new invert setting
    */
   public void setInvertSelection(boolean invert) {
      m_Invert = invert;
   }

   /** 
    * Returns true if selection attribute is nominal.
    *
    * @return true if selection attribute is nominal
    */
   public boolean isNominal() {
      if (getInputFormat() == null) {
         return false;
      } else {
         return getInputFormat().attribute(m_AttIndex.getIndex()).isNominal();
      }
   }
   
   /** 
    * Returns true if selection attribute is ranking.
    *
    * @return true if selection attribute is ranking
    */
   public boolean isRanking() {
      if (getInputFormat() == null) {
         return false;
      } else {
         return getInputFormat().attribute(m_AttIndex.getIndex()).isRanking();
      }
   }
   
   /**
    * determines the values to retain, it is always at least 1
    * and up to the maximum number of distinct values
    * 
    * @param inst the Instances to determine the values from which are kept  
    */
   public void determineValues(Instances inst) {
      int					i;
      AttributeStats		stats;
      int					attIdx;
      int					min;
      int					max;
      int					count;

      m_AttIndex.setUpper(inst.numAttributes() - 1);
      attIdx = m_AttIndex.getIndex();
      
      // init names
      m_Values = new HashSet();
      
      if (inst == null)
         return;
      
      // number of values to retain
      stats = inst.attributeStats(attIdx);
      if (m_Invert)
         count = stats.nominalCounts.length - m_NumValues;
      else
         count = m_NumValues;
      // out of bounds? -> fix
      if (count < 1)
         count = 1;  // at least one value!
      if (count > stats.nominalCounts.length)
         count = stats.nominalCounts.length;  // at max the existing values
      
      // determine min/max occurences
      Arrays.sort(stats.nominalCounts);
      if (m_LeastValues) {
         min = stats.nominalCounts[0];
         max = stats.nominalCounts[count - 1];
      }
      else {
         min = stats.nominalCounts[(stats.nominalCounts.length - 1) - count + 1];
         max = stats.nominalCounts[stats.nominalCounts.length - 1];
      }
      
      // add values if they are inside min/max (incl. borders) and not more than count
      stats = inst.attributeStats(attIdx);
      for (i = 0; i < stats.nominalCounts.length; i++) {
         if ( (stats.nominalCounts[i] >= min) && (stats.nominalCounts[i] <= max) && (m_Values.size() < count) )
            m_Values.add(inst.attribute(attIdx).value(i));
      }
   }
   
   /**
    * modifies the header of the Instances and returns the format w/o 
    * any instances
    * 
    * @param instanceInfo the instances structure to modify
    * @return the new structure (w/o instances!)
    */
   protected Instances modifyHeader(Instances instanceInfo) {
      instanceInfo = new Instances(getInputFormat(), 0); // copy before modifying
      Attribute oldAtt = instanceInfo.attribute(m_AttIndex.getIndex());
      int [] selection = new int[m_Values.size()];
      Iterator iter = m_Values.iterator();
      int i = 0;
      while (iter.hasNext()) {
         selection[i] = oldAtt.indexOfValue(iter.next().toString());
         i++;
      }
      FastVector newVals = new FastVector();
      for (i = 0; i < selection.length; i++) {
         newVals.addElement(oldAtt.value(selection[i]));
      }
      instanceInfo.deleteAttributeAt(m_AttIndex.getIndex());
      instanceInfo.insertAttributeAt(new Attribute(oldAtt.name(), newVals),
            m_AttIndex.getIndex());
      m_NominalMapping = new int [oldAtt.numValues()];
      for (i = 0; i < m_NominalMapping.length; i++) {
         boolean found = false;
         for (int j = 0; j < selection.length; j++) {
            if (selection[j] == i) {
               m_NominalMapping[i] = j;
               found = true;
               break;
            }
         }
         if (!found) {
            m_NominalMapping[i] = -1;
         }
      }
      
      return instanceInfo;
   }

   /** 
    * Returns the Capabilities of this filter.
    *
    * @return            the capabilities of this object
    * @see               Capabilities
    */
   public Capabilities getCapabilities() {
     Capabilities result = super.getCapabilities();
     result.disableAll();

     // attributes
     result.enableAllAttributes();
     result.enable(Capability.MISSING_VALUES);
     
     // class
     result.enableAllClasses();
     result.enable(Capability.MISSING_CLASS_VALUES);
     result.enable(Capability.NO_CLASS);
     
     return result;
   }

   /**
    * Sets the format of the input instances.
    *
    * @param instanceInfo an Instances object containing the input instance
    * structure (any instances contained in the object are ignored - only the
    * structure is required).
    * @return true if the outputFormat can be collected immediately
    * @throws UnsupportedAttributeTypeException if the specified attribute
    *         is not nominal.
    */
   public boolean setInputFormat(Instances instanceInfo) throws Exception {
      super.setInputFormat(instanceInfo);
      
      m_AttIndex.setUpper(instanceInfo.numAttributes() - 1);

      if (!isNominal() && !isRanking())
         throw new UnsupportedAttributeTypeException("Can only handle nominal attributes.");
      
      m_Values = null;
         
      return false;
   }
   
   /**
    * Set the output format. Takes the currently defined Values to retain and 
    * m_InputFormat and calls setOutputFormat(Instances) appropriately. 
    * Those instances that have a value to retain are "push"ed to the output.
    */
   protected void setOutputFormat() {
      Instances      instances;
      int            i;
      Instance       instance;
      
      if (m_Values == null) {
         setOutputFormat(null);
         return;
      }
      
      // get structure
      if (getModifyHeader())
         instances = modifyHeader(getInputFormat());
      else
         instances = new Instances(getInputFormat(), 0);
      setOutputFormat(instances);
      
      // remove instances with unwanted values, for the others change the values
      // value if m_ModifyHeader is set
      for (i = 0; i < getInputFormat().numInstances(); i++) {
         instance = getInputFormat().instance(i);
      
         if (m_Values.contains(instance.stringValue(m_AttIndex.getIndex()))) {
            if (getModifyHeader()) {
               instance.setValue(m_AttIndex.getIndex(),
                     m_NominalMapping[(int)instance.value(m_AttIndex.getIndex())]);
            }
            push(instance);
         }
      }
   }
   
   /**
    * Input an instance for filtering. Ordinarily the instance is processed
    * and made available for output immediately. Some filters require all
    * instances be read before producing output.
    *
    * @param instance the input instance
    * @return true if the filtered instance may now be
    * collected with output().
    * @throws IllegalStateException if no input format has been set.
    */
   public boolean input(Instance instance) {
      if (getInputFormat() == null) {
         throw new IllegalStateException("No input instance format defined");
      }
      
      if (m_NewBatch) {
         resetQueue();
         m_NewBatch = false;
      }

      if (isFirstBatchDone()) {
	push(instance);
	return true;
      } 
      else {
	bufferInput(instance);
	return false;
      }
   }

   /**
    * Signifies that this batch of input to the filter is finished. If the 
    * filter requires all instances prior to filtering, output() may now 
    * be called to retrieve the filtered instances.
    *
    * @return true if there are instances pending output
    * @throws IllegalStateException if no input structure has been defined
    */
   public boolean batchFinished() {
      if (getInputFormat() == null) {
         throw new IllegalStateException("No input instance format defined");
      }

      // process input
      if (m_Values == null) {
         determineValues(getInputFormat());
         setOutputFormat();
      } 
      flushInput();
      
      m_NewBatch = true;
      m_FirstBatchDone = true;

      return (numPendingOutput() != 0);
   }
   
   /**
    * Returns the revision string.
    * 
    * @return		the revision
    */
   public String getRevision() {
     return RevisionUtils.extract("$Revision: 5499 $");
   }
   
   /**
    * Main method for testing this class.
    *
    * @param argv should contain arguments to the filter: 
    * use -h for help
    */
   public static void main(String[] argv) {
      runFilter(new RemoveFrequentValues(), argv);
   }
}

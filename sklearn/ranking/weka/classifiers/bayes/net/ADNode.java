/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * ADNode.java
 * Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 * 
 */

package weka.classifiers.bayes.net;

import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformationHandler;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;

import java.io.FileReader;
import java.io.Serializable;

/**
 * The ADNode class implements the ADTree datastructure which increases
 * the speed with which sub-contingency tables can be constructed from
 * a data set in an Instances object. For details, see: <p/>
 *
 <!-- technical-plaintext-start -->
 * Andrew W. Moore, Mary S. Lee (1998). Cached Sufficient Statistics for Efficient Machine Learning with Large Datasets. Journal of Artificial Intelligence Research. 8:67-91.
 <!-- technical-plaintext-end -->
 * <p/>
 *
 <!-- technical-bibtex-start -->
 * BibTeX:
 * <pre>
 * &#64;article{Moore1998,
 *    author = {Andrew W. Moore and Mary S. Lee},
 *    journal = {Journal of Artificial Intelligence Research},
 *    pages = {67-91},
 *    title = {Cached Sufficient Statistics for Efficient Machine Learning with Large Datasets},
 *    volume = {8},
 *    year = {1998}
 * }
 * </pre>
 * <p/>
 <!-- technical-bibtex-end -->
 *
 * @author Remco Bouckaert (rrb@xm.co.nz)
 * @version $Revision: 1.7 $
 */
public class ADNode 
	implements Serializable, TechnicalInformationHandler, RevisionHandler {
  
  	/** for serialization */
  	static final long serialVersionUID = 397409728366910204L;
  
        final static int MIN_RECORD_SIZE = 0;
	
	/** list of VaryNode children **/
	public VaryNode [] m_VaryNodes;
	/** list of Instance children (either m_Instances or m_VaryNodes is instantiated) **/
	public Instance [] m_Instances;

        /** count **/
	public int m_nCount;

        /** first node in VaryNode array **/
        public int m_nStartNode;

        /** Creates new ADNode */
        public ADNode() {
        }

        /**
         * Returns an instance of a TechnicalInformation object, containing 
         * detailed information about the technical background of this class,
         * e.g., paper reference or book this class is based on.
         * 
         * @return the technical information about this class
         */
        public TechnicalInformation getTechnicalInformation() {
          TechnicalInformation 	result;
          
          result = new TechnicalInformation(Type.ARTICLE);
          result.setValue(Field.AUTHOR, "Andrew W. Moore and Mary S. Lee");
          result.setValue(Field.YEAR, "1998");
          result.setValue(Field.TITLE, "Cached Sufficient Statistics for Efficient Machine Learning with Large Datasets");
          result.setValue(Field.JOURNAL, "Journal of Artificial Intelligence Research");
          result.setValue(Field.VOLUME, "8");
          result.setValue(Field.PAGES, "67-91");
          
          return result;
        }

	/** create sub tree
	 * @param iNode index of the lowest node in the tree
	 * @param nRecords set of records in instances to be considered
	 * @param instances data set
         * @return VaryNode representing part of an ADTree
 	 **/
	public static VaryNode makeVaryNode(int iNode, FastVector nRecords, Instances instances) {
		VaryNode _VaryNode = new VaryNode(iNode);
		int nValues = instances.attribute(iNode).numValues();
                

		// reserve memory and initialize
		FastVector [] nChildRecords = new FastVector[nValues];
		for (int iChild = 0; iChild < nValues; iChild++) {
			nChildRecords[iChild] = new FastVector();
		}
		// divide the records among children
		for (int iRecord = 0; iRecord < nRecords.size(); iRecord++) {
			int iInstance = ((Integer) nRecords.elementAt(iRecord)).intValue();
			nChildRecords[(int) instances.instance(iInstance).value(iNode)].addElement(new Integer(iInstance));
		}

		// find most common value
		int nCount = nChildRecords[0].size();
		int nMCV = 0; 
		for (int iChild = 1; iChild < nValues; iChild++) {
			if (nChildRecords[iChild].size() > nCount) {
				nCount = nChildRecords[iChild].size();
				nMCV = iChild;
			}
		}
                _VaryNode.m_nMCV = nMCV;

                // determine child nodes
                _VaryNode.m_ADNodes = new ADNode[nValues];
		for (int iChild = 0; iChild < nValues; iChild++) {
			if (iChild == nMCV || nChildRecords[iChild].size() == 0) {
				_VaryNode.m_ADNodes[iChild] = null;
			} else {
				_VaryNode.m_ADNodes[iChild] = makeADTree(iNode + 1, nChildRecords[iChild], instances);
			}
		}
		return _VaryNode;
	} // MakeVaryNode

	/** 
	 * create sub tree
	 * 
	 * @param iNode index of the lowest node in the tree
	 * @param nRecords set of records in instances to be considered
	 * @param instances data set
         * @return ADNode representing an ADTree
	 */
	public static ADNode makeADTree(int iNode, FastVector nRecords, Instances instances) {
		ADNode _ADNode = new ADNode();
                _ADNode.m_nCount = nRecords.size();
                _ADNode.m_nStartNode = iNode;
                if (nRecords.size() < MIN_RECORD_SIZE) {
                  _ADNode.m_Instances = new Instance[nRecords.size()];
                  for (int iInstance = 0; iInstance < nRecords.size(); iInstance++) {
                    _ADNode.m_Instances[iInstance] = instances.instance(((Integer) nRecords.elementAt(iInstance)).intValue());
                  }
                } else {
                  _ADNode.m_VaryNodes = new VaryNode[instances.numAttributes() - iNode];
                  for (int iNode2 = iNode; iNode2 < instances.numAttributes(); iNode2++) {
                          _ADNode.m_VaryNodes[iNode2 - iNode] = makeVaryNode(iNode2, nRecords, instances);
                  }
                }
		return _ADNode;
	} // MakeADTree

	/** 
	 * create AD tree from set of instances
	 * 
	 * @param instances data set
         * @return ADNode representing an ADTree
	 */
	public static ADNode makeADTree(Instances instances) {
          FastVector nRecords = new FastVector(instances.numInstances());
          for (int iRecord = 0; iRecord < instances.numInstances(); iRecord++) {
            nRecords.addElement(new Integer(iRecord));
          }
          return makeADTree(0, nRecords, instances);
        } // MakeADTree
        
          /** 
           * get counts for specific instantiation of a set of nodes
           * 
           * @param nCounts - array for storing counts
           * @param nNodes - array of node indexes 
           * @param nOffsets - offset for nodes in nNodes in nCounts
           * @param iNode - index into nNode indicating current node
           * @param iOffset - Offset into nCounts due to nodes below iNode
           * @param bSubstract - indicate whether counts should be added or substracted
           */
        public void getCounts(
              int [] nCounts, 
              int [] nNodes, 
              int [] nOffsets, 
              int iNode, 
              int iOffset,
              boolean bSubstract
        ) {
//for (int iNode2 = 0; iNode2 < nCounts.length; iNode2++) {
//   System.out.print(nCounts[iNode2] + " ");
//}
//System.out.println();
          if (iNode >= nNodes.length) {
            if (bSubstract) {
              nCounts[iOffset] -= m_nCount;
            } else {
              nCounts[iOffset] += m_nCount;
            }
            return;
          } else {
            if (m_VaryNodes != null) {
              m_VaryNodes[nNodes[iNode] - m_nStartNode].getCounts(nCounts, nNodes, nOffsets, iNode, iOffset, this, bSubstract);
            } else {
              for (int iInstance = 0; iInstance < m_Instances.length; iInstance++) {
                int iOffset2 = iOffset;
                Instance instance = m_Instances[iInstance];
                for (int iNode2 = iNode; iNode2 < nNodes.length; iNode2++) {
                  iOffset2 = iOffset2 + nOffsets[iNode2] * (int) instance.value(nNodes[iNode2]);
                }
                if (bSubstract) {
	                nCounts[iOffset2]--;
                } else {
                	nCounts[iOffset2]++;
                }
              }
            }
          }
        } // getCounts


        /** 
         * print is used for debugging only and shows the ADTree in ASCII graphics
         */
        public void print() {
          String sTab = new String();for (int i = 0; i < m_nStartNode; i++) {
              sTab = sTab + "  ";
          }
          System.out.println(sTab + "Count = " + m_nCount);
          if (m_VaryNodes != null) {
	          for (int iNode = 0; iNode < m_VaryNodes.length; iNode++) {
	            System.out.println(sTab + "Node " + (iNode + m_nStartNode));
	            m_VaryNodes[iNode].print(sTab);
	          }
          } else {
              System.out.println(m_Instances);
          }
        }
        
        /**
         * for testing only
         * 
         * @param argv the commandline options
         */
        public static void main(String [] argv) {
            try {
                Instances instances = new Instances(new FileReader("\\iris.2.arff"));
                ADNode ADTree = ADNode.makeADTree(instances);
                int [] nCounts = new int[12];
                int [] nNodes = new int[3];
                int [] nOffsets = new int[3];
                nNodes[0] = 0;
                nNodes[1] = 3;
                nNodes[2] = 4;
                nOffsets[0] = 2;
                nOffsets[1] = 1;
                nOffsets[2] = 4;
                ADTree.print();
                ADTree.getCounts(nCounts, nNodes, nOffsets,0, 0, false); 
                
            } catch (Throwable t) {
                t.printStackTrace();
            }
        } // main
        
        /**
         * Returns the revision string.
         * 
         * @return		the revision
         */
        public String getRevision() {
          return RevisionUtils.extract("$Revision: 1.7 $");
        }
} // class ADNode

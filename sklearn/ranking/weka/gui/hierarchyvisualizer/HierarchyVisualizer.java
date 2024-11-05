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
 * HierarchicalClusterer.java
 * Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
*/

package weka.gui.hierarchyvisualizer;
/**
 * Shows cluster trees represented in Newick format as dendrograms.
 *
 * @author Remco Bouckaert (rrb@xm.co.nz, remco@cs.waikato.ac.nz)
 * @version $Revision: 5996 $
 */
import java.awt.Color;
import java.awt.Container;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.Vector;

import javax.swing.JFrame;

import weka.gui.visualize.PrintablePanel;

public class HierarchyVisualizer extends PrintablePanel implements ComponentListener {
	private static final long serialVersionUID = 1L;

	String m_sNewick;
	Node m_tree;
	int m_nLeafs;
	double m_fHeight;
	double m_fScaleX = 10;
	double m_fScaleY = 10;

	public HierarchyVisualizer(String sNewick) {
		try {
			parseNewick(sNewick);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		addComponentListener(this);
	} // c'tor

	int positionLeafs(Node node, int nPosX) {
		if (node.isLeaf()) {
			node.m_fPosX = nPosX + 0.5;
			nPosX++;
			return nPosX;
		} else {
			for (int i = 0; i < node.m_children.length; i++) {
				nPosX = positionLeafs(node.m_children[i], nPosX);
			}
		}
		return nPosX;
	}
	double positionRest(Node node) {
		if (node.isLeaf()) {
			return node.m_fPosX;
		} else {
			double fPosX = 0;
			for (int i = 0; i < node.m_children.length; i++) {
				fPosX += positionRest(node.m_children[i]);
			}
			fPosX /= node.m_children.length;
			node.m_fPosX = fPosX;
			return fPosX;
		}
	}
	double positionHeight(Node node, double fOffSet) {
		if (node.isLeaf()) {
			node.m_fPosY = fOffSet + node.m_fLength;
			return node.m_fPosY;
		} else {
			double fPosY = fOffSet + node.m_fLength;
			double fYMax = 0;
			for (int i = 0; i < node.m_children.length; i++) {
				fYMax = Math.max(fYMax, positionHeight(node.m_children[i], fPosY));
			}
			node.m_fPosY = fPosY;
			return fYMax;
		}
	}

	//Vector<String> m_sMetaDataNames;
	/** class for nodes in building tree data structure **/
	class Node {
		double m_fLength = -1;
		double m_fPosX = 0;
		double m_fPosY = 0;
		String m_sLabel;
		//Vector<String> m_sMetaDataValues;
		String m_sMetaData;


		/** list of children of this node **/
		Node[] m_children;
		/** parent node in the tree, null if root **/
		Node m_Parent = null;

		Node getParent() {
			return m_Parent;
		}

		void setParent(Node parent) {
			m_Parent = parent;
		}

		/** check if current node is root node **/
		boolean isRoot() {
			return m_Parent == null;
		}
		boolean isLeaf() {
			return m_children == null;
		}

		/** return nr of children **/
		int getChildCount() {
//			}
			if (m_children == null) {return 0;}
			return m_children.length;
		}

		Node getChild(int iChild) {
			return m_children[iChild];
		}


		/** count number of nodes in tree, starting with current node **/
		int getNodeCount() {
			if (m_children == null) {
				return 1;
			}
			int n = 1;
			for (int i = 0; i < m_children.length; i++) {
				n += m_children[i].getNodeCount();
			}
			return n;
		}
		public String toString() {
			StringBuffer buf = new StringBuffer();
			if (m_children != null) {
				buf.append("(");
				for (int i = 0; i < m_children.length-1; i++) {
					buf.append(m_children[i].toString());
					buf.append(',');
				}
				buf.append(m_children[m_children.length - 1].toString());
				buf.append(")");
			} else {
				buf.append(m_sLabel);
			}
			if (m_sMetaData != null) {
				buf.append('[');
				buf.append(m_sMetaData);
				buf.append(']');
			}
			buf.append(":" + m_fLength);
			return buf.toString();
		}

		double draw(Graphics g) {
			if (isLeaf()) {
				int x = (int)(m_fPosX * m_fScaleX);
				int y = (int)(m_fPosY * m_fScaleY);
				g.drawString(m_sLabel, x, y);
				g.drawLine((int)(m_fPosX * m_fScaleX), (int)(m_fPosY * m_fScaleY), (int)(m_fPosX * m_fScaleX), (int)((m_fPosY - m_fLength) * m_fScaleY));
			} else {
				double fPosX1 = Double.MAX_VALUE;
				double fPosX2 = -Double.MAX_VALUE;
				for (int i = 0; i < m_children.length; i++) {
					double f = m_children[i].draw(g);
					if (f < fPosX1) {fPosX1 = f;}
					if (f > fPosX2) {fPosX2 = f;}
				}
				g.drawLine((int)(m_fPosX * m_fScaleX), (int)(m_fPosY * m_fScaleY), (int)(m_fPosX * m_fScaleX), (int)((m_fPosY - m_fLength) * m_fScaleY));
				g.drawLine((int)(fPosX1 * m_fScaleX), (int)(m_fPosY * m_fScaleY), (int)(fPosX2 * m_fScaleX), (int)(m_fPosY* m_fScaleY));
			}
			return m_fPosX;
		}
	} // class Node

	/** helper method for parsing Newick tree **/
	int nextNode(String sStr, int i) {
		int nBraces = 0;
		char c = sStr.charAt(i);
		do {
			i++;
			if (i < sStr.length()) {
				c = sStr.charAt(i);
				// skip meta data block
				if (c == '[') {
					while (i < sStr.length() && sStr.charAt(i)!=']') {
						i++;
					}
					i++;
					if(i < sStr.length()) {
						c = sStr.charAt(i);
					}
				}

				switch (c) {
				case '(':
					nBraces++;
					break;
				case ')':
					nBraces--;
					break;
				default:
					break;
				}
			}
		} while (i < sStr.length()
				&& (nBraces > 0 || (c != ','&&c != ')'&&c != '(')));
		if (i >= sStr.length() || nBraces < 0) {
			return -1;
		} else if (sStr.charAt(i) == ')') {
			i++;
			if (sStr.charAt(i) == '[') {
				while (i < sStr.length() && sStr.charAt(i)!=']') {
					i++;
				}
				i++;
				if (i >= sStr.length()) {
					return -1;
				}
			}
			if (sStr.charAt(i) == ':') {
				i++;
				c = sStr.charAt(i);
				while (i < sStr.length() && (c=='.' || Character.isDigit(c))) {
					i++;
					if (i < sStr.length()) {
						c = sStr.charAt(i);
					}
				}
			}
		}
		return i;
	}

	/**
	 * convert string containing Newick tree into tree datastructure but only in
	 * the limited format as contained in m_sTrees
	 *
	 * @param sStr
	 * @return tree consisting of a Node
	 */
	void parseNewick(String sNewick) throws Exception {
		m_sNewick = sNewick;
		int i = m_sNewick.indexOf('(');
		if (i > 0) {
			m_sNewick = m_sNewick.substring(i);
		}
		System.err.println(m_sNewick);
		m_tree = parseNewick2(m_sNewick);
		System.err.println(m_tree.toString());
		m_nLeafs = positionLeafs(m_tree, 0);
		positionRest(m_tree);
		m_fHeight = positionHeight(m_tree, 0);
	}

	Node parseNewick2(String sStr) throws Exception {
		// System.out.println(sStr);
		if (sStr == null || sStr.length() == 0) {
			return null;
		}
		Node node = new Node();
		if (sStr.startsWith("(")) {
			int i1 = nextNode(sStr, 0);
			int i2 = nextNode(sStr, i1);
			node.m_children = new Node[2];
			node.m_children[0] = parseNewick2(sStr.substring(1, i1));
			node.m_children[0].m_Parent = node;
			String sStr2 = sStr.substring(i1+1, (i2>0?i2:sStr.length()));
			node.m_children[1] = parseNewick2(sStr2);
			node.m_children[1].m_Parent = node;
			if (sStr.lastIndexOf('[') > sStr.lastIndexOf(')')) {
				sStr = sStr.substring(sStr.lastIndexOf('['));
				i2 = sStr.indexOf(']');
				if (i2 < 0) {
					throw new Exception("unbalanced square bracket found:" + sStr);
				}
				sStr2 = sStr.substring(1, i2);
				node.m_sMetaData = sStr2;
			}
			if (sStr.lastIndexOf(':') > sStr.lastIndexOf(')')) {
				sStr = sStr.substring(sStr.lastIndexOf(':'));
				sStr = sStr.replaceAll("[,\\):]", "");
				node.m_fLength = new Double(sStr);
			} else {
				node.m_fLength = 1;
			}
		} else {
			// it is a leaf
			if (sStr.contains("[")) {
				// grab metadata
				int i1 = sStr.indexOf('[');
				int i2 = sStr.indexOf(']');
				if (i2 < 0) {
					throw new Exception("unbalanced square bracket found:" + sStr);
				}
				String sStr2 = sStr.substring(i1+1, i2);
				sStr = sStr.substring(0, i1) +sStr.substring(i2+1);
				node.m_sMetaData = sStr2;
			}
			if (sStr.indexOf(')') >= 0) {
				sStr = sStr.substring(0, sStr.indexOf(')'));
			}
			sStr = sStr.replaceFirst("[,\\)]", "");
			// System.out.println("parsing <<"+sStr+">>");
			if (sStr.length() > 0) {
				if (sStr.indexOf(':') >= 0) {
					int iColon = sStr.indexOf(':');
					node.m_sLabel = sStr.substring(0, iColon);
					if (sStr.indexOf(':', iColon+1) >= 0) {
						int iColon2 = sStr.indexOf(':', iColon+1);
						node.m_fLength = new Double(sStr.substring(iColon+1, iColon2));
						m_fTmpLength = new Double(sStr.substring(iColon2+1));
					} else {
						node.m_fLength = new Double(sStr.substring(iColon+1));
					}
				} else {
					node.m_sLabel = sStr;
					node.m_fLength = 1;
				}
			} else {
				return null;
			}
		}
		return node;
	}
	double m_fTmpLength;

	/**
	 * Fits the tree to the current screen size. Call this after window has been
	 * created to get the entire tree to be in view upon launch.
	 */
	public void fitToScreen() {
		m_fScaleX = 10;
		int nW = getWidth();
		if (m_nLeafs > 0) {
			m_fScaleX = nW / m_nLeafs;
		}
		m_fScaleY = 10;
		int nH = getHeight();
		if (m_fHeight > 0) {
			m_fScaleY = (nH - 10) / m_fHeight;
		}
		repaint();
	}

	/**
	 * Updates the screen contents.
	 *
	 * @param g
	 *            the drawing surface.
	 */
	public void paintComponent(Graphics g) {
		Color oldBackground = ((Graphics2D) g).getBackground();
		// if (m_BackgroundColor != null)
		((Graphics2D) g).setBackground(Color.WHITE);
		g.clearRect(0, 0, getSize().width, getSize().height);
		((Graphics2D) g).setBackground(oldBackground);
		g.setClip(3, 7, getWidth() - 6, getHeight() - 10);
		m_tree.draw(g);
		g.setClip(0, 0, getWidth(), getHeight());
	}


	public void componentHidden(ComponentEvent e) {}

	public void componentMoved(ComponentEvent e) {}

	public void componentResized(ComponentEvent e) {
		fitToScreen();
	}

	public void componentShown(ComponentEvent e) {}

	/**
	 * Main method for testing this class.
	 */
	public static void main(String[] args) {
	      //HierarchyVisualizer a = new HierarchyVisualizer("((((human:2.0,(chimp:1.0,bonobo:1.0):1.0):1.0,gorilla:3.0):1.0,siamang:4.0):1.0,orangutan:5.0)");
	      //HierarchyVisualizer a = new HierarchyVisualizer("(((human:2.0,(chimp:1.0,bonobo:1.0):1.0):1.0,gorilla:3.0):1.0,siamang:4.0)");
	      HierarchyVisualizer a = new HierarchyVisualizer(" (((5[theta=0.121335,lxg=0.122437]:0.00742795,3[theta=0.0972485,lxg=0.152762]:0.00742795)[theta=0.490359,lxg=0.0746703]:0.0183076,((2[theta=0.0866056,lxg=0.2295]:0.00993801,4[theta=0.135512,lxg=0.146674]:0.00993801)[theta=0.897783,lxg=0.0200762]:0.00901206,1[theta=0.200265,lxg=0.18925]:0.0189501)[theta=0.0946195,lxg=0.143427]:0.00678551)[theta=0.185562,lxg=0.139681]:0.0129598,(7[theta=0.176022,lxg=0.364039]:0.0320395,((0[theta=0.224286,lxg=0.156485]:0.0175487,8[theta=0.223313,lxg=0.157166]:0.0175487)[theta=0.631287,lxg=0.024042]:0.00758871,6[theta=0.337871,lxg=0.148799]:0.0251374)[theta=0.33847,lxg=0.040784]:0.00690208)[theta=0.209238,lxg=0.0636202]:0.00665587)[theta=0.560453,lxg=-0.138086]:0.01");
		  //HierarchyVisualizer a = new HierarchyVisualizer(" ((5[theta=0.121335,lxg=0.122437]:0.00742795,3[theta=0.0972485,lxg=0.152762]:0.00742795)[theta=0.490359,lxg=0.0746703]:0.0183076,2[theta=0.0866056,lxg=0.2295]:0.00993801)[theta=0.897783,lxg=0.0200762]:0.00901206");
	      a.setSize(800 ,600);
	      JFrame f;
	      f = new JFrame();
	      Container contentPane = f.getContentPane();
	      contentPane.add(a);
	      f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	      f.setSize(800,600);
	      f.setVisible(true);
	      a.fitToScreen();
	  }

} // class HierarchyVisualizer

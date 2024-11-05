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
 *    HierarchyPropertyParser.java
 *    Copyright (C) 2001 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.io.Serializable;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * This class implements a parser to read properties that have
 * a hierarchy(i.e. tree) structure.  Conceptually it's similar to 
 * the XML DOM/SAX parser but of course is much simpler and 
 * uses dot as the seperator of levels instead of back-slash.<br>
 * It provides interfaces to both build a parser tree and traverse
 * the tree. <br>
 * Note that this implementation does not lock the tree when different
 * threads are traversing it simultaneously, i.e. it's NOT synchronized
 * and multi-thread safe.  It is recommended that later implementation
 * extending this class provide a locking scheme and override the 
 * functions with the "synchronized" modifier (most of them are 
 * goToXXX() and information accessing functions).<p>
 *
 * @author Xin Xu (xx5@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class HierarchyPropertyParser
    implements Serializable {

    /** for serialization */
    private static final long serialVersionUID = -4151103338506077544L;

    /** Keep track of the root of the tree */
    private TreeNode m_Root;
    
    /** Keep track of the current node when traversing the tree */
    private TreeNode m_Current;

    /** The level separate in the path */
    private String m_Seperator = ".";
    
    /** The depth of the tree */
    private int m_Depth = 0;
	
    /** 
     * The inner class implementing a single tree node.
     * All fields are made public simply for convenient access,
     * Although a severe violation of OO Design principle.
     */
    private class TreeNode{
	/** The parent of this node */
	public TreeNode parent = null;
	
	/** The value of this node.  Always String  */
	public String value = null;

	/** The children of this node */
	public Vector children = null;	

	/** The level of this node */
	public int level = 0;

	/** The context of this node */
	public String context = null;
    }
    
    /** Default constructor */
    public HierarchyPropertyParser(){
	m_Root = new TreeNode();
	m_Root.parent = null;
	m_Root.children = new Vector();
	goToRoot();
    }
    
    /**
     * Constructor that builds a tree from the given property with
     * the given delimitor
     *
     * @param p the given property string
     * @param delim the given dilimitor
     */
    public HierarchyPropertyParser(String p, String delim) throws Exception {
	this();		
	build(p, delim);
    }
    
    /**
     * Set the seperator between levels.  Default is dot.
     *
     * @param s the seperator symbol
     */
    public void setSeperator(String s){
	m_Seperator = s;
    }

    /**
     * Get the seperator between levels.  Default is dot.
     *
     * @return the seperator symbol
     */
    public String getSeperator(){ return m_Seperator; }
   
    /** 
     * Build a tree from the given property with the given delimitor
     *
     * @param p the given property
     * @param delim the given delimitor
     */
    public void build(String p, String delim) throws Exception {
	StringTokenizer st = new StringTokenizer(p, delim);
	//System.err.println("delim: "+delim);
	while(st.hasMoreTokens()){
	    String property = st.nextToken().trim();
	    if(!isHierachic(property))
		throw new Exception("The given property is not in"+
				    "hierachy structure with seperators!");
	    add(property);	   
	}
	goToRoot();
    }

    /**
     * Add the given item of property to the tree
     *
     * @param property the given item
     */ 
    public synchronized void add(String property){
	String[] values = tokenize(property);
	if(m_Root.value == null)
	    m_Root.value = values[0];
	
	buildBranch(m_Root, values, 1);
    }
    
    /** 
     * Private function to build one branch of the tree
     * based on one property 
     * 
     * @param parent the parent of the node to be built
     * @param values the value of one property
     * @param lvl the level of the node to be built in the tree
     */
    private void buildBranch(TreeNode parent, String[] values, int lvl){
	// Precondition: children is not null
	
	if(lvl == values.length){  // Parent is leaf
	    parent.children = null;
	    return;
	}
	
	if(lvl > (m_Depth-1))
	    m_Depth = lvl+1;  // Depth starts from 1
	
	Vector kids = parent.children;
	int index = search(kids, values[lvl]);
	if(index != -1){
	    TreeNode newParent = (TreeNode)kids.elementAt(index);
	    if(newParent.children == null)
		newParent.children = new Vector();
	    buildBranch(newParent, values, lvl+1);
	}
	else{
	    TreeNode added = new TreeNode();
	    added.parent = parent;
	    added.value = values[lvl];
	    added.children = new Vector();
	    added.level = lvl;
	    if(parent != m_Root)
		added.context = parent.context + m_Seperator + parent.value;
	    else
		added.context = parent.value;
	    
	    kids.addElement(added);
	    buildBranch(added, values, lvl+1);
	}
    }
    
    /**
     * Tokenize the given string based on the seperator and
     * put the tokens into an array of strings
     *
     * @param rawString the given string
     * @return an array of strings
     */
    public String[] tokenize(String rawString){
	Vector result = new Vector();
	StringTokenizer tk = new StringTokenizer(rawString, m_Seperator);
	while(tk.hasMoreTokens())
	    result.addElement(tk.nextToken());
	
	String[] newStrings = new String[result.size()];
	for(int i=0; i < result.size(); i++)
	    newStrings[i] = (String)result.elementAt(i);
	
	return newStrings;
    }
    
    /**
     * Whether the HierarchyPropertyParser contains the given
     * string
     *
     * @param string the given string
     * @return whether contains
     */
    public boolean contains(String string){
	String[] item = tokenize(string);
	if(!item[0].equals(m_Root.value))
	    return false;
	
	return isContained(m_Root, item, 1);
    }
    
    /** 
     * Private function to decide whether one level of one branch 
     * contains the relevant values
     * 
     * @param parent the parent of the node to be searched
     * @param values the value of one property
     * @param lvl the level of the node in question
     * @return whether this branch contains the corresponding values
     */
    private boolean isContained(TreeNode parent, String[] values, int lvl){
	if(lvl == values.length)  // Parent is leaf
	    return true;
	else if(lvl > values.length)
	    return false;
	else{
	    Vector kids = parent.children;
	    int index = search(kids, values[lvl]);
	    if(index != -1){
		TreeNode newParent = (TreeNode)kids.elementAt(index);
		return isContained(newParent, values, lvl+1);
	    }
	    else
		return false;
	}
    }
    
    /** 
     * Whether the given string has a hierachy structure with
     * the seperators
     *
     * @param string the given string
     */
    public boolean isHierachic(String string){
	int index = string.indexOf(m_Seperator); 
	// Seperator not occur or first occurance at the end
	if((index == (string.length()-1)) || (index == -1))
	    return false;
	
	return true;
    }

    /**
     * Helper function to search for the given target string in a 
     * given vector in which the elements' value may hopefully is equal
     * to the target.  If such elements are found the first index
     * is returned, otherwise -1
     *
     * @param vct the given vector
     * @param target the given target string
     * @return the index of the found element, -1 if not found
     */
    public int search(Vector vct, String target){
	if(vct == null)
	    return -1;
	
	for(int i=0; i < vct.size(); i++)
	    if(target.equals(((TreeNode)vct.elementAt(i)).value))
		return i;
	
	return -1;
    }

    /**
     * Go to a certain node of the tree according to the specified path
     * Note that the path must be absolute path from the root.  <br>
     * For relative path, see goDown(String path).
     *
     * @param path the given absolute path
     * @return whether the path exists, if false the current position does
     * not move
     */
    public synchronized boolean goTo(String path){
	if(!isHierachic(path)){
	    if(m_Root.value.equals(path)){
		goToRoot();
		return true;
	    }
	    else
		return false;
	}
	
	TreeNode old = m_Current;
	m_Current = new TreeNode(); 
	goToRoot();
	String[] nodes = tokenize(path);
	if(!m_Current.value.equals(nodes[0]))
	    return false;

	for(int i=1; i < nodes.length; i++){
	    int pos = search(m_Current.children, nodes[i]);
	    if(pos == -1){
		m_Current = old;
		return false;
	    }
	    m_Current = (TreeNode)m_Current.children.elementAt(pos);
	}
	
	return true;
    }

    /**
     * Go to a certain node of the tree down from the current node 
     * according to the specified relative path.  The path does not
     * contain the value of current node
     * 
     * @param path the given relative path
     * @return whether the path exists, if false the current position does
     * not move
     */
    public synchronized boolean goDown(String path){
	if(!isHierachic(path))
	    return goToChild(path);

	TreeNode old = m_Current;
	m_Current = new TreeNode(); 
	String[] nodes = tokenize(path);
	int pos = search(old.children, nodes[0]);
	if(pos == -1){
	    m_Current = old;
	    return false;
	}
      
	m_Current = (TreeNode)old.children.elementAt(pos);
	for(int i=1; i < nodes.length; i++){
	    pos = search(m_Current.children, nodes[i]);
	    if(pos == -1){
		m_Current = old;
		return false;
	    }
	    
	    m_Current = (TreeNode)m_Current.children.elementAt(pos);
	}
	
	return true;
    }

    /**
     * Go to the root of the tree
     */    
    public synchronized void goToRoot(){
	m_Current=m_Root;
    }

    /**
     * Go to the parent from the current position in the tree
     * If the current position is the root, it stays there and does 
     * not move
     */
    public synchronized void goToParent(){
	if(m_Current.parent != null)  // Not root
	    m_Current = m_Current.parent;
    }
    
    /**
     * Go to one child node from the current position in the tree
     * according to the given value <br>
     * If the child node with the given value cannot be found it
     * returns false, true otherwise.  If false, the current position
     * does not change
     * 
     * @param value the value of the given child
     * @return whether the child can be found
     */
    public synchronized boolean goToChild(String value){
	if(m_Current.children == null) // Leaf
	    return false;
	
	int pos = search(m_Current.children, value);
	if(pos == -1)
	    return false;
	
	m_Current = (TreeNode)m_Current.children.elementAt(pos);	
	return true;
    }
    
    /**
     * Go to one child node from the current position in the tree
     * according to the given position <br>
     *
     * @param pos the position of the given child
     * @exception Exception if the position is out of range or leaf is reached
     */
    public synchronized void goToChild(int pos) throws Exception {
	if((m_Current.children == null) || 
	   (pos < 0) || (pos >= m_Current.children.size()))
	    throw new Exception("Position out of range or leaf reached");
	
	m_Current = (TreeNode)m_Current.children.elementAt(pos);
    }

    /**
     * The number of the children nodes. If current node is leaf, it
     * returns 0.
     * 
     * @return the number of the children nodes of the current position
     */    
    public synchronized int numChildren(){
	if(m_Current.children == null) // Leaf
	    return 0;
	
	return m_Current.children.size(); 
    }
    
    /**
     * The value in the children nodes. If current node is leaf, it
     * returns null.
     * 
     * @return the value in the children nodes
     */    
    public synchronized String[] childrenValues(){
	if(m_Current.children == null) // Leaf
	    return null;
	else{
	    Vector kids = m_Current.children;
	    String[] values = new String[kids.size()];
	    for(int i=0; i < kids.size(); i++)
		values[i] = ((TreeNode)kids.elementAt(i)).value;
	    return values;
	}
    }

    /**
     * The value in the parent node. If current node is root, it
     * returns null.
     * 
     * @return the value in the parent node
     */    
    public synchronized String parentValue(){
	if(m_Current.parent != null)  // Not root
	    return m_Current.parent.value;
	else return null;
    }

    /**
     * Whether the current position is a leaf
     *
     * @return whether the current position is a leaf
     */
    public synchronized boolean isLeafReached(){
	return (m_Current.children == null);
    }
    
    /**
     * Whether the current position is the root
     *
     * @return whether the current position is the root
     */
    public synchronized boolean isRootReached(){
	return (m_Current.parent == null);
    }
    
    /** 
     * Get the value of current node
     *
     * @return value level
     */ 
    public synchronized String getValue(){ return m_Current.value; }
    
    /** 
     * Get the level of current node.  Note the level starts from 0
     *
     * @return the level
     */ 
    public synchronized int getLevel(){ return m_Current.level; }

    /** 
     * Get the depth of the tree, i.e. (the largest level)+1
     *
     * @return the depth of the tree
     */
    public int depth(){ return m_Depth; }
    
    /**
     * The context of the current node, i.e. the path from the
     * root to the parent node of the current node, seperated by
     * the seperator.  If root, it returns null
     *
     * @return the context path
     */
    public synchronized String context(){
	return m_Current.context;
    }
    
    /**
     * The full value of the current node, i.e. its context + seperator
     * + its value.  For root, only its value.
     *
     * @return the context path
     */
    public synchronized String fullValue(){
	if(m_Current == m_Root)
	    return m_Root.value;
	else    
	    return (m_Current.context + m_Seperator + m_Current.value);
    }


    /**
     * Show the whole tree in text format
     *
     * @return the whole tree in text format
     */
    public String showTree(){
	return showNode(m_Root, null);
    }
    
    /**
     * Show one node of the tree in text format
     *
     * @param node the node in question
     * @return the node in text format
     */
    private String showNode(TreeNode node, boolean[] hasBar){
	StringBuffer text =  new StringBuffer();	

	for(int i=0; i < (node.level - 1); i++)
	    if(hasBar[i])
		text.append("  |       ");
	    else
		text.append("          ");

	if(node.level != 0)
	    text.append("  |------ ");
	text.append(node.value+"("+node.level+")"+"["+node.context+"]\n");

	if(node.children != null)
	    for(int i=0; i < node.children.size(); i++){
		boolean[] newBar = new boolean[node.level+1];
		int lvl = node.level;

		if(hasBar != null)
		    for(int j=0; j < lvl; j++)
			newBar[j] = hasBar[j];
		
		if((i == (node.children.size()-1)))
		    newBar[lvl] = false;
		else
		    newBar[lvl] = true;

		text.append(showNode((TreeNode)node.children.elementAt(i), newBar));		
	    }
	
	return text.toString();
    }
    
    /**
     * Tests out the parser.
     *
     * @param args should contain nothing
     */
    public static void main(String args[]){
	StringBuffer sb = new StringBuffer();
	sb.append("node1.node1_1.node1_1_1.node1_1_1_1, ");
	sb.append("node1.node1_1.node1_1_1.node1_1_1_2, ");
	sb.append("node1.node1_1.node1_1_1.node1_1_1_3, ");
	sb.append("node1.node1_1.node1_1_2.node1_1_2_1, ");
	sb.append("node1.node1_1.node1_1_3.node1_1_3_1, ");
	sb.append("node1.node1_2.node1_2_1.node1_2_1_1, ");
	sb.append("node1.node1_2.node1_2_3.node1_2_3_1, ");
	sb.append("node1.node1_3.node1_3_3.node1_3_3_1, ");
	sb.append("node1.node1_3.node1_3_3.node1_3_3_2, ");

	String p = sb.toString();
	try{
	    HierarchyPropertyParser hpp = new HierarchyPropertyParser(p, ", ");
	    System.out.println("seperator: "+hpp.getSeperator());
	    System.out.println("depth: "+hpp.depth());
	    System.out.println("The tree:\n\n"+hpp.showTree());
	    hpp.goToRoot();
	    System.out.println("goto: "+hpp.goTo("node1.node1_2.node1_2_1")
			       +": "+hpp.getValue()+" | "+hpp.fullValue()+
			       " leaf? "+ hpp.isLeafReached());
	    System.out.println("go down(wrong): "+hpp.goDown("node1"));
	    System.out.println("Stay still? "+hpp.getValue());
	    System.out.println("go to child: "+hpp.goToChild("node1_2_1_1")
			       +": "+hpp.getValue()+" | "+hpp.fullValue()+
			       " leaf? "+ hpp.isLeafReached()
			       +" root? "+ hpp.isRootReached());
	    System.out.println("parent: "+hpp.parentValue());
	    System.out.println("level: "+hpp.getLevel());
	    System.out.println("context: "+hpp.context());
	    hpp.goToRoot();
	    System.out.println("After gotoRoot. leaf? "+ hpp.isLeafReached()
			       +" root? "+ hpp.isRootReached());
	    System.out.println("Go down(correct): "+
			       hpp.goDown("node1_1.node1_1_1")+
			       " value: "+hpp.getValue()+" | "+hpp.fullValue()
			       +" level: "+hpp.getLevel()
			       +" leaf? "+ hpp.isLeafReached()
			       +" root? "+ hpp.isRootReached());
	    hpp.goToParent();
	    System.out.println("value: "+hpp.getValue()+" | "+hpp.fullValue());
	    System.out.println("level: "+hpp.getLevel());
	    
	    String[] chd = hpp.childrenValues();
	    for(int i=0; i < chd.length; i++){
		System.out.print("children "+i+": "+chd[i]);
		hpp.goDown(chd[i]);
		System.out.println("real value: "+hpp.getValue()+" | "+
				   hpp.fullValue()+
				   "(level: "+hpp.getLevel()+")");
		hpp.goToParent();
	    }
	    
	    System.out.println("Another way to go to root:"+hpp.goTo("node1")
			       +": "+hpp.getValue()+" | "+hpp.fullValue()); 
   	}catch(Exception e){
	    System.out.println(e.getMessage());
	    e.printStackTrace();
	}
    }
}
